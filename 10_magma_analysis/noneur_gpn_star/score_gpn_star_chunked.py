#!/usr/bin/env python
"""
Optional, step 2 of 3: GPN-STAR variant effect scores for one chromosome's SNVs, computed in chunks with a
progress line per chunk. Each chunk runs the GPN-STAR command line tool
    python -m gpn.star.inference vep <input> <MSA> <WINDOW> <MODEL> <output> --is_file ...
(https://github.com/songlab-cal/gpn) and is checkpointed, so a rerun resumes after the last finished
chunk. Needs a GPU, the gpn package, the GPN-STAR model and its multiple sequence alignment.

usage:   python score_gpn_star_chunked.py IN.parquet OUT.parquet MSA WINDOW MODEL
             [--chunk 100000] [--batch 16] [--workers 8] [--workdir DIR]
         IN:  variant_id, chrom (without "chr"), pos, ref, alt (aggregate_noneur_variants.py)
         OUT: variant_id, gpn_star_llr, in the input order
"""
import argparse
import math
import os
import subprocess
import sys
import time

import pandas as pd


def fmt(s):
    s = int(s)
    h, r = divmod(s, 3600)
    m, sec = divmod(r, 60)
    return f"{h}h{m:02d}m" if h else (f"{m}m{sec:02d}s" if m else f"{sec}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("msa")
    ap.add_argument("win", type=int)
    ap.add_argument("model")
    ap.add_argument("--chunk", type=int, default=100_000)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--workdir", default=None)
    a = ap.parse_args()

    tag = os.path.basename(a.inp).replace(".parquet", "")
    cols = ["variant_id", "chrom", "pos", "ref", "alt"]
    df = pd.read_parquet(a.inp, columns=cols)
    N = len(df)
    workdir = a.workdir or os.path.join(os.path.dirname(a.out) or ".", "work", tag)
    os.makedirs(workdir, exist_ok=True)
    print(f"[{tag}] {N:,} SNVs | chunk={a.chunk:,} | workdir={workdir}", flush=True)
    if N == 0:
        pd.DataFrame({"variant_id": [], "gpn_star_llr": []}).to_parquet(a.out, index=False)
        return

    nch = math.ceil(N / a.chunk)
    t_run, durs, done = time.monotonic(), [], 0
    for k in range(nch):
        lo, hi = k * a.chunk, min((k + 1) * a.chunk, N)
        ckpt = f"{workdir}/chunk_{k:04d}.parquet"
        if os.path.exists(ckpt):
            try:
                if len(pd.read_parquet(ckpt, columns=["variant_id"])) == hi - lo:
                    done += hi - lo
                    print(f"[{tag}] chunk {k + 1}/{nch} already scored", flush=True)
                    continue
            except Exception:
                pass
        cin, cout = f"{workdir}/_in_{k:04d}.parquet", f"{workdir}/_sc_{k:04d}.parquet"
        df.iloc[lo:hi][cols].reset_index(drop=True).to_parquet(cin, index=False)
        t0 = time.monotonic()
        subprocess.run([sys.executable, "-m", "gpn.star.inference", "vep", cin, a.msa, str(a.win),
                        a.model, cout, "--is_file", "--per_device_batch_size", str(a.batch),
                        "--dataloader_num_workers", str(a.workers)], check=True)
        dt = time.monotonic() - t0
        durs.append(dt)
        s = pd.read_parquet(cout)
        sc = s["score"] if "score" in s.columns else s.iloc[:, 0]
        assert len(sc) == hi - lo, f"row mismatch chunk {k}: {len(sc)} vs {hi - lo}"
        pd.DataFrame({"variant_id": df["variant_id"].values[lo:hi],
                      "gpn_star_llr": sc.values}).to_parquet(ckpt, index=False)
        os.remove(cin)
        os.remove(cout)
        done += hi - lo
        eta = (nch - (k + 1)) * sum(durs) / len(durs)
        print(f"[{tag}] chunk {k + 1}/{nch} | {done:,}/{N:,} ({100 * done / N:.1f}%) | {dt:.0f}s "
              f"({(hi - lo) / dt:.0f} SNV/s) | elapsed {fmt(time.monotonic() - t_run)} | ETA {fmt(eta)}", flush=True)

    res = pd.concat([pd.read_parquet(f"{workdir}/chunk_{k:04d}.parquet") for k in range(nch)], ignore_index=True)
    assert len(res) == N, f"final row mismatch: {len(res)} vs {N}"
    res.to_parquet(a.out, index=False)
    print(f"[{tag}] {N:,} SNVs -> {a.out} in {fmt(time.monotonic() - t_run)}", flush=True)


if __name__ == "__main__":
    main()
