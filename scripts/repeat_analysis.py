#!/usr/bin/env python
# coding: utf-8

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def window_density(track: np.ndarray, center: int, size: int) -> float:
    start = max(0, center - size // 2)
    end = min(len(track), center + size // 2)
    return float(track[start:end].mean())

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run repeat density analysis and ztf-6 motif distribution plotting."
    )

    parser.add_argument(
        "--rmsk-path",
        type=Path,
        required=True,
        help="Path to RepeatMasker table (rmsk.txt).",
    )
    parser.add_argument(
        "--repeat-summary-out",
        type=Path,
        required=True,
        help="Output TSV path for repeat density summary.",
    )
    parser.add_argument(
        "--repeat-fig-dir",
        type=Path,
        required=True,
        help="Output directory for repeat density PDFs.",
    )

    parser.add_argument(
        "--fimo-path",
        type=Path,
        required=True,
        help="Path to FIMO TSV file for ztf-6 motif hits.",
    )
    parser.add_argument(
        "--motif-fig-out",
        type=Path,
        required=True,
        help="Output PDF path for ztf-6 motif distribution plot.",
    )

    parser.add_argument(
        "--roll-window",
        type=int,
        default=150_000,
        help="Rolling window size in bp for repeat density.",
    )
    parser.add_argument(
        "--roi-centers",
        type=int,
        nargs="+",
        default=[18_100_000, 19_000_000, 19_750_000],
        help="ROI centers in bp for repeat density percentile calculation.",
    )
    parser.add_argument(
        "--random-n",
        type=int,
        default=10_000,
        help="Number of random windows for percentile background.",
    )
    parser.add_argument(
        "--subsampling-param",
        type=int,
        default=500,
        help="Downsampling step for plotting rolling density.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed.",
    )

    parser.add_argument(
        "--motif-b1",
        type=int,
        default=19_000_000,
        help="First motif ROI center in bp.",
    )
    parser.add_argument(
        "--motif-b2",
        type=int,
        default=19_750_000,
        help="Second motif ROI center in bp.",
    )
    parser.add_argument(
        "--motif-window",
        type=int,
        default=75_000,
        help="Half-window size for motif ROI shading in bp.",
    )
    parser.add_argument(
        "--motif-bin-size",
        type=int,
        default=100_000,
        help="Histogram bin size in bp for motif plot.",
    )
    parser.add_argument(
        "--chr-length",
        type=int,
        default=20_924_180,
        help="Chromosome length in bp for motif plot axis.",
    )

    return parser.parse_args()

def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    args.repeat_fig_dir.mkdir(parents=True, exist_ok=True)
    args.repeat_summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.motif_fig_out.parent.mkdir(parents=True, exist_ok=True)

    cols = [
        "bin",
        "swScore",
        "milliDiv",
        "milliDel",
        "milliIns",
        "chrom",
        "start",
        "end",
        "genoLeft",
        "strand",
        "repName",
        "repClass",
        "repFamily",
        "repStart",
        "repEnd",
        "repLeft",
        "id",
    ]

    rmsk = pd.read_csv(
        args.rmsk_path,
        sep="\t",
        names=cols,
        usecols=["chrom", "start", "end", "repClass"],
    )

    chr_v = rmsk[rmsk.chrom == "chrV"].copy()
    chr_v_len = int(chr_v.end.max())

    classes = ["TOTAL"] + sorted(chr_v.repClass.unique())
    repclasses_to_save = ["TOTAL", "RC", "Satellite", "Simple_repeat", "Unknown"]

    results = []

    for repclass in classes:
        if repclass == "TOTAL":
            df = chr_v
        else:
            df = chr_v[chr_v.repClass == repclass]

        n_repeats = len(df)
        mean_length = (df["end"] - df["start"]).mean() if n_repeats > 0 else 0.0

        print(
            f"\\n=== {repclass} ({n_repeats} intervals, mean length: {mean_length:.1f} bp) ==="
        )

        track = np.zeros(chr_v_len, dtype=np.uint8)
        for start, end in df[["start", "end"]].values:
            track[start:end] = 1

        csum = np.cumsum(track, dtype=np.int64)
        rolling = (csum[args.roll_window:] - csum[:-args.roll_window]) / args.roll_window

        edge_trim = args.roll_window // 2
        valid_start = edge_trim
        valid_end = chr_v_len - edge_trim

        indices = np.arange(0, len(rolling), args.subsampling_param)
        rolling_ds = rolling[indices]
        xvals = (valid_start + indices) / 1e6
        mean_density = float(rolling.mean())

        plt.figure(figsize=(15, 5))
        plt.plot(xvals, rolling_ds, color="blue", linewidth=1, label="Rolling density")

        for roi in args.roi_centers:
            plt.axvline(
                roi / 1e6,
                color="green",
                linestyle="--",
                linewidth=2,
                label="ROI" if roi == args.roi_centers[0] else "",
            )

        plt.axhline(mean_density, color="red", linestyle="--", label="Mean density")
        plt.ylabel("Repeat density", fontsize=12)
        plt.xlabel("chrV position (Mb)", fontsize=12)
        plt.title(f"{repclass} repeat density ({args.roll_window // 1000} kb rolling)", fontsize=13)
        plt.xlim(valid_start / 1e6, valid_end / 1e6)
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()

        if repclass in repclasses_to_save:
            out_pdf = args.repeat_fig_dir / f"{repclass}_repeat_density.pdf"
            plt.savefig(out_pdf)

        plt.show()

        roi_vals = [window_density(track, center, args.roll_window) for center in args.roi_centers]
        random_pos = rng.integers(valid_start, valid_end, args.random_n)
        rand_density = np.array(
            [window_density(track, pos, args.roll_window) for pos in random_pos]
        )

        print(f"  Genome mean: {rand_density.mean():.6f}")
        for i, val in enumerate(roi_vals):
            perc = (rand_density < val).mean()
            print(f"  ROI{i + 1} density: {val:.6f}, percentile: {perc:.4f}")

        row = {
            "class": repclass,
            "n_repeats": n_repeats,
            "mean_length": mean_length,
            "genome_mean": float(rand_density.mean()),
        }

        for i, val in enumerate(roi_vals):
            row[f"ROI{i + 1}_density"] = float(val)
            row[f"ROI{i + 1}_percentile"] = float((rand_density < val).mean())

        results.append(row)

    res_table = pd.DataFrame(results)

    print("\n" + "=" * 80)
    print("Repeat density summary:")
    print("=" * 80 + "\n")
    print(res_table.to_string(index=False))

    res_table.to_csv(args.repeat_summary_out, sep="\t", index=False)
    print(f"\nSaved to: {args.repeat_summary_out}")

    print("\n\nGenerating ztf-6 motif distribution plot...")

    df_motif = pd.read_csv(args.fimo_path, sep="\t")
    pos_motif = df_motif["start"]

    n_bins = args.chr_length // args.motif_bin_size + 1

    plt.figure(figsize=(14, 4))
    plt.hist(pos_motif, bins=n_bins, histtype="step", color="green", label="ztf-6 motif hits")

    plt.axvspan(
        args.motif_b1 - args.motif_window,
        args.motif_b1 + args.motif_window,
        facecolor="red",
        alpha=0.2,
        label="ROI around 19.0 Mb",
    )
    plt.axvspan(
        args.motif_b2 - args.motif_window,
        args.motif_b2 + args.motif_window,
        facecolor="blue",
        alpha=0.2,
        label="ROI around 19.75 Mb",
    )

    minor_ticks = np.arange(0, args.chr_length + 1, 200_000)
    major_ticks = np.arange(0, args.chr_length + 1, 1_000_000)
    plt.xticks(major_ticks, labels=[str(int(x / 1e6)) for x in major_ticks])
    plt.gca().set_xticks(minor_ticks, minor=True)

    plt.xlim(-0.1, args.chr_length)
    plt.xlabel("chrV genomic position (Mb)")
    plt.ylabel(f"ztf-6 motif hits in a {args.motif_bin_size // 1000} kb window")
    plt.title("Distribution of ztf-6 motif matches across chrV")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.motif_fig_out)
    plt.show()

    print(f"Saved ztf-6 plot to: {args.motif_fig_out}")

if __name__ == "__main__":
    main()
