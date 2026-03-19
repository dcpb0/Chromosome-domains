#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib import font_manager
import matplotlib as mpl
import matplotlib.transforms as transforms
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot condensin/TF/histone tracks and expression heatmaps."
    )
    parser.add_argument(
        "--condensin-dir",
        type=Path,
        required=True,
        help="Directory with condensin wig/bed files.",
    )
    parser.add_argument(
        "--tf-dir",
        type=Path,
        required=True,
        help="Directory with TF bedgraph and peak files.",
    )
    parser.add_argument(
        "--histone-dir",
        type=Path,
        required=True,
        help="Directory with histone bedgraph files.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Root directory for supporting input data.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Root directory for analysis result tables.",
    )
    parser.add_argument(
        "--font-ttf",
        type=Path,
        default=None,
        help="Optional TTF font path.",
    )
    parser.add_argument(
        "--hot-sites-bed",
        type=Path,
        required=True,
        help="BED file with HOT site midpoints.",
    )
    parser.add_argument(
        "--active-domains-bed",
        type=Path,
        required=True,
        help="BED file with active domains.",
    )
    parser.add_argument(
        "--all-domains-bed",
        type=Path,
        required=True,
        help="BED-like file with all chromatin domains.",
    )
    return parser.parse_args()

def chrom_match(c_file, chrom_query):
    if c_file == chrom_query:
        return True
    if c_file == chrom_query.replace("chr", ""):
        return True
    if "chr" + c_file == chrom_query:
        return True
    return False

def smooth_signal(x, y, window_bp):
    step = int(np.median(np.diff(x)))
    window = max(1, int(window_bp / step))
    kernel = np.ones(window) / window
    y_smooth = np.convolve(y, kernel, mode="same")
    return y_smooth

def load_bedgraph_region(path, chrom, region_start=None, region_end=None):
    x_list = []
    y_list = []
    with open(path) as f:
        for line in f:
            if line.startswith(("track", "#")):
                continue
            c, start, end, val = line.strip().split()
            start = int(start)
            end = int(end)
            val = float(val)
            if not chrom_match(c, chrom):
                continue
            mid = (start + end) // 2
            if region_start is not None and mid < region_start:
                continue
            if region_end is not None and mid > region_end:
                continue
            x_list.append(mid)
            y_list.append(val)
    return np.array(x_list), np.array(y_list)

def load_wig_region(filename, chrom, start, end):
    positions = []
    values = []
    keep = False

    with open(filename) as f:
        for line in f:
            if line.startswith("variableStep"):
                keep = f"chrom={chrom}" in line
                continue

            if keep:
                pos, val = line.split()
                pos = int(pos)

                if pos > end:
                    break

                if start <= pos <= end:
                    positions.append(pos)
                    values.append(float(val))

    return np.array(positions), np.array(values)

def load_bed_region(filename, chrom, start, end):
    peaks = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if not chrom_match(fields[0], chrom):
                continue

            peak_start = int(fields[1])
            peak_end = int(fields[2])

            if peak_end < start:
                continue
            if peak_start > end:
                continue

            peaks.append((peak_start, peak_end))

    return peaks

def main():
    args = parse_args()
    if args.font_ttf is not None and args.font_ttf.exists():
        font_manager.fontManager.addfont(str(args.font_ttf))
    base = args.condensin_dir
    tf_base = args.tf_dir
    
    files = [
        "GSE45678_CAPG1_N2_MxEmb",
        "GSE45678_DPY26_N2_MxEmb",
        "GSE45678_DPY28_N2_MxEmb",
        "GSE45678_SMC4_N2_MxEmb",
        "GSE45678_MIX1_N2_MxEmb"
    ]
    
    tf_signal_files = {
        "CEH-30": tf_base
        / "CEH-30 Combined -GFP ChIP- recalled peaks w peaks-V-1..20924180.bedgraph",
        "LIN-13": tf_base
        / "LIN-13 Combined -GFP ChIP-- LIN13_GFP_emb w peaks-V-1..20924180.bedgraph",
        "MEP-1": tf_base / "MEP-1 Combined -GFP ChIP- w peaks-V-1..20924180.bedgraph",
        "CEH-26": tf_base / "CEH-26 Combined -GFP ChIP- w peaks-V-1..20924180.bedgraph",
        "HLH-1": tf_base / "HLH-1 Combined -GFP ChIP- w peaks-V-1..20924180.bedgraph",
    }
    
    tf_peak_files = {
        "CEH-30": tf_base / "CEH-30_peaks.bed",
        "LIN-13": tf_base / "LIN-13_peaks.bed",
        "MEP-1": tf_base / "MEP-1_peaks.bed",
        "CEH-26": tf_base / "CEH-26_peaks.bed",
        "HLH-1": tf_base / "HLH-1_peaks.bed",
    }

    

    # ==============================
    # Parameters
    # ==============================
    chrom = "chrV"
    region_start = 19_500_000
    region_end = 19_800_000
    
    # ==============================
    # Preload Condensin data
    # ==============================
    condensin_data = {}
    for prefix in files:
        wig = base / f"{prefix}_average_var_chr.wig"
        bed = base / f"{prefix}_peaks_chr.bed"
    
        x, y = load_wig_region(wig, chrom, region_start, region_end)
        peaks = load_bed_region(bed, chrom, region_start, region_end)
    
        condensin_data[prefix] = {"x": x, "y": y, "peaks": peaks}
    
    # ==============================
    # Preload TF data
    # ==============================
    tf_data = {}
    for name in tf_signal_files.keys():
        sig_path = tf_signal_files[name]
        peak_path = tf_peak_files[name]
    
        x, y = load_bedgraph_region(sig_path, chrom, region_start, region_end)
        peaks = load_bed_region(peak_path, chrom, region_start, region_end)
    
        tf_data[name] = {"x": x, "y": y, "peaks": peaks}


    # ==============================
    # Parameters
    # ==============================
    condensin_color = "#08306b"
    tf_color = "#4292c6"
    condensin_peak_edge = "crimson"
    tf_peak_edge = "darkgreen"
    
    LABEL_X_POS = -0.05
    
    fig_width, fig_height = 23.38, 16.54
    
    plt.rcParams.update({"font.family": "Arial", "font.size": 12})
    
    all_tracks = list(condensin_data.keys()) + list(tf_data.keys())
    n_tracks = len(all_tracks)
    
    fig, axes = plt.subplots(
        n_tracks, 1, figsize=(fig_width, fig_height), dpi=200, sharex=True
    )
    
    # ------------------------------
    # Condensin tracks
    # ------------------------------
    for i, prefix in enumerate(condensin_data.keys()):
        ax = axes[i]
        data = condensin_data[prefix]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        # y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y, color=condensin_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(
                Rectangle(
                    (start, stripe_bottom),
                    end - start,
                    stripe_height,
                    facecolor="white",
                    edgecolor=condensin_peak_edge,
                    linewidth=1.0,
                    clip_on=False,
                    zorder=10,
                )
            )
    
        ax.set_ylabel(prefix.split("_")[1], rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    # ------------------------------
    # TF tracks
    # ------------------------------
    offset = len(condensin_data)
    for j, name in enumerate(tf_data.keys()):
        ax = axes[offset + j]
        data = tf_data[name]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        # y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y, color=tf_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(
                Rectangle(
                    (start, stripe_bottom),
                    end - start,
                    stripe_height,
                    facecolor="white",
                    edgecolor=tf_peak_edge,
                    linewidth=0.8,
                    clip_on=False,
                    zorder=10,
                )
            )
    
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    axes[0].xaxis.set_label_position("top")
    axes[0].xaxis.tick_top()
    axes[0].tick_params(axis="x", which="both", top=True, labeltop=True, bottom=False)
    
    tick_step = 20_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    axes[0].set_xticks(xticks)
    axes[0].set_xticklabels([f"{x/1e6:.2f}" for x in xticks])
    axes[0].set_title(
        f"{chrom}, {region_start/1000000}-{region_end/1000000} Mb, non-smoothed signal",
        pad=35,
    )
    
    for ax in axes[1:]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    
    plt.xlim(region_start, region_end)
    
    plt.subplots_adjust(hspace=0.5, left=0.15, right=0.95)
    
    legend_elements = [
        Line2D([0], [0], color=condensin_color, lw=2, label="Condensin signal"),
        Rectangle(
            (0, 0),
            1,
            1,
            facecolor="white",
            edgecolor=condensin_peak_edge,
            label="Condensin peaks",
        ),
        Line2D([0], [0], color=tf_color, lw=2, label="TF signal"),
        Rectangle(
            (0, 0), 1, 1, facecolor="white", edgecolor=tf_peak_edge, label="TF peaks"
        ),
    ]
    
    fig.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=False,
        bbox_to_anchor=(0.98, 0.98),
    )
    plt.savefig("chip_tracks_condensin_tfs_peaks_zoomed_19mb_right.pdf")
    plt.show()


    # ==============================
    # Parameters
    # ==============================
    condensin_peak_edge = "crimson"
    tf_color = "#4292c6"
    tf_peak_edge = "darkgreen"
    LABEL_X_POS = -0.05
    
    fig_width, fig_height = 23.38, 16.54
    plt.rcParams.update({"font.family": "Arial", "font.size": 12})
    
    n_condensin = len(condensin_data)
    n_tf = len(tf_data)
    
    height_ratios = [0.14]*n_condensin + [1]*n_tf
    
    fig, axes = plt.subplots(
        n_condensin + n_tf,
        1,
        figsize=(fig_width, fig_height),
        dpi=200,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios, 'hspace': 0.4}
    )
    
    # ===============================
    # Condensin peak tracks (short)
    # ===============================
    for i, prefix in enumerate(condensin_data.keys()):
        ax = axes[i]
        peaks = condensin_data[prefix]["peaks"]
    
        y_bottom = 0.1
        y_height = 0.8
        for start, end in peaks:
            ax.add_patch(Rectangle((start, y_bottom), end - start, y_height,
                                   facecolor="white", edgecolor=condensin_peak_edge,
                                   linewidth=1.0, clip_on=False, zorder=10))
    
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel(prefix.split("_")[1], rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
    # ===============================
    # TF tracks
    # ===============================
    offset = n_condensin
    for j, name in enumerate(tf_data.keys()):
        ax = axes[offset + j]
        data = tf_data[name]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        ax.plot(x, y, color=tf_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(Rectangle((start, stripe_bottom), end - start, stripe_height,
                                   facecolor="white", edgecolor=tf_peak_edge,
                                   linewidth=0.8, clip_on=False, zorder=10))
    
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # ===============================
    # X-axis formatting
    # ===============================
    
    tick_step = 20_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels([f"{x/1e6:.2f}" for x in xticks])
    axes[-1].tick_params(axis="x", which="both", top=False, bottom=True, labelbottom=True, pad=20)
    axes[-1].set_xlabel("Position (Mb)", labelpad=15)
    
    fig.suptitle(f"{chrom}, {region_start/1e6}-{region_end/1e6} Mb, non-smoothed signal",
                 fontsize=14, y=0.995)
    
    plt.subplots_adjust(hspace=0.5, left=0.15, right=0.95, bottom=0.18)
    plt.xlim(region_start, region_end)
    
    # ===============================
    # Legend
    # ===============================
    legend_elements = [
        Rectangle((0,0),1,1,facecolor='white', edgecolor=condensin_peak_edge, label='Condensin peaks'),
        Line2D([0],[0], color=tf_color, lw=2, label='TF signal'),
        Rectangle((0,0),1,1,facecolor='white', edgecolor=tf_peak_edge, label='TF peaks'),
    ]
    
    fig.legend(handles=legend_elements, loc='upper right', frameon=False, bbox_to_anchor=(0.98,0.98))
    
    plt.savefig("chip_tracks_condensin_tfs_peaks_cut_zoomed_19mb_right.pdf")
    plt.show()

    

    # ==============================
    # Parameters
    # ==============================
    region_start = 18_900_000
    region_end = 19_100_000
    
    # ==============================
    # Preload Condensin data
    # ==============================
    condensin_data = {}
    for prefix in files:
        wig = base / f"{prefix}_average_var_chr.wig"
        bed = base / f"{prefix}_peaks_chr.bed"
    
        x, y = load_wig_region(wig, chrom, region_start, region_end)
        peaks = load_bed_region(bed, chrom, region_start, region_end)
    
        condensin_data[prefix] = {"x": x, "y": y, "peaks": peaks}
    
    # ==============================
    # Preload TF data
    # ==============================
    tf_data = {}
    for name in tf_signal_files.keys():
        sig_path = tf_signal_files[name]
        peak_path = tf_peak_files[name]
    
        x, y = load_bedgraph_region(sig_path, chrom, region_start, region_end)
        peaks = load_bed_region(peak_path, chrom, region_start, region_end)
    
        # optionally: chromosome-wide z-score for TFs
        # if len(y) > 0:
        #     mean = np.mean(y)
        #     std = np.std(y)
        #     y = (y - mean) / std if std > 0 else y
    
        tf_data[name] = {"x": x, "y": y, "peaks": peaks}


    # ==============================
    # Parameters
    # ==============================
    condensin_color = "#08306b"
    tf_color = "#4292c6"
    condensin_peak_edge = "crimson"
    tf_peak_edge = "darkgreen"
    
    LABEL_X_POS = -0.05
    
    fig_width, fig_height = 23.38, 16.54
    
    plt.rcParams.update({"font.family": "Arial", "font.size": 12})
    
    all_tracks = list(condensin_data.keys()) + list(tf_data.keys())
    n_tracks = len(all_tracks)
    
    fig, axes = plt.subplots(
        n_tracks, 1, figsize=(fig_width, fig_height), dpi=200, sharex=True
    )
    
    # ------------------------------
    # Condensin tracks
    # ------------------------------
    for i, prefix in enumerate(condensin_data.keys()):
        ax = axes[i]
        data = condensin_data[prefix]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        # y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y, color=condensin_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(
                Rectangle(
                    (start, stripe_bottom),
                    end - start,
                    stripe_height,
                    facecolor="white",
                    edgecolor=condensin_peak_edge,
                    linewidth=1.0,
                    clip_on=False,
                    zorder=10,
                )
            )
    
        ax.set_ylabel(prefix.split("_")[1], rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    # ------------------------------
    # TF tracks
    # ------------------------------
    offset = len(condensin_data)
    for j, name in enumerate(tf_data.keys()):
        ax = axes[offset + j]
        data = tf_data[name]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        # y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y, color=tf_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(
                Rectangle(
                    (start, stripe_bottom),
                    end - start,
                    stripe_height,
                    facecolor="white",
                    edgecolor=tf_peak_edge,
                    linewidth=0.8,
                    clip_on=False,
                    zorder=10,
                )
            )
    
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    axes[0].xaxis.set_label_position("top")
    axes[0].xaxis.tick_top()
    axes[0].tick_params(axis="x", which="both", top=True, labeltop=True, bottom=False)
    
    tick_step = 20_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    axes[0].set_xticks(xticks)
    axes[0].set_xticklabels([f"{x/1e6:.2f}" for x in xticks])
    axes[0].set_title(
        f"{chrom}, {region_start/1000000}-{region_end/1000000} Mb, non-smoothed signal",
        pad=35,
    )
    
    for ax in axes[1:]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    
    plt.xlim(region_start, region_end)
    
    plt.subplots_adjust(hspace=0.5, left=0.15, right=0.95)
    
    legend_elements = [
        Line2D([0], [0], color=condensin_color, lw=2, label="Condensin signal"),
        Rectangle(
            (0, 0),
            1,
            1,
            facecolor="white",
            edgecolor=condensin_peak_edge,
            label="Condensin peaks",
        ),
        Line2D([0], [0], color=tf_color, lw=2, label="TF signal"),
        Rectangle(
            (0, 0), 1, 1, facecolor="white", edgecolor=tf_peak_edge, label="TF peaks"
        ),
    ]
    
    fig.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=False,
        bbox_to_anchor=(0.98, 0.98),
    )
    plt.savefig("chip_tracks_condensin_tfs_peaks_zoomed_19mb_left.pdf")
    plt.show()


    # ==============================
    # Parameters
    # ==============================
    condensin_peak_edge = "crimson"
    tf_color = "#4292c6"
    tf_peak_edge = "darkgreen"
    LABEL_X_POS = -0.05
    
    fig_width, fig_height = 23.38, 16.54
    plt.rcParams.update({"font.family": "Arial", "font.size": 12})
    
    n_condensin = len(condensin_data)
    n_tf = len(tf_data)
    
    height_ratios = [0.14]*n_condensin + [1]*n_tf
    
    fig, axes = plt.subplots(
        n_condensin + n_tf,
        1,
        figsize=(fig_width, fig_height),
        dpi=200,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios, 'hspace': 0.4}
    )
    
    # ===============================
    # Condensin peak tracks
    # ===============================
    for i, prefix in enumerate(condensin_data.keys()):
        ax = axes[i]
        peaks = condensin_data[prefix]["peaks"]
    
        # Use fixed small y values for peaks (0–1)
        y_bottom = 0.1
        y_height = 0.8
        for start, end in peaks:
            ax.add_patch(Rectangle((start, y_bottom), end - start, y_height,
                                   facecolor="white", edgecolor=condensin_peak_edge,
                                   linewidth=1.0, clip_on=False, zorder=10))
    
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel(prefix.split("_")[1], rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
    # ===============================
    # TF tracks
    # ===============================
    offset = n_condensin
    for j, name in enumerate(tf_data.keys()):
        ax = axes[offset + j]
        data = tf_data[name]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        ax.plot(x, y, color=tf_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(Rectangle((start, stripe_bottom), end - start, stripe_height,
                                   facecolor="white", edgecolor=tf_peak_edge,
                                   linewidth=0.8, clip_on=False, zorder=10))
    
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # ===============================
    # X-axis formatting
    # ===============================
    
    tick_step = 20_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels([f"{x/1e6:.2f}" for x in xticks])
    axes[-1].tick_params(axis="x", which="both", top=False, bottom=True, labelbottom=True, pad=20)
    axes[-1].set_xlabel("Position (Mb)", labelpad=15)
    
    fig.suptitle(f"{chrom}, {region_start/1e6}-{region_end/1e6} Mb, non-smoothed signal",
                 fontsize=14, y=0.995)
    
    plt.subplots_adjust(hspace=0.5, left=0.15, right=0.95, bottom=0.18)
    plt.xlim(region_start, region_end)
    
    # ===============================
    # Legend
    # ===============================
    legend_elements = [
        Rectangle((0,0),1,1,facecolor='white', edgecolor=condensin_peak_edge, label='Condensin peaks'),
        Line2D([0],[0], color=tf_color, lw=2, label='TF signal'),
        Rectangle((0,0),1,1,facecolor='white', edgecolor=tf_peak_edge, label='TF peaks'),
    ]
    
    fig.legend(handles=legend_elements, loc='upper right', frameon=False, bbox_to_anchor=(0.98,0.98))
    
    plt.savefig("chip_tracks_condensin_tfs_peaks_cut_zoomed_19mb_left.pdf")
    plt.show()

    

    # ==============================
    # Parameters
    # ==============================
    region_start = 16_660_000
    region_end = 20_910_000
    
    # ==============================
    # Preload Condensin data
    # ==============================
    condensin_data = {}
    for prefix in files:
        wig = base / f"{prefix}_average_var_chr.wig"
        bed = base / f"{prefix}_peaks_chr.bed"
    
        x, y = load_wig_region(wig, chrom, region_start, region_end)
        peaks = load_bed_region(bed, chrom, region_start, region_end)
    
        condensin_data[prefix] = {"x": x, "y": y, "peaks": peaks}
    
    # ==============================
    # Preload TF data
    # ==============================
    tf_data = {}
    for name in tf_signal_files.keys():
        sig_path = tf_signal_files[name]
        peak_path = tf_peak_files[name]
    
        x, y = load_bedgraph_region(sig_path, chrom, region_start, region_end)
        peaks = load_bed_region(peak_path, chrom, region_start, region_end)
    
        # optionally: chromosome-wide z-score for TFs
        # if len(y) > 0:
        #     mean = np.mean(y)
        #     std = np.std(y)
        #     y = (y - mean) / std if std > 0 else y
    
        tf_data[name] = {"x": x, "y": y, "peaks": peaks}


    # ==============================
    # Parameters
    # ==============================
    window_bp = 2_000
    condensin_color = "#08306b"
    tf_color = "#4292c6"
    condensin_peak_edge = "crimson"
    tf_peak_edge = "darkgreen"
    
    LABEL_X_POS = -0.05
    
    fig_width, fig_height = 23.38, 16.54
    
    plt.rcParams.update({"font.family": "Arial", "font.size": 12})
    
    all_tracks = list(condensin_data.keys()) + list(tf_data.keys())
    n_tracks = len(all_tracks)
    
    fig, axes = plt.subplots(
        n_tracks, 1, figsize=(fig_width, fig_height), dpi=200, sharex=True
    )
    
    # ------------------------------
    # Condensin tracks
    # ------------------------------
    for i, prefix in enumerate(condensin_data.keys()):
        ax = axes[i]
        data = condensin_data[prefix]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y_smooth, color=condensin_color, linewidth=1.8)
    
        # Calculate peak positions
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(
                Rectangle(
                    (start, stripe_bottom),
                    end - start,
                    stripe_height,
                    facecolor="white",
                    edgecolor=condensin_peak_edge,
                    linewidth=1.0,
                    clip_on=False,
                    zorder=10,
                )
            )
    
        ax.set_ylabel(prefix.split("_")[1], rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    # ------------------------------
    # TF tracks
    # ------------------------------
    offset = len(condensin_data)
    for j, name in enumerate(tf_data.keys()):
        ax = axes[offset + j]
        data = tf_data[name]
        x, y, peaks = data["x"], data["y"], data["peaks"]
    
        y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y_smooth, color=tf_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(
                Rectangle(
                    (start, stripe_bottom),
                    end - start,
                    stripe_height,
                    facecolor="white",
                    edgecolor=tf_peak_edge,
                    linewidth=0.8,
                    clip_on=False,
                    zorder=10,
                )
            )
    
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    # ------------------------------
    # Top axis formatting
    # ------------------------------
    axes[0].xaxis.set_label_position("top")
    axes[0].xaxis.tick_top()
    axes[0].tick_params(axis="x", which="both", top=True, labeltop=True, bottom=False)
    
    tick_step = 200_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    axes[0].set_xticks(xticks)
    axes[0].set_xticklabels([f"{x/1e6:.1f}" for x in xticks])
    axes[0].set_title(
        f"{chrom}, {region_start/1000000}-{region_end/1000000} Mb | rolling window = {window_bp/1000:.0f} kb",
        pad=35,
    )
    
    for ax in axes[1:]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    
    plt.xlim(region_start, region_end)
    
    plt.subplots_adjust(hspace=0.5, left=0.15, right=0.95)
    
    legend_elements = [
        Line2D([0], [0], color=condensin_color, lw=2, label="Condensin signal"),
        Rectangle(
            (0, 0),
            1,
            1,
            facecolor="white",
            edgecolor=condensin_peak_edge,
            label="Condensin peaks",
        ),
        Line2D([0], [0], color=tf_color, lw=2, label="TF signal"),
        Rectangle(
            (0, 0), 1, 1, facecolor="white", edgecolor=tf_peak_edge, label="TF peaks"
        ),
    ]
    
    fig.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=False,
        bbox_to_anchor=(0.98, 0.98),
    )
    
    plt.savefig("chip_tracks_condensin_tfs_peaks_16_21mb.pdf")
    plt.show()


    # ==============================
    # Parameters
    # ==============================
    condensin_peak_edge = "crimson"
    tf_color = "#4292c6"
    tf_peak_edge = "darkgreen"
    LABEL_X_POS = -0.05
    window_bp = 2_000

    fig_width, fig_height = 23.38, 16.54
    plt.rcParams.update({"font.family": "Arial", "font.size": 12})
    
    n_condensin = len(condensin_data)
    n_tf = len(tf_data)
    
    height_ratios = [0.14]*n_condensin + [1]*n_tf
    
    fig, axes = plt.subplots(
        n_condensin + n_tf,
        1,
        figsize=(fig_width, fig_height),
        dpi=200,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios, 'hspace': 0.4}
    )
    
    # ===============================
    # Condensin peak tracks
    # ===============================
    for i, prefix in enumerate(condensin_data.keys()):
        ax = axes[i]
        peaks = condensin_data[prefix]["peaks"]
    
        y_bottom = 0.1
        y_height = 0.8
        for start, end in peaks:
            ax.add_patch(Rectangle((start, y_bottom), end - start, y_height,
                                   facecolor="white", edgecolor=condensin_peak_edge,
                                   linewidth=1.0, clip_on=False, zorder=10))
    
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel(prefix.split("_")[1], rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
    # ===============================
    # TF tracks
    # ===============================
    offset = n_condensin
    for j, name in enumerate(tf_data.keys()):
        ax = axes[offset + j]
        data = tf_data[name]
        x, y, peaks = data["x"], data["y"], data["peaks"]
        y_smooth = smooth_signal(x, y, window_bp=window_bp)
    
        ax.plot(x, y_smooth, color=tf_color, linewidth=1.8)
    
        ymin, ymax = ax.get_ylim()
        stripe_height = (ymax - ymin) * 0.12
        stripe_bottom = ymin - (stripe_height * 1.5)
    
        for start, end in peaks:
            ax.add_patch(Rectangle((start, stripe_bottom), end - start, stripe_height,
                                   facecolor="white", edgecolor=tf_peak_edge,
                                   linewidth=0.8, clip_on=False, zorder=10))
    
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # ===============================
    # X-axis formatting
    # ===============================
    
    tick_step = 200_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels([f"{x/1e6:.2f}" for x in xticks])
    axes[-1].tick_params(axis="x", which="both", top=False, bottom=True, labelbottom=True, pad=20)
    axes[-1].set_xlabel("Position (Mb)", labelpad=15)
    
    fig.suptitle(f"{chrom}, {region_start/1000000}-{region_end/1000000} Mb | rolling window = {window_bp/1000:.0f} kb",
                 fontsize=14, y=0.995)
    
    plt.subplots_adjust(hspace=0.5, left=0.15, right=0.95, bottom=0.18)
    plt.xlim(region_start, region_end)
    
    # ===============================
    # Legend
    # ===============================
    legend_elements = [
        Rectangle((0,0),1,1,facecolor='white', edgecolor=condensin_peak_edge, label='Condensin peaks'),
        Line2D([0],[0], color=tf_color, lw=2, label='TF signal'),
        Rectangle((0,0),1,1,facecolor='white', edgecolor=tf_peak_edge, label='TF peaks'),
    ]
    
    fig.legend(handles=legend_elements, loc='upper right', frameon=False, bbox_to_anchor=(0.98,0.98))
    
    plt.savefig("chip_tracks_condensin_tfs_peaks_cut_16_21mb.pdf")
    plt.show()

    

    histone_files = {
        "H3K9me3": args.histone_dir / "H3K9me3_2.bedgraph",
        "H3K27me3": args.histone_dir / "H3K27me3.bedgraph",
        "H3K36me3": args.histone_dir / "H3K36me3.bedgraph",
    }
    
    histone_data = {}
    
    for name, path in histone_files.items():
        x, y = load_bedgraph_region(path, chrom, region_start, region_end)
        histone_data[name] = {"x": x, "y": y}


    base_folder = args.data_dir
    gen_folder = base_folder / "rnaseq"
    res_folder = args.results_dir / "DE"


    r_de_results = pd.read_csv(
        res_folder / "DEA_results_renamed_shrunk_270126.tsv", sep="\t", index_col=0
    )


    gtf_ann = pd.read_csv(
        gen_folder / "caenorhabditis_elegans.PRJNA13758.WBPS18.canonical_geneset.gtf",
        sep="\t",
        header=None,
        names=[
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    ).drop(0)


    gtf_ann_tr = gtf_ann.loc[gtf_ann.feature == "gene"]
    gtf_ann_tr.index = gtf_ann_tr.attribute.apply(lambda x: x.split()[1][1:-2])


    start_region = 16_600_000
    end_region = 20_900_000
    
    all_tr_list = gtf_ann_tr.attribute.apply(lambda x: x.split()[1][1:-2]).unique()
    gtf_ann_sub = gtf_ann_tr.loc[gtf_ann_tr.chr == "V"]
    gtf_ann_sub = gtf_ann_sub.loc[
        (gtf_ann_sub.start >= start_region) & (gtf_ann_sub.start <= end_region)
    ]
    tr_list = gtf_ann_sub.attribute.apply(lambda x: x.split()[1][1:-2]).unique()


    r_de_results_sub_reg = r_de_results.loc[
        list(set(r_de_results.index).intersection(set(tr_list)))
    ]
    r_de_results_sub_reg = r_de_results_sub_reg.loc[
        gtf_ann_tr.loc[r_de_results_sub_reg.index].sort_values("start").index
    ]
    r_de_results_sub_reg = r_de_results_sub_reg.loc[
        gtf_ann_tr.loc[r_de_results_sub_reg.index].start < end_region
    ]
    log2fc_columns = [
        col for col in r_de_results_sub_reg.columns if "log2FoldChange" in col
    ]
    heatmap_data = r_de_results_sub_reg[log2fc_columns]
    heatmap_data = heatmap_data.iloc[:, :-1]  # exclude 141+
    heatmap_data = heatmap_data[heatmap_data.notna().all(axis=1)]
    genomic_positions = gtf_ann_tr.loc[heatmap_data.index, "start"]
    genomic_positions_Mb = (genomic_positions / 1e6).round(1)
    
    expr_columns = [col for col in r_de_results_sub_reg.columns if "_meanExpr" in col]
    heatmap_data = r_de_results_sub_reg[expr_columns].iloc[:, :-1]
    
    n = 20
    coef_multi = 100
    
    gene_distances = np.diff(
        np.concatenate(
            ([region_start], gtf_ann_tr.loc[heatmap_data.index, "start"], [region_end])
        )
    )
    
    proportional_widths = gene_distances / gene_distances.sum() * heatmap_data.shape[0]
    proportional_widths = (proportional_widths * coef_multi).astype(int)
    
    interpolated_heatmap_data = np.zeros(
        (heatmap_data.shape[1], int(sum(proportional_widths)))
    )
    
    start = 0
    for idx, width in enumerate(proportional_widths[:-1]):
        if width > 0:
            gene_data = heatmap_data.iloc[idx].values
            for row in range(gene_data.shape[0]):
                interpolated_heatmap_data[row, start : start + width] = gene_data[row]
        start += width


    hot_path = args.hot_sites_bed

    hot_df = pd.read_csv(
        hot_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "id", "score", "strand"],
    )


    hot_reg = hot_df[
        (hot_df["chrom"] == chrom)
        & (hot_df["start"] >= region_start)
        & (hot_df["start"] <= region_end)
    ]


    active_path = args.active_domains_bed
    all_domains_path = args.all_domains_bed
    
    active_df = pd.read_csv(
        active_path, sep="\t", header=None, names=["chrom", "start", "end"]
    )
    
    active_reg = active_df[
        (active_df["chrom"] == chrom)
        & (active_df["end"] >= region_start)
        & (active_df["start"] <= region_end)
    ]


    all_domains_df = pd.read_csv(
        all_domains_path,
        sep="\t",
        skiprows=2,
        names=["chrom", "start", "end", "feature", "dd", "d", "color"],
    )
    all_domains_df = all_domains_df.drop(["dd", "d", "color"], axis=1)
    reg_domains_df = all_domains_df.loc[all_domains_df.feature == "Regulated"]
    reg_domains_region = reg_domains_df[
        (reg_domains_df["chrom"] == chrom)
        & (reg_domains_df["end"] >= region_start)
        & (reg_domains_df["start"] <= region_end)
    ]


    # ==============================
    # Parameters
    # ==============================
    window_bp = 2_000
    histone_color = "#6a51a3"
    LABEL_X_POS = -0.08
    
    n_histone = len(histone_data)

    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 12
    
    fig = plt.figure(figsize=(23.38, 16.54), dpi=300)
    
    gs = GridSpec(
        n_histone + 1,
        2,
        height_ratios=[1] * n_histone + [3],
        width_ratios=[50, 1],
        hspace=0.25,
        wspace=0.05,
    )
    
    # ==============================
    # HISTONE TRACKS
    # ==============================
    axes_histone = []
    
    for i, name in enumerate(histone_data.keys()):
        ax = fig.add_subplot(gs[i, 0])
        axes_histone.append(ax)
    
        x = histone_data[name]["x"]
        y = histone_data[name]["y"]
    
        if i == 2:
            for _, row in active_reg.iterrows():
                ax.axvspan(
                    row["start"],
                    row["end"],
                    color="#ffcccc",
                    alpha=0.25,
                    zorder=0,
                )
        else:
            for _, row in reg_domains_region.iterrows():
                ax.axvspan(
                    row["start"],
                    row["end"],
                    color="#ADD8E6",
                    alpha=0.25,
                    zorder=0,
                )
    
        y_smooth = smooth_signal(x, y, window_bp=window_bp)
        ax.plot(x, y_smooth, color=histone_color, linewidth=1.8)
    
        ax.set_xlim(region_start, region_end)
        ax.set_ylabel(name, rotation=0, ha="left", va="center")
        ax.yaxis.set_label_coords(LABEL_X_POS, 0.5)
    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    
    ax_overlay = axes_histone[2]
    
    trans = transforms.blended_transform_factory(ax_overlay.transData, ax_overlay.transAxes)
    
    rect_height = 0.12
    gap = 0.05
    
    for _, row in hot_reg.iterrows():
        x0 = (row["start"] + row["end"]) / 2
        width = row["end"] - x0
    
        rect = Rectangle(
            (x0, -(rect_height + gap)),
            max(width, 5000),
            rect_height,
            transform=trans,
            facecolor="red",
            edgecolor="none",
            alpha=0.9,
            clip_on=False,
            zorder=10,
        )
        ax_overlay.add_patch(rect)
    
    hot_patch = mpatches.Patch(color="red", fill=False, label="Hot Regions")
    active_domain_patch = mpatches.Patch(color="#ffcccc", alpha=0.5, label="Active Domains")
    reg_domain_patch = mpatches.Patch(color="#ADD8E6", alpha=0.5, label="Regulated Domains")
    
    ax_overlay.legend(
        handles=[hot_patch, active_domain_patch, reg_domain_patch],
        loc="upper left",
        bbox_to_anchor=(1, 1),
        frameon=False,
        fontsize=12,
    )
    
    # ==============================
    # HEATMAP
    # ==============================
    ax_heatmap = fig.add_subplot(gs[n_histone, 0], sharex=axes_histone[0])
    
    ax_heatmap.patch.set_visible(False)
    
    hm = ax_heatmap.imshow(
        interpolated_heatmap_data + 1,
        aspect="auto",
        cmap=sns.color_palette("crest", as_cmap=True),
        norm=LogNorm(vmin=1, vmax=interpolated_heatmap_data.max() + 1),
        extent=[region_start, region_end, 0, interpolated_heatmap_data.shape[0]],
        interpolation="nearest",
        zorder=1,
    )
    
    ax_heatmap.set_xlabel("Genomic position (Mb)")
    
    tick_step = 200_000
    xticks = np.arange(region_start, region_end + 1, tick_step)
    
    ax_heatmap.set_xticks(xticks)
    ax_heatmap.set_xticklabels([f"{x/1e6:.1f}" for x in xticks])
    
    axes_histone[0].xaxis.set_label_position("top")
    axes_histone[0].tick_params(axis="x", which="both", top=True, labeltop=True)
    axes_histone[0].set_xticks(xticks)
    axes_histone[0].set_xticklabels([f"{x/1e6:.1f}" for x in xticks])
    
    ax_heatmap.set_yticks(np.arange(0.5, interpolated_heatmap_data.shape[0] + 0.5))
    ax_heatmap.set_yticklabels(
        [label.split("_")[0] for label in heatmap_data.columns], rotation=0
    )
    
    cax = fig.add_subplot(gs[n_histone, 1])
    colorbar = plt.colorbar(hm, cax=cax)
    colorbar.set_label("log10(RPKM + 1)", fontsize=12)
    
    fig.suptitle(
        f"{chrom}, {region_start/1000000}-{region_end/1000000} Mb | rolling window = {window_bp/1000:.0f} kb",
        fontsize=12,
        y=0.94,
    )
    
    plt.subplots_adjust(top=0.90, bottom=0.1, left=0.1, right=0.9)
    plt.savefig("histone_expression_plot_A4.pdf", bbox_inches="tight")
    plt.show()


    

if __name__ == "__main__":
    main()
