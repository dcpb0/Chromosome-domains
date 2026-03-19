#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate scRNA-seq figures from DE and annotation tables."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Root directory with used input datasets.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Root directory with results tables.",
    )
    return parser.parse_args()

def main():
    args = parse_args()
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

    all_tr_list = gtf_ann_tr.attribute.apply(lambda x: x.split()[1][1:-2]).unique()

    gtf_ann_sub = gtf_ann_tr.loc[gtf_ann_tr.chr == "V"]

    start_region = 16_600_000
    end_region = 20_900_000
    
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

    display_interval = 50
    mask = np.zeros(heatmap_data.shape[0], dtype=bool)
    mask[::display_interval] = True
    
    selected_positions_labels = genomic_positions_Mb[mask].astype(str)
    
    plt.figure(figsize=(25, 5))
    ax = sns.heatmap(heatmap_data.T, cmap="vlag", center=0, annot=False)
    plt.title("LFC heatmap of genes from chrV:16.6-20.9 Mb", y=1.05)
    
    ax.set_xticks(np.arange(heatmap_data.shape[0])[mask])
    ax.set_xticklabels(selected_positions_labels, rotation=45, ha="right")
    
    original_labels = heatmap_data.columns
    clean_labels = np.array([label.split("_")[0] for label in original_labels])
    ax.set_yticklabels(clean_labels, rotation=0)
    
    plt.xlabel("Genomic Position (Mb)")
    plt.savefig(res_folder / "whole_region_lfc_without_insulation.pdf", bbox_inches="tight")
    plt.show()

    plt.figure(figsize=(15, 8))
    
    gs = plt.GridSpec(
        2, 2, height_ratios=[1, 3], width_ratios=[50, 1], wspace=0.05, hspace=0.05
    )
    
    n = 20
    coef_multi = 100
    chrom = "V"
    
    insulation_data = pd.read_csv(base_folder / "hic/aggregated_insul_score.csv")
    
    filtered_insulation_data = insulation_data[
        (insulation_data["chrom"] == chrom)
        & (insulation_data["chromStart"] >= start_region)
        & (insulation_data["chromEnd"] <= end_region)
    ]
    
    # Calculate proportional widths based on genomic distances
    gene_start_positions_mb = gtf_ann_tr.loc[heatmap_data.index, "start"] / 1e6
    gene_distances = np.diff(
        np.concatenate(
            ([start_region], gtf_ann_tr.loc[heatmap_data.index, "start"], [end_region])
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
    
    # Insulation scores plot (top-left)
    ax_insulation = plt.subplot(gs[0, 0])
    normalized_positions = (
        (filtered_insulation_data["chromStart"] - start_region)
        / (end_region - start_region)
        * interpolated_heatmap_data.shape[1]
    )
    ax_insulation.plot(
        normalized_positions,
        filtered_insulation_data["insulationScore"],
        label="Insulation Score",
        color="black",
    )
    ax_insulation.set_ylabel("Insulation Score")
    ax_insulation.set_xlim([0, interpolated_heatmap_data.shape[1]])
    ax_insulation.tick_params(
        axis="x", which="both", bottom=False, top=False, labelbottom=False
    )
    ax_insulation.set_title(
        f"Insulation score & Log2FC heatmap of genes from chr{chrom}:{start_region/1e6:.1f}-{end_region/1e6:.1f} Mb"
    )
    
    # Heatmap (bottom-left)
    ax_heatmap = plt.subplot(gs[1, 0], sharex=ax_insulation)
    hm = sns.heatmap(
        interpolated_heatmap_data,
        cmap="vlag",
        center=0,
        annot=False,
        cbar=False,
        ax=ax_heatmap,
    )
    ax_heatmap.set_xlabel("Genomic Position (Mb)")
    
    # Homogeneous ticks
    start_mb = start_region / 1_000_000
    end_mb = end_region / 1_000_000
    tick_step = 0.1
    ticks = np.arange(start_mb, end_mb + tick_step, tick_step)
    tick_positions = (
        (ticks - start_mb) / (end_mb - start_mb) * interpolated_heatmap_data.shape[1]
    )
    tick_labels = [f"{tick:.1f}" for tick in ticks]
    ax_heatmap.set_xticks(tick_positions)
    ax_heatmap.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax_heatmap.set_yticklabels(
        np.array([label.split("_")[0] for label in heatmap_data.columns]), rotation=0
    )
    
    # Colorbar (bottom-right, dedicated axis)
    cax = plt.subplot(gs[1, 1])
    colorbar = plt.colorbar(hm.collections[0], cax=cax)
    colorbar.set_label("L2FC")
    
    plt.tight_layout()
    plt.savefig(
        res_folder / "whole_region_lfc_with_insulation_plot.pdf", bbox_inches="tight"
    )
    plt.show()

    plt.figure(figsize=(15, 8))
    
    gs = GridSpec(
        2, 2, height_ratios=[1, 3], width_ratios=[50, 1], wspace=0.05, hspace=0.05
    )
    
    n = 20
    coef_multi = 100
    
    # Insulation data
    insulation_data = pd.read_csv(base_folder / "hic/aggregated_insul_score.csv")
    
    # Heatmap data
    expr_columns = [col for col in r_de_results_sub_reg.columns if "_meanExpr" in col]
    heatmap_data = r_de_results_sub_reg[expr_columns].iloc[:, :-1]
    
    filtered_insulation_data = insulation_data[
        (insulation_data["chrom"] == chrom)
        & (insulation_data["chromStart"] >= start_region)
        & (insulation_data["chromEnd"] <= end_region)
    ]
    
    # Calculate proportional widths based on genomic distances
    gene_distances = np.diff(
        np.concatenate(
            ([start_region], gtf_ann_tr.loc[heatmap_data.index, "start"], [end_region])
        )
    )
    proportional_widths = gene_distances / gene_distances.sum() * heatmap_data.shape[0]
    proportional_widths = (proportional_widths * coef_multi).astype(int)
    
    # Build interpolated heatmap
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
    
    # Top: insulation scores
    ax_insulation = plt.subplot(gs[0, 0])
    normalized_positions = (
        (filtered_insulation_data["chromStart"] - start_region)
        / (end_region - start_region)
        * interpolated_heatmap_data.shape[1]
    )
    ax_insulation.plot(
        normalized_positions, filtered_insulation_data["insulationScore"], color="black"
    )
    ax_insulation.set_ylabel("Insulation Score")
    ax_insulation.set_xlim([0, interpolated_heatmap_data.shape[1]])
    ax_insulation.tick_params(
        axis="x", which="both", bottom=False, top=False, labelbottom=False
    )
    ax_insulation.set_title(
        f"Insulation score & expression heatmap of genes from chr{chrom}:{start_region/1e6:.1f}-{end_region/1e6:.1f} Mb"
    )
    
    # Bottom: heatmap
    ax_heatmap = plt.subplot(gs[1, 0], sharex=ax_insulation)
    hm = sns.heatmap(
        interpolated_heatmap_data + 1,
        cmap=sns.color_palette("crest", as_cmap=True),
        annot=False,
        cbar=False,
        ax=ax_heatmap,
        norm=LogNorm(vmin=1, vmax=interpolated_heatmap_data.max() + 1),
    )
    ax_heatmap.set_xlabel("Genomic Position (Mb)")
    
    # Homogeneous ticks
    start_mb = start_region / 1_000_000
    end_mb = end_region / 1_000_000
    tick_step = 0.1
    ticks = np.arange(start_mb, end_mb + tick_step, tick_step)
    tick_positions = (
        (ticks - start_mb) / (end_mb - start_mb) * interpolated_heatmap_data.shape[1]
    )
    tick_labels = [f"{tick:.1f}" for tick in ticks]
    ax_heatmap.set_xticks(tick_positions)
    ax_heatmap.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax_heatmap.set_yticklabels(
        np.array([label.split("_")[0] for label in heatmap_data.columns]), rotation=0
    )
    
    # Colorbar (bottom-right)
    cax = plt.subplot(gs[1, 1])
    colorbar = plt.colorbar(hm.collections[0], cax=cax)
    colorbar.set_label("log10(RPKM + 1)")
    
    plt.tight_layout()
    plt.savefig(
        res_folder / "whole_region_expression_with_insulation_plot.pdf", bbox_inches="tight"
    )
    plt.show()

if __name__ == "__main__":
    main()
