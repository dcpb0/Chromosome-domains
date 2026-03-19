#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from pathlib import Path
import collections
import argparse

import cooler
from pycircos import pycircos

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

try:
    from IPython.display import display
except Exception:
    def display(obj):
        print(obj)


# ============================================================================
# Module-level constants
# ============================================================================

RNAPOL_DICT = {
    "protein_coding": "RNAPol II",
    "piRNA": "RNAPol II",
    "ncRNA": "Varies",
    "pseudogene": "RNAPol II",
    "tRNA": "RNAPol III",
    "snoRNA": "RNAPol II/III",
    "miRNA": "RNAPol II",
    "lincRNA": "RNAPol II",
    "snRNA": "RNAPol II/III",
    "antisense_RNA": "RNAPol II",
    "rRNA": "RNAPol I/III",
}

POLYA_DICT = {
    "protein_coding": "Poly(A) Tailed",
    "piRNA": "Non-Poly(A) Tailed",
    "ncRNA": "Varies",
    "pseudogene": "Varies",
    "tRNA": "Non-Poly(A) Tailed",
    "snoRNA": "Non-Poly(A) Tailed",
    "miRNA": "Non-Poly(A) Tailed",
    "lincRNA": "Poly(A) Tailed",
    "snRNA": "Varies",
    "antisense_RNA": "Varies",
    "rRNA": "Non-Poly(A) Tailed",
}

RNAPOL_COLOR_DICT = {
    "RNAPol I": "#8B4513",
    "RNAPol II": "#1E90FF",
    "RNAPol III": "#32CD32",
    "Varies": "#DAA520",
    "RNAPol II/III": "#6A5ACD",
    "RNAPol I/III": "#2F4F4F",
}

GT_COLOR_DICT = {
    "protein_coding": "#FF4500",
    "piRNA": "#00008B",
    "ncRNA": "#008000",
    "pseudogene": "#FFD700",
    "tRNA": "#800080",
    "snoRNA": "#00CED1",
    "miRNA": "#C71585",
    "lincRNA": "#A52A2A",
    "snRNA": "#FFA07A",
    "antisense_RNA": "#D2691E",
    "rRNA": "#006400",
}

POLYA_COLOR_DICT = {
    "Poly(A) Tailed": "#1E90FF",
    "Non-Poly(A) Tailed": "#228B22",
    "Varies": "#DAA520",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate FitHiC-derived tables and circos plots for chr V."
    )
    parser.add_argument(
        "--hic-dir",
        type=Path,
        required=True,
        help="Directory with HIC helper tables (chr sizes, insulation, gene windows).",
    )
    parser.add_argument(
        "--gse168803-dir",
        type=Path,
        required=True,
        help="Directory with GSE168803 mcool/FitHiC files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gen_folder = args.hic_dir
    gse168803_folder = args.gse168803_dir
    file_path = gse168803_folder / "N2_emb_DJ73_DJ74_30.mcool"
    
    resolution = 10000
    uri = f"{file_path}::resolutions/{resolution}"
    
    clr = cooler.Cooler(uri)
    mat = clr.matrix(balance=False, sparse=True)[:]
    
    # row sums
    marginals = np.array(mat.sum(axis=1)).flatten()
    bins = clr.bins()[:]
    bins["mid"] = (bins["start"] + bins["end"]) // 2

    pixels = clr.pixels()[:]
    
    interactions = pd.DataFrame(
        {
            "chr1": bins.loc[pixels["bin1_id"], "chrom"].values,
            "mid1": bins.loc[pixels["bin1_id"], "mid"].values,
            "chr2": bins.loc[pixels["bin2_id"], "chrom"].values,
            "mid2": bins.loc[pixels["bin2_id"], "mid"].values,
            "count": pixels["count"].values,
        }
    )
    interactions = interactions[interactions["count"] > 0]

    low = np.quantile(marginals, 0.01)
    high = np.quantile(marginals, 0.999)
    mask = (marginals > low) & (marginals < high)

    weights = bins["weight"].copy()
    weights = weights.replace([np.inf, -np.inf], np.nan)
    
    bias_scaled = pd.DataFrame({"chr": bins["chrom"], "mid": bins["mid"], "bias": weights})
    bias_scaled["bias"] /= bias_scaled["bias"].mean()

    mask = mask * bias_scaled["bias"].notna().values

    bias_scaled["bias"].fillna(1, inplace=True)

    frags = pd.DataFrame(
        {
            "chr": bins["chrom"],
            "extra": 0,
            "mid": bins["mid"],
            "marginal": marginals.astype(int),
            "mappable": mask.astype(int),
        }
    )

    df = pd.read_csv(
        gse168803_folder / "fithic_all/FitHiC.spline_pass1.res10000.significances.txt.gz",
        sep="\t",
    )
    df = df[df["chr1"] != df["chr2"]]
    print("Number of tested intrachromosomal interactions:", len(df))
    print("Number of interactions with count > 0:", (df["contactCount"] > 0).sum())
    print("Q-value stats:")
    print(df["q-value"].describe())

    mask = df["chr2"] == "V"
    df.loc[
        mask, ["chr1", "fragmentMid1", "chr2", "fragmentMid2", "bias1", "bias2"]
    ] = df.loc[
        mask, ["chr2", "fragmentMid2", "chr1", "fragmentMid1", "bias2", "bias1"]
    ].values

    roi_start = 16_600_000
    roi_chr = 'V'
    
    roi = df[(df["chr1"] == roi_chr) & (df["fragmentMid1"] >= roi_start)]

    for thr in [0.05, 0.01, 1e-3, 1e-4, 1e-5]:
        roi_sig = roi[roi["q-value"] < thr]
        roi_sig = roi_sig[roi_sig["ExpCC"] > 0.1]
        print(len(roi_sig))
        display(roi_sig.chr2.value_counts())

    sign_fdr_thr = 1e-4

    roi_sig = roi[roi["q-value"] < sign_fdr_thr]

    ax = sns.scatterplot(
        x=roi["ExpCC"] + 1, y=roi["contactCount"] + 1, hue=roi["q-value"], alpha=0.3
    )
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.xlabel("Expected")
    plt.ylabel("Observed")
    plt.title("Interchromosomal contacts: observed vs expected")
    plt.show()

    plt.scatter(x=roi_sig["ExpCC"], y=roi_sig["contactCount"], alpha=0.3)
    plt.xlabel("Expected")
    plt.ylabel("Observed")
    plt.title("Interchromosomal contacts: observed vs expected")
    plt.show()

    sign_fdr_thr = 1e-3
    roi_start = 16_000_000
    roi_end = 19_000_000
    
    roi = df[
        (df["chr1"] == roi_chr)
        & (df["fragmentMid1"] >= roi_start)
        & (df["fragmentMid1"] <= roi_end)
    ]
    roi_sig = roi[roi["q-value"] < sign_fdr_thr]

    resolution = 10000
    half_bin = resolution // 2
    
    df_bedpe = pd.DataFrame(
        {
            "chr1": roi_sig["chr1"],
            "start1": roi_sig["fragmentMid1"] - half_bin,
            "end1": roi_sig["fragmentMid1"] + half_bin,
            "chr2": roi_sig["chr2"],
            "start2": roi_sig["fragmentMid2"] - half_bin,
            "end2": roi_sig["fragmentMid2"] + half_bin,
        }
    )
    # df_bedpe.to_csv(gse168803_folder/'chr_V_16.csv', index=False)

    whole_region = (16680000, 20923630)
    corg = (16_680_000, 17_980_000)
    region_5s = (17180000, 17280000)
    
    regions_dict = {'whole_region': whole_region,
                   'corg': corg,
                   "region_5s":region_5s}

    sign_fdr_thr = 1e-3
    resolution = 10000
    half_bin = resolution // 2
    
    for name, borders in regions_dict.items():
        roi = df[
            (df["chr1"] == "V")
            & (df["fragmentMid1"] >= borders[0])
            & (df["fragmentMid1"] <= borders[1])
        ]
        roi_sig = roi[roi["q-value"] < sign_fdr_thr]
        df_bedpe = pd.DataFrame(
            {
                "chr1": roi_sig["chr1"],
                "start1": roi_sig["fragmentMid1"] - half_bin,
                "end1": roi_sig["fragmentMid1"] + half_bin,
                "chr2": roi_sig["chr2"],
                "start2": roi_sig["fragmentMid2"] - half_bin,
                "end2": roi_sig["fragmentMid2"] + half_bin,
            }
        )
        print(df_bedpe.shape)
        df_bedpe.to_csv(gse168803_folder/f'{name}_chr_V_contacts.csv', index=False)

    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle

    for region_name in regions_dict.keys():
        Garc = pycircos.Garc
        Gcircle = pycircos.Gcircle
        pal_chr = sns.color_palette("hls", 6)
        circle = Gcircle(figsize=(8, 8))
        with open(gen_folder / "chrsizes_celegans_general.csv") as f:
            f.readline()
            for i, line in enumerate(f):
                line = line.rstrip().split(",")
                name = line[0]
                length = int(line[-1])
                arc = Garc(
                    arc_id=name,
                    size=length,
                    interspace=2,
                    raxis_range=(935, 985),
                    labelposition=100,
                    label_visible=True,
                    facecolor=pal_chr[i],
                )
                circle.add_garc(arc)
        circle.set_garcs(-65, 245)
        for arc_id in circle.garc_dict:
            if arc_id == "V":
                ticklabelschr = list(np.arange(0, circle.garc_dict[arc_id].size // 1000000 + 1))
            else:
                ticklabelschr = None
            print(ticklabelschr)
            circle.tickplot(
                arc_id, raxis_range=(985, 1000), tickinterval=1000000, ticklabels=ticklabelschr
            )
        
        values_all = []
        arcdata_dict = collections.defaultdict(dict)
        with open(gse168803_folder / f'{region_name}_chr_V_contacts.csv') as f:
            f.readline()
            for line in f:
                line = line.rstrip().split(",")
                name1 = line[0]
                start1 = int(line[1]) - 1
                end1 = int(line[2])
                name2 = line[3]
                start2 = int(line[4]) - 1
                end2 = int(line[5])
                source = (name1, start1, end1, 835)
                destination = (name2, start2, end2, 835)
                circle.chord_plot(
                    source, destination, facecolor=circle.garc_dict[name2].facecolor
                )
        
        arcdata_dict = collections.defaultdict(dict)
        with open(gen_folder / "aggregated_gene_types_per_window.csv") as f:
            f.readline()
            for line in f:
                line = line.rstrip().split(",")
                name = line[0]
                start = int(line[1])
                width = int(line[2]) - (int(line[1])) + 1
                if name not in arcdata_dict:
                    arcdata_dict[name]["positions"] = []
                    arcdata_dict[name]["widths"] = []
                    arcdata_dict[name]["colors"] = []
                arcdata_dict[name]["positions"].append(start)
                arcdata_dict[name]["widths"].append(width)
                arcdata_dict[name]["colors"].append(gt_color_dict[line[-1]])
        
        for key in arcdata_dict:
            circle.barplot(
                key,
                data=[1] * len(arcdata_dict[key]["positions"]),
                positions=arcdata_dict[key]["positions"],
                width=arcdata_dict[key]["widths"],
                raxis_range=[840, 855],
                facecolor=arcdata_dict[key]["colors"],
            )
        # line plot
        values_all = []
        arcdata_dict = collections.defaultdict(dict)
        with open(gen_folder / "aggregated_insul_score.csv") as f:
            f.readline()
            for line in f:
                line = line.rstrip().split(",")
                name = line[0]
                start = int(line[1]) - 1
                end = int(line[2])
                mid = (start + end) / 2
                value = float(line[-1])
                values_all.append(value)
                if name not in arcdata_dict:
                    arcdata_dict[name]["positions"] = []
                    arcdata_dict[name]["values"] = []
                arcdata_dict[name]["positions"].append(mid)
                arcdata_dict[name]["values"].append(value)
        
        vmin, vmax = min(values_all), max(values_all)
        for key in arcdata_dict:
            circle.lineplot(
                key,
                data=arcdata_dict[key]["values"],
                positions=arcdata_dict[key]["positions"],
                rlim=[vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)],
                raxis_range=(860, 925),
                linecolor="royalblue",
                spine=False,
            )
        
        # Create a list of patches for the legend
        legend_patches = [
            mpatches.Patch(color=color, label=gene_type)
            for gene_type, color in gt_color_dict.items()
        ]
        
        # Add the legend to the plot
        plt.legend(
            handles=legend_patches,
            loc="upper right",
            bbox_to_anchor=(1.25, 1),
            title="Gene Types",
        )
        
        # Show the plot
        plt.savefig(f"{region_name}_gtypes_circplot.pdf", bbox_inches="tight")
        plt.show()

if __name__ == "__main__":
    main()
