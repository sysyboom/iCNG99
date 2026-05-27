import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


data = {
    "model": [
        "iUL909",
        "iUma22",
        "papla-GEM",
        "iPrub22",
        "iUmbe1",
        "iYali5",
        "iCNG99",
        "iRV890",
        "iYL_w29",
        "iNR1684",
        "AltGEM_iDD1552"
    ],

    "biomass": [67, 41, 48, 48, 48, 48, 103, 66, 47, 43, 46],

    "lipids": [23, 8, 12, 4, 15, 16, 55, 23, 14, 11, 14]
}

df = pd.DataFrame(data)

df = df.sort_values(
    "biomass",
    ascending=True
).reset_index(drop=True)

out_dir = "/mnt/NFS/fengch/new_score_models"
os.makedirs(out_dir, exist_ok=True)

highlight_model = "iCNG99"

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

fig, ax = plt.subplots(figsize=(8.2, 5.8))

y = np.arange(len(df))
offset = 0.16

biomass_color = "#43aa8b"
lipid_color = "#f9c74f"

for i in range(len(df)):

    model = df.loc[i, "model"]

    is_highlight = model == highlight_model

    alpha = 1.0 if is_highlight else 0.7

    biomass_value = df.loc[i, "biomass"]
    lipid_value = df.loc[i, "lipids"]

    ax.hlines(
        y=y[i] + offset,
        xmin=0.8,
        xmax=max(0.8, biomass_value - 1.0),

        linewidth=2.0 if is_highlight else 2.0,

        color=biomass_color,

        linestyle="-",

        alpha=alpha
    )

    ax.scatter(
        biomass_value,
        y[i] + offset,

        s=60 if is_highlight else 20,

        facecolors="white",

        edgecolors=biomass_color,

        linewidth=1.2 if is_highlight else 0.6,

        alpha=alpha,

        zorder=3
    )

    ax.hlines(
        y=y[i] - offset,
        xmin=0.8,
        xmax=max(0.8, lipid_value - 0.8),

        linewidth=2.2 if is_highlight else 1.3,

        color=lipid_color,

        linestyle="--",

        alpha=alpha
    )

    ax.scatter(
        lipid_value,
        y[i] - offset,

        s=65 if is_highlight else 22,

        color=lipid_color,

        edgecolors="black" if is_highlight else lipid_color,

        linewidth=0.5 if is_highlight else 0.2,

        alpha=alpha,

        zorder=3
    )

for i in range(len(df)):

    model = df.loc[i, "model"]

    is_highlight = model == highlight_model

    biomass_value = df.loc[i, "biomass"]
    lipid_value = df.loc[i, "lipids"]

    # biomass
    ax.text(
        biomass_value + 1.2,
        y[i] + offset,

        str(biomass_value),

        va="center",

        fontsize=8.5 if is_highlight else 8.5,

        color=biomass_color,

        fontweight="bold" if is_highlight else "normal",

        alpha=1.0 if is_highlight else 0.7
    )

    # lipid
    ax.text(
        lipid_value + 1.2,
        y[i] - offset,

        str(lipid_value),

        va="center",

        fontsize=8.5 if is_highlight else 7.2,

        color=lipid_color,

        fontweight="bold" if is_highlight else "normal",

        alpha=1.0 if is_highlight else 0.7
    )

ax.set_yticks(y)
ax.set_yticklabels(
    df["model"],
    fontsize=9
)

ax.set_xlabel(
    "Number of Components",
    fontsize=11,
    labelpad=10
)

ax.set_title(
    "Comparison of biomass formulations across GEMs",
    fontsize=13,
    fontweight="bold",
    pad=14
)

legend_elements = [

    Line2D(
        [0], [0],

        color=biomass_color,

        marker="o",

        markerfacecolor="white",

        markeredgecolor=biomass_color,

        linestyle="-",

        linewidth=1.5,

        markersize=5,

        label="Biomass components"
    ),

    Line2D(
        [0], [0],

        color=lipid_color,

        marker="o",

        markerfacecolor=lipid_color,

        markeredgecolor=lipid_color,

        linestyle="-",

        linewidth=1.5,

        markersize=5,

        label="Lipid components"
    )
]

ax.legend(
    handles=legend_elements,

    frameon=False,

    fontsize=8.5,

    loc="lower right"
)

ax.set_xlim(0, 115)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.grid(
    axis="x",

    linestyle="--",

    linewidth=0.5,

    alpha=0.22
)

ax.tick_params(labelsize=9)

plt.tight_layout()

pdf_path = os.path.join(
    out_dir,
    "Biomass_Lipid_Lollipop_Final.pdf"
)


plt.savefig(
    pdf_path,
    bbox_inches="tight"
)

plt.savefig(
    png_path,
    dpi=600,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print(pdf_path)
