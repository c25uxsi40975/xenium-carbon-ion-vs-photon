from pathlib import Path

import matplotlib as mpl
import numpy as np
import scanpy as sc

MM = 1 / 25.4
COL1 = 88 * MM
COL2 = 180 * MM
DPI = 300

RC = {
    "pdf.fonttype": 42,
    "font.sans-serif": ["Arial"],
    "mathtext.fontset": "custom", "mathtext.sf": "Arial:italic",
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "figure.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "legend.title_fontsize": 6,
    "axes.linewidth": 0.5, "lines.linewidth": 0.75, "lines.markersize": 3, "patch.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "axes.titlepad": 3, "axes.labelpad": 2,
    "legend.frameon": False, "legend.handletextpad": 0.3, "legend.borderaxespad": 0.2,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
}

BLACK, ORANGE, SKY, GREEN, YELLOW, BLUE, VERMILLION, PURPLE = (
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
GREY, LIGHTGREY = "#888888", "#CCCCCC"
UP, DOWN = VERMILLION, BLUE
CION, GAMMA, XRAY = BLUE, VERMILLION, VERMILLION
BOTH = PURPLE

def _tint(color, f):
    rgb = np.array(mpl.colors.to_rgb(color))
    return mpl.colors.to_hex(rgb + (1 - rgb) * f)

ROI_COLORS = {
    "SF1": LIGHTGREY, "SF2": "#AAAAAA", "Xray_Under": _tint(XRAY, 0.5), "Xray_Late": XRAY, "Cion_Late": CION,
    "C-ion_Sham": LIGHTGREY, "C-ion_24_h": _tint(CION, 0.6), "C-ion_48_h": _tint(CION, 0.3), "C-ion_3.5_d": CION,
    "gamma-ray_Sham": "#AAAAAA", "gamma-ray_24_h": _tint(GAMMA, 0.6), "gamma-ray_48_h": _tint(GAMMA, 0.3),
    "gamma-ray_3.5_d": GAMMA,
}

LAYER_COLORS = {"Villus": ORANGE, "Crypt": BLUE, "Muscularis": GREEN, "Lymphoid aggregate": BLACK}

CELL_COLORS = {
    "Epithelial cell": "#C39800", "Enterocyte": "#C59F72", "Progenitor cell": "#BCB7FF", "Paneth cell": "#D70000",
    "Goblet cell": "#EEB9B9", "Enteroendocrine cell": "#028800",
    "Fibroblast": "#00ACC7", "Telocyte": "#FF7FD1", "Trophocyte": "#FFA530", "Adamdec1+ fibroblast": "#00ACC7",
    "Villus fibroblast": "#8C3CFF",
    "Endothelial cell": "#C55746", "Vascular EC": "#A1756A", "Lymphatic EC": "#95B578",
    "Muscularis SMC": "#C004B9", "Villus SMC": "#0774D8",
    "Immune cell": "#6A68B0", "Macrophage": "#04C600", "Eosinophil": "#5E7E66", "Neutrophil": "#CBC400",
    "ILC": "#EC0077", "T cell": "#8F7A00", "Plasma cell": "#FF7266",
    "Enteric neuron": "#A67BB9", "Mesothelial cell": "#9E4B00", "Unassigned": "#C7C7C7",
}

HUMAN_CELL_COLORS = {
    "Myofiber": "#E5C494", "Fibroblast": "#1F77B4", "MuSC": "#FF7F0E", "Perineurial cell": "#9EDAE5",
    "Schwann cell": "#8C564B", "Adipocyte": "#FFBB78",
    "Arterial EC": "#E377C2", "EC": "#AA40FC", "Venous EC": "#17BECF", "Lymphatic EC": "#C49C94",
    "Arterial SMC": "#279E68", "Venous SMC": "#B5BD61", "Pericyte": "#F7B6D2",
    "Macrophage": "#C5B0D5", "Neutrophil": "#8C6D31", "T cell": "#FFD92F", "Plasma cell": "#98DF8A",
    "Mast cell": "#E7969C", "Unassigned": "#C7C7C7",
    "Endothelial cell": "#E377C2", "Smooth muscle cell": "#279E68", "Immune cell": "#C5B0D5",
}

def palette(categories):
    return [CELL_COLORS.get(c, LIGHTGREY) for c in categories]

CONDITION = {"Xray_Late": "Xray Late", "Cion_Late": "C-ion Late", "Xray_Under": "Xray Under", "gamma": "γ-ray", "Cion": "C-ion"}

def label(x):
    return CONDITION.get(x, x)

def gene(name):
    return "$\\mathsf{" + name + "}$"

def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{path}.{ext}", dpi=DPI)

sc.set_figure_params(vector_friendly=True, frameon=False)
mpl.rcParams.update(RC)
