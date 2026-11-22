# Xenium analysis comparing carbon-ion irradiation with gamma-ray and X-ray irradiation

```text
mouse_analysis/   mouse jejunal samples, 14 sections on two slides (No1, No2)
human_analysis/   human hypopharyngeal samples, 5 sections on one slide (SF1, SF2, Xray_Under, Xray_Late, Cion_Late)
xenium_data_*/    merge-rois-*.ipynb merges the ROI zarrs into integrated/integrated_*.zarr
figstyle.py       plot settings
```

## Running

00 made `data/integrated_*_preprocessed.zarr` and is kept as a record. Harmony does not always give the same result on a different machine, and the clusters and annotation in 01 depend on it, so start from 01.

Before running, put the deposited files in these folders (relative to this repository):

| File | Folder |
| --- | --- |
| `integrated_mouse_preprocessed.zarr` | `mouse_analysis/data/` |
| `qc_per_cell_before.parquet` (mouse) | `mouse_analysis/results/00_preprocess_and_cluster/` |
| `roi_numbers.csv` | `mouse_analysis/results/00_preprocess_and_cluster/` |
| `integrated_mouse.zarr` | `xenium_data_mouse/integrated/` |
| `integrated_human_preprocessed.zarr` | `human_analysis/data/` |
| `qc_per_cell_before.parquet` (human) | `human_analysis/results/00_preprocess_and_cluster/` |

- mouse: 01 to 06 in order
- human: 01, 02, 03, `Rscript aux_downsample_fibroblast.R`, 04

Notebooks write to `data/` and `results/<notebook>/`.

## Environment

Python 3.12 (`requirements.txt`). R 4.4.2 with Seurat 5.5.0, SeuratObject 5.4.0, DropletUtils 1.26.0 and tidyverse 2.0.0.
