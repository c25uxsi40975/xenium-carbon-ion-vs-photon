# Xenium analysis comparing carbon-ion irradiation with gamma-ray and X-ray irradiation

```text
mouse_analysis/       mouse jejunal samples, 14 sections on two slides (No1, No2)
human_analysis/       human hypopharyngeal samples, 5 sections on one slide (SF1, SF2, Xray_Under, Xray_Late, Cion_Late)
xenium_data_*/        merge-rois-*.ipynb merges the ROI zarrs into integrated/integrated_*.zarr
dystroglycan_image_analysis/  β-dystroglycan-positive area of vessel walls measured with MIPAR (one row per ROI; not the Xenium sections)
build_from_geo.ipynb  rebuilds the inputs of the notebooks from the GEO files
export_h5ad.py        wrote the h5ad files deposited in GEO
figstyle.py           plot settings
```

## Running

The data are in GEO (<GEO accession>). Put all files of the series and its samples in `geo/` in this repository and run `build_from_geo.ipynb` from the repository root. It writes:

| File | From |
| --- | --- |
| `*_analysis/data/integrated_*_annotated.zarr` | `integrated_*_annotated.h5ad` |
| `*_analysis/data/integrated_*_preprocessed.zarr` | `integrated_*_annotated.h5ad`, without the cell types and other results of 01 |
| `*_analysis/results/00_preprocess_and_cluster/qc_per_cell_before.parquet` | `human_qc_per_cell_before.parquet`, `mouse_qc_per_cell_before.parquet` |
| `mouse_analysis/results/00_preprocess_and_cluster/roi_numbers.csv` | `mouse_roi_numbers.csv` |
| `xenium_data_*/integrated/integrated_*.zarr` | `integrated_*.h5ad` and the files of each sample (images, cell and nucleus outlines, transcripts) |

00 and 01 are kept as a record of how the clusters and cell types were made. Harmony does not always give the same result on a different machine, so rerunning them may not give the same annotation. `integrated_*_annotated.h5ad` in GEO holds the result of 01, so the analysis can be reproduced from 02:

- mouse: 02 to 06 in order (06 also reads the results of 03 and 04)
- human: 02, 03, `Rscript aux_downsample_fibroblast.R`, 04

`dystroglycan_image_analysis/dystroglycan_violin.ipynb` does not use the GEO files. It reads `dystroglycan_vascular_area_fraction.csv` in the same folder, which has one row per ROI (one vessel):

- `group`: SF, C-ion or X-ray
- `roi`: number of the ROI within the group (numbers do not correspond between groups)
- `area_fraction_pct`: β-dystroglycan-positive area as a percentage of the vessel wall area, measured with MIPAR

Notebooks write to `data/` and `results/<notebook>/`.

## Environment

Python 3.12 (`requirements.txt`). R 4.4.2 with Seurat 5.5.0, SeuratObject 5.4.0, DropletUtils 1.26.0 and tidyverse 2.0.0.
