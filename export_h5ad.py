import os

import anndata as ad

os.makedirs("geo_h5ad", exist_ok=True)
for s in ("human", "mouse"):
    sources = {
        f"integrated_{s}": f"xenium_data_{s}/integrated/integrated_{s}.zarr/tables/table",
        f"integrated_{s}_annotated": f"{s}_analysis/data/integrated_{s}_annotated.zarr",
    }
    for name, src in sources.items():
        ad.read_zarr(src).write_h5ad(f"geo_h5ad/{name}.h5ad", compression="gzip", convert_strings_to_categoricals=False)
