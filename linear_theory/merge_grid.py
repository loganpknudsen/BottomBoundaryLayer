"""Combine the chunk files from stability_grid.py into one dataset.

Usage: python merge_grid.py GRID OUTDIR  ->  OUTDIR/GRID.nc
"""
import glob
import os
import sys

import xarray as xr

name, outdir = sys.argv[1:3]
files = sorted(glob.glob(os.path.join(outdir, f"{name}_[0-9][0-9].nc")))
ds = xr.concat([xr.open_dataset(f) for f in files], dim="slope_burger_number")
ds = ds.sortby("slope_burger_number")
ds.to_netcdf(os.path.join(outdir, f"{name}.nc"))
print(f"merged {len(files)} files -> {os.path.join(outdir, name + '.nc')}: {dict(ds.sizes)}")
