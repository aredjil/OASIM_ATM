#!/usr/bin/env python3
import time
import numpy as np
from netCDF4 import Dataset
from oasim import oasim_lib, calc_unit  # assuming your wrapper is in oasim.py


def benchmark():
    # --- Parameters ---
    n_iter = 300
    n_points = 412300
    n_waves = 33

    # --- Physical constants ---
    a1 = 611.21      # Pa
    a3 = 17.502
    a4 = 32.19       # K
    To = 273.16      # K
    b1 = 0.14e-2     # cm/Pa
    b2 = 0.21        # cm

    # --- Load NetCDF data ---
    home = os.environ["HOME"]
    FILE_NAME = f"{home}/projects/OASIM_ATM/data000.nc"

    print(f"Reading NetCDF data from: {FILE_NAME}")
    nc = Dataset(FILE_NAME, "r")

    # Read all required variables
    iyr = nc.variables["iyr"][:]
    iday = nc.variables["iday"][:]
    sec_b = nc.variables["sec_b"][:]
    sec_e = nc.variables["sec_e"][:]
    sp = nc.variables["sp"][:, :]
    msl = nc.variables["msl"][:, :]
    ws10 = nc.variables["ws10"][:, :]
    tco3 = nc.variables["tco3"][:, :]
    t2m = nc.variables["t2m"][:, :]
    d2m = nc.variables["d2m"][:, :]
    tcc = nc.variables["tcc"][:, :]
    tclw = nc.variables["tclw"][:, :]
    cdrem = nc.variables["cdrem"][:, :]
    taua = nc.variables["taua"][:, :, :]
    asymp = nc.variables["asymp"][:, :, :]
    ssalb = nc.variables["ssalb"][:, :, :]
    lat = nc.variables["lat"][:]
    lon = nc.variables["lon"][:]
    nc.close()

    print("✅ NetCDF data loaded successfully.")

    # --- Initialize OASIM ---
    olib = oasim_lib("../../OASIMlib/liboasim-py.so", "config.yaml", lat, lon)
    cunit = calc_unit(n_points, olib)

    # --- Preallocate arrays ---
    points = np.arange(1, n_points + 1, dtype=np.int32)
    edout = np.zeros((n_points, n_waves), dtype=np.float64, order="F")
    esout = np.zeros((n_points, n_waves), dtype=np.float64, order="F")

    # --- Benchmark Loop ---
    print("🚀 Starting benchmark loop...")
    start = time.time()

    for iter in range(n_iter):
        if iter % 10 == 0:
            print(f"  Iteration {iter+1}/{n_iter}")

        # Compute derived physical quantities
        T = t2m[:, iter]
        Td = d2m[:, iter]
        es_Td = a1 * np.exp(a3 * (Td - To) / (Td - a4))
        es_T = a1 * np.exp(a3 * (T - To) / (T - a4))
        rh = 100.0 * es_Td / es_T
        wv = b1 * es_Td * (sp[:, iter] / 100.0) / (msl[:, iter] / 100.0) + b2

        # Extract optical properties
        taua_slice = taua[:, :, iter]
        asymp_slice = asymp[:, :, iter]
        ssalb_slice = ssalb[:, :, iter]

        # Perform OASIM monrad computation
        ed, es = cunit.monrad(
            points,
            int(iyr[iter]),
            int(iday[iter]),
            float(sec_b[iter]),
            float(sec_e[iter]),
            sp[:, iter],
            msl[:, iter],
            ws10[:, iter],
            tco3[:, iter],
            t2m[:, iter],
            d2m[:, iter],
            tcc[:, iter],
            tclw[:, iter],
            cdrem[:, iter],
            taua_slice,
            asymp_slice,
            ssalb_slice,
        )

        edout += ed
        esout += es

    elapsed = time.time() - start
    print(f"\n✅ Benchmark completed in {elapsed:.2f} seconds.")

    # --- Cleanup ---
    del cunit
    del olib
    print("All done! Results are ready.")

    return edout, esout


if __name__ == "__main__":
    import os
    edout, esout = benchmark()
