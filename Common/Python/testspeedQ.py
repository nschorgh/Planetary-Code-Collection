#! /usr/bin/env python3
"""
Measuring the speed of conductionQ.py using tqdm library.
author: Cyril Mergny
On Intel CPU i7 10750
Before optimization ~ 6500 it/s 
After optimization ~ 16500 it/s
Close to 3 times faster.
"""

import numpy as np
from tqdm import tqdm

import grids
from conductionQ import conductionQ
from flux_noatm import flux_noatm


def check_nochanges(T):
    """
    Test that after the changes, the temperature outputed is still the same.
    """
    Tref = np.load("Tspeedtest.npy")
    error = abs(T - Tref)
    er_thresh = 1e-10
    if np.max(error) < er_thresh:
        print("Test passed !")
    else:
        print("Test failed !")
        print(f"Max error = {np.max(error)}")
    return np.max(error)


if __name__ == "__main__":
    sigSB = 5.6704e-8
    Period = 88775.244 * 670  # [seconds]
    NSTEPS = 50000
    emiss = 1.0
    Fgeo = 0.2  # [W/m^2]

    STEPSPERSOL = 120
    dt = Period / STEPSPERSOL
    thIn = 120.0  # thermal inertia
    albedo = 0.2
    latitude = 5.0  # [degree]

    nz = 60
    zmax = 2.5
    zfac = 1.05

    rhocv = np.full(nz + 1, 1200.0 * 800.0)  # (density) * (heat capacity)
    delta = thIn / rhocv[1] * np.sqrt(Period / np.pi)
    ti = np.full(nz + 1, thIn)
    Rau = 1.52
    Decl = 0.0
    print("Skin depth= ", delta)
    print("Time step=", dt)
    print("zmax=", zmax)
    print("Thermal inertia=", thIn, " Period=", Period)
    print("Heat Flux at bottom boundary=", -Fgeo)

    T = np.full(nz + 1, 210.0)
    Tmean = np.zeros_like(T)
    z = grids.setgrid(nz, zmax, zfac)
    latitude = np.deg2rad(latitude)
    time = 0.0

    Fmean = 0.0
    HA = 0.0
    Qn = (1 - albedo) * flux_noatm(Rau, Decl, latitude, HA, 0.0, 0.0)
    Fsurf = 0.0

    Qns = np.zeros(NSTEPS + 1)  # storing all fluxs
    Qns[0] = Qn

    # Precalculating flux for all times
    # To not influence the speed of conductionQ
    for n in range(0, NSTEPS + 1):
        time = (n + 1) * dt  #   time at n+1;
        HA = 2 * np.pi * (time / Period % 1.0)  #  hour angle
        Qns[n] = (1 - albedo) * flux_noatm(Rau, Decl, latitude, HA, 0.0, 0.0)

    # SPEED MEASUREMENT
    for n in tqdm(range(0, NSTEPS + 1)):
        conductionQ(nz, z, dt, Qn, Qns[n], T, ti, rhocv, emiss, Fgeo, Fsurf)

    error = check_nochanges(T)
