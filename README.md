# Earth-Moon Patched-Conic Trajectory Solver
# Required Functions : kepler_E.m , sv_from_coe.m
# Analysis.m performs trade study
## Overview
This project is a custom MATLAB-based flight dynamics and trajectory optimization tool. It simulates Earth-to-Moon transit missions using a 3D patched-conic approximation. The primary purpose of this tool is to analyze and compare deep space mission architectures—specifically comparing classic Apollo-style Equatorial Low Lunar Orbits (LLO) against modern Artemis/Gateway-style Near Rectilinear Halo Orbits (NRHO).

## What It Does
The solver takes initial Earth departure conditions and lunar target parameters to compute the full mission trajectory from launch to Lunar Orbit Insertion (LOI). 

Key capabilities include:
- **Highly Elliptical Earth Staging:** Models departure from an elliptical phasing orbit to calculate fuel savings via the Oberth Effect.
- **Arrival Angle Optimization:** Uses a B-Plane targeting approximation to iteratively find the optimal lunar arrival angle for a specified perilune radius.
- **3D Keplerian Propagation:** Propagates orbits using Kepler's equations and rotates them into 3D inertial space using 3-1-3 Euler rotations to account for target Right Ascension of the Ascending Node (RAAN) and Inclination.
- **Delta-V (ΔV) Budgeting:** Calculates the impulsive velocity changes required for Trans-Lunar Injection (TLI) and Lunar Orbit Insertion (LOI).
- **Trade Study Analysis:** Exports mission data to conduct spacecraft mass-scaling sweeps using the Tsiolkovsky Rocket Equation.
- **Interactive Dashboard:** Generates a 4-panel 3D visualization of the Earth departure plane, geocentric transfer, lunar hyperbolic arrival, and final parking orbit.

## Why It Is Used
In professional spacecraft mission design, high-fidelity N-body numerical propagators (like GMAT, FreeFlyer, or STK) are highly sensitive. If provided with a random launch date and velocity, their differential correctors will fail to converge, and the simulated spacecraft will get lost in deep space.

This analytical solver bridges that gap. It provides the **mathematically sound "initial guess"** required to seed high-fidelity software. By rapidly generating baseline ΔV budgets, true anomalies, and required launch window phase angles, mission designers can evaluate the feasibility and mass limits of an architecture (Phase A Preliminary Design) before committing to computationally expensive numerical integration.

## Core Files
1. **`LunarPatchedConic.m`**: The primary solver. Handles the patched-conic math, optimization loops, 3D visualization dashboard, and `.mat` data export.
2. **`Analysis.m`**: The trade study script. Loads exported mission data (e.g., LLO vs. NRHO) to generate comparative bar charts and propellant vs. initial mass line graphs.
3. **Helper Functions**: Requires `kepler_E.m` (solves Kepler's equation for eccentric anomaly) and `sv_from_coe.m` (converts Classical Orbital Elements to Cartesian state vectors).

## Example Use Case: The NRHO Advantage
A standard run of this tool demonstrates why the Lunar Gateway targets an NRHO. By setting the target to a polar orbit (90° inclination) with a 70,000 km apolune, the tool reveals a ~50% reduction in LOI ΔV compared to a 100 km circular LLO. This massive fuel saving enables the delivery of heavy station modules to lunar orbit.
