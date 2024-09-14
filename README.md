# Multithreaded Eulerian Fluid Simulation

## Overview
This project implements an Eulerian fluid solver in C++, capable of simulating fluid flow around objects such as a circle or an airfoil. The solver discretizes the fluid domain using a grid-based approach, solving the Euler equations numerically to capture the fluid's motion.


https://github.com/user-attachments/assets/f3e7131d-6818-44fa-ac38-f218c1b43f15


## Principles
The simulation is based on solving the Euler equations for inviscid flow using the finite difference method. A structured grid is used to discretize the domain, with each cell representing velocity and pressure values. The simulation progresses over time, solving for changes in these values across the grid.

## Optimizations
Writing this in C++ already offers excellent performance, but I plan to implement additional optimizaitons in the future, and this section will be updated accordingly.

## Results
The solver successfully simulates fluid flow in various setups, such as a wind tunnel around a circle or an airfoil, showing vortex shedding and flow separation where applicable. The results match expected physical behaviors for incompressible flows.

## Current Next Steps
- Clean-up code in main function, make everything more modular
- Allow importing of arbitrary boundary models (like from a PNG image)
- Basic optimizations (limiting memcpy, optimzing loops, etc)
- Add more visualizations (streamlines, pressure, vorticity, velocity)
- Deeper optimzation (utilizing AVX512 instructions for parellelization, and finally multithread everything)

### Resources Used
- https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)
- https://matthias-research.github.io/pages/tenMinutePhysics/17-fluidSim.pdf
- https://quangduong.me/notes/eulerian_fluid_sim_p1/
- https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
