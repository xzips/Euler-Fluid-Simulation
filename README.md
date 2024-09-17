# Multithreaded, Cross-platform, Optimized Eulerian Fluid Simulation

## Overview
This project implements an Eulerian fluid solver in C++ and SFML, capable of simulating fluid flow around objects such as a circle or an airfoil. The solver discretizes the fluid domain using a grid-based approach, solving the Euler equations numerically to capture the fluid's motion. Several optimizations are implemented to speed up execution, including Red-Black Gauss-Seidel method, with all performance-critical functions running on all available CPU cores using OpenMP.

## Realtime Demo
https://github.com/user-attachments/assets/a9e8d0e7-eed9-4988-8715-9b5f4afdeb4a




## Principles
The simulation is based on solving the Euler equations for inviscid flow using the finite difference method. Each iteration of the simulation looks like this:
- Add or modify boundary conditions, i.e. dragging the model, setting inlet velocity, adding dye to see streamlines, etc.
- Solve the incompressibility condition (constant density everywhere) using 50-200 iterations of the Gauss-Seidel method over the entire simulation grid. This step takes around 90% of the compute time, and hence much time was spent optimizing it.
- Solve boundaries conditions at edges of simulation and for the body in the flow.
- Advect velocity from each cell to the next in the appropriate directions.
- Advect dye (visualization).
- Render results, handle user interface, etc.

## Optimizations
Fluid simulations are computationally expensive, and optimization was always a goal for this project. Currently, the following techniques are implemented:
- Writing everything in C++ and compiling with maximum auto-optimization offers good baseline performance
- OpenMP is used for multicore processing of all simulation aspects, greatly accelerating solve time. Benchmarks showed that the difference between OpenMP on/off on my 12-core Ryzen 3.8GHz system was on the order of **~10x performance gain** (451ms solve time to 44.5ms solve time).
- Using a slight variation on the traditional method, the **Gauss-Seidel Red-Black** approach to solving the incompressibility condition allows for parallelizing an otherwise purely sequential problem, which contributed to the 10x performance gain quoted above.
- Rendering functions operate on raw pixel and data arrays as much as possible, and only use expensive draw and display buffer operations at the very end of rendering where it's necessary. This improved frame render time from 25ms to 0.36ms (70x improvement) compared to the original function.

## Results
The solver successfully simulates fluid flow in various setups, such as a wind tunnel around a circle or an airfoil, showing vortex shedding and flow separation where applicable. The results match expected physical behaviors for incompressible flows. Additionally, since C++ and SFML is multiplatform, **I was able to compile and run the simulation in realtime on my smartphone device with only minor tweaks.**

### Cross-Platform Demo

https://github.com/user-attachments/assets/654c3dff-21dc-403d-966d-42f1c26822f5

The performance on mobile is degraded a bit since I didn't bother to install OpenMP for this demo, but it would run significantly faster (at a higher simulation resolution more specifically) if it could utilize all the ARM cores in my phone.

## Next Steps
- I hope to explore more advanced algorithms for solving incompressibility condition that are more mathematically efficient, as I believe I've currently maxed out the performance from programmatic optimizations such as multithreading. Two possible candidates include multigrid methods, and conjugate gradient methods, both of which will require significant time investment to understand and implement correctly.
- Adding more visualization would be nice, such as computed streamlines, pressure, vorticity, and velocity. These are relatively simple to add, but would just take a bit of time get the colors and graphics right.
- Finally, adding multiple bodies at the same time, or perhaps runtime-editable boundary conditions or fluid sources could be interesting.

- 
### Resources for Initial Implementation
- https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)
- https://matthias-research.github.io/pages/tenMinutePhysics/17-fluidSim.pdf
- https://quangduong.me/notes/eulerian_fluid_sim_p1/
- https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

### Resources for Optimization
- https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
- https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
- https://www.wias-berlin.de/people/john/LEHRE/MULTIGRID/multigrid.pdf
- https://math.mit.edu/classes/18.086/2006/am63.pdf
