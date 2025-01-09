# Fluid Simulation

Showcase video: https://youtu.be/HiwcsvRoU-Q
</br>Project created in Unity 2022.3, run on MacBook Air M1 2020

</br>Screenshot:

![SPH screenshot](https://github.com/ilialek/Resources/blob/main/SPH_GitHub.png)

This project implements a **Smoothed Particle Hydrodynamics (SPH)** fluid simulation using **Unity DOTS** for high-performance, multi-threaded computation. It includes features such as surface tension, viscosity, and collision handling with obstacles, powered by spatial hashing for efficient neighbor searches.

## Features
- **SPH Fluid Simulation**:
  - Particle-based simulation using the SPH method.
  - Calculates density, pressure, and forces based on fluid dynamics principles.
- **Physics Features**:
  - **Surface Tension**: Simulates cohesive forces between particles at the fluid's surface.
  - **Viscosity**: Models internal fluid friction for smooth and realistic behavior.
  - **Obstacle Collisions**: Particles interact with and are constrained by static obstacles.
- **Performance**:
  - **Unity DOTS** (ECS, Jobs, Burst) for multi-threaded, high-performance CPU calculations.
  - **Spatial Hashing** for efficient neighbor search and collision detection.

