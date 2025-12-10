# GrainSolute-KMC  
### A Kinetic Monte Carlo Solver for Grain Growth, Solute Segregation, and Full Nanocrystalline Stabilization

GrainSolute-KMC implements a fast Kinetic Monte Carlo (KMC) model that couples a q-state Potts model for crystallographic orientations with a lattice-gas model for solute thermodynamics and diffusion.  
This model is directly based on the formulation developed in the Acta Materialia publication:

O. Hussein & Y. Mishin, “A model of full thermodynamic stabilization of nanocrystalline alloys,” Acta Materialia 301 (2025) 121545.

---

## 1. Model Overview
This implementation reproduces the behavior of the full Potts + lattice-gas KMC model described in Hussein & Mishin (2025). The model uses:

### 1.1 Degrees of Freedom
- **Site orientation**: σₖ ∈ {1,…,q}  
- **Solute occupancy**: ξₖ ∈ {0,1}

---

## 2. Energy Function
E = E_cryst + E_sg + E_ss

### 2.1 Crystallographic (Potts) energy
E_cryst = Σₖ J_gg nₖ

### 2.2 Solute–GB interaction
E_sg = Σₖ ξₖ J_s + Σₖ ξₖ J_sg φ(nₖ)

### 2.3 Solute–solute interactions
E_ss = ½ Σ ξₖ ξₗ J_ss + ½ Σ ξₖ ξₗ J_ssg φ(nₖₗ)

---

## 3. Kinetic Monte Carlo Algorithm
Two event types:

1. **Orientation flips**  
2. **Solute jumps**

Transition rates follow harmonic TST:

ν_ij = ν₀ exp(−ε_ij / kT)

With nonlinear barrier law:

ε_ij =
• ε₀ exp(E_ij / 2ε₀), for E_ij ≤ 0  
• E_ij + ε₀ exp(−E_ij / 2ε₀), for E_ij > 0  

---

## 4. Physical Capabilities
Includes:

- Grain growth/shrinkage  
- Solute segregation  
- GB migration  
- GB-enhanced diffusion  
- Fully stabilized nanocrystalline states  
- Phase transitions and dynamic equilibrium states  

---

## 5. Repository Layout
(main.c, grain_init.c, solute_init.c, transition_rates.c, etc.)

---

## 6. Build Instructions
gcc -std=c11 -O3 -Wall -Wextra ...
Produces `potts`.

---

## 7. Running
./potts ./run_directory  
SEED=42 ./potts ./run1

---

## 8. Parameter File
Maps exactly to the model terms: J_gg, J_sg, J_s, J_ss, J_sgs, T, Ed_oo, Eg_oo, Ecorr, etc.

---

## 9. Visualization
Use visualizer.py or OvitoConverter.py.

---

## 10. Batch Runs
file_gen.sh automates parameter sweeps.

---

## 11. Citation
O. Hussein & Y. Mishin, Acta Materialia 301 (2025) 121545.
