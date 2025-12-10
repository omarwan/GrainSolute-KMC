# GrainSolute-KMC

Kinetic Monte Carlo implementation of a 2D q-state Potts model with mobile solutes. Transition events (grain-state flips and solute swaps) are sampled with Fenwick trees for O(log N) selection from the rate lists. Outputs are grid snapshots and timing data for post-processing.

## Repository layout
- `main.c` — simulation driver, time loop, and I/O.
- `grain_init.c`, `solute_init.c` — grid initialization (bulk/random/file-based).
- `transition_rates.c`, `energy_calculation.c`, `local_energy.c` — event construction and energetics.
- `fenwick.c` — Fenwick tree utilities for rate sampling.
- `csv_saver.c` — dumps grid states to CSV.
- `struc_gen.py` — generates sample grain/solute structure CSVs.
- `visualizer.py`, `OvitoConverter.py` — Python helpers for inspecting outputs.
- `file_gen.sh` — example batch runner that builds parameter sweeps and launches jobs.
- `params.txt` — example parameter file; `generated_data_256x1024/` holds matching structure CSVs.

## Requirements
- C11 compiler (`gcc` or `clang`) and standard math library.
- Python 3 (optional) with `numpy`, `matplotlib`, `pandas`, `ipywidgets` for structure generation and visualization.
- `screen` (optional) if you want to use `file_gen.sh`.

## Build
From the repository root:
```bash
gcc -std=c11 -O3 -Wall -Wextra \
  main.c affected_sites.c allocation.c boundary_sites.c csv_saver.c \
  energy_calculation.c fenwick.c grain_init.c local_energy.c periodic_neighbor.c \
  printer.c read_params.c rng.c solute_init.c transition_rates.c update_boundary.c \
  -lm -o potts
```
This produces the `potts` binary in the current directory.

## Running a simulation
1) Prepare a run directory containing `params.txt` (see below). The example `params.txt` here points to the provided `generated_data_256x1024` structures.  
2) Ensure the output folder referenced by `output_dir` exists (e.g., `mkdir -p output`).  
3) Run the solver, passing the directory that holds `params.txt`:
```bash
# use the current directory's params.txt
./potts ./

# or with a fixed seed for reproducibility
SEED=42 ./potts ./run1
```
The environment variable `SEED` overrides the RNG seed (default is `12345` if unset or invalid).

During the run, grids are written every `e_step` iterations to `output_dir/grains_grid_step_<step>.csv` and `output_dir/solutes_grid_step_<step>.csv`. A running log of simulated time is appended to `output_dir/simulated_time.txt` with columns:
```
Step, Flip Iteration, Flip Time, Swap Iteration, Swap Time, Simulated Time
```

## Parameter file (`params.txt`)
Key fields consumed by `read_params.c`:
- `grid_size_x`, `grid_size_y` — lattice dimensions.
- `init_type` — 1: all grains set to state 1; 2: random grain states in `[1, q_no]`; 3: read grain and solute CSVs from the paths below.
- `q_no` — total number of grain states; `q_selected` — preferred state used in the force/n4 terms for flips.
- `density` — solute fraction when `init_type` ≠ 3 (randomly placed with unique labels).
- `Jgg`, `Js`, `Jsg`, `Jss`, `Jsgs`, `F` — interaction/field parameters used in energy and rate calculations.
- `n_step` — total KMC iterations; `e_step` — interval for saving CSV snapshots.
- `T` — temperature in the Arrhenius transition probability.
- `Ed_oo`, `Eg_oo`, `Ecorr` — activation/energetic terms in the rate expressions.
- `output_dir` — folder where CSVs and timing are written (must exist).
- `grain_structure_file`, `solute_structure_file` — CSV paths used when `init_type=3`. Each CSV must have `grid_size_x` rows and `grid_size_y` comma-separated integers.

The provided `params.txt` is a working starting point for the sample 256×1024 case.

## Generating initial structures
`struc_gen.py` creates `generated_data_256x1024/structure_grains.csv` and a series of `structure_solutes_c<density>.csv` files. Adjust grid sizes/densities at the top of the script, then run:
```bash
python3 struc_gen.py
```

## Visualization and post-processing
- `visualizer.py` shows grain states with solute positions for saved steps. Update `directory`/`steps` in the script and run with Python.
- `OvitoConverter.py` converts saved grids to LAMMPS-style dump files for OVITO; adjust the `directory` and output paths as needed.

## Batch runs
`file_gen.sh` demonstrates how to sweep parameter combinations: it builds run directories with templated `params.txt` files and launches each case in a detached `screen` session via `./potts ./`. Edit the parameter arrays at the top before using it.

## Tips
- Keep `q_selected` within `1..q_no`; the code exits otherwise.
- Absolute paths in `grain_structure_file` and `solute_structure_file` avoid confusion when launching from other directories.
- Large grids and small `e_step` values create many CSVs; monitor disk usage.
