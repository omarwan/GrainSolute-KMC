#!/bin/bash

# Define arrays for each parameter that can have multiple values
GRID_SIZE_X_VALUES=(64)
GRID_SIZE_Y_VALUES=(64)
INIT_TYPE_VALUES=(1)
Q_NO_VALUES=(50)
DENSITY_VALUES=(0.2)
JGG_VALUES=(1.0)
JS_VALUES=(0.0)
JSG_VALUES=(-1.5 -1.6)
JSS_VALUES=(0.25)
JSGS_VALUES=(0.05)
N_STEPS_VALUES=(1000000000)
E_STEPS_VALUES=(10000000)
T_VALUES=(0.10)
ED_OO_VALUES=(1)
EG_OO_VALUES=(1)

# Create main output directory
MAIN_OUTPUT_DIR="./Runs"
mkdir -p "$MAIN_OUTPUT_DIR"

# Nested loops to iterate over each parameter value
for GRID_SIZE_X in "${GRID_SIZE_X_VALUES[@]}"; do
for GRID_SIZE_Y in "${GRID_SIZE_Y_VALUES[@]}"; do
for INIT_TYPE in "${INIT_TYPE_VALUES[@]}"; do
for Q_NO in "${Q_NO_VALUES[@]}"; do
for DENSITY in "${DENSITY_VALUES[@]}"; do
for JGG in "${JGG_VALUES[@]}"; do
for JS in "${JS_VALUES[@]}"; do
for JSG in "${JSG_VALUES[@]}"; do
for JSS in "${JSS_VALUES[@]}"; do
for JSGS in "${JSGS_VALUES[@]}"; do
for N_STEPS in "${N_STEPS_VALUES[@]}"; do
for E_STEPS in "${E_STEPS_VALUES[@]}"; do
for T in "${T_VALUES[@]}"; do
for ED_OO in "${ED_OO_VALUES[@]}"; do
for EG_OO in "${EG_OO_VALUES[@]}"; do

    # Create unique directory for each combination
    UNIQUE_DIR="$MAIN_OUTPUT_DIR/Grid${GRID_SIZE_X}x${GRID_SIZE_Y}_Q${Q_NO_VALUES}_Init${INIT_TYPE}_Density${DENSITY}_JGG${JGG}_JS${JS}_JSG${JSG}_JSS${JSS}_JSGS${JSGS}_Steps${N_STEPS}_${E_STEPS}_Temp${T}_EdOO${ED_OO}_EgOO${EG_OO}"
    mkdir -p "$UNIQUE_DIR"
    mkdir -p "$UNIQUE_DIR/output"

    # Create params.txt file in the unique folder
    PARAMS_FILE="$UNIQUE_DIR/params.txt"
    cat > "$PARAMS_FILE" <<EOF
# Simulation Grid Configuration
grid_size_x=$GRID_SIZE_X
grid_size_y=$GRID_SIZE_Y

# Simulation Initialization Parameters
init_type=$INIT_TYPE
q_no=$Q_NO
q_selected=1
density=$DENSITY

# Interaction Energy Parameters
Jgg=$JGG
Js=$JS
Jsg=$JSG
Jss=$JSS
Jsgs=$JSGS
F=0.0

# Simulation Steps Configuration
n_step=$N_STEPS
e_step=$E_STEPS

# Temperature Values
T=$T

# Diffusion Parameters
Ed_oo=$ED_OO
Eg_oo=$EG_OO
Ecorr=0.5
# Output Configuration
output_dir=output

# Files for structural information, applicable if INIT_TYPE=3
grain_structure_file=./structure_grains.csv
solute_structure_file=./structure_solutes_c$DENSITY.csv
EOF

screen -dmS "run_${GRID_SIZE_X}x${GRID_SIZE_Y}_Q${Q_NO_VALUES}_Init${INIT_TYPE}_Density${DENSITY}_JGG${JGG}_JS${JS}_JSG${JSG}_JSS${JSS}_JSGS${JSGS}_Steps${N_STEPS}_${E_STEPS}_Temp${T}_EdOO${ED_OO}_EgOO${EG_OO}" bash -c "cd $UNIQUE_DIR; ./potts ./"

cd ../../

done; done; done; done; done; done; done; done; done; done; done; done; done; done; done

# Additional information
echo "Simulation jobs have been submitted for all parameter combinations."