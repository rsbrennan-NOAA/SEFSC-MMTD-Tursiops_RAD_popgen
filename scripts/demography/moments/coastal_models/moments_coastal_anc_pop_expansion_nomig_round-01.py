#!/usr/bin/env python3
"""
Moments demographic inference for Tursiops populations - Coastal only
"""

import moments
import numpy as np
import matplotlib.pylab as plt
import demes
import demesdraw
import os
import matplotlib.cm as cm
import time
from datetime import datetime
import csv
import argparse
import yaml

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run demographic inference for coastal split only. Ancestral expansion with migration")
parser.add_argument("rep_number", help="Repetition number for the model run.")  # REPNUMBER as string
parser.add_argument("model_name", help="Model name (e.g., 'coastal_anc_expansion_mig').")

args = parser.parse_args()

REPNUMBER = args.rep_number 
MODEL_NAME = args.model_name

# define file paths
PERTURB = 0
MAXITER = 5
METHOD = "fmin"

BASE_DEME_GRAPH_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models/demographic_models/fourpop_{MODEL_NAME}.yaml"
REP_DEME_GRAPH_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models/{MODEL_NAME}/fourpop_{MODEL_NAME}_rep{REPNUMBER}.yaml"
OPTIONS_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models/option_files/options_fourpop_{MODEL_NAME}.yaml"
OUTPUT_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/{MODEL_NAME}/output_yaml/{MODEL_NAME}_rep{REPNUMBER}_{METHOD}.yaml"

OUTPUT_RESULTS_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/{MODEL_NAME}/optimization_{MODEL_NAME}_{REPNUMBER}_{METHOD}_results.txt"
PLOT_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/{MODEL_NAME}/optimization_{MODEL_NAME}_{REPNUMBER}_{METHOD}.png"
SUMMARY_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/{MODEL_NAME}/{MODEL_NAME}_{METHOD}_summary.csv"

# Create output directories if they don't exist
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/{MODEL_NAME}", exist_ok=True)
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/{MODEL_NAME}/output_yaml", exist_ok=True)
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/{MODEL_NAME}", exist_ok=True)
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models/{MODEL_NAME}", exist_ok=True)


# Print file paths for verification
print(f"Deme graph file: {BASE_DEME_GRAPH_PATH}")
print(f"rep deme file: {REP_DEME_GRAPH_PATH}")
print(f"options file: {OPTIONS_PATH}")
print(f"Output file: {OUTPUT_PATH}")

# record start time
start_time = time.time()
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"[{current_time}] Starting moments analysis for coastal, rep {REPNUMBER}")

# Generate random starting values

# Bottleneck start
bottleneck_time_val = int(np.random.uniform(5000, 50000))
# Split time (must be < bottleneck_time)
split_time_val = int(np.random.uniform(500, min(20000, bottleneck_time_val)))
# Expansion start: when descendant populations start expanding (must be < split_time)
expansion_start_val = int(np.random.uniform(300, min(15000, split_time_val)))

print(f"Bottleneck time: {bottleneck_time_val}")
print(f"Split time: {split_time_val}")
print(f"Expansion start time: {expansion_start_val}")

n_ancestor = int(np.random.uniform(500, 40000)) # Ancestor initial size
n_bottleneck = int(np.random.uniform(100, 5000)) # Bottleneck size
n_cgulf_start = int(np.random.uniform(50, 5000))  
n_catlantic_start = int(np.random.uniform(50, 5000))
n_catlantic_end = int(np.random.uniform(n_bottleneck, 40000)) # CAtlantic final size
n_cgulf_end = int(np.random.uniform(n_bottleneck, 40000)) # CGulf final size

print(f"N_Ancestor: {n_ancestor}")
print(f"N_bottleneck: {n_bottleneck}")
print(f"N_CGulf_start: {n_cgulf_start}")
print(f"N_CAtlantic_start: {n_catlantic_start}")
print(f"N_CAtlantic_end: {n_catlantic_end}")
print(f"N_CGulf_end: {n_cgulf_end}")

# Load and modify the base YAML file
print(f"\nCreating replicate-specific YAML file...")
with open(BASE_DEME_GRAPH_PATH, 'r') as f:
    yaml_data = yaml.safe_load(f)

# Update the YAML with random starting values
for deme in yaml_data['demes']:
    if deme['name'] == 'Ancestor':
        # Epoch 0: Initial population before bottleneck
        deme['epochs'][0]['end_time'] = bottleneck_time_val
        deme['epochs'][0]['start_size'] = n_ancestor   
        # Epoch 1: Bottleneck continues until split
        deme['epochs'][1]['end_time'] = split_time_val 
        deme['epochs'][1]['start_size'] = n_bottleneck 
    elif deme['name'] == 'CAtlantic':
        # Epoch 0: Constant size after split until expansion
        deme['epochs'][0]['end_time'] = expansion_start_val
        deme['epochs'][0]['start_size'] = n_catlantic_start 
        # Epoch 1: Expansion to final size
        deme['epochs'][1]['end_size'] = n_catlantic_end 
    elif deme['name'] == 'CGulf':
        # Epoch 0: Constant size after split until expansion
        deme['epochs'][0]['end_time'] = expansion_start_val 
        deme['epochs'][0]['start_size'] = n_cgulf_start
        # Epoch 1: Expansion to final size
        deme['epochs'][1]['end_size'] = n_cgulf_end   


# Save the modified YAML for this replicate
with open(REP_DEME_GRAPH_PATH, 'w') as f:
    yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False, width=float('inf'))

print(f"Created replicate-specific YAML: {REP_DEME_GRAPH_PATH}")

# load the sfs 
print(f"Loading SFS")
fs = moments.Spectrum.from_file("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf.sfs")
    
# pop ids are in wrong order- fix them
new_pop_ids = ['CAtlantic', 'CGulf']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
print(f"Population IDs: {fs.pop_ids}")

# Run optimization
print(f"Starting demographic model optimization for {MODEL_NAME}, rep {REPNUMBER}")
ret = moments.Demes.Inference.optimize(
    REP_DEME_GRAPH_PATH,
    OPTIONS_PATH,
    fs,
    method=METHOD,
    output=OUTPUT_PATH,
    overwrite=True,
    verbose=1,
    perturb=PERTURB,
    maxiter=MAXITER
)

# Process results
param_names, opt_params, LL = ret
print(f"Log-likelihood: {-LL}")
print("Best fit parameters:")

# save to file
#with open(OUTPUT_RESULTS_FILE, "w") as f:
#    f.write(f"Log-likelihood: {-LL}\n")
#    f.write("Best fit parameters\n")
for n, p in zip(param_names, opt_params):
    line = f"{n}\t{p}"
#        f.write(line + "\n")
    print(line)

# Create a dictionary of results to save to the CSV summary file
results_dict = {
    'Run': REPNUMBER,
    'LogLikelihood': -LL
}

# Add all parameters to the dictionary
for n, p in zip(param_names, opt_params):
    results_dict[n] = p

# Check if summary file exists and create/append appropriately
file_exists = os.path.isfile(SUMMARY_FILE)

with open(SUMMARY_FILE, 'a', newline='') as csvfile:
    fieldnames = ['Run', 'LogLikelihood'] + list(param_names)
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    # Write header only if the file doesn't exist yet
    if not file_exists:
        writer.writeheader()
    
    writer.writerow(results_dict)

print(f"Results appended to summary file: {SUMMARY_FILE}")

end_time = time.time()
end_time_print = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

elapsed_time = end_time - start_time
print(f"Analysis for {MODEL_NAME}, rep {REPNUMBER} ended at {end_time_print}")
print(f"Analysis completed for {MODEL_NAME}, rep {REPNUMBER} in {elapsed_time / 3600:.2f} hours")


# plot
try:
    graph = demes.load(OUTPUT_PATH)
    fig = plt.figure(figsize=(10, 8))
    demesdraw.tubes(graph, ax=fig.add_subplot(111))
    plt.savefig(PLOT_FILE, dpi=300)
    print(f"Saved plot to {PLOT_FILE}")
except Exception as e:
    print(f"Failed to generate plot for {MODEL_NAME} rep {REPNUMBER}: {e}")

