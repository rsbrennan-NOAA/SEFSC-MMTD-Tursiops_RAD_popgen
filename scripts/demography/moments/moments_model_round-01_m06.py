#!/usr/bin/env python3
"""
Moments demographic inference for Tursiops populations - Model 06
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
parser = argparse.ArgumentParser(description="Run demographic inference for Tursiops populations.")
parser.add_argument("model_number", help="Model number to specify the model (e.g., '01', '02').")
parser.add_argument("rep_number", help="Repetition number for the model run.")  # REPNUMBER as string

args = parser.parse_args()

MODEL_NUMBER = args.model_number 
REPNUMBER = args.rep_number 

# define file paths
PERTURB = 0
MAXITER = 15
METHOD = "fmin"

BASE_DEME_GRAPH_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/fourpop_{MODEL_NUMBER}.yaml"
REP_DEME_GRAPH_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/round_01/fourpop_{MODEL_NUMBER}_rep{REPNUMBER}.yaml"
OPTIONS_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/options_fourpop_{MODEL_NUMBER}_Simple.yaml"
OUTPUT_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/initialRuns/output_yaml/model_{MODEL_NUMBER}_rep{REPNUMBER}_{METHOD}_downSample.yaml"

OUTPUT_RESULTS_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/initialRuns/optimization_model_{MODEL_NUMBER}_{REPNUMBER}_{METHOD}_downSample_results.txt"
PLOT_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/optimization_model_{MODEL_NUMBER}_{REPNUMBER}_{METHOD}_downSample.png"
SUMMARY_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/initialRuns/model_{MODEL_NUMBER}_{METHOD}_downSample_summary.csv"

# Print file paths for verification
print(f"Deme graph file: {BASE_DEME_GRAPH_PATH}")
print(f"rep deme file: {REP_DEME_GRAPH_PATH}")
print(f"options file: {OPTIONS_PATH}")
print(f"Output file: {OUTPUT_PATH}")

# record start time
start_time = time.time()
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"[{current_time}] Starting moments analysis for model {MODEL_NUMBER}, rep {REPNUMBER}")

# Generate random starting values

# Generate hierarchical time parameters for Model 02
# Split1 (Ancestor): randomly choose from wide range
split1_val = int(np.random.uniform(10000, 100000))
print(f"Split1 chosen: {split1_val}")

# Split2 (Coastal): randomly choose between 500 and Admix_Time
split2_val = int(np.random.uniform(500, split1_val))
print(f"Split2 chosen: {split2_val}")

# Admix_Time: randomly choose between Split2 and Split1
admix_time_val = int(np.random.uniform(split2_val, split1_val))
print(f"Admix_Time chosen: {admix_time_val} ")

# Admixture proportion (Offshore contribution)
admix_prop_val = np.random.uniform(0.1, 0.9)
print(f"Admix_Prop chosen: {admix_prop_val:.3f}")

# Generate random population sizes
n_ancestor = int(np.random.uniform(100, 30000))
n_offshore = int(np.random.uniform(100, 30000))
n_coastal = int(np.random.uniform(100, 30000))
n_intermediate = int(np.random.uniform(100, 30000))
n_catlantic = int(np.random.uniform(100, 30000))
n_cgulf = int(np.random.uniform(100, 30000))

print(f"\nPopulation sizes chosen:")
print(f"N_Ancestor: {n_ancestor}")
print(f"N_Offshore: {n_offshore}")
print(f"N_Coastal: {n_coastal}")
print(f"N_Intermediate: {n_intermediate}")
print(f"N_CAtlantic: {n_catlantic}")
print(f"N_CGulf: {n_cgulf}")


# Generate random migration rates (8 migration parameters)
migration_rates = []
migration_names = [
    "m_CGulf_to_Int", "m_Int_to_CGulf", "m_offshore_to_int", "m_int_to_offshore",
    "m_offshore_to_CGulf", "m_CGulf_to_offshore", "m_intermediate_to_coastal", "m_coastal_to_intermediate"
]

print(f"\nMigration rates chosen:")
for i, name in enumerate(migration_names):
    log_rate = np.random.uniform(np.log10(1e-10), np.log10(0.001))
    rate = 10**log_rate
    migration_rates.append(rate)
    print(f"{name}: {rate:.2e}")

# Load and modify the base YAML file
print(f"\nCreating replicate-specific YAML file...")
with open(BASE_DEME_GRAPH_PATH, 'r') as f:
    yaml_data = yaml.safe_load(f)

# Update the YAML with random starting values for Model 02
for deme in yaml_data['demes']:
    if deme['name'] == 'Ancestor':
        deme['epochs'][0]['end_time'] = split1_val
        deme['epochs'][0]['start_size'] = n_ancestor
    elif deme['name'] == 'Offshore':
        deme['epochs'][0]['start_size'] = n_offshore
    elif deme['name'] == 'Coastal':
        deme['epochs'][0]['end_time'] = split2_val
        deme['epochs'][0]['start_size'] = n_coastal
    elif deme['name'] == 'CAtlantic':
        deme['epochs'][0]['start_size'] = n_catlantic
    elif deme['name'] == 'CGulf':
        deme['epochs'][0]['start_size'] = n_cgulf
    elif deme['name'] == 'Intermediate':
        deme['start_time'] = admix_time_val
        deme['proportions'] = [admix_prop_val, 1-admix_prop_val]
        deme['epochs'][0]['start_size'] = n_intermediate



# Update migration rates - convert to Python float to avoid numpy serialization
for i, migration in enumerate(yaml_data['migrations']):
    migration['rate'] = float(migration_rates[i])  # Add float() conversion


# Save the modified YAML for this replicate
with open(REP_DEME_GRAPH_PATH, 'w') as f:
    yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False, width=float('inf'))

print(f"Created replicate-specific YAML: {REP_DEME_GRAPH_PATH}")

# load the sfs 
print(f"Loading SFS")
fs = moments.Spectrum.from_file("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_downsample/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs")
    
# pop ids are in wrong order- fix them
new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
print(f"Population IDs: {fs.pop_ids}")

# Run optimization
print(f"Starting demographic model optimization for model {MODEL_NUMBER}, rep {REPNUMBER}")
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
print(f"Analysis for model {MODEL_NUMBER}, rep {REPNUMBER} ended at {end_time_print}")
print(f"Analysis completed for model {MODEL_NUMBER}, rep {REPNUMBER} in {elapsed_time / 3600:.2f} hours")


# plot
try:
    graph = demes.load(OUTPUT_PATH)
    fig = plt.figure(figsize=(10, 8))
    demesdraw.tubes(graph, ax=fig.add_subplot(111))
    plt.savefig(PLOT_FILE, dpi=300)
    print(f"Saved plot to {PLOT_FILE}")
except Exception as e:
    print(f"Failed to generate plot for model {MODEL_NUMBER}, rep {REPNUMBER}: {e}")

