#!/usr/bin/env python3
"""
Moments demographic inference
Uses best-fit parameters from previous round as starting points with perturbation
"""

import moments
import numpy as np
import pandas as pd
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
parser = argparse.ArgumentParser(description="Run demographic inference with configurable parameters")
parser.add_argument("model_number", help="Model number to specify the model (e.g., '01', '02').")
parser.add_argument("rep_number", help="Repetition number for the model run.")
parser.add_argument("round_number", type=int, help="Current round number.")
parser.add_argument("perturb", type=int, help="Perturbation parameter.")
parser.add_argument("maxiter", type=int, help="Maximum iterations.")
parser.add_argument("method", help="Optimization method (e.g., 'fmin').")
parser.add_argument("sfs_file", help="Path to the SFS file.")

args = parser.parse_args()

MODEL_NUMBER = args.model_number 
REPNUMBER = args.rep_number 
CURRENT_ROUND = args.round_number
PERTURB = args.perturb
MAXITER = args.maxiter
METHOD = args.method
SFS_FILE = args.sfs_file

# Calculate previous round after getting current round from arguments
PREVIOUS_ROUND = CURRENT_ROUND - 1


def get_best_rep_for_model(model_number):
    """find the best rep from previous round"""
    summary_file = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/round{PREVIOUS_ROUND}/model_{MODEL_NUMBER}_fmin_downSample_summary.csv"
    df = pd.read_csv(summary_file)
    
    # Find the run with highest log-likelihood
    best_row = df.loc[df['LogLikelihood'].idxmax()]
    best_rep = best_row['Run']
    
    # Convert to integer first, then to zero-padded string- 01, 02, etc.
    best_rep_int = int(float(best_rep))  #both int and float
    best_rep_str = str(best_rep_int).zfill(2)
    
    print(f"Best rep for model {model_number} from round {PREVIOUS_ROUND}: {best_rep_int} (LL: {best_row['LogLikelihood']:.1f})")
    return best_rep_str


BEST_REP = get_best_rep_for_model(MODEL_NUMBER)


BEST_YAML_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/round{PREVIOUS_ROUND}/output_yaml/model_{MODEL_NUMBER}_rep{BEST_REP}_{METHOD}_downSample-Round{PREVIOUS_ROUND}.yaml"
OPTIONS_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/options_fourpop_{MODEL_NUMBER}_Simple.yaml"
OUTPUT_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/round{CURRENT_ROUND}/output_yaml/model_{MODEL_NUMBER}_rep{REPNUMBER}_{METHOD}_downSample-Round{CURRENT_ROUND}.yaml"

PLOT_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/mod7_simple/round{CURRENT_ROUND}/round{CURRENT_ROUND}_optimization_model_{MODEL_NUMBER}_{REPNUMBER}_{METHOD}_downSample.png"
SUMMARY_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/round{CURRENT_ROUND}/model_{MODEL_NUMBER}_{METHOD}_downSample_summary.csv"

# Create output directories
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/round{CURRENT_ROUND}", exist_ok=True)
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/round{CURRENT_ROUND}/output_yaml", exist_ok=True)
os.makedirs(f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/mod7_simple/round{CURRENT_ROUND}/", exist_ok=True)

# Print file paths for verification
print(f"Best YAML file: {BEST_YAML_PATH}")
print(f"Options file: {OPTIONS_PATH}")
print(f"Output file: {OUTPUT_PATH}")

# Record start time
start_time = time.time()
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"[{current_time}] Starting Round {CURRENT_ROUND} run for model {MODEL_NUMBER}, rep {REPNUMBER}")

# Load best YAML file from previous round
print(f"Using best fit from Round {PREVIOUS_ROUND} model {MODEL_NUMBER}, rep {BEST_REP}")

# Load the SFS
print(f"\nLoading SFS...")
fs = moments.Spectrum.from_file(SFS_FILE)

# Fix population IDs order
new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
print(f"Population IDs: {fs.pop_ids}")

# Run optimization
print(f"\nStarting Round {CURRENT_ROUND} model optimization for model {MODEL_NUMBER}, rep {REPNUMBER}")
print(f"Using {MAXITER} iterations with {PERTURB} perturbation")

ret = moments.Demes.Inference.optimize(
    BEST_YAML_PATH,
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
print(f"\nOptimization completed!")
print(f"Log-likelihood: {-LL}")
print("Best fit parameters:")

for n, p in zip(param_names, opt_params):
    line = f"{n}\t{p}"
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

print(f"\nResults appended to summary file: {SUMMARY_FILE}")

# Calculate and print timing
end_time = time.time()
end_time_print = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
elapsed_time = end_time - start_time

print(f"\nRound {CURRENT_ROUND} analysis for model {MODEL_NUMBER}, rep {REPNUMBER} ended at {end_time_print}")
print(f"Analysis completed in {elapsed_time / 3600:.2f} hours")

# Generate plot
print(f"\nGenerating demographic plot...")
graph = demes.load(OUTPUT_PATH)
fig = plt.figure(figsize=(10, 8))
demesdraw.tubes(graph, ax=fig.add_subplot(111))
plt.title(f"Round {CURRENT_ROUND} - Model {MODEL_NUMBER}, Rep {REPNUMBER}\nLog-likelihood: {-LL:.1f}")
plt.savefig(PLOT_FILE, dpi=300, bbox_inches='tight')
print(f"Saved plot to {PLOT_FILE}")
plt.close()


print(f"\nRound {CURRENT_ROUND} complete for model {MODEL_NUMBER}, rep {REPNUMBER}")