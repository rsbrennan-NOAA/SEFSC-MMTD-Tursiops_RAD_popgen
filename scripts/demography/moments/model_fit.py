#!/usr/bin/env python3
"""
Moments demographic inference for Tursiops populations
"""
# mamba activate moments

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


MODEL_NUMBER = "07"
CURRENT_ROUND = 4
METHOD = "fmin"
SFS_FILE = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_downsample/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs"
RUN= "run2"

def get_best_rep_for_model(model_number):
    summary_file = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_expansion_{RUN}/round{CURRENT_ROUND}/model_{MODEL_NUMBER}_{METHOD}_downSample_summary.csv"
    df = pd.read_csv(summary_file)
    
    # Find the run with highest log-likelihood
    best_row = df.loc[df['LogLikelihood'].idxmax()]
    best_rep = best_row['Run']
    
    # Convert to integer first, then to zero-padded string- 01, 02, etc.
    best_rep_int = int(float(best_rep))  #both int and float
    best_rep_str = str(best_rep_int).zfill(2)
    
    print(f"Best rep for model {model_number} from round {CURRENT_ROUND}: {best_rep_int} (LL: {best_row['LogLikelihood']:.1f})")
    return best_rep_str

BEST_REP = get_best_rep_for_model(MODEL_NUMBER)

fs = moments.Spectrum.from_file(SFS_FILE)
new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)

#fs.mask[1, :, :, :] = True  # CAtlantic singletons
#fs.mask[:, 1, :, :] = True  # CGulf singletons  
#fs.mask[:, :, 1, :] = True  # Intermediate singletons
#fs.mask[:, :, :, 1] = True  # Offshore singletons


# look at model fits:
BEST_YAML_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments//mod7_expansion_{RUN}/round{CURRENT_ROUND}/output_yaml/model_{MODEL_NUMBER}_rep{BEST_REP}_{METHOD}_downSample-Round{CURRENT_ROUND}.yaml"
fs2 = moments.Spectrum.from_demes(BEST_YAML_PATH, sampled_demes= new_pop_ids,sample_sizes=[30,30,30,30])
#fs2.mask[1, :, :, :] = True  # CAtlantic singletons
#fs2.mask[:, 1, :, :] = True  # CGulf singletons  
#fs2.mask[:, :, 1, :] = True  # Intermediate singletons
#fs2.mask[:, :, :, 1] = True  # Offshore singletons



#plt.figure(figsize=(12, 8))  # width, height in inches
#moments.Plotting.plot_4d_comp_Poisson(fs2, fs)
#plt.savefig("/home/rbrennan/Tursiops-RAD-popgen/model_comparison.png",dpi=300, bbox_inches='tight')
#plt.show()

plt.close('all')

pop_names = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']

for i, pop_name in enumerate(pop_names):
    # Get marginals
    other_pops = [j for j in range(4) if j != i]
    fs2_marg = fs2.marginalize(other_pops)
    fs_marg = fs.marginalize(other_pops)
    
    moments.Plotting.plot_1d_comp_multinom(fs2_marg, fs_marg)
    plt.title(f'{pop_name} - Model vs Data')
    
    plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/model_{pop_name}_comparison.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Saved {pop_name}_comparison.png")

#fs2_pair = fs2.marginalize([2,3])  # CAtlantic vs CGulf
#fs_pair = fs_downsampled.marginalize([2,3])
#moments.Plotting.plot_2d_comp_multinom(fs2_pair, fs_pair)



# All possible pairs of populations
pairs = [
    ([2,3], 'CAtlantic vs CGulf'),
    ([1,3], 'CAtlantic vs Intermediate'), 
    ([1,2], 'CAtlantic vs Offshore'),
    ([0,3], 'CGulf vs Intermediate'),
    ([0,2], 'CGulf vs Offshore'),
    ([0,1], 'Intermediate vs Offshore')
]


# Create each plot separately and save individually
for i, (marginalize_pops, title) in enumerate(pairs):
    plt.figure(figsize=(8, 6))
    
    # Get 2D marginal SFS for this pair
    fs2_pair = fs2.marginalize(marginalize_pops)
    fs_pair = fs.marginalize(marginalize_pops)
    
    moments.Plotting.plot_2d_comp_multinom(fs2_pair, fs_pair, resid_range=3, vmin=1)
    plt.title(title, fontsize=14)
    
    # Save each plot separately
    filename = f"/home/rbrennan/Tursiops-RAD-popgen/model_2d_{title.replace(' vs ', '_vs_').replace(' ', '_')}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Saved: {filename}")


















### ---------------------------------
# # for this, I ran a bunch of independent reps for each model. need to go and grab the best one, then plot.

## Coastal, with expansion


#!/usr/bin/env python3
"""
Plot best models from coastal run directories
"""
import moments
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import os
import glob
from pathlib import Path

BASE_DIR = "analysis/moments"
ROUND = 4
SFS_FILE = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf.sfs"

def find_run_directories(base_path, prefix="coastal"):
    """Find all run directories with the specified prefix"""
    pattern = f"{prefix}_run*"
    search_pattern = os.path.join(base_path, pattern)
    run_dirs = glob.glob(search_pattern)
    run_dirs.sort()  # Sort for consistent ordering
    print(f"Found run directories: {run_dirs}")
    return run_dirs

def get_best_run_for_directory(dir_path):
    """Find the best run from a coastal directory"""
    summary_file = os.path.join(dir_path, f"round{ROUND}", "coastal_fmin_summary.csv")
    
    if not os.path.exists(summary_file):
        print(f"Warning: Summary file not found: {summary_file}")
        return None
        
    df = pd.read_csv(summary_file)
    
    # Find the run with highest log-likelihood
    best_row = df.loc[df['LogLikelihood'].idxmax()]
    best_run = best_row['Run']
    
    print(f"Best run for {os.path.basename(dir_path)}: {best_run} (LL: {best_row['LogLikelihood']:.1f})")
    return best_run, best_row['LogLikelihood']

# test that they work
run_dirs = find_run_directories("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/", "coastal")
for dir_path in run_dirs:
    result = get_best_run_for_directory(dir_path)
    if result:
        run, ll = result
        print(f"  -> Run {run}, LogLikelihood: {ll}")


# look at model fits:

RUN="run15"
BEST_REP="15"
INPUT_DIR="/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/"
MODEL="coastal_pop_bottle_expansion_mig_run15"
SUM_NAME="coastal_pop_bottle_expansion_mig_run15_fmin_summary.csv"

fs = moments.Spectrum.from_file(SFS_FILE)
new_pop_ids = ['CAtlantic', 'CGulf']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)

BEST_YAML_PATH = f"{INPUT_DIR}/{MODEL}_{RUN}/round4/output_yaml/{MODEL}_rep{BEST_REP}_fmin-Round4.yaml"
fs2 = moments.Spectrum.from_demes(BEST_YAML_PATH, sampled_demes= new_pop_ids,sample_sizes=[66,64])

# Create each marginal plot individually using the built-in method
pop_names = ['CAtlantic', 'CGulf']

for i, pop_name in enumerate(pop_names):
    # Get marginals
    other_pops = [j for j in range(2) if j != i]
    fs2_marg = fs2.marginalize(other_pops)
    fs_marg = fs.marginalize(other_pops)
    
    # Use the built-in plotting function (creates its own figure)
    moments.Plotting.plot_1d_comp_multinom(fs2_marg, fs_marg)
    plt.title(f'{pop_name} - Model vs Data')
    
    # Save each plot individually
    plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_coastal_{pop_name}_comparison.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Saved {pop_name}_comparison.png")

# Create each plot separately and save individually
plt.figure(figsize=(8, 6))

moments.Plotting.plot_2d_comp_multinom(fs2, fs, resid_range=3, vmin=1)
plt.title("coastal", fontsize=14)
    
# Save each plot separately
plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_coastal_2d_comparison.png", dpi=300, bbox_inches='tight')

    
# --------------------------------------------
# coastal no expansion

BASE_DIR = "analysis/moments"
ROUND = 4
SFS_FILE = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf.sfs"

def find_run_directories(base_path, prefix="coastal"):
    """Find all run directories with the specified prefix"""
    pattern = f"{prefix}_run*"
    search_pattern = os.path.join(base_path, pattern)
    run_dirs = glob.glob(search_pattern)
    run_dirs.sort()  # Sort for consistent ordering
    print(f"Found run directories: {run_dirs}")
    return run_dirs

def get_best_run_for_directory(dir_path, prefix="coastal"):
    """Find the best run from a coastal directory"""
    summary_file = os.path.join(dir_path, f"round{ROUND}", prefix + "_fmin_summary.csv")
    
    if not os.path.exists(summary_file):
        print(f"Warning: Summary file not found: {summary_file}")
        return None
        
    df = pd.read_csv(summary_file)
    
    # Find the run with highest log-likelihood
    best_row = df.loc[df['LogLikelihood'].idxmax()]
    best_run = best_row['Run']
    
    print(f"Best run for {os.path.basename(dir_path)}: {best_run} (LL: {best_row['LogLikelihood']:.1f})")
    return best_run, best_row['LogLikelihood']

# test that they work
run_dirs = find_run_directories("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/", "coastal_no_expansion")
for dir_path in run_dirs:
    result = get_best_run_for_directory(dir_path, "coastal_no_expansion")
    if result:
        run, ll = result
        print(f"  -> Run {run}, LogLikelihood: {ll}")

# all good. 
# look at model fits:

RUN="run7"
BEST_REP="15"
INPUT_DIR="/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/"
MODEL="coastal_no_expansion"
SUM_NAME="coastal_no_expansion_fmin_summary.csv"


fs = moments.Spectrum.from_file(SFS_FILE)
new_pop_ids = ['CAtlantic', 'CGulf']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)


BEST_YAML_PATH = f"{INPUT_DIR}/{MODEL}_{RUN}/round4/output_yaml/{MODEL}_rep{BEST_REP}_fmin-Round4.yaml"
fs2 = moments.Spectrum.from_demes(BEST_YAML_PATH, sampled_demes= new_pop_ids,sample_sizes=[66,64])

# Create each marginal plot individually using the built-in method
pop_names = ['CAtlantic', 'CGulf']

for i, pop_name in enumerate(pop_names):
    # Get marginals
    other_pops = [j for j in range(2) if j != i]
    fs2_marg = fs2.marginalize(other_pops)
    fs_marg = fs.marginalize(other_pops)
    
    # Use the built-in plotting function (creates its own figure)
    moments.Plotting.plot_1d_comp_multinom(fs2_marg, fs_marg)
    plt.title(f'{pop_name} - Model vs Data')
    
    # Save each plot individually
    plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_coastal_no_expansion_{pop_name}_comparison.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Saved {pop_name}_comparison.png")

# Create each plot separately and save individually
plt.figure(figsize=(8, 6))

moments.Plotting.plot_2d_comp_multinom(fs2, fs, resid_range=3, vmin=1)
plt.title("coastal", fontsize=14)
    
# Save each plot separately
plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_coastal_no_expansion_2d_comparison.png", dpi=300, bbox_inches='tight')

    







# --------------------------------------------
# coastal pop specific expansion

BASE_DIR = "analysis/moments"
ROUND = 4
SFS_FILE = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf.sfs"

def find_run_directories(base_path, prefix="coastal_pop_expansion"):
    """Find all run directories with the specified prefix"""
    pattern = f"{prefix}_run*"
    search_pattern = os.path.join(base_path, pattern)
    run_dirs = glob.glob(search_pattern)
    run_dirs.sort()  # Sort for consistent ordering
    print(f"Found run directories: {run_dirs}")
    return run_dirs

def get_best_run_for_directory(dir_path, prefix="coastal"):
    """Find the best run from a coastal directory"""
    summary_file = os.path.join(dir_path, f"round{ROUND}", prefix + "_fmin_summary.csv")
    
    if not os.path.exists(summary_file):
        print(f"Warning: Summary file not found: {summary_file}")
        return None
        
    df = pd.read_csv(summary_file)
    
    # Find the run with highest log-likelihood
    best_row = df.loc[df['LogLikelihood'].idxmax()]
    best_run = best_row['Run']
    
    print(f"Best run for {os.path.basename(dir_path)}: {best_run} (LL: {best_row['LogLikelihood']:.1f})")
    return best_run, best_row['LogLikelihood']

# test that they work
run_dirs = find_run_directories("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/", "coastal_pop_expansion")
for dir_path in run_dirs:
    result = get_best_run_for_directory(dir_path, "coastal_pop_pop_expansion")
    if result:
        run, ll = result
        print(f"  -> Run {run}, LogLikelihood: {ll}")

# all good. 
# look at model fits:

RUN="run15"
BEST_REP="15"
INPUT_DIR="/home/rbrennan/Tursiops-RAD-popgen/analysis/moments"
MODEL="coastal_pop_bottle_expansion_mig"
SUM_NAME="coastal_pop_bottle_expansion_mig_fmin_summary.csv"


fs = moments.Spectrum.from_file(SFS_FILE)
new_pop_ids = ['CAtlantic', 'CGulf']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)

BEST_YAML_PATH = f"{INPUT_DIR}/{MODEL}_{RUN}/round4/output_yaml/{MODEL}_rep{BEST_REP}_fmin-Round4.yaml"
fs2 = moments.Spectrum.from_demes(BEST_YAML_PATH, sampled_demes= new_pop_ids,sample_sizes=[66,64])

# Create each marginal plot individually using the built-in method
pop_names = ['CAtlantic', 'CGulf']

for i, pop_name in enumerate(pop_names):
    # Get marginals
    other_pops = [j for j in range(2) if j != i]
    fs2_marg = fs2.marginalize(other_pops)
    fs_marg = fs.marginalize(other_pops)
    
    # Use the built-in plotting function (creates its own figure)
    moments.Plotting.plot_1d_comp_multinom(fs2_marg, fs_marg)
    plt.title(f'{pop_name} - Model vs Data')
    
    # Save each plot individually
    plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_test_{pop_name}_comparison.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Saved {pop_name}_comparison.png")

# Create each plot separately and save individually
plt.figure(figsize=(8, 6))

moments.Plotting.plot_2d_comp_multinom(fs2, fs, resid_range=3, vmin=1)
plt.title("coastal", fontsize=14)
    
# Save each plot separately
plt.savefig(f"/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_coastal_pop_expansion_2d_comparison.png", dpi=300, bbox_inches='tight')

# plot
try:
graph = demes.load(BEST_YAML_PATH)
fig = plt.figure(figsize=(10, 8))
demesdraw.tubes(graph, ax=fig.add_subplot(111))
plt.savefig("/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/test.png", dpi=300)
print(f"Saved plot ")
except Exception as e:
    print(f"Failed to generate plot ")


