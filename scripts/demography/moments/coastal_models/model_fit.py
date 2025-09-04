#!/usr/bin/env python3
"""
Moments demographic inference for Tursiops populations
"""
# mamba activate moments

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# edit to make more efficient below here. 

# I have 14 models. 

# 1. Bottleneck in ancestor, expansion in ancestor- Migration
# - `coastal_anc_expansion_mig`
# 2. Bottleneck in ancestor, expansion in ancestor- no Migration
# - `coastal_anc_expansion_nomig`
# 3. Bottleneck in ancestor, instant expansion in ancestor- Migration
# - `coastal_anc_inst_expansion_mig`
# 4. Bottleneck in ancestor, instant expansion in ancestor- no Migration
# - `coastal_anc_inst_expansion_nomig`
# 5.  Bottleneck in ancestor, expansion in populations- Migration
# - `coastal_anc_pop_expansion_mig`
# - `submit_snake_$MODEL.sh` - DONE
# 6.  Bottleneck in ancestor, expansion in populations- no Migration
# - `coastal_anc_pop_expansion_nomig`
# 7. Bottleneck in ancestor, instant expansion in populations- Migration
# - `coastal_anc_pop_inst_expansion_mig`
# 8. Bottleneck in ancestor, instant expansion in populations- no Migration
# - `coastal_anc_pop_inst_expa
# 9. Bottleneck in populations, expansion in populations- Migration
# - `coastal_pop_bottle_expans
# 10. Bottleneck in populations, expansion in populations- no Migration
# - `coastal_pop_bottle_expans
# 11. Bottleneck in populations, instant expansion in populations- Migration
# - `coastal_pop_bottle_inst_e
# 12. Bottleneck in populations, instant expansion in populations- no Migration
# - `coastal_pop_bottle_inst_e
# 13. No expansion, with migration
# - `coastal_no_expansion_mig`
# 14. No expansion, no migration. 
# - `coastal_no_expansion_nomig`

import os
import glob
import pandas as pd
import moments
import matplotlib.pyplot as plt

BASE_DIR = "analysis/moments"
ROUND = 4
SFS_FILE = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf.sfs"
INPUT_DIR = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments"

MODELS = [
    "coastal_anc_expansion_mig",
    "coastal_anc_expansion_nomig", 
    "coastal_anc_inst_expansion_mig",
    "coastal_anc_inst_expansion_nomig",
    "coastal_anc_pop_expansion_mig",
    "coastal_anc_pop_expansion_nomig",
    "coastal_anc_pop_inst_expansion_mig",
    "coastal_anc_pop_inst_expansion_nomig",
    "coastal_pop_bottle_expansion_mig",
    "coastal_pop_bottle_expansion_nomig",
    "coastal_pop_bottle_inst_expansion_mig",
    "coastal_pop_bottle_inst_expansion_nomig",
    "coastal_no_expansion_mig",
    "coastal_no_expansion_nomig"
]

def find_run_directories(base_path, prefix):
    """Find all run directories with the specified prefix"""
    pattern = f"{prefix}_run*"
    search_pattern = os.path.join(base_path, pattern)
    run_dirs = glob.glob(search_pattern)
    run_dirs.sort()
    print(f"Found {len(run_dirs)} run directories for {prefix}")
    return run_dirs



def get_best_run_for_model(model_name, base_path, round_num=4):
    """Find the best run across all run directories for a given model"""
    run_dirs = find_run_directories(base_path, model_name)
    
    best_overall_ll = float('-inf')
    best_overall_info = None
    
    for dir_path in run_dirs:
        result = get_best_run_for_directory(dir_path, model_name)
        if result:
            run, ll = result
            if ll > best_overall_ll:
                best_overall_ll = ll
                best_overall_info = {
                    'model': model_name,
                    'directory': dir_path,
                    'run': run,
                    'log_likelihood': ll
                }
    
    if best_overall_info:
        print(f"Best for {model_name}: Run {best_overall_info['run']} (LL: {best_overall_info['log_likelihood']:.1f})")
    
    return best_overall_info

def get_best_run_for_directory(dir_path, model_name):
    """Find the best run from a model directory"""
    summary_file = os.path.join(dir_path, f"round{ROUND}", f"{model_name}_fmin_summary.csv")
    
    if not os.path.exists(summary_file):
        print(f"Warning: Summary file not found: {summary_file}")
        return None
        
    df = pd.read_csv(summary_file)
    
    # Convert LogLikelihood to numeric, handling any non-numeric values
    df['LogLikelihood'] = pd.to_numeric(df['LogLikelihood'], errors='coerce')
        
    # Find the run with highest log-likelihood
    best_row = df.loc[df['LogLikelihood'].idxmax()]
    best_run = best_row['Run']
    
    return best_run, best_rep, best_row['LogLikelihood']

# test
best_runs_per_model = {}

for model in MODELS:
    best_info = get_best_run_for_model(model, INPUT_DIR)
    if best_info:
        best_runs_per_model[model] = best_info


# now find overall best model:

def find_overall_best_model(best_runs_dict):
    """Find the model with the highest log-likelihood across all models"""
    if not best_runs_dict:
        print("No models to compare")
        return None
    
    best_overall = None
    best_ll = float('-inf')
    
    print("\n" + "="*80)
    print("MODEL COMPARISON - All Best Runs:")
    print("="*80)
    
    # Sort models by log-likelihood for nice display
    sorted_models = sorted(best_runs_dict.items(), 
                          key=lambda x: x[1]['log_likelihood'], 
                          reverse=True)
    
    for i, (model_name, model_info) in enumerate(sorted_models, 1):
        ll = model_info['log_likelihood']
        run = model_info['run']
        
        print(f"{i:2d}. {model_name:<40} LL: {ll:8.1f} (Run {run})")
        
        if ll > best_ll:
            best_ll = ll
            best_overall = model_info
    
    print("="*80)
    if best_overall:
        print(f"OVERALL BEST MODEL: {best_overall['model']}")
        print(f"Log-Likelihood: {best_overall['log_likelihood']:.1f}")
        print(f"Run: {best_overall['run']}")
        print(f"Directory: {best_overall['directory']}")
    print("="*80)
    
    return best_overall

# Find the overall best model
overall_best = find_overall_best_model(best_runs_per_model)


def create_model_plots(model_info, fs_observed, output_dir="/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit"):
    """Create 1D marginal and 2D comparison plots for a single model"""
    model_name = model_info['model']
    run = model_info['run']
    ll = model_info['log_likelihood']
    directory = model_info['directory']
    
    print(f"Creating plots for {model_name}...")
    
    # Extract the run directory name (e.g., "run14" from the directory path)
    run_dir = os.path.basename(directory)  # Gets something like "coastal_pop_bottle_inst_expansion_mig_run14"
    run_number = run_dir.split('_')[-1]  # Gets "run14"
    
    # Construct path with zero-padded rep number
    yaml_path = f"{INPUT_DIR}/{model_name}_{run_number}/round4/output_yaml/{model_name}_rep{int(run):02d}_fmin-Round4.yaml"
    
    print(f"Looking for YAML file: {yaml_path}")
    
    if not os.path.exists(yaml_path):
        print(f"Warning: YAML file not found: {yaml_path}")
        return False
    
    try:
        #  demographic plot
        print(f"\nGenerating demographic plot...")
        graph = demes.load(yaml_path)
        fig = plt.figure(figsize=(10, 8))
        demesdraw.tubes(graph, ax=fig.add_subplot(111))
        plt.title(f"{model_name} - Rep {int(run)}\nLog-likelihood: {ll:.1f}")
        
        filename = f"model_{model_name}_demography_LL{ll:.1f}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved demographic plot to {filename}")
        plt.show()
        plt.close()
        
        # Load model prediction for SFS comparison
        fs_model = moments.Spectrum.from_demes(yaml_path, 
                                             sampled_demes=['CAtlantic', 'CGulf'],
                                             sample_sizes=[66, 64])
       # Create 1D marginal plots
        pop_names = ['CAtlantic', 'CGulf']
        for i, pop_name in enumerate(pop_names):
            other_pops = [j for j in range(2) if j != i]
            fs_model_marg = fs_model.marginalize(other_pops)
            fs_obs_marg = fs_observed.marginalize(other_pops)
            
            moments.Plotting.plot_1d_comp_multinom(fs_model_marg, fs_obs_marg)
            plt.title(f'{model_name} - {pop_name}\nLL: {ll:.1f}', fontsize=12)
            
            filename = f"model_{model_name}_{pop_name}_1D_LL{ll:.1f}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.show()
            print(f"  Saved {filename}")
        
        # Create 2D comparison plot
        plt.figure(figsize=(8, 6))
        moments.Plotting.plot_2d_comp_multinom(fs_model, fs_observed, resid_range=3, vmin=1)
        plt.title(f'{model_name}\nLL: {ll:.1f}', fontsize=14)
        
        filename = f"model_{model_name}_2D_LL{ll:.1f}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"  Saved {filename}")
        
        return True
        
    except Exception as e:
        print(f"Error creating plots for {model_name}: {str(e)}")
        return False

# Load observed data
fs_observed = moments.Spectrum.from_file(SFS_FILE)
new_pop_ids = ['CAtlantic', 'CGulf']
fs_observed = moments.Spectrum(fs_observed, pop_ids=new_pop_ids)


output_dir = "/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit"
create_model_plots(overall_best, fs_observed, output_dir)
    
for model_name in MODELS:
    print(f"\nPlotting {model_name}...")
    create_model_plots(best_runs_per_model[model_name], fs_observed, output_dir)




def create_model_summary(best_runs_dict, overall_best, output_file="/home/rbrennan/Tursiops-RAD-popgen/figures/model_fit/model_comparison_summary2.csv"):
    """Create a summary CSV file with all model information"""
    
    summary_data = []
    
    for model_name, model_info in best_runs_dict.items():
        # Extract run directory info
        run_dir = os.path.basename(model_info['directory'])
        run_number = run_dir.split('_')[-1]
        
        # Construct YAML path
        yaml_path = f"{INPUT_DIR}/{model_name}_{run_number}/round4/output_yaml/{model_name}_rep{int(model_info['run']):02d}_fmin-Round4.yaml"
        
        # Mark if this is the overall best
        is_best = (model_info == overall_best)
        
        summary_data.append({
            'Model': model_name,
            'Run': int(model_info['run']),
            'LogLikelihood': model_info['log_likelihood'],
            'Directory': model_info['directory'],
            'YAML_Path': yaml_path,
            'Run_Directory': run_number,
            'Overall_Best': is_best,
            'Rank': 0  # Will fill this in next
        })
    
    # Convert to DataFrame and sort by log-likelihood
    df = pd.DataFrame(summary_data)
    df = df.sort_values('LogLikelihood', ascending=False).reset_index(drop=True)
    
    # Add rank column
    df['Rank'] = range(1, len(df) + 1)
    
    # Reorder columns for better readability
    df = df[['Rank', 'Model', 'LogLikelihood', 'Run', 'Overall_Best', 'Directory', 'Run_Directory', 'YAML_Path']]
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"Model summary saved to: {output_file}")
    
    # Also display the summary
    print("\nMODEL COMPARISON SUMMARY:")
    print("="*100)
    print(df.to_string(index=False))
    print("="*100)
    
    return df

# Create the summary
create_model_summary(best_runs_per_model, overall_best)