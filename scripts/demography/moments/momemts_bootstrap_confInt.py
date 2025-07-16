"""
Moments confidence intervals
Uses bootstrapped sfs from moments_bootstrap_sfs.nf
"""

import moments
import os
import demes
import demesdraw
import pandas as pd
import glob
import os
from pathlib import Path
import matplotlib.pylab as plt


def read_bootstrap_sfs(base_dir):
   """
   read bootstrap SFS into a list for moments.
   
   Args:
       base_dir: path to the bootstrap directory containing bootstrap_* subdirectories
   
   Returns:
       list of moments.Spectrum objects to use with uncerts
   """
   bootstraps = []
   
   # Get all bootstrap directories
   bootstrap_dirs = [d for d in os.listdir(base_dir) if d.startswith('bootstrap_')]
   #bootstrap_dirs.sort(key=lambda x: int(x.split('_')[1]))  # Sort numerically
   
   sfs_filename = "Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs"
   
   for bootstrap_dir in bootstrap_dirs:
       sfs_path = os.path.join(base_dir, bootstrap_dir, "dadi", sfs_filename)
       
       if os.path.exists(sfs_path):
           # Load the SFS file
           fs = moments.Spectrum.from_file(sfs_path)
           
           # Fix population IDs order
           new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
           fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
           
           bootstraps.append(fs)
       else:
           print(f"SFS not found in {bootstrap_dir}")
   
   return bootstraps

# read in files
base_dir = "~/Tursiops-RAD-popgen/analysis/moments/bootstrap_sfs/"
base_dir = os.path.expanduser(base_dir)  # Expand ~ to full path
bootstraps = load_bootstrap_sfs(base_dir)
print(f"Loaded {len(bootstraps)} bootstrap SFS files")


def find_best(base_path="~/Tursiops-RAD-popgen/analysis/moments"):
    """find the best run in mod7_expansion_run directories"""
    base_path = os.path.expanduser(base_path)
    run_dirs = glob.glob(os.path.join(base_path, "mod7_expansion_run*"))
    
    best_likelihood = float('-inf')
    best_info = None
    
    for run_dir in run_dirs:
        csv_path = os.path.join(run_dir, "round4", "model_07_fmin_downSample_summary.csv")
        
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            max_likelihood = df['LogLikelihood'].max()
            max_row = df.loc[df['LogLikelihood'].idxmax()]
            
            if max_likelihood > best_likelihood:
                best_likelihood = max_likelihood
                best_info = {
                    'run_dir': run_dir,
                    'run_number': int(max_row['Run']),  # convert to int
                    'likelihood': max_likelihood,
                    'full_row': max_row
                }
    
    return best_info

# Find the best model
best_model = find_best()
print(f"Best LogLikelihood: {best_model['likelihood']}")
print(f"Found in: {best_model['run_dir']}")
print(f"Run number: {best_model['run_number']}")

def get_best_yaml(best_model_info):
    """Get the path to the YAML file corresponding to the best model"""
    run_dir = best_model_info['run_dir']
    run_number = int(best_model_info['run_number'])  # Convert to int
    
    yaml_filename = f"model_07_rep{run_number:02d}_fmin_downSample-Round4.yaml"
    yaml_path = os.path.join(run_dir, "round4", "output_yaml", yaml_filename)
    
    if os.path.exists(yaml_path):
        return yaml_path
    else:
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")

output = get_best_yaml(best_model)

OPTIONS_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/options_fourpop_07_expansion.yaml"
fs = moments.Spectrum.from_file("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_downsample/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs")
# pop ids are in wrong order- fix them
new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
print(f"Population IDs: {fs.pop_ids}")

std_err = moments.Demes.Inference.uncerts(
    output,
    OPTIONS_PATH,
    fs,
    bootstraps=bootstraps,
    method="FIM",
    log=False,
    fit_ancestral_misid=False,
    output="/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/uncerts_mod7.txt",
    overwite=True
)

# need to get opt_params as well

print("Standard errors:")
print("param\t\topt\t\tstderr")
for n, p, e in zip(param_names, opt_params, std_err):
    print(f"{n}\t{p:-11g}\t{e:-14g}")


# draw figure


PLOT_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/bestmodel.pdf"

graph = demes.load(output)
fig = plt.figure(figsize=(10, 8))
demesdraw.tubes(graph, ax=fig.add_subplot(111))
plt.savefig(PLOT_FILE, dpi=300)
print(f"Saved plot to {PLOT_FILE}")
