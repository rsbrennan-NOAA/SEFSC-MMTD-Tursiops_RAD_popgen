
#!/usr/bin/env python3
"""
Moments demographic inference for Tursiops populations
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


# record start time
start_time = time.time()
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"[{current_time}] Starting moments analysis")


# load the sfs 
print("Loading SFS")
fs = moments.Spectrum.from_file("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs")

    
# pop ids are in wrong order- fix them
new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
print(f"Population IDs: {fs.pop_ids}")
 
# define file paths
deme_graph = "/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/fourpop_02_Admix_AncCoast.yaml"
options = "/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/options_fourpop_02_Simple.yaml"
output = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/model_02-bestfit_fmin.yaml"

print(f"deme_graph file: {deme_graph}")
print(f"options file: {options}")
print(f"output file: {output}")

# Run optimization
print("Starting demographic model optimization")
ret = moments.Demes.Inference.optimize(
    deme_graph,
    options,
    fs,
    method="fmin",
    output=output,
    overwrite=True,
    verbose=1,
    perturb = 0.7
)


# Process results
param_names, opt_params, LL = ret
print(f"Log-likelihood: {-LL}")
print("Best fit parameters:")

# save to file
with open("optimization_model_02_results.txt", "w") as f:
    f.write(f"Log-likelihood: {-LL}\n")
    f.write("Best fit parameters\n")
    for n, p in zip(param_names, opt_params):
        line = f"{n}\t{p:.3f}"
        f.write(line + "\n")
        print(line)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Analysis ended at {end_time}")
print(f"Analysis completed in {elapsed_time/3600:.2f} hours")

# plot
try:
     graph = demes.load(output)
     fig = plt.figure(figsize=(10, 8))
     demesdraw.tubes(graph, ax=fig.add_subplot(111))
     plt.savefig("demes_model_02.png", dpi=300)
     print("Saved demographic model plot")
except Exception as e:
     print(f"Failed to generate plot: {e}")

