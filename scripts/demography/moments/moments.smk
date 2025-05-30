# snakemake pipeline of moments.

#mamba activate snakemake_9.5.1
# --dry-run


# Snakemake pipeline for moments demographic modeling
# Converted from SLURM array job

# Configuration
MODELS = [f"{i:02d}" for i in range(1, 9)]  # Models 01-08
REPS = [f"{i:02d}" for i in range(1, 25)]   # Reps 01-24 (192 total jobs / 8 models = 24 reps per model)
METHOD = "fmin"  # Optimization method
ROUND4_REPS = [f"{i:02d}" for i in range(1, 11)]   # Reps 01-10


# Paths
ANALYSIS_DIR = "~/Tursiops-RAD-popgen/analysis/moments/"
SCRIPT_DIR = "~/Tursiops-RAD-popgen/scripts/demography/moments/"
LOG_DIR = "~/Tursiops-RAD-popgen/logout/"
rule all:
    input:
        # Round 3 outputs (Round 1 and 2 will run automatically as dependencies)
        expand("~/Tursiops-RAD-popgen/analysis/moments/round4/output_yaml/model_{model}_rep{rep}_{method}_downSample-Round4.yaml",
               model=MODELS, rep=ROUND4_REPS, method=METHOD),
        expand("~/Tursiops-RAD-popgen/analysis/moments/round4/model_{model}_{method}_downSample_summary.csv",
               model=MODELS, method=METHOD)

rule run_moments_model_round1:
    output:
        yaml_out = f"~/Tursiops-RAD-popgen/analysis/moments/round_1/output_yaml/model_{{model}}_rep{{rep}}_{METHOD}_downSample.yaml",
        plot = f"~/Tursiops-RAD-popgen/figures/moments/optimization_model_{{model}}_{{rep}}_{METHOD}_downSample.png",
        rep_yaml = "~/Tursiops-RAD-popgen/scripts/demography/moments/round_1/fourpop_{model}_rep{rep}.yaml"
        # with f strings, will replace the {} immediately. So the doubles are escaped, remain for wildcard to be replaced later
        # with regulat, they are not replaced immediately. 
    params:
        model = "{model}", # will replace in the script below, not here. 
        rep = "{rep}",
        analysis_dir = ANALYSIS_DIR, # this is a python variable, so use it directly.
        script_path = f"{SCRIPT_DIR}moments_model_round-01_m{{model}}.py" # this is a mix of python vairable and snakemake wildcard
    log:
        "~/Tursiops-RAD-popgen/logout/moments_round_1_m{model}_r{rep}.log"
    conda:
        "moments"
    resources:
        mem_mb = 6000,
        runtime = 4320,  # 3 days in minutes
        cpus_per_task = 1
    shell:
        """
        cd {params.analysis_dir}
        echo "model number: {params.model}" > {log}
        echo "rep number: {params.rep}" >> {log}
        echo "Starting analysis at $(date)" >> {log}
        
        python -u {params.script_path} {params.model} {params.rep} >> {log} 2>&1
        
        echo "Finished analysis at $(date)" >> {log}
        """


# Rule to ensure Round 1 CSV summary files are complete
rule complete_round1_summary:
    input:
        # This rule triggers after all Round 1 reps for a model are done
        yaml_files = expand("~/Tursiops-RAD-popgen/analysis/moments/initialRuns/output_yaml/model_{{model}}_rep{rep}_{method}_downSample.yaml", 
                           rep=REPS, method=METHOD)
    output:
        summary = f"~/Tursiops-RAD-popgen/analysis/moments/round_1/model_{{model}}_{METHOD}_downSample_summary.csv"
    params:
        expected_reps = len(REPS)
    shell:
        """
        # csv is made by the python scripts
        # This rule makes sure all reps are completed before considering the step done
        # check if the file exists and has the expected number of rows (header + reps)
        if [ -f {output.summary} ]; then
            EXPECTED_ROWS=$(($(({params.expected_reps})) + 1))  # number of reps + header
            ACTUAL_ROWS=$(wc -l < {output.summary})
            if [ $ACTUAL_ROWS -eq $EXPECTED_ROWS ]; then
                echo "Round 1 summary file {output.summary} is complete with $ACTUAL_ROWS rows"
                touch {output.summary} # to mark the file as edited.
            else
                echo "Warning: Round 1 summary file {output.summary} has only $ACTUAL_ROWS rows, expected $EXPECTED_ROWS"
                exit 1
            fi
        else
            echo "Error: Round 1 summary file {output.summary} does not exist"
            exit 1
        fi
        """

# Rule for running round 2
rule run_moments_model_round2:
    input:
        # Round 2 depends on the Round 1 summary file to find best parameters
        round1_summary = f"~/Tursiops-RAD-popgen/analysis/moments/initialRuns/model_{{model}}_{METHOD}_downSample_summary.csv"
    output:
        yaml_out = f"~/Tursiops-RAD-popgen/analysis/moments/round2/output_yaml/model_{{model}}_rep{{rep}}_{METHOD}_downSample-Round2.yaml",
        plot = f"~/Tursiops-RAD-popgen/figures/moments/round2/round2_optimization_model_{{model}}_{{rep}}_{METHOD}_downSample.png"
    params:
        model = "{model}",
        rep = "{rep}",
        analysis_dir = ANALYSIS_DIR,
        script_path = f"{SCRIPT_DIR}moments_model_round-02.py"
    log:
        "~/Tursiops-RAD-popgen/logout/moments_round2_m{model}_r{rep}.log"
    conda:
        "moments"
    resources:
        mem_mb = 8000,
        runtime = 4320,  # 3 days in minutes
        cpus_per_task = 1 
    shell:
        """
        cd {params.analysis_dir}
        echo "Round 2 - model number: {params.model}" > {log}
        echo "Round 2 - rep number: {params.rep}" >> {log}
        echo "Starting Round 2 analysis at $(date)" >> {log}
        
        python -u {params.script_path} {params.model} {params.rep} >> {log} 2>&1
        
        echo "Finished Round 2 analysis at $(date)" >> {log}
        """


# make sure csv done:
rule complete_round2_summary:
    input:
        yaml_files = expand("~/Tursiops-RAD-popgen/analysis/moments/round2/output_yaml/model_{{model}}_rep{rep}_{method}_downSample-Round2.yaml", 
                           rep=REPS, method=METHOD)# this makes sure all the yaml files exist
    output:
        summary = f"~/Tursiops-RAD-popgen/analysis/moments/round2/model_{{model}}_{METHOD}_downSample_summary.csv"
    params:
        expected_reps = len(REPS)
    shell:
        """
        if [ -f {output.summary} ]; then
            EXPECTED_ROWS=$(($(({params.expected_reps})) + 1))  # number of reps + header
            ACTUAL_ROWS=$(wc -l < {output.summary})
            if [ $ACTUAL_ROWS -eq $EXPECTED_ROWS ]; then
                echo "Round 2 summary file {output.summary} is complete with $ACTUAL_ROWS rows"
                touch {output.summary}
            else
                echo "Warning: Round 2 summary file {output.summary} has only $ACTUAL_ROWS rows, expected $EXPECTED_ROWS"
                exit 1
            fi
        else
            echo "Error: Round 2 summary file {output.summary} does not exist"
            exit 1
        fi
        """


# Rule for running round 3
rule run_moments_model_round3:
    input:
        # Round 2 depends on the Round 1 summary file to find best parameters
        round2_summary = f"~/Tursiops-RAD-popgen/analysis/moments/round2/model_{{model}}_{METHOD}_downSample_summary.csv"
    output:
        yaml_out = f"~/Tursiops-RAD-popgen/analysis/moments/round3/output_yaml/model_{{model}}_rep{{rep}}_{METHOD}_downSample-Round3.yaml",
        plot = f"~/Tursiops-RAD-popgen/figures/moments/round3/round3_optimization_model_{{model}}_{{rep}}_{METHOD}_downSample.png"
    params:
        model = "{model}",
        rep = "{rep}",
        analysis_dir = ANALYSIS_DIR,
        script_path = f"{SCRIPT_DIR}moments_model_round-03.py"
    log:
        "~/Tursiops-RAD-popgen/logout/moments_round3_m{model}_r{rep}.log"
    conda:
        "moments"
    resources:
        mem_mb = 8000,  # Increased memory
        runtime = 4320,
        cpus_per_task  = 1
    shell:
        """
        cd {params.analysis_dir}
        echo "Round 3 - model number: {params.model}" > {log}
        echo "Round 3 - rep number: {params.rep}" >> {log}
        echo "Starting Round 3 analysis at $(date)" >> {log}
        
        python -u {params.script_path} {params.model} {params.rep} >> {log} 2>&1
        
        echo "Finished Round 3 analysis at $(date)" >> {log}
        """


# make sure csv done:
rule complete_round3_summary:
    input:
        yaml_files = expand("~/Tursiops-RAD-popgen/analysis/moments/round3/output_yaml/model_{{model}}_rep{rep}_{method}_downSample-Round3.yaml", 
                           rep=REPS, method=METHOD)# this makes sure all the yaml files exist
    output:
        summary = f"~/Tursiops-RAD-popgen/analysis/moments/round3/model_{{model}}_{METHOD}_downSample_summary.csv"
    params:
        expected_reps = len(REPS)
    shell:
        """
        if [ -f {output.summary} ]; then
            EXPECTED_ROWS=$(($(({params.expected_reps})) + 1))  # number of reps + header
            ACTUAL_ROWS=$(wc -l < {output.summary})
            if [ $ACTUAL_ROWS -eq $EXPECTED_ROWS ]; then
                echo "Round 3 summary file {output.summary} is complete with $ACTUAL_ROWS rows"
                touch {output.summary}
            else
                echo "Warning: Round 3 summary file {output.summary} has only $ACTUAL_ROWS rows, expected $EXPECTED_ROWS"
                exit 1
            fi
        else
            echo "Error: Round 3 summary file {output.summary} does not exist"
            exit 1
        fi
        """


# Rule for running round 4
rule run_moments_model_round4:
    input:
        # Round 4 depends on the Round 3 summary file to find best parameters
        round3_summary = f"~/Tursiops-RAD-popgen/analysis/moments/round3/model_{{model}}_{METHOD}_downSample_summary.csv"
    output:
        yaml_out = f"~/Tursiops-RAD-popgen/analysis/moments/round4/output_yaml/model_{{model}}_rep{{rep}}_{METHOD}_downSample-Round4.yaml",
        plot = f"~/Tursiops-RAD-popgen/figures/moments/round4/round4_optimization_model_{{model}}_{{rep}}_{METHOD}_downSample.png"
    params:
        model = "{model}",
        rep = "{rep}",
        analysis_dir = ANALYSIS_DIR,
        script_path = f"{SCRIPT_DIR}moments_model_round-04-allSFS.py"
    log:
        "~/Tursiops-RAD-popgen/logout/moments_round4_m{model}_r{rep}.log"
    conda:
        "moments"
    resources:
        mem_mb = 10000,  # Increased memory
        runtime = 7200,
        cpus_per_task  = 1,
        round_4_runs = 1
    shell:
        """
        cd {params.analysis_dir}
        echo "Round 4 - model number: {params.model}" > {log}
        echo "Round 4 - rep number: {params.rep}" >> {log}
        echo "Starting Round 4 analysis at $(date)" >> {log}
        
        python -u {params.script_path} {params.model} {params.rep} >> {log} 2>&1
        
        echo "Finished Round 4 analysis at $(date)" >> {log}
        """


# make sure csv done:
rule complete_round4_summary:
    input:
        yaml_files = expand("~/Tursiops-RAD-popgen/analysis/moments/round4/output_yaml/model_{{model}}_rep{rep}_{method}_downSample-Round4.yaml", 
                           rep=ROUND4_REPS, method=METHOD)# this makes sure all the yaml files exist
    output:
        summary = f"~/Tursiops-RAD-popgen/analysis/moments/round4/model_{{model}}_{METHOD}_downSample_summary.csv"
    params:
        expected_reps = len(ROUND4_REPS)
    shell:
        """
        if [ -f {output.summary} ]; then
            EXPECTED_ROWS=$(($(({params.expected_reps})) + 1))  # number of reps + header
            ACTUAL_ROWS=$(wc -l < {output.summary})
            if [ $ACTUAL_ROWS -eq $EXPECTED_ROWS ]; then
                echo "Round 4 summary file {output.summary} is complete with $ACTUAL_ROWS rows"
                touch {output.summary}
            else
                echo "Warning: Round 4 summary file {output.summary} has only $ACTUAL_ROWS rows, expected $EXPECTED_ROWS"
                exit 1
            fi
        else
            echo "Error: Round 4 summary file {output.summary} does not exist"
            exit 1
        fi
        """