# moments analysis - Coastal model
# Each round waits for the previous round to complete

# Configuration
METHOD = "fmin"
RUN_NAME = config.get("run_name", "run_default")
MODEL_NAME = config.get("model_name")  # Get from command line

SFS_FILE_COASTAL = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf.sfs"
SCRIPT_DIR = "/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models"
ANALYSIS_DIR = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments"

# Define rep counts for each round
ROUND1_REPS = [f"{i:02d}" for i in range(1, 61)]   # 60 reps
ROUND2_REPS = [f"{i:02d}" for i in range(1, 21)]   # 20 reps  
ROUND3_REPS = [f"{i:02d}" for i in range(1, 21)]   # 20 reps
ROUND4_REPS = [f"{i:02d}" for i in range(1, 21)]   # 20 reps

# Target rule - run all rounds
rule all:
    input:
        f"{ANALYSIS_DIR}/{MODEL_NAME}_{RUN_NAME}"

# Round 1: 40 reps with random starting values
rule round1:
    output:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/output_yaml/{MODEL_NAME}_rep{{rep}}_{METHOD}.yaml"
    log:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/logs/round1_{MODEL_NAME}_rep{{rep}}.log"
    name: f"r1-{MODEL_NAME}"
    params:
        rep = "{rep}",
        script = f"{SCRIPT_DIR}/moments_{MODEL_NAME}_round-01.py",
    conda:
        "moments"
    resources:
        mem_mb = 6000,
        runtime = "24h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = f"--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round1_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round1_%j.err"
    retries: 2
    shell:
        """
        cd {ANALYSIS_DIR}
        python -u {params.script} {params.rep} {MODEL_NAME}
        """

# Round 2: 20 reps (waits for all round 1 to finish)
rule round2:
    input:
        expand(f"{ANALYSIS_DIR}/{MODEL_NAME}/output_yaml/{MODEL_NAME}_rep{{r1_rep}}_{METHOD}.yaml",
               r1_rep=ROUND1_REPS)
    output:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/round2/output_yaml/{MODEL_NAME}_rep{{rep}}_{METHOD}-Round2.yaml"
    name: f"r2-{MODEL_NAME}"
    log:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/logs/round2_{MODEL_NAME}_rep{{rep}}.log"
    params:
        rep = "{rep}",
        script = f"{SCRIPT_DIR}/moments_coastal_generic_round-02.py",
    conda:
        "moments"
    resources:
        mem_mb = 8000,
        runtime = "24h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = f"--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round2_rep%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round2_rep%j.err"
    retries: 2
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        
        cd {ANALYSIS_DIR}
        python -u {params.script} {params.rep} {MODEL_NAME}
        """

# Round 3: 10 reps (waits for all round 2 to finish)
rule round3:
    input:
        expand(f"{ANALYSIS_DIR}/{MODEL_NAME}/round2/output_yaml/{MODEL_NAME}_rep{{r2_rep}}_{METHOD}-Round2.yaml",
               r2_rep=ROUND2_REPS)
    output:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/round3/output_yaml/{MODEL_NAME}_rep{{rep}}_{METHOD}-Round3.yaml"
    name: f"r3-{MODEL_NAME}"
    log:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/logs/round3_{MODEL_NAME}_rep{{rep}}.log"
    conda:
        "moments"  # environment name
    params:
        rep = "{rep}",
        script = f"{SCRIPT_DIR}/moments_coastal_generic_round-03_plus.py",
    resources:
        mem_mb = 8000,
        runtime = "30h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = f"--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round3_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round3_%j.err"
    retries: 2
    shell:
        """
        python -u {params.script} {params.rep} 3 2 5 {METHOD} {SFS_FILE_COASTAL} {MODEL_NAME}
        """

# Round 4: 20 reps (waits for all round 3 to finish)
rule round4:
    input:
        expand(f"{ANALYSIS_DIR}/{MODEL_NAME}/round3/output_yaml/{MODEL_NAME}_rep{{r3_rep}}_{METHOD}-Round3.yaml",
               r3_rep=ROUND3_REPS)
    output:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/round4/output_yaml/{MODEL_NAME}_rep{{rep}}_{METHOD}-Round4.yaml"
    name: f"r4-{MODEL_NAME}"
    log:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/logs/round4_{MODEL_NAME}_rep{{rep}}.log"
    conda:
        "moments"  
    params:
        rep = "{rep}",
        script = f"{SCRIPT_DIR}/moments_coastal_generic_round-03_plus.py",
    resources:
        mem_mb = 8000,
        runtime = "30h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = f"--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round4_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_round4_%j.err"
    retries: 2
    shell:
        """
        python -u {params.script} {params.rep} 4 1 20 {METHOD} {SFS_FILE_COASTAL} {MODEL_NAME}
        """

# Final step: rename the directory
rule rename_final:
    input:
        expand(f"{ANALYSIS_DIR}/{MODEL_NAME}/round4/output_yaml/{MODEL_NAME}_rep{{rep}}_{METHOD}-Round4.yaml",
               rep=ROUND4_REPS)
    output:
        directory(f"{ANALYSIS_DIR}/{MODEL_NAME}_{RUN_NAME}")
    log:
        f"{ANALYSIS_DIR}/{MODEL_NAME}/logs/rename_{MODEL_NAME}_{RUN_NAME}.log"
    resources:
        mem_mb = 6000,
        runtime = "1h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = f"--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_rename_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_{MODEL_NAME}_rename_%j.err"
    retries: 2
    shell:
        """
        sleep 60

        if [ -d {ANALYSIS_DIR}/{MODEL_NAME} ]; then
            mv {ANALYSIS_DIR}/{MODEL_NAME} {ANALYSIS_DIR}/{MODEL_NAME}_{RUN_NAME}
        fi
        """