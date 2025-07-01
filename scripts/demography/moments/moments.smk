# moments analysis
# Each round waits for the previous round to complete

# Configuration
MODEL_NUMBER = "07"
METHOD = "fmin"
RUN_NAME = "run2"  

SFS_FILE_DOWNSAMPLE = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_downsample/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs"
SFS_FILE_EQUALSS = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs"
SCRIPT_DIR = "/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments"
ANALYSIS_DIR = "/home/rbrennan/Tursiops-RAD-popgen/analysis/moments"

# Define rep counts for each round
ROUND1_REPS = [f"{i:02d}" for i in range(1, 81)]   # 80 reps
ROUND2_REPS = [f"{i:02d}" for i in range(1, 41)]   # 40 reps  
ROUND3_REPS = [f"{i:02d}" for i in range(1, 11)]   # 10 reps
ROUND4_REPS = [f"{i:02d}" for i in range(1, 21)]   # 20 reps

# Target rule - run all rounds
rule all:
    input:
        f"{ANALYSIS_DIR}/mod7_expansion_{RUN_NAME}"

# Round 1: 80 reps with random starting values
rule round1:
    output:
        f"{ANALYSIS_DIR}/mod7_expansion/output_yaml/model_{MODEL_NUMBER}_rep{{rep}}_{METHOD}_downSample.yaml"
    log:
        f"{ANALYSIS_DIR}/mod7_expansion/logs/round1_model_{MODEL_NUMBER}_rep{{rep}}.log"
    name: "mts_r1"
    params:
        rep = "{rep}",
        model = MODEL_NUMBER,
        script = f"{SCRIPT_DIR}/moments_mod7_expansion_round-01.py",
    conda:
        "moments"
    resources:
        mem_mb = 6000,
        runtime = "24h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = "--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round1_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round1_%j.err"
    shell:
        """
        cd {ANALYSIS_DIR}
        python -u {params.script} {params.model} {params.rep}
        """

# Round 2: 40 reps (waits for all round 1 to finish)
rule round2:
    input:
        expand(f"{ANALYSIS_DIR}/mod7_expansion/output_yaml/model_{MODEL_NUMBER}_rep{{r1_rep}}_{METHOD}_downSample.yaml",
               r1_rep=ROUND1_REPS)
    output:
        f"{ANALYSIS_DIR}/mod7_expansion/round2/output_yaml/model_{MODEL_NUMBER}_rep{{rep}}_{METHOD}_downSample-Round2.yaml"
    name: "mts_r2"
    log:
        f"{ANALYSIS_DIR}/mod7_expansion/logs/round2_model_{MODEL_NUMBER}_rep{{rep}}.log"
    params:
        rep = "{rep}",
        model = MODEL_NUMBER,
        script = f"{SCRIPT_DIR}/moments_mod7_expansion_round-02.py",
    conda:
        "moments"
    resources:
        mem_mb = 8000,
        runtime = "24h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = "--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round2_rep%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round2_rep%j.err"
    shell:
        """
        cd {ANALYSIS_DIR}
        python -u {params.script} {params.model} {params.rep}
        """

# Round 3: 10 reps (waits for all round 2 to finish)
rule round3:
    input:
        expand(f"{ANALYSIS_DIR}/mod7_expansion/round2/output_yaml/model_{MODEL_NUMBER}_rep{{r2_rep}}_{METHOD}_downSample-Round2.yaml",
               r2_rep=ROUND2_REPS)
    output:
        f"{ANALYSIS_DIR}/mod7_expansion/round3/output_yaml/model_{MODEL_NUMBER}_rep{{rep}}_{METHOD}_downSample-Round3.yaml"
    name: "mts_r3"
    log:
        f"{ANALYSIS_DIR}/mod7_expansion/logs/round3_model_{MODEL_NUMBER}_rep{{rep}}.log"
    conda:
        "moments"  # environment name
    params:
        rep = "{rep}",
        model = MODEL_NUMBER,
        script = f"{SCRIPT_DIR}/moments_mod7_expansion_round-03_plus.py",
    resources:
        mem_mb = 8000,
        runtime = "24h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = "--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round3_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round3_%j.err"
    shell:
        """
        python -u {params.script} {params.model} {params.rep} 3 2 5 {METHOD} {SFS_FILE_EQUALSS}
        """

# Round 4: 20 reps (waits for all round 3 to finish)
rule round4:
    input:
        expand(f"{ANALYSIS_DIR}/mod7_expansion/round3/output_yaml/model_{MODEL_NUMBER}_rep{{r3_rep}}_{METHOD}_downSample-Round3.yaml",
               r3_rep=ROUND3_REPS)
    output:
        f"{ANALYSIS_DIR}/mod7_expansion/round4/output_yaml/model_{MODEL_NUMBER}_rep{{rep}}_{METHOD}_downSample-Round4.yaml"
    name: "mts_r4"
    log:
        f"{ANALYSIS_DIR}/mod7_expansion/logs/round4_model_{MODEL_NUMBER}_rep{{rep}}.log"
    conda:
        "moments"  
    params:
        rep = "{rep}",
        model = MODEL_NUMBER,
        script = f"{SCRIPT_DIR}/moments_mod7_expansion_round-03_plus.py",
    resources:
        mem_mb = 8000,
        runtime = "24h",
        cpus_per_task = 1,
        slurm_partition = "standard",
        slurm_extra = "--output=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round4_%j.out --error=/home/rbrennan/Tursiops-RAD-popgen/logout/moments_round4_%j.err"
    shell:
        """
        python -u {params.script} {params.model} {params.rep} 4 1 5 {METHOD} {SFS_FILE_EQUALSS}
        """
# Final step: rename the directory
rule rename_final:
    input:
        expand(f"{ANALYSIS_DIR}/mod7_expansion/round4/output_yaml/model_{MODEL_NUMBER}_rep{{rep}}_{METHOD}_downSample-Round4.yaml",
               rep=ROUND4_REPS)
    output:
        directory(f"{ANALYSIS_DIR}/mod7_expansion_{RUN_NAME}")
    shell:
        """
        if [ -d {ANALYSIS_DIR}/mod7_expansion ]; then
            mv {ANALYSIS_DIR}/mod7_expansion {ANALYSIS_DIR}/mod7_expansion_{RUN_NAME}
        fi
        """