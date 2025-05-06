#!/bin/bash
#SBATCH --job-name=gen_boots
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mem=28G
#SBATCH --partition=standard
#SBATCH --time=72:00:00

module load bio/bedtools/2.31.1
module load bio/vcftools/0.1.16

VCF_FILE=~/Tursiops-RAD-popgen/analysis/variants/filtered.final.noMAF.vcf.gz
NUM_BOOTSTRAPS=100
OUTPUT_DIR=~/Tursiops-RAD-popgen/analysis/moments/bootstrap_sfs
GENOME_FILE=~/reference_genomes/bottlenose_dolphin/GCF_011762595.1_mTurTru1.mat.Y_genomic_genomeFile.txt  # File with chromosome sizes

# Create output directories

cd ~/Tursiops-RAD-popgen/analysis/moments
mkdir -p chunks

# Step 1: get chr ids from the vcf. this is somewhat slow, so do it ahead of time (depending).
#zcat $VCF_FILE | cut -f1| grep -v "^#" | sort | uniq > ~/Tursiops-RAD-popgen/analysis/moments/chromosomes.txt
#CHROMOSOMES=($(cat ~/Tursiops-RAD-popgen/analysis/moments/chromosomes.txt))
CHROMOSOMES=($(zcat $VCF_FILE | grep -v "^#" | cut -f1 | sort | uniq)) # or do it here. uncomment this

echo "Found ${#CHROMOSOMES[@]} chromosomes"

# Step 3: 
echo "Splitting each chromosome into 5 parts"
for CHR in "${CHROMOSOMES[@]}"; do
  # Get chromosome length
  CHR_LENGTH=$(grep -w "$CHR" $GENOME_FILE | cut -f2)
  
  # Calculate chunk sizes (divide by 3)
  PART_SIZE=$((CHR_LENGTH / 5))
  
  # Create 5 chunks
  for PART in {0..4}; do
    START=$((PART * PART_SIZE))
    # For the last part, ensure we go to the end of the chromosome
    if [ $PART -eq 2 ]; then
      END=$CHR_LENGTH
    else
      END=$(((PART + 1) * PART_SIZE))
    fi
    
    # make filename
    MOD_CHR=$(echo $CHR | tr '.' '_')
    CHUNK_ID="${MOD_CHR}_part$((PART + 1))"
    echo "Creating $CHUNK_ID: $CHR:$START-$END"
    
    # Create BED file for this chunk
    echo -e "${CHR}\t${START}\t${END}" > chunks/${CHUNK_ID}.bed
    
    # Extract VCF data for this chunk
     bedtools intersect -a $VCF_FILE -b chunks/${CHUNK_ID}.bed -header > chunks/${CHUNK_ID}.vcf
  done
done

# Count total number of chunks (should be number of chromosomes x 5)
NUM_CHUNKS=$(ls chunks/*_part*.vcf | wc -l)
echo "Created $NUM_CHUNKS chunks (there should be 5 parts per each of the 22 chromosomes)"

# get header to use later
grep "^#" $VCF_FILE > header.txt

# Step 4: Bootstrap loop
for i in $(seq 1 $NUM_BOOTSTRAPS); do
  echo "Creating bootstrap replicate $i"
  
  # Create directory for this bootstrap
  mkdir -p $OUTPUT_DIR/bootstrap_$i
  
  # start with header, then will append variants after
  cat header.txt > $OUTPUT_DIR/bootstrap_$i/bootstrap_$i.vcf
  
  # sampling chunks with replacement
  # sample NUM_CHUNKS chunks to maintain the same genome size roughly
  for j in $(seq 1 $NUM_CHUNKS); do
    # select a chunk
    RANDOM_CHUNK=$(ls chunks/*_part*.vcf | shuf -n 1)
    echo "Selected chunk $RANDOM_CHUNK for position $j"
    
    # append to bootstrap file
    grep -v "^#" $RANDOM_CHUNK >> $OUTPUT_DIR/bootstrap_$i/bootstrap_$i.vcf
  done
  
  # Run easySFS on bootstrapped VCF
  #echo "Running easySFS on bootstrap $i..."
  #easySFS.py -i $OUTPUT_DIR/bootstrap_$i/bootstrap_$i.vcf -p  ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop_all.clust \
  #  --order Coastal_Atlantic,Coastal_Gulf,Intermediate,Offshore --proj 224,64,78,88 --outdir $OUTPUT_DIR/bootstrap_$i
  
  echo "Completed bootstrap replicate $i"
done
