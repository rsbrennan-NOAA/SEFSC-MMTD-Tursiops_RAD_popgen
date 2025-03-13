
```bash

mamba activate vcflib-1.0.9
module load bio/stacks/2.65 
module load bio/vcftools
module load lib64/gsl


cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix

# append aduncus samples onto these pop files:

#echo -e "SRR5357655\tAduncus" >> fourpop.clust
#echo -e "SRR5357656\tAduncus" >> fourpop.clust
#echo -e "SRR5357657\tAduncus" >> fourpop.clust

#echo -e "SRR5357655\tAduncus" >> sixpop.clust
#echo -e "SRR5357656\tAduncus" >> sixpop.clust
#echo -e "SRR5357657\tAduncus" >> sixpop.clust


cut -f 1 fourpop.clust > fourpop.keep 
cut -f 1 sixpop.clust  > sixpop.keep 

vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz --keep fourpop.keep --recode --recode-INFO-all --stdout |  bgzip > fourpop.vcf.gz
zcat  fourpop.vcf.gz | awk '{if ($1 == "#CHROM"){print NF-9; exit}}'
# there should be 248 indivs.


vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz --keep sixpop.keep --recode --recode-INFO-all --stdout |  bgzip > sixpop.vcf.gz
zcat  sixpop.vcf.gz | awk '{if ($1 == "#CHROM"){print NF-9; exit}}'
# After filtering, kept 213 out of 337 Individuals

# convert to treemix
populations --in-vcf fourpop.vcf.gz --treemix -O . -M fourpop.clust
# fourpop.p.treemix

populations --in-vcf sixpop.vcf.gz --treemix -O . -M sixpop.clust
#sixpop.p.treemix

tail -n +2 fourpop.p.treemix |  gzip  > fourpop.p.treemix.gz
tail -n +2 sixpop.p.treemix  | gzip > sixpop.p.treemix.gz

# run treemix

treemix -i fourpop.p.treemix.gz -o fourpop_nomigration -root Aduncus
treemix -i sixpop.p.treemix.gz -o sixpop_nomigration -root Aduncus

# first get consensus tree
# https://github.com/carolindahms/TreeMix/blob/main/Step1_TreeMix.sh

```


```bash
mamba activate vcflib-1.0.9
module load bio/stacks/2.65 
module load bio/vcftools
module load lib64/gsl

infile=fourpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=100		
pathP=~/bin/phylip-3.697/exe/consense	
outname=fourpop		
minM=1 
maxM=8
migrep=10

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/fourpop

sh ~/Tursiops-RAD-popgen/scripts/treemix_step1.sh $infile $ncore $blockk $outgroup $nboot $pathP $outname $minM $maxM $migrep


infile=sixpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=500		
pathP=~/bin/phylip-3.697/exe/consense	
outname=sixpop
minM=1 
maxM=10
migrep=10

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/sixpop

sh ~/Tursiops-RAD-popgen/scripts/treemix_step1.sh $infile $ncore $blockk $outgroup $nboot $pathP $outname $minM $maxM $migrep

```

run R step 1


```bash

infile=fourpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=500		
mig=2
outname=fourpop
runs=30
tree=fourpop_constree.newick
pathP=~/bin/phylip-3.697/exe/consense	

sh ~/Tursiops-RAD-popgen/scripts/treemix_step3.sh $infile $ncore $blockk $outgroup $nboot $mig $outname $runs $tree $pathP

```


test script from bite R. 

```bash

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/fourpop_new

infile=fourpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=100		
numk=2
outname=fourpop
pathP=~/bin/phylip-3.697/exe/consense	

sh ~/Tursiops-RAD-popgen/scripts/treemix_bootstraps.sh $infile $numk $ncore $blockk $outgroup $nboot $pathP $outname

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/sixpop_new
infile=sixpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=100		
numk=4
outname=sixpop
pathP=~/bin/phylip-3.697/exe/consense	

sh ~/Tursiops-RAD-popgen/scripts/treemix_bootstraps.sh $infile $numk $ncore $blockk $outgroup $nboot $pathP $outname


cd /mnt/c/Users/Reid.Brennan/Documents/projects/Tursiops_RAD_popgen/analysis/pop_structure/

```








































```bash

#--------------------------------------------------------------------------------------------
# four pop run

outname=fourpop

for i in {1..100}
do
  sleep 0.5

  seed=$((RANDOM + 1$(date +%N) % 1000000))

  treemix -i fourpop.p.treemix.gz -o bootstrap_fourpop/fourpop_constree_bootrep_${i} -root Aduncus -bootstrap -k 75 -seed $seed -se 
done

### Create a file with all the bootstrapped trees

rm fourpop_bootconstree.tre

for i in `seq 1 100`
do
 bootf="bootstrap_fourpop/fourpop_constree_bootrep_${i}.treeout.gz"
 gunzip -c $bootf | head -1 >> fourpop_bootconstree.tre
done


#### Run PHYLIP on the bootstrapped trees to obtain a consensus tree

# Create parameters file
outgroup=Aduncus
outname=fourpop

# Find the position of Outgroup population
posOutgroup=`head -1 $outname"_bootconstree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
# echo $posOutgroup
echo $outname"_bootconstree.tre" > $outname"_PhylipInputFile"
echo "O" >> $outname"_PhylipInputFile"
echo $posOutgroup >> $outname"_PhylipInputFile"
echo "Y" >> $outname"_PhylipInputFile"

rm outfile
rm outtree

# Run Phylip
~/bin/phylip-3.697/exe/consense < fourpop_PhylipInputFile

### The output from Phylip will be modified because:
### 1) TreeMix accepts trees with only one line
### 2) TreeMix accepts newick format file 

cat outtree | tr -d "\n" > $outname"_constree.newick"
echo >> $outname"_constree.newick"


#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
# 6 pop run

outname=sixpop

for i in {1..100}
do

  sleep 0.5
  seed=$((RANDOM + 1$(date +%N) % 1000000))

  treemix -i sixpop.p.treemix.gz -o bootstrap_sixpop/sixpop_constree_bootrep_${i} -root Aduncus -bootstrap -k 75 -seed $seed -se
done

### Create a file with all the bootstrapped trees
rm sixpop_bootconstree.tre
for i in `seq 1 100`
do
 bootf="bootstrap_sixpop/sixpop_constree_bootrep_${i}.treeout.gz"
 gunzip -c $bootf | head -1 >> sixpop_bootconstree.tre
done

# Create parameters file

outgroup=Aduncus
outname=sixpop

# Find the position of Outgroup population
posOutgroup=`head -1 $outname"_bootconstree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
# echo $posOutgroup
echo $outname"_bootconstree.tre" > $outname"_PhylipInputFile"
echo "O" >> $outname"_PhylipInputFile"
echo $posOutgroup >> $outname"_PhylipInputFile"
echo "Y" >> $outname"_PhylipInputFile"

rm outfile
rm outtree

# Run Phylip
~/bin/phylip-3.697/exe/consense < sixpop_PhylipInputFile

cat outtree | tr -d "\n" > $outname"_constree.newick"
echo >> $outname"_constree.newick"

##################################################
###### Run TreeMix with the consensus tree #######
##################################################

#mkdir -p test_migrations_fourpop
#mkdir -p test_migrations_sixpop


for pop in sixpop
  do
    for m in {1..10}
     do
        for i in {1..10}
        do
        sleep 0.2
          seed=$((RANDOM + 1$(date +%N) % 1000000))

        treemix \
             -i ${pop}.p.treemix.gz \
             -o ./"test_migrations_"${pop}/treemix.${i}.${m} \
             -global \
             -m ${m} \
             -k 75 \
             -seed $seed \
  	         -root Aduncus \
             -tf ${pop}"_constree.newick" 
        done 
  done
done

for pop in fourpop
  do
    for m in {1..5}
     do
        for i in {1..10}
        do
        sleep 0.2
          seed=$((RANDOM + 1$(date +%N) % 1000000))

        treemix \
             -i ${pop}.p.treemix.gz \
             -o ./"test_migrations_"${pop}/treemix.${i}.${m} \
             -global \
             -m ${m} \
             -k 75 \
             -seed $seed \
  	         -root Aduncus \
             -tf ${pop}"_constree.newick" 
        done 
  done
done


```

##### step 2 in R



##### step 3 treemix:

```bash
mamba activate vcflib-1.0.9
module load bio/stacks/2.65 
module load bio/vcftools
module load lib64/gsl


cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix
#####################################################
######### Step 3: Final runs with optimum m #########
#####################################################
# This script builds a consensus tree with the optimum number of migration events from a specified number of bootstraps and uses it to run multiple final independent runs. 
# You will need to have installed TreeMix, Parallel and PHYLIP Consense.

infile=fourpop.p.treemix.gz           # TreeMix input file
ncore=1 	    # maximum number of cores to use
blockk=75 	    # SNP block size
outgroup=Aduncus         # set outgroup, for an unrooted ML tree put here 'noRoot' (without quotes)
nboot=20	    # number of bootstrap replicates of tree with migration
mig=2              # number of migration events
outname=fourpop_final	    # name for output file
runs=10             # number of independent runs (N)
tree=fourpop_constree.newick             # name of consensus tree build without migration events (in newick format)
pathP=~/bin/phylip-3.697/exe/consense	    # path to Phylip consense program
outdir=final_runs_fourpop

############################
### Bootstrap procedure ####
############################

#mkdir final_runs_fourpop
#mkdir final_runs_fourpop/bootstrap

for i in {1..100}
do
  sleep 0.4
  seed=$((RANDOM + 1$(date +%N) % 1000000))

  treemix -i $infile -bootstrap -k $blockk -m $mig -se -root Aduncus -seed $seed -tf $tree -o $outdir"/bootstrap/"$outname"_constree_bootrep_"$i 
done

rm $outname"_boottree.tre"

### Create a file with all the bootstrapped trees
for a in `seq 1 100`
do
 bootf=$outdir"/bootstrap/"$outname"_constree_bootrep_"$a".treeout.gz"
 gunzip -c $bootf | head -1 >> $outname"_boottree.tre"
done

#########################################################################
#### Run PHYLIP on the bootstrapped trees to obtain a consensus tree ####
#########################################################################
echo "***** Phylip - consensus tree construction: START *****"
### Clean the environment
rm -rf outfile outtree screanout

# Create parameters file
# Find the position of Outgroup population
posOutgroup=`head -1 $outname"_boottree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
# echo $posOutgroup
echo $outname"_boottree.tre" > $outname"_PhylipInputFile"
echo "O" >> $outname"_PhylipInputFile"
echo $posOutgroup >> $outname"_PhylipInputFile"
echo "Y" >> $outname"_PhylipInputFile"


# Run Phylip
$pathP < $outname"_PhylipInputFile" > screanout

cat outtree | tr -d "\n" > $outname"_finalconstree.newick"
echo >> $outname"_finalconstree.newick"

##################################################
#### Run TreeMix with the new consensus tree #####
##################################################

for i in {1..100}
  do

  sleep 0.2  
  seed=$((RANDOM + 1$(date +%N) % 1000000))

  treemix \
  -i ${infile} \
  -o $outdir"/final/"${outname}_${mig}"mig_finalrun_"${i} \
  -global \
  -m ${mig} \
  -k 75 \
  -seed $seed \
  -root Aduncus \
  -se \
  -tf $outname"_finalconstree.newick"

done








#--------------------------------------------------------------------------------------------------------------------
# six pop

infile=sixpop.p.treemix.gz           # TreeMix input file
ncore=1 	    # maximum number of cores to use
blockk=75 	    # SNP block size
outgroup=Aduncus         # set outgroup, for an unrooted ML tree put here 'noRoot' (without quotes)
nboot=20	    # number of bootstrap replicates of tree with migration
mig=4              # number of migration events
outname=sixpop_final	    # name for output file
runs=10             # number of independent runs (N)
tree=sixpop_constree.newick             # name of consensus tree build without migration events (in newick format)
pathP=~/bin/phylip-3.697/exe/consense	    # path to Phylip consense program
outdir=final_runs_sixpop

############################
### Bootstrap procedure ####
############################

#mkdir final_runs_fourpop
#mkdir final_runs_fourpop/bootstrap

for i in {1..100}
do
  sleep 0.4
  seed=$((RANDOM + 1$(date +%N) % 1000000))

  treemix -i $infile -bootstrap -k $blockk -m $mig -se -root Aduncus -seed $seed -tf $tree -o $outdir"/bootstrap/"$outname"_constree_bootrep_"$i 
done

rm $outname"_boottree.tre"

### Create a file with all the bootstrapped trees
for a in `seq 1 100`
do
 bootf=$outdir"/bootstrap/"$outname"_constree_bootrep_"$a".treeout.gz"
 gunzip -c $bootf | head -1 >> $outname"_boottree.tre"
done

#########################################################################
#### Run PHYLIP on the bootstrapped trees to obtain a consensus tree ####
#########################################################################
rm -rf outfile outtree screanout

# Create parameters file
# Find the position of Outgroup population
posOutgroup=`head -1 $outname"_boottree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
# echo $posOutgroup
echo $outname"_boottree.tre" > $outname"_PhylipInputFile"
echo "O" >> $outname"_PhylipInputFile"
echo $posOutgroup >> $outname"_PhylipInputFile"
echo "Y" >> $outname"_PhylipInputFile"


# Run Phylip
$pathP < $outname"_PhylipInputFile" > screanout

cat outtree | tr -d "\n" > $outname"_finalconstree.newick"
echo >> $outname"_finalconstree.newick"

##################################################
#### Run TreeMix with the new consensus tree #####
##################################################

for i in {1..100}
  do

  sleep 0.4
  seed=$((RANDOM + 1$(date +%N) % 1000000))

  treemix \
  -i ${infile} \
  -o $outdir"/final/"${outname}${mig}"mig_finalrun_"${i} \
  -global \
  -m ${mig} \
  -k 75 \
  -seed $seed \
  -root Aduncus \
  -se \
  -tf $outname"_finalconstree.newick"

done


# just ran the above. now go to R and check results.






```










