# treemix analysis

First, prepare input files for treemix. 

remember that the clust files are generated in `subset_pops.R`, at the end. In this, I subset the NC population bc it is hugely overrepresented in the data. 


```bash

mamba activate vcflib-1.0.9
module load bio/stacks/2.65 
module load bio/vcftools
module load lib64/gsl


cd ~/Tursiops-RAD-popgen/analysis/pop_structure/

dos2unix fourpop.clust
dos2unix sixpop.clust

# append aduncus samples onto these pop files:

#echo -e "SRR5357655\tAduncus" >> fourpop.clust
#echo -e "SRR5357656\tAduncus" >> fourpop.clust
#echo -e "SRR5357657\tAduncus" >> fourpop.clust

#echo -e "SRR5357655\tAduncus" >> sixpop.clust
#echo -e "SRR5357656\tAduncus" >> sixpop.clust
#echo -e "SRR5357657\tAduncus" >> sixpop.clust

# replace names with shorter versions. it messes up treemix, the long names. 

sed -i 's/Coastal_Atlantic/CoastAtl/g' fourpop.clust
sed -i 's/Coastal_Gulf/CoastGulf/g' fourpop.clust

sed -i 's/Coastal_Atlantic/CoastAtl/g' sixpop.clust
sed -i 's/Coastal_Gulf/CoastGulf/g' sixpop.clust
sed -i 's/Intermediate_Atlantic/IntAtl/g' sixpop.clust
sed -i 's/Intermediate_Gulf/IntGulf/g' sixpop.clust
sed -i 's/Offshore_Atlantic/OffAtl/g' sixpop.clust
sed -i 's/Offshore_Gulf/OffGulf/g' sixpop.clust

cut -f 1 fourpop.clust > fourpop.keep 
cut -f 1 sixpop.clust  > sixpop.keep 

vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz --keep fourpop.keep --recode --recode-INFO-all --stdout |  bgzip > fourpop.vcf.gz
zcat  fourpop.vcf.gz | awk '{if ($1 == "#CHROM"){print NF-9; exit}}'
# there should be 250 indivs.


vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz --keep sixpop.keep --recode --recode-INFO-all --stdout |  bgzip > sixpop.vcf.gz
zcat  sixpop.vcf.gz | awk '{if ($1 == "#CHROM"){print NF-9; exit}}'
# 226

# convert to treemix
populations --in-vcf fourpop.vcf.gz --treemix -O treemix/ -M fourpop.clust
# fourpop.p.treemix

populations --in-vcf sixpop.vcf.gz --treemix -O treemix/ -M sixpop.clust
#sixpop.p.treemix

cd treemix
tail -n +2 fourpop.p.treemix |  gzip  > fourpop.p.treemix.gz
tail -n +2 sixpop.p.treemix  | gzip > sixpop.p.treemix.gz

# run treemix

treemix -i fourpop.p.treemix.gz -o fourpop_nomigration -root Aduncus
treemix -i sixpop.p.treemix.gz -o sixpop_nomigration -root Aduncus

# first get consensus tree
# https://github.com/carolindahms/TreeMix/blob/main/Step1_TreeMix.sh

```

#### step 1: determine number of migration events. 

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

cp ../fourpop.p.treemix.gz ./

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
cp ../sixpop.p.treemix.gz ./

sh ~/Tursiops-RAD-popgen/scripts/treemix_step1.sh $infile $ncore $blockk $outgroup $nboot $pathP $outname $minM $maxM $migrep

```

`treemix_step2.R` to determine how many migration events. 


#### step 3 run bootstraps

```bash

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/fourpop

infile=fourpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=100		
numk=2
outname=fourpop
pathP=~/bin/phylip-3.697/exe/consense	

sh ~/Tursiops-RAD-popgen/scripts/treemix_bootstraps.sh $infile $numk $ncore $blockk $outgroup $nboot $pathP $outname

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/sixpop
infile=sixpop.p.treemix.gz
ncore=4
blockk=75
outgroup=Aduncus
nboot=100		
numk=4
outname=sixpop
pathP=~/bin/phylip-3.697/exe/consense	

sh ~/Tursiops-RAD-popgen/scripts/treemix_bootstraps.sh $infile $numk $ncore $blockk $outgroup $nboot $pathP $outname


# cd /mnt/c/Users/Reid.Brennan/Documents/projects/Tursiops_RAD_popgen/analysis/pop_structure/

```



then analyze the bootstrap reps: `treemix_step4.R`



