
# first all individuals
library(dartR)

genl_all <- gl.read.vcf("analysis/variants/filtered.final_ids_LDthin_noX.vcf.gz")
genl_fourpop <- gl.read.vcf("analysis/variants/tursiops_subset_fourpop.vcf.gz")
genl_sixpop <- gl.read.vcf("analysis/variants/tursiops_subset_sixpop.vcf.gz")

genl_all
genl_fourpop
genl_sixpop

library(dartR.base)

genl <- genl_all
pops1 <- read.table("analysis/population_assignments_summary.txt", header=T)

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$sixpop)
head(pops)
pops$id == genl@ind.names

#reorder pops
pops <- pops[match(genl@ind.names, pops$id), ]

all(pops$id == genl@ind.names)

genl@pop <- as.factor(pops$pop)
table(genl@pop)

gl2bayesAss(genl, outpath="analysis/pop_structure/bayesass", 
            outfile="sixpop.BayesAss.txt")


pops <- data.frame(id = pops1$indiv,
                   pop = pops1$fourpop)
head(pops)
pops$id == genl@ind.names

#reorder pops
pops <- pops[match(genl@ind.names, pops$id), ]

all(pops$id == genl@ind.names)

genl@pop <- as.factor(pops$pop)
table(genl@pop)


gl2bayesAss(genl, outpath="analysis/pop_structure/bayesass", 
            outfile="fourpop.BayesAss.txt")

#---------------------------------------------------------------------------
# then the subset populations

genl_fourpop

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$fourpop)
head(pops)
pops$id == genl_fourpop@ind.names

#reorder pops
pops <- pops[match(genl_fourpop@ind.names, pops$id), ]
nrow(pops)
all(pops$id == genl_fourpop@ind.names)

genl_fourpop@pop <- as.factor(pops$pop)
table(genl_fourpop@pop)

gl2bayesAss(genl_fourpop, outpath="analysis/pop_structure/bayesass", 
            outfile="fourpop_subset.BayesAss.txt")

# sixpop

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$sixpop)
head(pops)
pops$id == genl_sixpop@ind.names

#reorder pops
pops <- pops[match(genl_sixpop@ind.names, pops$id), ]
nrow(pops)
all(pops$id == genl_sixpop@ind.names)

genl_sixpop@pop <- as.factor(pops$pop)
table(genl_sixpop@pop)

gl2bayesAss(genl_sixpop, outpath="analysis/pop_structure/bayesass", 
            outfile="sixpop_subset.BayesAss.txt")

