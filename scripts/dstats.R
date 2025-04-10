
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#d stats
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

dstat_four <- read.table("analysis/pop_structure/dstats/fourpop_all_BBAA.txt", header=T)
dstat_four

dstat_four$p.value_multTesting <- p.adjust(dstat_four$p.value,method="bonferroni")

write.table(dstat_four,"analysis/pop_structure/dstats/dstats_fourpop.txt",
            sep="\t", quote=F, row.names=F)
# remember, always arranged to be positive. 
# p1 and p2 can be flipped, would make negative.
# as is, indicates gene flow between p2 and p3. 
# so 
dstat_six <- read.table("analysis/pop_structure/dstats/sixpop_all_BBAA.txt", header=T)
dstat_six$p.value_multTesting <- p.adjust(dstat_six$p.value,method="bonferroni")

dstat_six[dstat_six$p.value_multTesting < 0.05,]
write.table(dstat_six,"analysis/pop_structure/dstats/dstats_sixpop.txt",
            sep="\t", quote=F, row.names=F)
#

