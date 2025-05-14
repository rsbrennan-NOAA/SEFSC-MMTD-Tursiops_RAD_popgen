fourpop_int <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_intermediate/fourpop_intermediate.final.summary",header=T)

plot((fourpop_int$year/1000),(fourpop_int$Ne_median), 
     type="n", xlab="Time (1k years ago)", 
     ylab="Effective Population Size",ylim=c(0,50000), xlim=c(0,300))

lines((fourpop_int$year/1000),(fourpop_int$Ne_median),type="s",col="#4DAF4A",lwd = 5)
lines((fourpop_int$year/1000),(fourpop_int$Ne_2.5.),type="s",col="#4DAF4A",lty=3)
lines((fourpop_int$year/1000),(fourpop_int$Ne_97.5.),type="s",col="#4DAF4A",lty=3)

fourpop_offshore <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_offshore/fourpop_offshore.final.summary",header=T)

lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_median),type="s",col="#E69F00",lwd = 5)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_2.5.),type="s",col="#E69F00",lty=3)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_97.5.),type="s",col="#E69F00",lty=3)


fourpop_offshore <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_coastal_atl/fourpop_Coastal_Atlantic.final.summary",header=T)

lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_median),type="s",col="#56B4E9",lwd = 5)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_2.5.),type="s",col="#56B4E9",lty=3)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_97.5.),type="s",col="#56B4E9",lty=3)

fourpop_offshore <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_coastal_gulf/fourpop_Coastal_Gulf.final.summary",header=T)

lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_median),type="s",col="#004488",lwd = 5)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_2.5.),type="s",col="#004488",lty=3)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_97.5.),type="s",col="#004488",lty=3)

# coastal
  #56B4E9
  #004488

#Offshore:
  #F0B800
  #B65A00

#intermediate:
  #1B9E77
  #66A61E

#------------------------------------------------------------------
# no singletons

fourpop_Coastal_Gulf <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_intermediate_noSingletons/fourpop_intermediate.final.summary",header=T)

plot((fourpop_Coastal_Gulf$year/1000),(fourpop_Coastal_Gulf$Ne_median), 
     type="n", xlab="Time (1k years ago)", 
     ylab="Effective Population Size (1k individuals)",ylim=c(0,50000), xlim=c(0,200))

lines((fourpop_Coastal_Gulf$year/1000),(fourpop_Coastal_Gulf$Ne_median),type="s",col="dodgerblue3",lwd = 5)
lines((fourpop_Coastal_Gulf$year/1000),(fourpop_Coastal_Gulf$Ne_2.5.),type="s",col="dodgerblue3",lty=3)
lines((fourpop_Coastal_Gulf$year/1000),(fourpop_Coastal_Gulf$Ne_97.5.),type="s",col="dodgerblue3",lty=3)


fourpop_offshore <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_offshore_noSingletons/fourpop_offshore.final.summary",header=T)

lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_median),type="s",col="orange",lwd = 5)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_2.5.),type="s",col="orange",lty=3)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_97.5.),type="s",col="orange",lty=3)


fourpop_offshore <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_coastal_atl_noSingletons/fourpop_Coastal_Atlantic.final.summary",header=T)

lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_median),type="s",col="red",lwd = 5)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_2.5.),type="s",col="red",lty=3)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_97.5.),type="s",col="red",lty=3)

fourpop_offshore <- read.table("analysis/moments/stairway_plot_v2.1.1/fourpop_coastal_gulf_noSingletons/fourpop_Coastal_Gulf.final.summary",header=T)

lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_median),type="s",col="green",lwd = 5)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_2.5.),type="s",col="green",lty=3)
lines((fourpop_offshore$year/1000),(fourpop_offshore$Ne_97.5.),type="s",col="green",lty=3)








#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------
# 6 pop
#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------


sixpop_coastal_atlantic <- read.table("analysis/moments/stairway_plot_v2.1.1/sixpop_coastal_atlantic/sixpop_coastal_atlantic.final.summary",header=T)
sixpop_coastal_gulf <- read.table("analysis/moments/stairway_plot_v2.1.1/sixpop_coastal_gulf/sixpop_coastal_gulf.final.summary",header=T)
sixpop_int_atlantic <- read.table("analysis/moments/stairway_plot_v2.1.1/sixpop_int_atlantic/sixpop_Intermediate_Atlantic.final.summary",header=T)
sixpop_int_gulf <- read.table("analysis/moments/stairway_plot_v2.1.1/sixpop_int_gulf/sixpop_Intermediate_Gulf.final.summary",header=T)
sixpop_off_atlantic <- read.table("analysis/moments/stairway_plot_v2.1.1/sixpop_off_atlantic/sixpop_Offshore_Atlantic.final.summary",header=T)
sixpop_off_gulf <- read.table("analysis/moments/stairway_plot_v2.1.1/sixpop_off_gulf/sixpop_Offshore_Gulf.final.summary",header=T)

pop <- sixpop_coastal_atlantic
col <- "orange"

plot((sixpop_coastal_atlantic$year/1000),(sixpop_coastal_atlantic$Ne_median), 
     type="n", xlab="Time (1k years ago)", 
     ylab="Effective Population Size (1k individuals)",ylim=c(0,50000), xlim=c(0,200))

lines((pop$year/1000),(pop$Ne_median),type="s",col=col,lwd = 5)
lines((pop$year/1000),(pop$Ne_2.5.),type="s",col=col,lty=3)
lines((pop$year/1000),(pop$Ne_97.5.),type="s",col=col,lty=3)

pop <- sixpop_coastal_gulf
col <- "blue"

lines((pop$year/1000),(pop$Ne_median),type="s",col=col,lwd = 5)
lines((pop$year/1000),(pop$Ne_2.5.),type="s",col=col,lty=3)
lines((pop$year/1000),(pop$Ne_97.5.),type="s",col=col,lty=3)

pop <- sixpop_int_atlantic
col <- "purple"

lines((pop$year/1000),(pop$Ne_median),type="s",col=col,lwd = 5)
lines((pop$year/1000),(pop$Ne_2.5.),type="s",col=col,lty=3)
lines((pop$year/1000),(pop$Ne_97.5.),type="s",col=col,lty=3)

pop <- sixpop_int_gulf
col <- "orange"

lines((pop$year/1000),(pop$Ne_median),type="s",col=col,lwd = 5)
lines((pop$year/1000),(pop$Ne_2.5.),type="s",col=col,lty=3)
lines((pop$year/1000),(pop$Ne_97.5.),type="s",col=col,lty=3)

pop <- sixpop_off_atlantic
col <- "red"

lines((pop$year/1000),(pop$Ne_median),type="s",col=col,lwd = 5)
lines((pop$year/1000),(pop$Ne_2.5.),type="s",col=col,lty=3)
lines((pop$year/1000),(pop$Ne_97.5.),type="s",col=col,lty=3)

pop <- sixpop_off_gulf
col <- "red"

lines((pop$year/1000),(pop$Ne_median),type="s",col=col,lwd = 5)
lines((pop$year/1000),(pop$Ne_2.5.),type="s",col=col,lty=3)
lines((pop$year/1000),(pop$Ne_97.5.),type="s",col=col,lty=3)


