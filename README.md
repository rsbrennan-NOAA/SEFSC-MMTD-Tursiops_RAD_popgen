# Tursiops popgen and hybridization

## Analysis

### alignment and variant calling

- trim: `fastp.nf`  
- align: `bwa.nf` and merge: `merge_bams.nf`  
- Remove low mapping samples: `filter_aligned.R`  
- Call Variants: `freebayes.sh`  
- Concat variants: `vcf_concat.sh`  
- Filter: `filter.1.sh`, `filter.2.sh`, `filter_missing.R`, `filter.3.sh`  
- ld thin: `ld_thin.sh`  

### Map

Initial map and calculation of depth, distance from shore: `map.R`

Map with population assignments, from population structure analyses below: `map_populations.R`

### population structure

- Admixture: `admixture.sh` + `admixture.R`
- DAPC: `dapc.R`
- Fst: `fst.R`

### outgroup

- download aduncus for outgroup: `sra_aduncus.sh`
- align: `bwa_aduncus.nf` and mark duplicates `mkdups_aduncus.sh`
- calculate alignment stats: `bam_stats_aduncus.sh`
- call variants: `freebayes_aduncus.sh` and merge `vcf_concat_aduncus.sh`
- merge aduncus with tursiops: `tursiops_aduncus_intersect.sh`

### all sites vcf

- Freebayes with all sites output: `freebayes_allsites.sh`
- concat and filter: `data_processing/filter_cat_allSites.nf`

- Genetic diversity:  
  - `pixy_fourpop.sh` +  
  - Heterozygosity and FIS: `Fis_het.sh`
  - plot and analyze: `diversity.R`

### Introgression

- f stats and admixture graphs: `admixtools.R`
- ABBA/BABA: `dstats.sh` and `dstats.R`
- Treemix: `treemix.md`. Note there are references to other scripts in here.  
- New Hybrids: `to_newhybrids.sh`, `newhybrids.sh`, `newhybrids.R`
- triangle plot: `triangle_plot.R`
- bgc: `bgc.R`

### Figures

Fig. 1: `Figure_01.R`  
Fig. 2: `Figure_02.R`  
Fig. 3: `Figure_03.R`  

#### Supplemental figures

Still need to add where/how these are generated.  

---

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
