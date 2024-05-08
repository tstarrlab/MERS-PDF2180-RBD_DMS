Compute per-barcode binding to mAb samples
================
Tyler Starr
4/24/2024

This notebook reads in per-barcode counts from `count_variants.ipynb`
for mAb-binding titration experiments, computes functional scores for
RBD binding values via delta-AUC metrics, and does some basic QC on
variant binding functional scores.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$mAb_EC50_dir)){
  dir.create(file.path(config$mAb_EC50_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.8 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/linux-rocky8-nehalem/gcc-8.5.0/intel-oneapi-mkl-2021.4.0-h43nkmwzvaltaa6ii5l7n6e7ruvjbmnv/mkl/2021.4.0/lib/intel64/libmkl_rt.so.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [5] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ##  [9] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.30        haven_2.4.3      colorspace_2.0-3
    ##  [5] vctrs_0.5.2      generics_0.1.2   htmltools_0.5.2  utf8_1.2.2      
    ##  [9] rlang_1.0.6      pillar_1.7.0     glue_1.6.2       withr_2.5.0     
    ## [13] DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1    
    ## [17] lifecycle_1.0.3  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_1.0.2      evaluate_0.15    knitr_1.37       tzdb_0.2.0      
    ## [25] fastmap_1.1.0    fansi_1.0.2      broom_0.7.12     Rcpp_1.0.11     
    ## [29] backports_1.4.1  scales_1.2.1     jsonlite_1.8.7   fs_1.5.2        
    ## [33] hms_1.1.1        digest_0.6.29    stringi_1.7.6    grid_4.1.3      
    ## [37] cli_3.6.0        tools_4.1.3      magrittr_2.0.2   crayon_1.5.0    
    ## [41] pkgconfig_2.0.3  ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1    
    ## [45] lubridate_1.8.0  rstudioapi_0.13  assertthat_0.2.1 rmarkdown_2.13  
    ## [49] httr_1.4.7       R6_2.5.1         compiler_4.1.3

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#eliminate rows from barcode_runs that are not from an tite-seq experiment
barcode_runs <- barcode_runs[barcode_runs$sample_type %in% c("LCA60","S41","4C2","MERS4"),]

#read file giving count of each barcode in each sort partition
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F))

#eliminate rows from counts that are not part of an titration bin sample
counts <- subset(counts, sample %in% barcode_runs[barcode_runs$sample_type %in% c("LCA60","S41","4C2","MERS4"),"sample"])

#read in barcode-variant lookup tables
dt <- data.table(read.csv(file=config$codon_variant_table_file_pools,stringsAsFactors=F))

dt <- merge(counts, dt, by=c("library","barcode"));rm(counts)

samples_LCA60 <- data.frame(sample=sort(unique(paste(rep("LCA60",6),formatC(barcode_runs[barcode_runs$sample_type=="LCA60","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))

samples_S41 <- data.frame(sample=sort(unique(paste(rep("S41",6),formatC(barcode_runs[barcode_runs$sample_type=="S41","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))

samples_4C2 <- data.frame(sample=sort(unique(paste(rep("4C2",6),formatC(barcode_runs[barcode_runs$sample_type=="4C2","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))

samples_MERS4 <- data.frame(sample=sort(unique(paste(rep("MERS4",6),formatC(barcode_runs[barcode_runs$sample_type=="LCA60","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))
```

Convert from Illumina read counts to estimates of the number of cells
that were sorted into a bin, and add some other useful information to
our data frame.

``` r
#for each bin, normalize the read counts to the observed ratio of cell recovery among bins
for(i in 1:nrow(barcode_runs)){
  lib <- as.character(barcode_runs$library[i])
  bin <- as.character(barcode_runs$sample[i])
  ratio <- sum(dt[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
  if(ratio<1){ #if there are fewer reads from a FACS bin than cells sorted
    dt[library==lib & sample==bin, count.norm := as.numeric(count)] #don't normalize cell counts, make count.norm the same as count
    print(paste("reads < cells for",lib,bin,", un-normalized (ratio",ratio,")")) #print to console to inform of undersampled bins
  }else{
    dt[library==lib & sample==bin, count.norm := as.numeric(count/ratio)] #normalize read counts by the average read:cell ratio, report in new "count.norm" column
    print(paste("read:cell ratio for",lib,bin,"is",ratio))
  }
}
```

    ## [1] "read:cell ratio for pool1 LCA60_01_bin1 is 3.77015541081399"
    ## [1] "read:cell ratio for pool1 LCA60_01_bin2 is 4.58799948081808"
    ## [1] "read:cell ratio for pool1 LCA60_01_bin3 is 3.62249833355553"
    ## [1] "read:cell ratio for pool1 LCA60_01_bin4 is 3.91952275787885"
    ## [1] "read:cell ratio for pool1 LCA60_02_bin1 is 2.84986281299026"
    ## [1] "read:cell ratio for pool1 LCA60_02_bin2 is 3.95152707139308"
    ## [1] "read:cell ratio for pool1 LCA60_02_bin3 is 2.97745819999743"
    ## [1] "read:cell ratio for pool1 LCA60_02_bin4 is 4.31000677112005"
    ## [1] "read:cell ratio for pool1 LCA60_03_bin1 is 4.55553321828747"
    ## [1] "read:cell ratio for pool1 LCA60_03_bin2 is 3.85768146939668"
    ## [1] "read:cell ratio for pool1 LCA60_03_bin3 is 3.18588031825172"
    ## [1] "read:cell ratio for pool1 LCA60_03_bin4 is 1.50531746716105"
    ## [1] "read:cell ratio for pool1 LCA60_04_bin1 is 4.76769112747233"
    ## [1] "read:cell ratio for pool1 LCA60_04_bin2 is 4.02765868858232"
    ## [1] "read:cell ratio for pool1 LCA60_04_bin3 is 5.5980047772938"
    ## [1] "read:cell ratio for pool1 LCA60_04_bin4 is 18.5"
    ## [1] "read:cell ratio for pool1 LCA60_05_bin1 is 4.16523109869017"
    ## [1] "read:cell ratio for pool1 LCA60_05_bin2 is 4.610827594581"
    ## [1] "read:cell ratio for pool1 LCA60_05_bin3 is 2"
    ## [1] "read:cell ratio for pool1 LCA60_05_bin4 is 2.21052631578947"
    ## [1] "read:cell ratio for pool1 LCA60_06_bin1 is 5.51551004386608"
    ## [1] "read:cell ratio for pool1 LCA60_06_bin2 is 5.73855111904362"
    ## [1] "read:cell ratio for pool1 LCA60_06_bin3 is 1.31137724550898"
    ## [1] "read:cell ratio for pool1 LCA60_06_bin4 is 2.21368948247078"
    ## [1] "read:cell ratio for pool1 S41_01_bin1 is 4.64436567665065"
    ## [1] "read:cell ratio for pool1 S41_01_bin2 is 4.65445419164548"
    ## [1] "read:cell ratio for pool1 S41_01_bin3 is 3.44066485635942"
    ## [1] "read:cell ratio for pool1 S41_01_bin4 is 4.66717458754908"
    ## [1] "read:cell ratio for pool1 S41_02_bin1 is 4.60872687829816"
    ## [1] "read:cell ratio for pool1 S41_02_bin2 is 4.14902206189109"
    ## [1] "read:cell ratio for pool1 S41_02_bin3 is 4.03591069751347"
    ## [1] "read:cell ratio for pool1 S41_02_bin4 is 4.13957506643311"
    ## [1] "read:cell ratio for pool1 S41_03_bin1 is 3.45874029208876"
    ## [1] "read:cell ratio for pool1 S41_03_bin2 is 3.46837327139065"
    ## [1] "read:cell ratio for pool1 S41_03_bin3 is 3.7162967371597"
    ## [1] "read:cell ratio for pool1 S41_03_bin4 is 6.48498012492902"
    ## [1] "read:cell ratio for pool1 S41_04_bin1 is 6.09435014515407"
    ## [1] "read:cell ratio for pool1 S41_04_bin2 is 5.35868824241185"
    ## [1] "reads < cells for pool1 S41_04_bin3 , un-normalized (ratio 0.869757599076568 )"
    ## [1] "reads < cells for pool1 S41_04_bin4 , un-normalized (ratio 0.8125 )"
    ## [1] "read:cell ratio for pool1 S41_05_bin1 is 4.35486890051068"
    ## [1] "read:cell ratio for pool1 S41_05_bin2 is 4.29209783631232"
    ## [1] "read:cell ratio for pool1 S41_05_bin3 is 2.91106194690266"
    ## [1] "read:cell ratio for pool1 S41_05_bin4 is 1.98213446168312"
    ## [1] "read:cell ratio for pool1 S41_06_bin1 is 5.51551004386608"
    ## [1] "read:cell ratio for pool1 S41_06_bin2 is 5.73855111904362"
    ## [1] "read:cell ratio for pool1 S41_06_bin3 is 1.31137724550898"
    ## [1] "read:cell ratio for pool1 S41_06_bin4 is 2.21368948247078"
    ## [1] "read:cell ratio for pool1 4C2_01_bin1 is 3.34665684940444"
    ## [1] "read:cell ratio for pool1 4C2_01_bin2 is 3.92576274681092"
    ## [1] "read:cell ratio for pool1 4C2_01_bin3 is 3.76618825226867"
    ## [1] "read:cell ratio for pool1 4C2_01_bin4 is 3.50577853459656"
    ## [1] "read:cell ratio for pool1 4C2_02_bin1 is 4.22405282707388"
    ## [1] "read:cell ratio for pool1 4C2_02_bin2 is 7.49815237450919"
    ## [1] "reads < cells for pool1 4C2_02_bin3 , un-normalized (ratio 4.47889639992706e-05 )"
    ## [1] "read:cell ratio for pool1 4C2_02_bin4 is 4.22631879449147"
    ## [1] "read:cell ratio for pool1 4C2_03_bin1 is 4.56609839988036"
    ## [1] "read:cell ratio for pool1 4C2_03_bin2 is 3.64583942348022"
    ## [1] "read:cell ratio for pool1 4C2_03_bin3 is 4.61955945121049"
    ## [1] "read:cell ratio for pool1 4C2_03_bin4 is 3.51061864208805"
    ## [1] "read:cell ratio for pool1 4C2_04_bin1 is 4.83989459142003"
    ## [1] "read:cell ratio for pool1 4C2_04_bin2 is 4.76298506618015"
    ## [1] "read:cell ratio for pool1 4C2_04_bin3 is 3.53496106750533"
    ## [1] "read:cell ratio for pool1 4C2_04_bin4 is 1.63157894736842"
    ## [1] "read:cell ratio for pool1 4C2_05_bin1 is 3.96224617000225"
    ## [1] "read:cell ratio for pool1 4C2_05_bin2 is 4.73221293117868"
    ## [1] "read:cell ratio for pool1 4C2_05_bin3 is 3.99208443271768"
    ## [1] "read:cell ratio for pool1 4C2_05_bin4 is 1.26086956521739"
    ## [1] "read:cell ratio for pool1 4C2_06_bin1 is 5.51551004386608"
    ## [1] "read:cell ratio for pool1 4C2_06_bin2 is 5.73855111904362"
    ## [1] "read:cell ratio for pool1 4C2_06_bin3 is 1.31137724550898"
    ## [1] "read:cell ratio for pool1 4C2_06_bin4 is 2.21368948247078"
    ## [1] "read:cell ratio for pool1 MERS4_01_bin1 is 5.79314803766131"
    ## [1] "read:cell ratio for pool1 MERS4_01_bin2 is 4.97573664798637"
    ## [1] "read:cell ratio for pool1 MERS4_01_bin3 is 6.0703982926285"
    ## [1] "read:cell ratio for pool1 MERS4_01_bin4 is 5.1728102794544"
    ## [1] "read:cell ratio for pool1 MERS4_02_bin1 is 5.33189792314439"
    ## [1] "read:cell ratio for pool1 MERS4_02_bin2 is 4.4887110387653"
    ## [1] "read:cell ratio for pool1 MERS4_02_bin3 is 5.28773963891581"
    ## [1] "read:cell ratio for pool1 MERS4_02_bin4 is 5.58049628926777"
    ## [1] "read:cell ratio for pool1 MERS4_03_bin1 is 5.12517993115677"
    ## [1] "read:cell ratio for pool1 MERS4_03_bin2 is 6.11450179409329"
    ## [1] "read:cell ratio for pool1 MERS4_03_bin3 is 5.34868313000737"
    ## [1] "read:cell ratio for pool1 MERS4_03_bin4 is 7.43693356667546"
    ## [1] "read:cell ratio for pool1 MERS4_04_bin1 is 5.67213537668535"
    ## [1] "read:cell ratio for pool1 MERS4_04_bin2 is 5.99532042698452"
    ## [1] "read:cell ratio for pool1 MERS4_04_bin3 is 5.24842452105958"
    ## [1] "read:cell ratio for pool1 MERS4_04_bin4 is 2.90506329113924"
    ## [1] "read:cell ratio for pool1 MERS4_05_bin1 is 5.90687141614486"
    ## [1] "read:cell ratio for pool1 MERS4_05_bin2 is 5.28934097733238"
    ## [1] "read:cell ratio for pool1 MERS4_05_bin3 is 1.02535211267606"
    ## [1] "reads < cells for pool1 MERS4_05_bin4 , un-normalized (ratio 0.512396694214876 )"
    ## [1] "read:cell ratio for pool1 MERS4_06_bin1 is 5.51551004386608"
    ## [1] "read:cell ratio for pool1 MERS4_06_bin2 is 5.73855111904362"
    ## [1] "read:cell ratio for pool1 MERS4_06_bin3 is 1.31137724550898"
    ## [1] "read:cell ratio for pool1 MERS4_06_bin4 is 2.21368948247078"
    ## [1] "read:cell ratio for pool2 MERS4_01_bin1 is 4.44889545001129"
    ## [1] "read:cell ratio for pool2 MERS4_01_bin2 is 4.17767535725974"
    ## [1] "read:cell ratio for pool2 MERS4_01_bin3 is 4.76705903090603"
    ## [1] "read:cell ratio for pool2 MERS4_01_bin4 is 5.71116749029079"
    ## [1] "read:cell ratio for pool2 MERS4_02_bin1 is 4.77060668058962"
    ## [1] "read:cell ratio for pool2 MERS4_02_bin2 is 5.33339280181035"
    ## [1] "read:cell ratio for pool2 MERS4_02_bin3 is 6.70501480393474"
    ## [1] "read:cell ratio for pool2 MERS4_02_bin4 is 4.61797943441304"
    ## [1] "read:cell ratio for pool2 MERS4_03_bin1 is 4.4332033919105"
    ## [1] "read:cell ratio for pool2 MERS4_03_bin2 is 5.09257701871542"
    ## [1] "read:cell ratio for pool2 MERS4_03_bin3 is 4.21683304099206"
    ## [1] "read:cell ratio for pool2 MERS4_03_bin4 is 4.69078278146955"
    ## [1] "read:cell ratio for pool2 MERS4_04_bin1 is 4.31566644201291"
    ## [1] "read:cell ratio for pool2 MERS4_04_bin2 is 4.41008194777639"
    ## [1] "read:cell ratio for pool2 MERS4_04_bin3 is 4.76515550688253"
    ## [1] "reads < cells for pool2 MERS4_04_bin4 , un-normalized (ratio 0.798507462686567 )"
    ## [1] "read:cell ratio for pool2 MERS4_05_bin1 is 4.29497990640783"
    ## [1] "read:cell ratio for pool2 MERS4_05_bin2 is 4.12011051796451"
    ## [1] "reads < cells for pool2 MERS4_05_bin3 , un-normalized (ratio 0.825153374233129 )"
    ## [1] "reads < cells for pool2 MERS4_05_bin4 , un-normalized (ratio 0.3 )"
    ## [1] "read:cell ratio for pool2 MERS4_06_bin1 is 4.44551094511568"
    ## [1] "read:cell ratio for pool2 MERS4_06_bin2 is 4.6397358415049"
    ## [1] "reads < cells for pool2 MERS4_06_bin3 , un-normalized (ratio 0.905781584582441 )"
    ## [1] "reads < cells for pool2 MERS4_06_bin4 , un-normalized (ratio 0.532051282051282 )"

``` r
#annotate each barcode as to whether it's a homolog variant, SARS-CoV-2 wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
dt[,variant_class:=as.character(NA)]
dt[sublibrary %in% c("lib53","lib54","lib84", "lib85") & n_codon_substitutions==0, variant_class := "wildtype"]
dt[sublibrary %in% c("lib53","lib54","lib84", "lib85") & n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[sublibrary %in% c("lib53","lib54","lib84", "lib85") & n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[sublibrary %in% c("lib53","lib54","lib84", "lib85") & n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[sublibrary %in% c("lib53","lib54","lib84", "lib85") & n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

dt[sublibrary %in% c("lib82", "lib83"), variant_class := target]

#cast the data frame into wide format
dt <- dcast(dt, library + sublibrary + barcode + target + variant_class + aa_substitutions + n_aa_substitutions ~ sample, value.var="count.norm")
```

## Calculating mean bin for each barcode at each sample concentration

Next, for each barcode at each of the sera concentrations, calculate the
“mean bin” response variable. This is calculated as a simple mean, where
the value of each bin is the integer value of the bin (bin1=unbound,
bin4=highly bound) – because of how bins are defined, the mean
fluorescence of cells in each bin are equally spaced on a log-normal
scale, so mean bin correlates with simple mean fluorescence.

We do not use the fluorescence boundaries of the FACS bins in our
calculations here, but we provide them for posterity’s sake below.

\`\`

``` r
#function that returns mean bin and sum of counts for four bins cell counts. Includes cutoffs for bimodal sample splits to filter out
calc.meanbin <- function(vec, split13filter=0.4, split24filter=0.4, split14filter=0.2){
  total <- sum(vec)
  if(is.na(total) | (vec[1] > split13filter*total & vec[3] > split13filter*total) | (vec[2] > split24filter*total & vec[4] > split24filter*total) | (vec[1] > split14filter*total & vec[4] > split14filter*total)){
    return(list(NA,NA))
  }else{
    return( list((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4]), total) )
  }
}
  

#iterate through titration samples, compute mean_bin and total_count for each barcode variant
#LCA60
for(i in 1:nrow(samples_LCA60)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_LCA60[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_LCA60[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_LCA60[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_LCA60[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_LCA60[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_LCA60[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#S41
for(i in 1:nrow(samples_S41)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_S41[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_S41[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_S41[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_S41[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_S41[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_S41[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#4C2
for(i in 1:nrow(samples_4C2)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_4C2[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_4C2[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_4C2[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_4C2[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_4C2[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_4C2[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#MERS4
for(i in 1:nrow(samples_MERS4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_MERS4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_MERS4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_MERS4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_MERS4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_MERS4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_MERS4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}
```

## Fit binding curves

We will calculate a simple AUC metric across each barcode’s titration
series. We will also include a minimum cell count that is required for a
meanbin estimate to be used in the titration fit, and a minimum number
of concentrations with determined meanbin that is required for a
titration to be reported.

``` r
#For QC and filtering, output columns giving the average number of cells that were sampled for a barcode across the 9 sample concentrations, and a value for the number of meanbin estimates that were removed for being below the # of cells cutoff
cutoff <- 2

dt[,`LCA60_avgcount` := mean(c(`LCA60_01_totalcount`,`LCA60_02_totalcount`,`LCA60_03_totalcount`,
                             `LCA60_04_totalcount`,`LCA60_05_totalcount`,`LCA60_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`LCA60_min_cell_filtered` := sum(c(c(`LCA60_01_totalcount`,`LCA60_02_totalcount`,`LCA60_03_totalcount`,
                                       `LCA60_04_totalcount`,`LCA60_05_totalcount`,`LCA60_06_totalcount`)<cutoff,
                                     is.na(c(`LCA60_01_totalcount`,`LCA60_02_totalcount`,`LCA60_03_totalcount`,
                                             `LCA60_04_totalcount`,`LCA60_05_totalcount`,`LCA60_06_totalcount`))),na.rm=T),by=c("library","barcode")]

dt[,`S41_avgcount` := mean(c(`S41_01_totalcount`,`S41_02_totalcount`,`S41_03_totalcount`,
                             `S41_04_totalcount`,`S41_05_totalcount`,`S41_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`S41_min_cell_filtered` := sum(c(c(`S41_01_totalcount`,`S41_02_totalcount`,`S41_03_totalcount`,
                                       `S41_04_totalcount`,`S41_05_totalcount`,`S41_06_totalcount`)<cutoff,
                                     is.na(c(`S41_01_totalcount`,`S41_02_totalcount`,`S41_03_totalcount`,
                                             `S41_04_totalcount`,`S41_05_totalcount`,`S41_06_totalcount`))),na.rm=T),by=c("library","barcode")]

dt[,`4C2_avgcount` := mean(c(`4C2_01_totalcount`,`4C2_02_totalcount`,`4C2_03_totalcount`,
                             `4C2_04_totalcount`,`4C2_05_totalcount`,`4C2_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`4C2_min_cell_filtered` := sum(c(c(`4C2_01_totalcount`,`4C2_02_totalcount`,`4C2_03_totalcount`,
                                       `4C2_04_totalcount`,`4C2_05_totalcount`,`4C2_06_totalcount`)<cutoff,
                                     is.na(c(`4C2_01_totalcount`,`4C2_02_totalcount`,`4C2_03_totalcount`,
                                             `4C2_04_totalcount`,`4C2_05_totalcount`,`4C2_06_totalcount`))),na.rm=T),by=c("library","barcode")]

dt[,`MERS4_avgcount` := mean(c(`MERS4_01_totalcount`,`MERS4_02_totalcount`,`MERS4_03_totalcount`,
                             `MERS4_04_totalcount`,`MERS4_05_totalcount`,`MERS4_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`MERS4_min_cell_filtered` := sum(c(c(`MERS4_01_totalcount`,`MERS4_02_totalcount`,`MERS4_03_totalcount`,
                                       `MERS4_04_totalcount`,`MERS4_05_totalcount`,`MERS4_06_totalcount`)<cutoff,
                                     is.na(c(`MERS4_01_totalcount`,`MERS4_02_totalcount`,`MERS4_03_totalcount`,
                                             `MERS4_04_totalcount`,`MERS4_05_totalcount`,`MERS4_06_totalcount`))),na.rm=T),by=c("library","barcode")]


#function that fits a nls regression to the titration series, including an option to filter below certain thresholds for average cells across all samples, and number of samples below a cutoff of cells
fit.titration <- function(y.vals,x.vals,count.vals,min.cfu=cutoff,
                          min.means=1,min.average=7,EC50.start=10,
                          a.start=3,a.lower=2,a.upper=3,
                          b.start=1,b.lower=1,b.upper=1.5,
                          n.start=1,n.lower=1/2,n.upper=2){
  indices <- count.vals>min.cfu & !is.na(y.vals)
  y <- y.vals[indices]
  x <- x.vals[indices]
  if((length(y) < min.means*length(y.vals)) | (mean(count.vals,na.rm=T) < min.average)){ #return NAs if < min.means fraction of concentrations have above min.cfu counts or if the average count across all concentrations is below min.average
    return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))
  }else{
    fit <- nls(y ~ (a/(1+(EC50/x)^n))+b,
               start=list(a=a.start,b=b.start,n=n.start,EC50=EC50.start),
               lower=list(a=a.lower,b=b.lower,n=n.lower,EC50=min(x.vals[x.vals>0])/1000), #constrain EC50 to be no lower than 1/100x the lowest concentration value
               upper=list(a=a.upper,b=b.upper,n=n.upper,EC50=1e5), #constrain EC50 to be no higher than the 10x highest concentration value
               algorithm="port")
    y.pred <- predict(fit,newdata=list(x=x))
    resid <- y - y.pred
    resid.norm <- resid/as.numeric(summary(fit)$coefficients["a","Estimate"])
    nMSR <- mean((resid.norm)^2,na.rm=T)
    return(list(as.numeric(summary(fit)$coefficients["EC50","Estimate"]),
                as.numeric(summary(fit)$coefficients["EC50","Std. Error"]),
                as.numeric(summary(fit)$coefficients["a","Estimate"]),
                as.numeric(summary(fit)$coefficients["b","Estimate"]),
                as.numeric(nMSR)))
  }
}

#fit titration to LCA60 data for each barcode
dt[,c("EC50_LCA60","EC50_SE_LCA60","response_LCA60","baseline_LCA60","nMSR_LCA60") :=
     tryCatch(fit.titration(y.vals=c(`LCA60_01_meanbin`,`LCA60_02_meanbin`,`LCA60_03_meanbin`,
                                     `LCA60_04_meanbin`,`LCA60_05_meanbin`,`LCA60_06_meanbin`),
                            x.vals=samples_LCA60$conc,
                            count.vals=c(`LCA60_01_totalcount`,`LCA60_02_totalcount`,`LCA60_03_totalcount`,
                                         `LCA60_04_totalcount`,`LCA60_05_totalcount`,`LCA60_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]


#fit titration to S41 data for each barcode
dt[,c("EC50_S41","EC50_SE_S41","response_S41","baseline_S41","nMSR_S41") :=
     tryCatch(fit.titration(y.vals=c(`S41_01_meanbin`,`S41_02_meanbin`,`S41_03_meanbin`,
                                     `S41_04_meanbin`,`S41_05_meanbin`,`S41_06_meanbin`),
                            x.vals=samples_S41$conc,
                            count.vals=c(`S41_01_totalcount`,`S41_02_totalcount`,`S41_03_totalcount`,
                                         `S41_04_totalcount`,`S41_05_totalcount`,`S41_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]


#fit titration to 4C2 data for each barcode
dt[,c("EC50_4C2","EC50_SE_4C2","response_4C2","baseline_4C2","nMSR_4C2") :=
     tryCatch(fit.titration(y.vals=c(`4C2_01_meanbin`,`4C2_02_meanbin`,`4C2_03_meanbin`,
                                     `4C2_04_meanbin`,`4C2_05_meanbin`,`4C2_06_meanbin`),
                            x.vals=samples_4C2$conc,
                            count.vals=c(`4C2_01_totalcount`,`4C2_02_totalcount`,`4C2_03_totalcount`,
                                         `4C2_04_totalcount`,`4C2_05_totalcount`,`4C2_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]


#fit titration to MERS4 data for each barcode
dt[,c("EC50_MERS4","EC50_SE_MERS4","response_MERS4","baseline_MERS4","nMSR_MERS4") :=
     tryCatch(fit.titration(y.vals=c(`MERS4_01_meanbin`,`MERS4_02_meanbin`,`MERS4_03_meanbin`,
                                     `MERS4_04_meanbin`,`MERS4_05_meanbin`,`MERS4_06_meanbin`),
                            x.vals=samples_MERS4$conc,
                            count.vals=c(`MERS4_01_totalcount`,`MERS4_02_totalcount`,`MERS4_03_totalcount`,
                                         `MERS4_04_totalcount`,`MERS4_05_totalcount`,`MERS4_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]

#save temp data for later if wnat to load post-fits to optimize script below
save(dt, file=paste(config$mAb_EC50_dir,"/EC50_temp_mAbs.Rda",sep=""))

#load(file=paste(config$mAb_EC50_dir,"/EC50_temp_mAbs.Rda",sep=""))
```

## QC and sanity checks

We will do some QC to make sure we got good titration curves for most of
our library barcodes. We will also spot check titration curves from
across our measurement range, and spot check curves whose fit parameters
hit the different boundary conditions of the fit variables.

We successfully generated EC50 estimates for 0.2507405 of our LCA60 mAb
measurements, 0.2591319 of the S41 mAb, 0.2259368 of 4C2, 0.473588 of
MERS4 mAb.

``` r
#make functions that allow me to plot a titration for any given row from the counts data frames, for spot checking curves
plot.titration <- function(row,mAb,output.text=F){
  y.vals <- c();for(sample in get(paste("samples_",mAb,sep=""))$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
  x.vals <- get(paste("samples_",mAb,sep=""))$conc
  count.vals <- c();for(sample in get(paste("samples_",mAb,sep=""))$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
  title <- paste(dt[row,target]," ",mAb,sep="")
  indices <- count.vals>cutoff & !is.na(count.vals)
  y.vals <- y.vals[indices]
  x.vals <- x.vals[indices]
  plot(x.vals,y.vals,xlab=paste("[",mAb,"] (ng/mL)",sep=""),
       ylab="mean bin",log="x",ylim=c(1,4),xlim=c(0.01,11000),pch=19,main=title)
  EC50_var <- paste("EC50_",mAb,sep="")
  fit <- nls(y.vals ~ (a/(1+(EC50/x.vals)^n))+b,
             start=list(a=3,b=1,EC50=dt[row,get(EC50_var)],n=1),
             lower=list(a=2,b=1,EC50=0.001,n=1/2),
             upper=list(a=3,b=1.25,EC50=100000,n=2),
             algorithm="port")
  if(!is.na(dt[row,get(EC50_var)])){
    lines(10^c(seq(-3,5,0.25)),predict(fit,newdata=list(x.vals=10^c(seq(-3,5,0.25)))))
    legend("topleft",bty="n",cex=1,legend=paste("EC50",format(dt[row,get(EC50_var)],digits=3),"ng/mL"))
  }
  if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
    vars <- c("barcode","variant_class","wildtype","position","mutant",as.character(paste(mAb,"_avgcount",sep="")),as.character(paste(mAb,"_min_cell_filtered",sep="")),as.character(paste("EC50_",mAb,sep="")),as.character(paste("EC50_SE_",mAb,sep="")),as.character(paste("baseline_",mAb,sep="")),as.character(paste("response_",mAb,sep="")),as.character(paste("nMSR_",mAb,sep="")))
    return(dt[row,..vars])
  }
}
```

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_LCA60>=10000)[1],"LCA60")
plot.titration(which(dt$EC50_LCA60>=10000)[2],"LCA60")
plot.titration(which(dt$EC50_LCA60>=10000)[3],"LCA60")
plot.titration(which(dt$EC50_LCA60>1 & dt$EC50_LCA60<10)[1],"LCA60")
plot.titration(which(dt$EC50_LCA60>1 & dt$EC50_LCA60<10)[2],"LCA60")
plot.titration(which(dt$EC50_LCA60>1 & dt$EC50_LCA60<10)[3],"LCA60")
plot.titration(which(dt$EC50_LCA60<1)[1],"LCA60")
plot.titration(which(dt$EC50_LCA60<1)[2],"LCA60")
plot.titration(which(dt$EC50_LCA60<1)[3],"LCA60")
```

<img src="compute_mAb_EC50_files/figure-gfm/EC50_LCA60-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_S41>=10000)[1],"S41")
plot.titration(which(dt$EC50_S41>=10000)[2],"S41")
plot.titration(which(dt$EC50_S41>=10000)[3],"S41")
plot.titration(which(dt$EC50_S41>1 & dt$EC50_S41<10)[1],"S41")
plot.titration(which(dt$EC50_S41>1 & dt$EC50_S41<10)[2],"S41")
plot.titration(which(dt$EC50_S41>1 & dt$EC50_S41<10)[3],"S41")
plot.titration(which(dt$EC50_S41<1)[1],"S41")
plot.titration(which(dt$EC50_S41<1)[2],"S41")
plot.titration(which(dt$EC50_S41<1)[3],"S41")
```

<img src="compute_mAb_EC50_files/figure-gfm/EC50_S41-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_4C2>=10000)[1],"4C2")
plot.titration(which(dt$EC50_4C2>=10000)[2],"4C2")
plot.titration(which(dt$EC50_4C2>=10000)[3],"4C2")
plot.titration(which(dt$EC50_4C2>1 & dt$EC50_4C2<10)[1],"4C2")
plot.titration(which(dt$EC50_4C2>1 & dt$EC50_4C2<10)[2],"4C2")
plot.titration(which(dt$EC50_4C2>1 & dt$EC50_4C2<10)[3],"4C2")
plot.titration(which(dt$EC50_4C2<0.5)[1],"4C2")
plot.titration(which(dt$EC50_4C2<0.5)[2],"4C2")
plot.titration(which(dt$EC50_4C2<0.5)[3],"4C2")
```

<img src="compute_mAb_EC50_files/figure-gfm/EC50_4C2-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_MERS4>=10000)[1],"MERS4")
plot.titration(which(dt$EC50_MERS4>=10000)[2],"MERS4")
plot.titration(which(dt$EC50_MERS4>=10000)[3],"MERS4")
plot.titration(which(dt$EC50_MERS4>1 & dt$EC50_MERS4<10)[1],"MERS4")
plot.titration(which(dt$EC50_MERS4>1 & dt$EC50_MERS4<10)[2],"MERS4")
plot.titration(which(dt$EC50_MERS4>1 & dt$EC50_MERS4<10)[3],"MERS4")
plot.titration(which(dt$EC50_MERS4<0.5)[1],"MERS4")
plot.titration(which(dt$EC50_MERS4<0.5)[2],"MERS4")
plot.titration(which(dt$EC50_MERS4<0.5)[3],"MERS4")
```

<img src="compute_mAb_EC50_files/figure-gfm/EC50_MERS4-1.png" style="display: block; margin: auto;" />

## Data filtering by fit quality

Next, let’s filter out poor fits using the value we previously computed,
the *normalized* mean square residual (nMSR). This metric computes the
residual between the observed response variable and that predicted from
the titration fit, normalizes this residual by the response range of the
titration fit (which is allowed to vary between 1.5 and 3 per the
titration fits above), and computes the mean-square of these normalized
residuals.

Distribution of the nMSR metric in each set of fits

``` r
par(mfrow=c(2,2))
hist(dt$nMSR_LCA60,main="LCA60",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
hist(dt$nMSR_S41,main="S41",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
hist(dt$nMSR_4C2,main="4C2",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
hist(dt$nMSR_MERS4,main="MERS4",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
```

<img src="compute_mAb_EC50_files/figure-gfm/nMSR_distribution-1.png" style="display: block; margin: auto;" />

As we would expect, the MSR stat decreases with cell count, indicating
that higher cell counts leads to better curve fits. Also show the cutoff
I’m proposing for nMSR (10x median across all fits), legend gives
percent of curve fits eliminated

``` r
median.nMSR <- median(c(dt$nMSR_LCA60,dt$nMSR_S41,dt$nMSR_4C2,dt$nMSR_MERS4),na.rm=T)

threshold <- 30

par(mfrow=c(3,3))
plot(log10(dt$`LCA60_avgcount`),dt$nMSR_LCA60,main="LCA60",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=threshold*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_LCA60 > 10*median.nMSR & !is.na(nMSR_LCA60),])/nrow(dt[!is.na(nMSR_LCA60),]),digits=3),"%"))

plot(log10(dt$`S41_avgcount`),dt$nMSR_S41,main="S41",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=threshold*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_S41 > 10*median.nMSR & !is.na(nMSR_S41),])/nrow(dt[!is.na(nMSR_S41),]),digits=3),"%"))

plot(log10(dt$`4C2_avgcount`),dt$nMSR_4C2,main="4C2",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=threshold*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_4C2 > 10*median.nMSR & !is.na(nMSR_4C2),])/nrow(dt[!is.na(nMSR_4C2),]),digits=3),"%"))

plot(log10(dt$`MERS4_avgcount`),dt$nMSR_MERS4,main="MERS4",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=threshold*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_MERS4 > 10*median.nMSR & !is.na(nMSR_MERS4),])/nrow(dt[!is.na(nMSR_MERS4),]),digits=3),"%"))
```

<img src="compute_mAb_EC50_files/figure-gfm/nMSR_v_cell_count-1.png" style="display: block; margin: auto;" />

Next, we will apply this filtering step on normalized MSR, removing
curves with nMSR \>10x the median across all experiments

``` r
dt[nMSR_LCA60 > threshold*median.nMSR,c("EC50_LCA60","EC50_SE_LCA60","response_LCA60","baseline_LCA60") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

dt[nMSR_S41 > threshold*median.nMSR,c("EC50_S41","EC50_SE_S41","response_S41","baseline_S41") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

dt[nMSR_4C2 > threshold*median.nMSR,c("EC50_4C2","EC50_SE_4C2","response_4C2","baseline_4C2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

dt[nMSR_MERS4 > threshold*median.nMSR,c("EC50_MERS4","EC50_SE_MERS4","response_MERS4","baseline_MERS4") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
```

## Final scores

Let’s visualize the EC50 binding measurements as violin plots for the
different wildtype targets, for each mAb.

LCA60

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib84","lib85") & !is.na(EC50_LCA60),],aes(x=variant_class,y=EC50_LCA60))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("LCA60 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_LCA60_DMS-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_LCA60_DMS.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib82","lib83") & !is.na(EC50_LCA60),],aes(x=target,y=EC50_LCA60))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("LCA60 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_LCA60_panmerbeco-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_LCA60_panmerbeco.pdf",sep="")))
```

S41

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib84","lib85") & !is.na(EC50_S41),],aes(x=variant_class,y=EC50_S41))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S41 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_S41_DMS-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_S41_DMS.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib82","lib83") & !is.na(EC50_S41),],aes(x=target,y=EC50_S41))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S41 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_S41_panmerbeco-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_S41_panmerbeco.pdf",sep="")))
```

4C2

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib84","lib85") & !is.na(EC50_4C2),],aes(x=variant_class,y=EC50_4C2))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("4C2 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_4C2_DMS-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_4C2_DMS.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib82","lib83") & !is.na(EC50_4C2),],aes(x=target,y=EC50_4C2))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("4C2 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_4C2_panmerbeco-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_4C2_panmerbeco.pdf",sep="")))
```

MERS4

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib84","lib85") & !is.na(EC50_MERS4),],aes(x=variant_class,y=EC50_MERS4))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("MERS4 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_MERS4_DMS-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_MERS4_DMS.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[sublibrary %in% c("lib82","lib83") & !is.na(EC50_MERS4),],aes(x=target,y=EC50_MERS4))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("MERS4 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,ncol=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_mAb_EC50_files/figure-gfm/binding_distribution_vioplot_MERS4_panmerbeco-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_MERS4_panmerbeco.pdf",sep="")))
```

## Save barcode-level metrics

In the next script, we will collapse bcs down to final
mutant/variant-level phenotypes, integrate things like expression
effects of variants, and visualize final phenotypes.

``` r
dt[,.(library,sublibrary,barcode,target,variant_class,aa_substitutions,n_aa_substitutions,
     `LCA60_avgcount`,EC50_LCA60,
     `S41_avgcount`,EC50_S41,
     `4C2_avgcount`,EC50_4C2,
     `MERS4_avgcount`,EC50_MERS4)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$mAb_EC50_file, row.names=F)
```
