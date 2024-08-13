Compute per-barcode binding to MERS serum sample
================
Tyler Starr
4/25/2024

This notebook reads in per-barcode counts from `count_variants.ipynb`
for sera-binding titration experiments, computes functional scores for
RBD binding values via delta-AUC metrics, and does some basic QC on
variant binding functional scores.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

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

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$sera_delta_AUC_dir)){
  dir.create(file.path(config$sera_delta_AUC_dir))
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
    ## [13] knitr_1.37       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.30        haven_2.4.3      colorspace_2.0-3
    ##  [5] vctrs_0.5.2      generics_0.1.2   htmltools_0.5.2  utf8_1.2.2      
    ##  [9] rlang_1.0.6      pillar_1.7.0     glue_1.6.2       withr_2.5.0     
    ## [13] DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1    
    ## [17] lifecycle_1.0.3  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_1.0.2      evaluate_0.15    tzdb_0.2.0       fastmap_1.1.0   
    ## [25] fansi_1.0.2      broom_0.7.12     Rcpp_1.0.11      backports_1.4.1 
    ## [29] scales_1.2.1     jsonlite_1.8.7   fs_1.5.2         hms_1.1.1       
    ## [33] digest_0.6.29    stringi_1.7.6    grid_4.1.3       cli_3.6.0       
    ## [37] tools_4.1.3      magrittr_2.0.2   crayon_1.5.0     pkgconfig_2.0.3 
    ## [41] ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0 
    ## [45] rstudioapi_0.13  assertthat_0.2.1 rmarkdown_2.13   httr_1.4.7      
    ## [49] R6_2.5.1         compiler_4.1.3

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#eliminate rows from barcode_runs that are not from an tite-seq experiment
barcode_runs <- barcode_runs[barcode_runs$sample_type == "serum",]

#read file giving count of each barcode in each sort partition
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F))

#eliminate rows from counts that are not part of an titration bin sample
counts <- subset(counts, sample %in% barcode_runs[barcode_runs$sample_type=="serum","sample"])

#read in barcode-variant lookup tables
dt <- data.table(read.csv(file=config$codon_variant_table_file_pools,stringsAsFactors=F))

dt <- merge(counts, dt, by=c("library","barcode"));rm(counts)

#make tables giving names of Titeseq samples and the corresponding serum incubation concentrations

samples_serum <- data.frame(sample=sort(unique(paste(rep("serum",5),formatC(barcode_runs[barcode_runs$sample_type=="serum","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/1000,1/10000,1/100000,0))
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

    ## [1] "read:cell ratio for pool1 serum_01_bin1 is 5.47109192170237"
    ## [1] "read:cell ratio for pool1 serum_01_bin2 is 4.64881856761146"
    ## [1] "read:cell ratio for pool1 serum_01_bin3 is 4.9089864744183"
    ## [1] "read:cell ratio for pool1 serum_01_bin4 is 4.78560181874053"
    ## [1] "read:cell ratio for pool1 serum_02_bin1 is 5.93467342641424"
    ## [1] "read:cell ratio for pool1 serum_02_bin2 is 4.86111104653731"
    ## [1] "read:cell ratio for pool1 serum_02_bin3 is 5.2010204185718"
    ## [1] "read:cell ratio for pool1 serum_02_bin4 is 6.65282145021973"
    ## [1] "read:cell ratio for pool1 serum_03_bin1 is 5.24823877929502"
    ## [1] "read:cell ratio for pool1 serum_03_bin2 is 5.27033347304287"
    ## [1] "read:cell ratio for pool1 serum_03_bin3 is 4.94938977166044"
    ## [1] "read:cell ratio for pool1 serum_03_bin4 is 9.46308724832215"
    ## [1] "read:cell ratio for pool1 serum_04_bin1 is 4.97256897910169"
    ## [1] "read:cell ratio for pool1 serum_04_bin2 is 5.12542338350269"
    ## [1] "reads < cells for pool1 serum_04_bin3 , un-normalized (ratio 0.600732600732601 )"
    ## [1] "reads < cells for pool1 serum_04_bin4 , un-normalized (ratio 0.521739130434783 )"
    ## [1] "read:cell ratio for pool1 serum_05_bin1 is 5.51551004386608"
    ## [1] "read:cell ratio for pool1 serum_05_bin2 is 5.73855111904362"
    ## [1] "read:cell ratio for pool1 serum_05_bin3 is 1.31137724550898"
    ## [1] "read:cell ratio for pool1 serum_05_bin4 is 2.21368948247078"
    ## [1] "read:cell ratio for pool2 serum_01_bin1 is 5.48320103212859"
    ## [1] "read:cell ratio for pool2 serum_01_bin2 is 5.17768314373178"
    ## [1] "read:cell ratio for pool2 serum_01_bin3 is 4.64636022588509"
    ## [1] "read:cell ratio for pool2 serum_01_bin4 is 5.25416894471761"
    ## [1] "read:cell ratio for pool2 serum_02_bin1 is 5.70586082941746"
    ## [1] "read:cell ratio for pool2 serum_02_bin2 is 5.47646290052734"
    ## [1] "read:cell ratio for pool2 serum_02_bin3 is 3.73659774334611"
    ## [1] "read:cell ratio for pool2 serum_02_bin4 is 5.55002494009923"
    ## [1] "read:cell ratio for pool2 serum_03_bin1 is 4.99246660232726"
    ## [1] "read:cell ratio for pool2 serum_03_bin2 is 3.38613510857589"
    ## [1] "read:cell ratio for pool2 serum_03_bin3 is 3.81561091508595"
    ## [1] "reads < cells for pool2 serum_03_bin4 , un-normalized (ratio 0.390986601705238 )"
    ## [1] "read:cell ratio for pool2 serum_04_bin1 is 4.543626618841"
    ## [1] "read:cell ratio for pool2 serum_04_bin2 is 5.56513459011993"
    ## [1] "reads < cells for pool2 serum_04_bin3 , un-normalized (ratio 0.791855203619909 )"
    ## [1] "reads < cells for pool2 serum_04_bin4 , un-normalized (ratio 0.608247422680412 )"
    ## [1] "read:cell ratio for pool2 serum_05_bin1 is 4.44551094511568"
    ## [1] "read:cell ratio for pool2 serum_05_bin2 is 4.6397358415049"
    ## [1] "reads < cells for pool2 serum_05_bin3 , un-normalized (ratio 0.905781584582441 )"
    ## [1] "reads < cells for pool2 serum_05_bin4 , un-normalized (ratio 0.532051282051282 )"

``` r
#annotate each barcode in MERS library as to whether it's a homolog variant, wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations, or just the name of the target for the pan-merbeco library
dt[,variant_class:=as.character(NA)]
dt[sublibrary %in% c("lib84", "lib85") & n_codon_substitutions==0, variant_class := "wildtype"]
dt[sublibrary %in% c("lib84", "lib85") & n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[sublibrary %in% c("lib84", "lib85") & n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[sublibrary %in% c("lib84", "lib85") & n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[sublibrary %in% c("lib84", "lib85") & n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

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
    return(list(as.numeric(NA),as.numeric(NA)))
  }else{
    return( list(as.numeric((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4])), as.numeric(total)) )
  }
}
  

#iterate through titration samples, compute mean_bin and total_count for each barcode variant
#serum
for(i in 1:nrow(samples_serum)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_serum[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_serum[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_serum[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_serum[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_serum[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_serum[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}
```

## Calculate per-bc AUC metrics

We will calculate a simple AUC metric across each barcode’s titration
series. We will also include a minimum cell count that is required for a
meanbin estimate to be used in the titration fit, and a minimum number
of concentrations with determined meanbin that is required for a
titration to be reported.

``` r
#For QC and filtering, output columns giving the average number of cells that were sampled for a barcode across the 9 sample concentrations, and a value for the number of meanbin estimates that were removed for being below the # of cells cutoff
cutoff <- 2

#serum
dt[,`serum_avgcount` := mean(c(`serum_01_totalcount`,`serum_02_totalcount`,`serum_03_totalcount`,
                             `serum_04_totalcount`,`serum_05_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`serum_min_cell_filtered` := sum(c(c(`serum_01_totalcount`,`serum_02_totalcount`,`serum_03_totalcount`,
                                       `serum_04_totalcount`,`serum_05_totalcount`)<cutoff,
                                     is.na(c(`serum_01_totalcount`,`serum_02_totalcount`,`serum_03_totalcount`,
                                       `serum_04_totalcount`,`serum_05_totalcount`))),na.rm=T),by=c("library","barcode")]

#function that calculates an AUC metric across four log-spaced points and a zero-point for substraction, including an option to filter below certain thresholds for average cells across all samples, and number of samples below a cutoff of cells
fit.auc <- function(x.vals,y.vals,zero.val,count.vals,zero.count.val,min.cfu=cutoff){
  if(sum(!is.na(y.vals))==4 & !is.na(zero.val)){
    if(sum(count.vals > min.cfu) == length(count.vals) & zero.count.val > min.cfu){
      y.bg <- y.vals - zero.val
      y.bg[y.bg<0] <- 0
      auc <- sum(diff(rev(x.vals)) * (head(rev(y.bg),-1)+tail(rev(y.bg),-1)))/2 #reverse order, I supply high to low but I want low to high
      return(auc)
    }else{
      return(as.numeric(NA))
    }
  }else{
    return(as.numeric(NA))
  }
}

#fit auc to serum serum sera data for each barcode
dt[,c("serum_AUC") := fit.auc(x.vals=log10(samples_serum$conc[1:4]),
                              y.vals=c(`serum_01_meanbin`,`serum_02_meanbin`,`serum_03_meanbin`,`serum_04_meanbin`),
                      zero.val=`serum_05_meanbin`,
                      count.vals=c(`serum_01_totalcount`,`serum_02_totalcount`,`serum_03_totalcount`,`serum_04_totalcount`),
                      zero.count.val=`serum_05_totalcount`),
   by=c("library","barcode")]
```

## QC and sanity checks

We will do some QC to make sure we got good titration curves for most of
our library barcodes.

Let’s visualize the AUC binding measurements as violin plots for the
different wildtype targets, for each serum metric. This is pretty
worthless!

``` r
p1 <- ggplot(dt[!is.na(`serum_AUC`),],aes(x=variant_class,y=`serum_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("serum sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~target,nrow=1)

grid.arrange(p1,ncol=1)
```

<img src="compute_serum_AUC_files/figure-gfm/binding_distribution_vioplot_serum-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_serum.pdf",sep="")))
```

## Save barcode-level metrics

Before we collapse bcs down to final mutant/variant-level phenotypes,
integrate things like expression effects of variants, and visualize
final phenotypes, save file of per-barcode measurements.

``` r
dt[,.(library,sublibrary,barcode,target,variant_class,aa_substitutions,n_aa_substitutions,
     `serum_avgcount`,`serum_AUC`)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$sera_delta_AUC_file, row.names=F)
```

## Plot representative binding curves

Want to illustrate representative binding curves from which AUCs were
measured. Will do, for e.g. wildtype MERS versus some other
representatives \* haven’t updated this code \*

``` r
# par(mfrow=c(1,6))
# 
# #SARS-CoV-2_WH1
# #set empty plot window
# plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="mouse1.4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
# #put in faint points per-replicate
# y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_05_meanbin`]
# points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
# y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_04_totalcount` > cutoff,`mouse1-4_04_meanbin`]
# points(rep(10^-5,length(y)),y,pch=16,col="#7f7f7f02")
# y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_03_totalcount` > cutoff,`mouse1-4_03_meanbin`]
# points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
# y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_02_totalcount` > cutoff,`mouse1-4_02_meanbin`]
# points(rep(10^-3,length(y)),y,pch=16,col="#7f7f7f02")
# y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_01_totalcount` > cutoff,`mouse1-4_01_meanbin`]
# points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
# #put in black line for the average
# y1 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_05_meanbin`],na.rm=T)
# points(10^-6,y1,pch=16,col="black")
# y2 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_04_totalcount` > cutoff,`mouse1-4_04_meanbin`],na.rm=T)
# points(10^-5,y2,pch=16,col="black")
# y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_03_totalcount` > cutoff,`mouse1-4_03_meanbin`],na.rm=T)
# points(10^-4,y3,pch=16,col="black")
# y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_02_totalcount` > cutoff,`mouse1-4_02_meanbin`],na.rm=T)
# points(10^-3,y4,pch=16,col="black")
# y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib47_SARSr-wts" & `mouse1-4_01_totalcount` > cutoff,`mouse1-4_01_meanbin`],na.rm=T)
# points(10^-2,y5,pch=16,col="black")
# #connect average points with lines
# lines(c(10^-5, 10^-4),c(y2, y3),lwd=1.5,col="black",lty=2)
# lines(c(10^-4, 10^-3),c(y3, y4),lwd=1.5,col="black",lty=2)
# lines(c(10^-3, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
# abline(h=y1,lty=2,col="gray50")
# legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
# 
# #save pdf
# invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_serum.pdf",sep="")))
```

## Plot all curves

``` r
par(mfrow=c(1,6))

# for(bg in c(unique(dt$target)){
#   #mouse1-4
#   #set empty plot window
#   plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main=paste(bg,", mouse1-4"),ylab="serum binding (mean FACS bin)",xlab="serum dilution (10^-6 = 0)")
#   #put in faint points per-barcode
#   y <- dt[target==bg & variant_class=="wildtype" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_05_meanbin`]
#   points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
#   y <- dt[target==bg & variant_class=="wildtype" & `mouse1-4_04_totalcount` > cutoff,`mouse1-4_04_meanbin`]
#   points(rep(10^-5,length(y)),y,pch=16,col="#7f7f7f02")
#   y <- dt[target==bg & variant_class=="wildtype" & `mouse1-4_03_totalcount` > cutoff,`mouse1-4_03_meanbin`]
#   points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
#   y <- dt[target==bg & variant_class=="wildtype" & `mouse1-4_02_totalcount` > cutoff,`mouse1-4_02_meanbin`]
#   points(rep(10^-3,length(y)),y,pch=16,col="#7f7f7f02")
#   y <- dt[target==bg & variant_class=="wildtype" & `mouse1-4_01_totalcount` > cutoff,`mouse1-4_01_meanbin`]
#   points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#   #put in black line for the average
#   y1 <- mean(dt[target==bg & variant_class=="wildtype" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_05_meanbin`],na.rm=T)
#   points(10^-6,y1,pch=16,col="black")
#   y2 <- mean(dt[target==bg & variant_class=="wildtype" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_04_meanbin`],na.rm=T)
#   points(10^-5,y2,pch=16,col="black")
#   y3 <- mean(dt[target==bg & variant_class=="wildtype" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_03_meanbin`],na.rm=T)
#   points(10^-4,y3,pch=16,col="black")
#   y4 <- mean(dt[target==bg & variant_class=="wildtype" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_02_meanbin`],na.rm=T)
#   points(10^-3,y4,pch=16,col="black")
#   y5 <- mean(dt[target==bg & variant_class=="wildtype" & `mouse1-4_05_totalcount` > cutoff,`mouse1-4_01_meanbin`],na.rm=T)
#   points(10^-2,y5,pch=16,col="black")
#   #connect average points with lines
#   lines(c(10^-5, 10^-4),c(y2, y3),lwd=1.5,col="black",lty=2)
#   lines(c(10^-4, 10^-3),c(y3, y4),lwd=1.5,col="black",lty=2)
#   lines(c(10^-3, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
#   abline(h=y1,lty=2,col="gray50")
#   
```
