Merge barcode-variant tables for libraries combined into different
experimental pools
================
Tyler Starr
02/20/2024

This notebook takes different pools as defined in the config.yaml file
and makes a merged barcode-variant lookup table that pools different
libraries according to experimental specifications.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","plotly","withr","htmlwidgets","knitr")
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
    ##  [1] knitr_1.37        htmlwidgets_1.5.4 withr_2.5.0       plotly_4.10.1    
    ##  [5] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [9] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ## [13] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2  xfun_0.30         haven_2.4.3       colorspace_2.0-3 
    ##  [5] vctrs_0.5.2       generics_0.1.2    viridisLite_0.4.0 htmltools_0.5.2  
    ##  [9] utf8_1.2.2        rlang_1.0.6       pillar_1.7.0      glue_1.6.2       
    ## [13] DBI_1.1.2         dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1     
    ## [17] lifecycle_1.0.3   munsell_0.5.0     gtable_0.3.0      cellranger_1.1.0 
    ## [21] rvest_1.0.2       evaluate_0.15     tzdb_0.2.0        fastmap_1.1.0    
    ## [25] fansi_1.0.2       broom_0.7.12      Rcpp_1.0.11       backports_1.4.1  
    ## [29] scales_1.2.1      jsonlite_1.8.7    fs_1.5.2          hms_1.1.1        
    ## [33] digest_0.6.29     stringi_1.7.6     grid_4.1.3        cli_3.6.0        
    ## [37] tools_4.1.3       magrittr_2.0.2    lazyeval_0.2.2    crayon_1.5.0     
    ## [41] pkgconfig_2.0.3   ellipsis_0.3.2    xml2_1.3.3        reprex_2.0.1     
    ## [45] lubridate_1.8.0   rstudioapi_0.13   assertthat_0.2.1  rmarkdown_2.13   
    ## [49] httr_1.4.7        R6_2.5.1          compiler_4.1.3

Read in tables of per-library barcode-variant lookups, and store the
underlying “library” in a new column called “sublibrary”

``` r
#read in barcode-variant lookup tables
dt_MERS <- data.table(read.csv(file=config$codon_variant_table_file_MERS,stringsAsFactors=F)); dt_MERS[,sublibrary:=library]
dt_PDF2180 <- data.table(read.csv(file=config$codon_variant_table_file_PDF2180,stringsAsFactors=F)); dt_PDF2180[,sublibrary:=library]

dt <- rbind(dt_MERS,dt_PDF2180)

head(dt)
```

    ##    target library          barcode variant_call_support codon_substitutions
    ## 1:   MERS   lib51 AAAAAAAAACATTCGT                    1           AAC145TTG
    ## 2:   MERS   lib51 AAAAAAAAAGACTTTC                    1           GTA158AAA
    ## 3:   MERS   lib51 AAAAAAAAAGCGATAG                    4           TCC128---
    ## 4:   MERS   lib51 AAAAAAAAATTGAGGT                    2           CCT139CAA
    ## 5:   MERS   lib51 AAAAAAAACAATCCCG                    4           ACA188TAT
    ## 6:   MERS   lib51 AAAAAAAACCCATGGA                    1            CCC54TGG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:            N145L                     1                  1      lib51
    ## 2:            V158K                     1                  1      lib51
    ## 3:            S128-                     1                  1      lib51
    ## 4:            P139Q                     1                  1      lib51
    ## 5:            T188Y                     1                  1      lib51
    ## 6:             P54W                     1                  1      lib51

Make new tables for each pool (better to do ‘by pool’ than ‘by library’,
since each pool is unique but any library can in theory be part of
multiple pools). Note, this code must be manually updated when new pools
are designated in the config.yaml file.

``` r
#pool1
dt_pool1 <- dt[sublibrary %in% config$pool1]
dt_pool1[,library:="pool1"]
head(dt_pool1)
```

    ##    target library          barcode variant_call_support codon_substitutions
    ## 1:   MERS   pool1 AAAAAAAAACATTCGT                    1           AAC145TTG
    ## 2:   MERS   pool1 AAAAAAAAAGACTTTC                    1           GTA158AAA
    ## 3:   MERS   pool1 AAAAAAAAAGCGATAG                    4           TCC128---
    ## 4:   MERS   pool1 AAAAAAAAATTGAGGT                    2           CCT139CAA
    ## 5:   MERS   pool1 AAAAAAAACAATCCCG                    4           ACA188TAT
    ## 6:   MERS   pool1 AAAAAAAACCCATGGA                    1            CCC54TGG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:            N145L                     1                  1      lib51
    ## 2:            V158K                     1                  1      lib51
    ## 3:            S128-                     1                  1      lib51
    ## 4:            P139Q                     1                  1      lib51
    ## 5:            T188Y                     1                  1      lib51
    ## 6:             P54W                     1                  1      lib51

``` r
#pool2
dt_pool2 <- dt[sublibrary %in% config$pool2]
dt_pool2[,library:="pool2"]
head(dt_pool2)
```

    ##    target library          barcode variant_call_support codon_substitutions
    ## 1:   MERS   pool2 AAAAAAAAAAACTTCC                    2  CCC54GAT ATC153TTT
    ## 2:   MERS   pool2 AAAAAAAAAACAAGTT                    2            ACG36ATT
    ## 3:   MERS   pool2 AAAAAAAAACATCGTA                    2   GAT8TAT ATT104TCT
    ## 4:   MERS   pool2 AAAAAAAAACTTATGT                    1                    
    ## 5:   MERS   pool2 AAAAAAAAAGTATCTT                    1            GGC15GAA
    ## 6:   MERS   pool2 AAAAAAAAAGTCTTAA                    1            GTA82ATT
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:       P54D I153F                     2                  2      lib52
    ## 2:             T36I                     1                  1      lib52
    ## 3:        D8Y I104S                     2                  2      lib52
    ## 4:                                      0                  0      lib52
    ## 5:             G15E                     1                  1      lib52
    ## 6:             V82I                     1                  1      lib52

Eliminate barcodes that are repeated between different variants within a
single pool.

``` r
#pool1
duplicates_pool1 <- dt_pool1[duplicated(dt_pool1,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool1[,duplicate:=FALSE]
if(nrow(duplicates_pool1) > 0){
  for(i in 1:nrow(duplicates_pool1)){
    dt_pool1[library==duplicates_pool1[i,library] & barcode==duplicates_pool1[i,barcode],duplicate:=TRUE]
  }
}
dt_pool1 <- dt_pool1[duplicate==FALSE,]; dt_pool1[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool1), "repeated barcodes from pool1"))
```

    ## [1] "Removed 3 repeated barcodes from pool1"

``` r
#pool2
duplicates_pool2 <- dt_pool2[duplicated(dt_pool2,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool2[,duplicate:=FALSE]
if(nrow(duplicates_pool2) > 0){
  for(i in 1:nrow(duplicates_pool2)){
    dt_pool2[library==duplicates_pool2[i,library] & barcode==duplicates_pool2[i,barcode],duplicate:=TRUE]
  }
}
dt_pool2 <- dt_pool2[duplicate==FALSE,]; dt_pool2[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool2), "repeated barcodes from pool2"))
```

    ## [1] "Removed 2 repeated barcodes from pool2"

Merge the per-pool tables back into one aggregate data table and save.

``` r
dt_final <- rbind(dt_pool1, dt_pool2)

dt_final %>%
  write.csv(file=config$codon_variant_table_file_pools, row.names=F)
```
