FORGENIUS-pipeline
================
SCGM & MW
2024-03-30

``` r
#IMPORTANT: don't change order of libraries as some mask others!
library(vcfR)
```

    ## Warning: package 'vcfR' was built under R version 4.3.3

    ## 
    ##    *****       ***   vcfR   ***       *****
    ##    This is vcfR 1.15.0 
    ##      browseVignettes('vcfR') # Documentation
    ##      citation('vcfR') # Citation
    ##    *****       *****      *****       *****

``` r
library(adegenet)
```

    ## Warning: package 'adegenet' was built under R version 4.3.3

    ## Loading required package: ade4

    ## Warning: package 'ade4' was built under R version 4.3.3

    ## 
    ##    /// adegenet 2.1.10 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

``` r
#library(genepop)
library(graph4lg)
```

    ## Warning: package 'graph4lg' was built under R version 4.3.3

    ## Welcome to 'graph4lg' package. Let's do landscape genetics analysis with graphs

``` r
library(hierfstat)
```

    ## Warning: package 'hierfstat' was built under R version 4.3.3

    ## 
    ## Attaching package: 'hierfstat'

    ## The following objects are masked from 'package:adegenet':
    ## 
    ##     Hs, read.fstat

``` r
#library(radiator)
#library(VariantAnnotation)
#library(genomation)
#library(pegas)
#library(GenomicRanges)
```

Read data and convert to different formats

``` r
data_vcf <- read.vcfR("data/Aalba_target.vcf.gz") #change species name
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7577
    ##   header_line: 7578
    ##   variant count: 3187
    ##   column count: 422
    ## Meta line 1000 read in.Meta line 2000 read in.Meta line 3000 read in.Meta line 4000 read in.Meta line 5000 read in.Meta line 6000 read in.Meta line 7000 read in.Meta line 7577 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 3187
    ##   Character matrix gt cols: 422
    ##   skip: 0
    ##   nrows: 3187
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant: 3187
    ## All variants processed

``` r
# subset VCF files based on bed file - do it in vcftools!

# adegenet formats (genind and genpop)
data_genind <- vcfR2genind(data_vcf)
pop(data_genind) <- substr(indNames(data_genind), 1,8) #takes pop name from the first 8 digits of sample name, e.g. AUT00215
data_genpop <- genind2genpop(data_genind)
```

    ## 
    ##  Converting data from a genind to a genpop object... 
    ## 
    ## ...done.

``` r
print(data_genpop)
```

    ## /// GENPOP OBJECT /////////
    ## 
    ##  // 18 populations; 3,187 loci; 6,374 alleles; size: 2.4 Mb
    ## 
    ##  // Basic content
    ##    @tab:  18 x 6374 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 2-2)
    ##    @loc.fac: locus factor for the 6374 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: genind2genpop(x = data_genind)
    ## 
    ##  // Optional content
    ##    - empty -

``` r
pop_names<-popNames(data_genpop)

# counting individuals per population
pop_labels <- data_genind@pop
individuals_per_pop <- table(pop_labels)

# hierfstat format
data_hierfstat <- genind2hierfstat(data_genind,pop=NULL)

# bayescan format to external file
write.bayescan(dat=data_hierfstat,diploid=TRUE,fn="data_bayescan")

# genepop format to external file
#genind_to_genepop(data_genind, output = "data_genepop.txt")
```

Compute number of polymorphic loci per population

``` r
# Estimate allele frequency in adegenet
allele_freq <- tab(data_genpop)

# Count monomorphic loci and missing loci per population
monomorphic_counts <- rowSums(allele_freq == 1)
nb_loci <- rowSums(allele_freq != "NA")/2
polymoprhic_loci <- 1-(monomorphic_counts/nb_loci)
```

Compute basic stats with hierfstat

``` r
Hs <- Hs(data_hierfstat)
Ho <- Ho(data_hierfstat)
Fis <- 1-(Ho/Hs)
```

Compute average pairwise genetic differentiation with hierfstat

``` r
# Average pairwise Fst
pairwise_WCfst <- pairwise.WCfst(data_hierfstat)
ave_pairwise_fst<-rowMeans(pairwise_WCfst, na.rm = TRUE)

# Population-specific Fst from bayescan

#If you are running RStudio in Windows 10
system2("powershell", arg = c("-file", "bayescan.ps1")) #check out running parameters in ps1 file
source ("bayescan/plot_R.r")
plot_bayescan("data_bayescan_fst.txt",add_text=T,pos=0.02,size=0.7,FDR=0.05)
```

![](FORGENIUS-pipeline_files/figure-gfm/genetic%20differentiation-1.png)<!-- -->

    ## $outliers
    ##   [1]   22   23   45   70  117  123  130  149  167  171  176  198  255  260  299
    ##  [16]  307  308  310  311  328  349  359  365  400  423  433  462  507  509  516
    ##  [31]  526  540  552  577  578  618  620  630  635  648  649  653  668  671  698
    ##  [46]  700  701  702  703  715  722  729  780  782  804  841  855  858  876  877
    ##  [61]  878  900  908  910  931  952  973  983  999 1005 1008 1014 1015 1028 1043
    ##  [76] 1044 1046 1047 1048 1049 1050 1051 1110 1113 1115 1120 1160 1211 1231 1234
    ##  [91] 1250 1283 1284 1296 1303 1307 1309 1316 1341 1364 1371 1379 1392 1401 1402
    ## [106] 1417 1434 1447 1449 1455 1457 1489 1499 1513 1521 1598 1651 1670 1728 1732
    ## [121] 1734 1782 1800 1801 1809 1812 1818 1820 1836 1846 1847 1903 1904 1905 1908
    ## [136] 1927 1930 1954 2030 2033 2035 2040 2073 2087 2097 2098 2103 2165 2178 2191
    ## [151] 2211 2221 2231 2273 2274 2275 2278 2388 2392 2394 2400 2412 2457 2467 2468
    ## [166] 2473 2480 2522 2535 2538 2567 2570 2581 2582 2585 2602 2604 2611 2623 2629
    ## [181] 2723 2782 2840 2841 2842 2843 2854 2861 2862 2863 2864 2871 2885 2909 2944
    ## [196] 2945 2956 2975 2977 2978 2985 3001 3002 3006 3007 3018 3028 3044 3045 3046
    ## [211] 3047 3048 3050 3087 3096 3098 3100 3125 3129 3184
    ## 
    ## $nb_outliers
    ## [1] 220

``` r
output_bayescan<-read.table("data_bayescan.sel")
output_bayescan$logL<-NULL
pop_specific_fst<-colMeans(output_bayescan, na.rm = TRUE)

#If you are running RStudio in Linux - not tested yet!
#scan.pops <- radiator::run_bayescan(
#    data = "data_bayescan",
#    n = 5000,
#    thin = 10,
#    nbp = 20,
#    pilot = 5000,
#    burn = 50000,
#    subsample = NULL,
#    iteration.subsample = 1,
#    parallel.core = parallel::detectCores() - 1,
#    pr_odds = 10,
#    bayescan.path = "bayescan/BayeScan2.1_linux64bits"
#    )
```

Print table

``` r
table_results<-data.frame(unlist(nb_loci), unlist(individuals_per_pop), unlist(polymoprhic_loci), unlist(Hs), unlist(Fis), unlist(ave_pairwise_fst), unlist(pop_specific_fst))
names (table_results) = c("nb_loci", "GCU", "N", "polymorphic_loci", "Hs", "Fis", "ave_pairwise_fst", "pop_specific_fst")
write.csv(table_results,"table_results_Aalba_target.csv") #change species name
```
