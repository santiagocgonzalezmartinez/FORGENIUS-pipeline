FORGENIUS-pipeline
================
SCGM & MW
2024-03-27

``` r
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
library(genepop)
```

    ## Warning: package 'genepop' was built under R version 4.3.3

``` r
library(graph4lg)
```

    ## Warning: package 'graph4lg' was built under R version 4.3.3

    ## Welcome to 'graph4lg' package. Let's do landscape genetics analysis with graphs

Read data and convert to different formats

``` r
data_vcf <- read.vcfR("data/Aalba_SNP_sampleFilt.vcf.gz")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7577
    ##   header_line: 7578
    ##   variant count: 66000
    ##   column count: 422
    ## Meta line 1000 read in.Meta line 2000 read in.Meta line 3000 read in.Meta line 4000 read in.Meta line 5000 read in.Meta line 6000 read in.Meta line 7000 read in.Meta line 7577 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 66000
    ##   Character matrix gt cols: 422
    ##   skip: 0
    ##   nrows: 66000
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant 15000Processed variant 16000Processed variant 17000Processed variant 18000Processed variant 19000Processed variant 20000Processed variant 21000Processed variant 22000Processed variant 23000Processed variant 24000Processed variant 25000Processed variant 26000Processed variant 27000Processed variant 28000Processed variant 29000Processed variant 30000Processed variant 31000Processed variant 32000Processed variant 33000Processed variant 34000Processed variant 35000Processed variant 36000Processed variant 37000Processed variant 38000Processed variant 39000Processed variant 40000Processed variant 41000Processed variant 42000Processed variant 43000Processed variant 44000Processed variant 45000Processed variant 46000Processed variant 47000Processed variant 48000Processed variant 49000Processed variant 50000Processed variant 51000Processed variant 52000Processed variant 53000Processed variant 54000Processed variant 55000Processed variant 56000Processed variant 57000Processed variant 58000Processed variant 59000Processed variant 60000Processed variant 61000Processed variant 62000Processed variant 63000Processed variant 64000Processed variant 65000Processed variant 66000Processed variant: 66000
    ## All variants processed

``` r
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
    ##  // 18 populations; 66,000 loci; 132,000 alleles; size: 49.1 Mb
    ## 
    ##  // Basic content
    ##    @tab:  18 x 132000 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 2-2)
    ##    @loc.fac: locus factor for the 132000 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: genind2genpop(x = data_genind)
    ## 
    ##  // Optional content
    ##    - empty -

``` r
pop_names<-popNames(data_genpop)

# genepop format
# data_genepop <- suppressWarnings((genind_to_genepop(data_genind, output = "data.frame")))
```

Compute number of polymorphic loci per population

``` r
# Estimate allele frequency in adegenet
allele_freq <- tab(data_genpop)

# Count monomorphic loci and missing loci per population
monomorphic_counts <- rowSums(allele_freq == 1)
nb_loci <- rowSums(allele_freq != "NA")/2

# Print % polymorphic loci per population
for (i in 1:length(monomorphic_counts)) {
  cat(pop_names[i], ":", 1-(monomorphic_counts[i]/nb_loci[i]), "\n")
}
```

    ## AUT00179 : 0.9339848 
    ## ITA00271 : 0.8660758 
    ## SVN00025 : 0.9252121 
    ## SVN00023 : 0.9186212 
    ## ITA00260 : 0.9337424 
    ## ROU00358 : 0.8841212 
    ## ROU00104 : 0.8807879 
    ## ROU00389 : 0.8859394 
    ## ROU00477 : 0.8843939 
    ## AUT00215 : 0.9374091 
    ## DEU00114 : 0.9281061 
    ## FRA00006 : 0.9395909 
    ## ITA00069 : 0.9244545 
    ## FRA00019 : 0.9345758 
    ## FRA00004 : 0.9188333 
    ## ESP00339 : 0.9408939 
    ## ITA00217 : 0.9175303 
    ## ITA00029 : 0.8441667

Compute genetic diversity - expected heterozygosity with genepop

``` r
#genedivFis(data_genepop, outputFile = 'genetic_diversity.txt')
```

Produce a nice report
