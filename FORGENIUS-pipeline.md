FORGENIUS-pipeline, test run 20 SNPs
================
SCGM & MW
2024-03-31

Load packages

``` r
#IMPORTANT: don't change order of libraries as some mask others!
library(vcfR)
library(adegenet)
#library(genepop)
library(graph4lg)
library(hierfstat)
#library(radiator)
#library(VariantAnnotation)
#library(genomation)
#library(pegas)
#library(GenomicRanges)
library(LEA)
#library(rnaturalearth)z
library(dplyr)
```

Read data and convert to different formats

``` r
data_vcf <- read.vcfR("data/Aalba_20SNPs.vcf.gz") #change species name
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7577
    ##   header_line: 7578
    ##   variant count: 20
    ##   column count: 422
    ## Meta line 1000 read in.Meta line 2000 read in.Meta line 3000 read in.Meta line 4000 read in.Meta line 5000 read in.Meta line 6000 read in.Meta line 7000 read in.Meta line 7577 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 20
    ##   Character matrix gt cols: 422
    ##   skip: 0
    ##   nrows: 20
    ##   row_num: 0
    ## Processed variant: 20
    ## All variants processed

``` r
# subset VCF files based on bed file - do it in vcftools!

# adegenet formats (genind and genpop)
data_genind <- vcfR2genind(data_vcf)
pop(data_genind) <- substr(indNames(data_genind), 1,8) #takes pop name from the first 8 digits of sample name, e.g. AUT00215
individual_names<-indNames(data_genind) #extract individual names
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
    ##  // 18 populations; 20 loci; 40 alleles; size: 19.7 Kb
    ## 
    ##  // Basic content
    ##    @tab:  18 x 40 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 2-2)
    ##    @loc.fac: locus factor for the 40 columns of @tab
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
write.bayescan(dat=data_hierfstat,diploid=TRUE,fn="data_bayescan") #writting this file can take about 6 hours for 90k SNPs; better to store these files once generated!

# structure format to external file
write.struct(dat=data_hierfstat,ilab=individual_names,pop=NULL,MARKERNAMES=FALSE,MISSING=-9,fname="data_structure") #writting this file can take about xxx hours for 90k SNPs; better to store these files once generated!

# geno format to external file
struct2geno (input.file="data_structure", ploidy=2, FORMAT=2, extra.row=0, extra.column=1)
```

    ## Input file in the STRUCTURE format. The genotypic matrix has 413 individuals and 20 markers. 
    ## The number of extra rows is 0 and the number of extra columns is 1 .
    ## Missing alleles are encoded as -9 , converted as 9.
    ## Output files: data_structure.geno  .lfmm.

``` r
# genepop format to external file
#genind_to_genepop(data_genind, output = "data_genepop.txt")
```

Run admixture analysis with LEA

``` r
# run SNMF
data_snmf <- snmf(input.file = "data_structure.geno",
                   K = 1:10, #indicate max number of K
                   entropy = TRUE,
                   repetitions = 5,
                   project = "new",
                   CPU = 10, #indicate the number of CPUs to use
                   alpha = 100
                  )
```

    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] 842197358
    ## [1] "*************************************"
    ## [1] "*          create.dataset            *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)                 413
    ##         -L (number of loci)                        20
    ##         -s (seed random init)                      842197358
    ##         -r (percentage of masked data)             0.05
    ##         -x (genotype file in .geno format)         C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -o (output file in .geno format)           C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ## 
    ##  Write genotype file with masked data, C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 1  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          1
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run1/data_structure_r1.1.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run1/data_structure_r1.1.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ## 
    ## Least-square error: 3364.135597
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run1/data_structure_r1.1.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run1/data_structure_r1.1.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      1
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run1/data_structure_r1.1.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run1/data_structure_r1.1.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.677097
    ## Cross-Entropy (masked data):  0.634293
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 2  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          2
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run1/data_structure_r1.2.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run1/data_structure_r1.2.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [========]
    ## Number of iterations: 21
    ## 
    ## Least-square error: 2936.240981
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run1/data_structure_r1.2.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run1/data_structure_r1.2.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      2
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run1/data_structure_r1.2.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run1/data_structure_r1.2.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.583018
    ## Cross-Entropy (masked data):  0.523973
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 3  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          3
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run1/data_structure_r1.3.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run1/data_structure_r1.3.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2689.119318
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run1/data_structure_r1.3.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run1/data_structure_r1.3.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      3
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run1/data_structure_r1.3.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run1/data_structure_r1.3.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.529672
    ## Cross-Entropy (masked data):  0.516561
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 4  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          4
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run1/data_structure_r1.4.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run1/data_structure_r1.4.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2558.900755
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run1/data_structure_r1.4.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run1/data_structure_r1.4.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      4
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run1/data_structure_r1.4.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run1/data_structure_r1.4.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.505152
    ## Cross-Entropy (masked data):  0.525259
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 5  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          5
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run1/data_structure_r1.5.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run1/data_structure_r1.5.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2427.049662
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run1/data_structure_r1.5.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run1/data_structure_r1.5.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      5
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run1/data_structure_r1.5.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run1/data_structure_r1.5.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.47274
    ## Cross-Entropy (masked data):  0.536019
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 6  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          6
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run1/data_structure_r1.6.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run1/data_structure_r1.6.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [============================]
    ## Number of iterations: 75
    ## 
    ## Least-square error: 2377.987568
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run1/data_structure_r1.6.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run1/data_structure_r1.6.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      6
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run1/data_structure_r1.6.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run1/data_structure_r1.6.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.458725
    ## Cross-Entropy (masked data):  0.519415
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 7  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          7
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run1/data_structure_r1.7.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run1/data_structure_r1.7.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=================]
    ## Number of iterations: 45
    ## 
    ## Least-square error: 2265.306926
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run1/data_structure_r1.7.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run1/data_structure_r1.7.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      7
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run1/data_structure_r1.7.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run1/data_structure_r1.7.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.439048
    ## Cross-Entropy (masked data):  0.491224
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 8  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          8
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run1/data_structure_r1.8.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run1/data_structure_r1.8.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [===============]
    ## Number of iterations: 39
    ## 
    ## Least-square error: 2311.855205
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run1/data_structure_r1.8.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run1/data_structure_r1.8.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      8
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run1/data_structure_r1.8.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run1/data_structure_r1.8.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.444871
    ## Cross-Entropy (masked data):  0.633399
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 9  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          9
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run1/data_structure_r1.9.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run1/data_structure_r1.9.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=======]
    ## Number of iterations: 20
    ## 
    ## Least-square error: 2196.299506
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run1/data_structure_r1.9.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run1/data_structure_r1.9.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      9
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run1/data_structure_r1.9.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run1/data_structure_r1.9.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.41706
    ## Cross-Entropy (masked data):  0.601451
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 10  repetition 1      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          10
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run1/data_structure_r1.10.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run1/data_structure_r1.10.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  842197358
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [===========================================]
    ## Number of iterations: 115
    ## 
    ## Least-square error: 2163.615891
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run1/data_structure_r1.10.Q:       OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run1/data_structure_r1.10.G:    OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      10
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run1/data_structure_r1.10.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run1/data_structure_r1.10.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.409355
    ## Cross-Entropy (masked data):  0.638742
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] 1679156594
    ## [1] "*************************************"
    ## [1] "*          create.dataset            *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)                 413
    ##         -L (number of loci)                        20
    ##         -s (seed random init)                      1679156594
    ##         -r (percentage of masked data)             0.05
    ##         -x (genotype file in .geno format)         C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -o (output file in .geno format)           C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ## 
    ##  Write genotype file with masked data, C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 1  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          1
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run2/data_structure_r2.1.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run2/data_structure_r2.1.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ## 
    ## Least-square error: 3372.290561
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run2/data_structure_r2.1.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run2/data_structure_r2.1.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      1
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run2/data_structure_r2.1.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run2/data_structure_r2.1.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.674249
    ## Cross-Entropy (masked data):  0.687993
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 2  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          2
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run2/data_structure_r2.2.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run2/data_structure_r2.2.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2898.316959
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run2/data_structure_r2.2.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run2/data_structure_r2.2.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      2
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run2/data_structure_r2.2.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run2/data_structure_r2.2.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.5735
    ## Cross-Entropy (masked data):  0.620373
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 3  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          3
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run2/data_structure_r2.3.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run2/data_structure_r2.3.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=======]
    ## Number of iterations: 20
    ## 
    ## Least-square error: 2610.521134
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run2/data_structure_r2.3.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run2/data_structure_r2.3.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      3
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run2/data_structure_r2.3.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run2/data_structure_r2.3.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.520623
    ## Cross-Entropy (masked data):  0.603099
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 4  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          4
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run2/data_structure_r2.4.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run2/data_structure_r2.4.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [============]
    ## Number of iterations: 33
    ## 
    ## Least-square error: 2563.797406
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run2/data_structure_r2.4.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run2/data_structure_r2.4.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      4
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run2/data_structure_r2.4.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run2/data_structure_r2.4.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.503521
    ## Cross-Entropy (masked data):  0.598619
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 5  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          5
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run2/data_structure_r2.5.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run2/data_structure_r2.5.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==============]
    ## Number of iterations: 38
    ## 
    ## Least-square error: 2525.528910
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run2/data_structure_r2.5.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run2/data_structure_r2.5.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      5
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run2/data_structure_r2.5.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run2/data_structure_r2.5.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.491977
    ## Cross-Entropy (masked data):  0.563826
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 6  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          6
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run2/data_structure_r2.6.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run2/data_structure_r2.6.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=======]
    ## Number of iterations: 18
    ## 
    ## Least-square error: 2287.572055
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run2/data_structure_r2.6.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run2/data_structure_r2.6.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      6
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run2/data_structure_r2.6.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run2/data_structure_r2.6.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.439225
    ## Cross-Entropy (masked data):  0.645944
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 7  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          7
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run2/data_structure_r2.7.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run2/data_structure_r2.7.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 25
    ## 
    ## Least-square error: 2210.496393
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run2/data_structure_r2.7.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run2/data_structure_r2.7.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      7
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run2/data_structure_r2.7.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run2/data_structure_r2.7.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.427691
    ## Cross-Entropy (masked data):  0.610556
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 8  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          8
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run2/data_structure_r2.8.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run2/data_structure_r2.8.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [===========]
    ## Number of iterations: 29
    ## 
    ## Least-square error: 2173.331010
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run2/data_structure_r2.8.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run2/data_structure_r2.8.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      8
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run2/data_structure_r2.8.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run2/data_structure_r2.8.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.421371
    ## Cross-Entropy (masked data):  0.73173
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 9  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          9
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run2/data_structure_r2.9.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run2/data_structure_r2.9.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==========]
    ## Number of iterations: 28
    ## 
    ## Least-square error: 2237.402117
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run2/data_structure_r2.9.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run2/data_structure_r2.9.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      9
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run2/data_structure_r2.9.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run2/data_structure_r2.9.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.419161
    ## Cross-Entropy (masked data):  0.682586
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 10  repetition 2      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          10
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run2/data_structure_r2.10.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run2/data_structure_r2.10.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1679156594
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 23
    ## 
    ## Least-square error: 2244.693532
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run2/data_structure_r2.10.Q:       OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run2/data_structure_r2.10.G:    OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      10
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run2/data_structure_r2.10.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run2/data_structure_r2.10.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.426771
    ## Cross-Entropy (masked data):  0.703344
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] 1677343043
    ## [1] "*************************************"
    ## [1] "*          create.dataset            *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)                 413
    ##         -L (number of loci)                        20
    ##         -s (seed random init)                      1677343043
    ##         -r (percentage of masked data)             0.05
    ##         -x (genotype file in .geno format)         C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -o (output file in .geno format)           C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ## 
    ##  Write genotype file with masked data, C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 1  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          1
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run3/data_structure_r3.1.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run3/data_structure_r3.1.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ## 
    ## Least-square error: 3352.251820
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run3/data_structure_r3.1.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run3/data_structure_r3.1.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      1
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run3/data_structure_r3.1.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run3/data_structure_r3.1.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.671765
    ## Cross-Entropy (masked data):  0.744507
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 2  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          2
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run3/data_structure_r3.2.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run3/data_structure_r3.2.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2904.378121
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run3/data_structure_r3.2.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run3/data_structure_r3.2.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      2
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run3/data_structure_r3.2.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run3/data_structure_r3.2.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.576942
    ## Cross-Entropy (masked data):  0.678082
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 3  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          3
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run3/data_structure_r3.3.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run3/data_structure_r3.3.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2638.649411
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run3/data_structure_r3.3.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run3/data_structure_r3.3.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      3
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run3/data_structure_r3.3.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run3/data_structure_r3.3.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.520422
    ## Cross-Entropy (masked data):  0.615716
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 4  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          4
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run3/data_structure_r3.4.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run3/data_structure_r3.4.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=======]
    ## Number of iterations: 18
    ## 
    ## Least-square error: 2560.167167
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run3/data_structure_r3.4.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run3/data_structure_r3.4.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      4
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run3/data_structure_r3.4.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run3/data_structure_r3.4.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.500577
    ## Cross-Entropy (masked data):  0.628617
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 5  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          5
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run3/data_structure_r3.5.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run3/data_structure_r3.5.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 23
    ## 
    ## Least-square error: 2353.801863
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run3/data_structure_r3.5.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run3/data_structure_r3.5.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      5
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run3/data_structure_r3.5.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run3/data_structure_r3.5.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.459681
    ## Cross-Entropy (masked data):  0.696075
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 6  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          6
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run3/data_structure_r3.6.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run3/data_structure_r3.6.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [========]
    ## Number of iterations: 21
    ## 
    ## Least-square error: 2367.270347
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run3/data_structure_r3.6.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run3/data_structure_r3.6.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      6
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run3/data_structure_r3.6.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run3/data_structure_r3.6.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.456115
    ## Cross-Entropy (masked data):  0.645343
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 7  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          7
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run3/data_structure_r3.7.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run3/data_structure_r3.7.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=============]
    ## Number of iterations: 34
    ## 
    ## Least-square error: 2238.835331
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run3/data_structure_r3.7.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run3/data_structure_r3.7.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      7
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run3/data_structure_r3.7.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run3/data_structure_r3.7.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.427595
    ## Cross-Entropy (masked data):  0.628173
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 8  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          8
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run3/data_structure_r3.8.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run3/data_structure_r3.8.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==========]
    ## Number of iterations: 27
    ## 
    ## Least-square error: 2218.500068
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run3/data_structure_r3.8.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run3/data_structure_r3.8.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      8
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run3/data_structure_r3.8.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run3/data_structure_r3.8.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.424783
    ## Cross-Entropy (masked data):  0.768654
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 9  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          9
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run3/data_structure_r3.9.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run3/data_structure_r3.9.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 25
    ## 
    ## Least-square error: 2117.303687
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run3/data_structure_r3.9.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run3/data_structure_r3.9.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      9
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run3/data_structure_r3.9.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run3/data_structure_r3.9.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.407749
    ## Cross-Entropy (masked data):  0.764352
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 10  repetition 3      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          10
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run3/data_structure_r3.10.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run3/data_structure_r3.10.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  1677343043
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==========]
    ## Number of iterations: 28
    ## 
    ## Least-square error: 2153.510037
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run3/data_structure_r3.10.Q:       OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run3/data_structure_r3.10.G:    OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      10
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run3/data_structure_r3.10.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run3/data_structure_r3.10.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.415239
    ## Cross-Entropy (masked data):  0.71189
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] 940919775
    ## [1] "*************************************"
    ## [1] "*          create.dataset            *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)                 413
    ##         -L (number of loci)                        20
    ##         -s (seed random init)                      940919775
    ##         -r (percentage of masked data)             0.05
    ##         -x (genotype file in .geno format)         C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -o (output file in .geno format)           C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ## 
    ##  Write genotype file with masked data, C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 1  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          1
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run4/data_structure_r4.1.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run4/data_structure_r4.1.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ## 
    ## Least-square error: 3376.130755
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run4/data_structure_r4.1.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run4/data_structure_r4.1.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      1
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run4/data_structure_r4.1.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run4/data_structure_r4.1.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.673313
    ## Cross-Entropy (masked data):  0.707482
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 2  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          2
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run4/data_structure_r4.2.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run4/data_structure_r4.2.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2923.627262
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run4/data_structure_r4.2.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run4/data_structure_r4.2.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      2
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run4/data_structure_r4.2.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run4/data_structure_r4.2.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.57966
    ## Cross-Entropy (masked data):  0.630144
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 3  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          3
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run4/data_structure_r4.3.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run4/data_structure_r4.3.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==========]
    ## Number of iterations: 28
    ## 
    ## Least-square error: 2689.769926
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run4/data_structure_r4.3.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run4/data_structure_r4.3.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      3
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run4/data_structure_r4.3.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run4/data_structure_r4.3.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.526898
    ## Cross-Entropy (masked data):  0.609011
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 4  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          4
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run4/data_structure_r4.4.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run4/data_structure_r4.4.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2563.981741
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run4/data_structure_r4.4.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run4/data_structure_r4.4.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      4
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run4/data_structure_r4.4.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run4/data_structure_r4.4.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.500755
    ## Cross-Entropy (masked data):  0.574162
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 5  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          5
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run4/data_structure_r4.5.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run4/data_structure_r4.5.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 23
    ## 
    ## Least-square error: 2474.684241
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run4/data_structure_r4.5.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run4/data_structure_r4.5.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      5
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run4/data_structure_r4.5.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run4/data_structure_r4.5.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.480091
    ## Cross-Entropy (masked data):  0.632101
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 6  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          6
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run4/data_structure_r4.6.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run4/data_structure_r4.6.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 24
    ## 
    ## Least-square error: 2298.113232
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run4/data_structure_r4.6.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run4/data_structure_r4.6.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      6
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run4/data_structure_r4.6.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run4/data_structure_r4.6.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.443376
    ## Cross-Entropy (masked data):  0.621337
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 7  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          7
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run4/data_structure_r4.7.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run4/data_structure_r4.7.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==========]
    ## Number of iterations: 28
    ## 
    ## Least-square error: 2255.444285
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run4/data_structure_r4.7.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run4/data_structure_r4.7.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      7
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run4/data_structure_r4.7.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run4/data_structure_r4.7.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.435388
    ## Cross-Entropy (masked data):  0.634454
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 8  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          8
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run4/data_structure_r4.8.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run4/data_structure_r4.8.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [===============]
    ## Number of iterations: 40
    ## 
    ## Least-square error: 2218.056050
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run4/data_structure_r4.8.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run4/data_structure_r4.8.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      8
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run4/data_structure_r4.8.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run4/data_structure_r4.8.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.427334
    ## Cross-Entropy (masked data):  0.657444
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 9  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          9
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run4/data_structure_r4.9.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run4/data_structure_r4.9.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [===========]
    ## Number of iterations: 30
    ## 
    ## Least-square error: 2217.703269
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run4/data_structure_r4.9.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run4/data_structure_r4.9.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      9
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run4/data_structure_r4.9.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run4/data_structure_r4.9.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.42048
    ## Cross-Entropy (masked data):  0.621633
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 10  repetition 4      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          10
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run4/data_structure_r4.10.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run4/data_structure_r4.10.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  940919775
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 23
    ## 
    ## Least-square error: 2164.757399
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run4/data_structure_r4.10.Q:       OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run4/data_structure_r4.10.G:    OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      10
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run4/data_structure_r4.10.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run4/data_structure_r4.10.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.408237
    ## Cross-Entropy (masked data):  0.708516
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] 38606813
    ## [1] "*************************************"
    ## [1] "*          create.dataset            *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)                 413
    ##         -L (number of loci)                        20
    ##         -s (seed random init)                      38606813
    ##         -r (percentage of masked data)             0.05
    ##         -x (genotype file in .geno format)         C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -o (output file in .geno format)           C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ## 
    ##  Write genotype file with masked data, C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 1  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          1
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run5/data_structure_r5.1.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run5/data_structure_r5.1.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ## 
    ## Least-square error: 3365.075065
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run5/data_structure_r5.1.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run5/data_structure_r5.1.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      1
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run5/data_structure_r5.1.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K1/run5/data_structure_r5.1.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.673771
    ## Cross-Entropy (masked data):  0.698163
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 2  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          2
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run5/data_structure_r5.2.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run5/data_structure_r5.2.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2912.573913
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run5/data_structure_r5.2.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run5/data_structure_r5.2.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      2
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run5/data_structure_r5.2.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K2/run5/data_structure_r5.2.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.57957
    ## Cross-Entropy (masked data):  0.621984
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 3  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          3
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run5/data_structure_r5.3.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run5/data_structure_r5.3.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2659.751783
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run5/data_structure_r5.3.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run5/data_structure_r5.3.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      3
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run5/data_structure_r5.3.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K3/run5/data_structure_r5.3.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.525515
    ## Cross-Entropy (masked data):  0.622782
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 4  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          4
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run5/data_structure_r5.4.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run5/data_structure_r5.4.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [======]
    ## Number of iterations: 16
    ## 
    ## Least-square error: 2510.924456
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run5/data_structure_r5.4.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run5/data_structure_r5.4.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      4
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run5/data_structure_r5.4.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K4/run5/data_structure_r5.4.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.488164
    ## Cross-Entropy (masked data):  0.590997
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 5  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          5
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run5/data_structure_r5.5.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run5/data_structure_r5.5.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=========]
    ## Number of iterations: 24
    ## 
    ## Least-square error: 2436.591077
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run5/data_structure_r5.5.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run5/data_structure_r5.5.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      5
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run5/data_structure_r5.5.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K5/run5/data_structure_r5.5.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.468014
    ## Cross-Entropy (masked data):  0.542277
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 6  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          6
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run5/data_structure_r5.6.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run5/data_structure_r5.6.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [========]
    ## Number of iterations: 21
    ## 
    ## Least-square error: 2346.588792
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run5/data_structure_r5.6.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run5/data_structure_r5.6.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      6
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run5/data_structure_r5.6.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K6/run5/data_structure_r5.6.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.455443
    ## Cross-Entropy (masked data):  0.541423
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 7  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          7
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run5/data_structure_r5.7.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run5/data_structure_r5.7.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [================]
    ## Number of iterations: 44
    ## 
    ## Least-square error: 2303.088874
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run5/data_structure_r5.7.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run5/data_structure_r5.7.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      7
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run5/data_structure_r5.7.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K7/run5/data_structure_r5.7.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.440064
    ## Cross-Entropy (masked data):  0.530408
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 8  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          8
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run5/data_structure_r5.8.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run5/data_structure_r5.8.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [===========]
    ## Number of iterations: 29
    ## 
    ## Least-square error: 2193.825810
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run5/data_structure_r5.8.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run5/data_structure_r5.8.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      8
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run5/data_structure_r5.8.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K8/run5/data_structure_r5.8.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.424931
    ## Cross-Entropy (masked data):  0.524056
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 9  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          9
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run5/data_structure_r5.9.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run5/data_structure_r5.9.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [==========]
    ## Number of iterations: 28
    ## 
    ## Least-square error: 2176.904248
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run5/data_structure_r5.9.Q:     OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run5/data_structure_r5.9.G:  OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      9
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run5/data_structure_r5.9.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K9/run5/data_structure_r5.9.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.41267
    ## Cross-Entropy (masked data):  0.581007
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")
    ## 
    ## [1] "*************************************"
    ## [1] "* sNMF K = 10  repetition 5      *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)             413
    ##         -L (number of loci)                    20
    ##         -K (number of ancestral pops)          10
    ##         -x (input file)                        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         -q (individual admixture file)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run5/data_structure_r5.10.Q
    ##         -g (ancestral frequencies file)        C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run5/data_structure_r5.10.G
    ##         -i (number max of iterations)          200
    ##         -a (regularization parameter)          100
    ##         -s (seed random init)                  38606813
    ##         -e (tolerance error)                   1E-05
    ##         -p (number of processes)               1
    ##         - diploid
    ## 
    ## Read genotype file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno:      OK.
    ## 
    ## 
    ## Main algorithm:
    ##  [                                                                           ]
    ##  [=============]
    ## Number of iterations: 36
    ## 
    ## Least-square error: 2135.876153
    ## Write individual ancestry coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run5/data_structure_r5.10.Q:       OK.
    ## Write ancestral allele frequency coefficient file C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run5/data_structure_r5.10.G:    OK.
    ## 
    ## [1] "*************************************"
    ## [1] "*    cross-entropy estimation       *"
    ## [1] "*************************************"
    ## summary of the options:
    ## 
    ##         -n (number of individuals)         413
    ##         -L (number of loci)                20
    ##         -K (number of ancestral pops)      10
    ##         -x (genotype file)                 C:\Users\Santi\Documents\GitHub\FORGENIUS-pipeline\data_structure.geno
    ##         -q (individual admixture)          C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run5/data_structure_r5.10.Q
    ##         -g (ancestral frequencies)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/K10/run5/data_structure_r5.10.G
    ##         -i (with masked genotypes)         C:/Users/Santi/Documents/GitHub/FORGENIUS-pipeline/data_structure.snmf/masked/data_structure_I.geno
    ##         - diploid
    ## 
    ## Cross-Entropy (all data):     0.413357
    ## Cross-Entropy (masked data):  0.634473
    ## The project is saved into :
    ##  data_structure.snmfProject 
    ## 
    ## To load the project, use:
    ##  project = load.snmfProject("data_structure.snmfProject")
    ## 
    ## To remove the project, use:
    ##  remove.snmfProject("data_structure.snmfProject")

Estimate admixture coefficients

``` r
# processing sNMF results
plot(data_snmf, cex = 1.2, col = "lightblue", pch = 19)
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-1.png)<!-- -->

``` r
best_K = 4  #select best K based of graph
ce <-  cross.entropy(data_snmf, K = best_K)
best_run <- which.min(ce)

# plot barchart for best run of best K
barchart (data_snmf, best_K, best_run, sort.by.Q = T, col = rainbow(best_K),border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients") #use "sort.by.Q = F, names.arg= individual_names, cex.names = 0.4, las=2, " to get ind labels
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-2.png)<!-- -->

    ## $order
    ##   [1]  24  30  32  34  54  59  60  62  63 114 119 130 136 138 142 154 164 165
    ##  [19] 178 180 188 189 195 196 197 198 202 214 240 241 249 251 252 254 255 257
    ##  [37] 258 260 263 264 265 266 286 293 301 302 307 311 313 316 322 323 326 327
    ##  [55] 329 330 332 336 338 341 342 344 348 351 354 360 364 365 369 376 377 380
    ##  [73] 385 387 391 394 396 397 399 400 403 405 406 410 411 413  14  16  79  83
    ##  [91] 143 206 208 211 212 269 273 274 275 279 281 288 305 310 315 318 320 324
    ## [109] 325 350 355 357 370 381 384   4   5  10  13  22  27  28  33  40  46  47
    ## [127]  52  61  64  68  71  73  76  86  99 100 101 102 103 105 122 132 167 170
    ## [145] 172 174 176 181 182 187 191 192 193 207 209 213 218 229 230 271 276 278
    ## [163] 280 282 285 287 291 292 296 319 331 340 363 366 371 372 373 382 392 401
    ## [181] 404 412   1   2   3   6   7   8   9  11  12  15  17  18  19  20  21  23
    ## [199]  25  26  29  31  35  36  37  38  39  41  42  43  44  45  48  49  50  51
    ## [217]  53  55  56  57  58  65  66  67  69  70  72  74  75  77  78  80  81  82
    ## [235]  84  85  87  88  89  90  91  92  93  94  95  96  97  98 104 106 107 108
    ## [253] 109 110 111 112 113 115 116 117 118 120 121 123 124 125 126 127 128 129
    ## [271] 131 133 134 135 137 139 140 141 144 145 146 147 148 149 150 151 152 153
    ## [289] 155 156 157 158 159 160 161 162 163 166 168 169 171 173 175 177 179 183
    ## [307] 184 185 186 190 194 199 200 201 203 204 205 210 215 216 217 219 220 221
    ## [325] 222 223 224 225 226 227 228 231 232 233 234 235 236 237 238 239 242 243
    ## [343] 244 245 246 247 248 250 253 256 259 261 262 267 268 270 272 277 283 284
    ## [361] 289 290 294 295 297 298 299 300 303 304 306 308 309 312 314 317 321 328
    ## [379] 333 334 335 337 339 343 345 346 347 349 352 353 356 358 359 361 362 367
    ## [397] 368 374 375 378 379 383 386 388 389 390 393 395 398 402 407 408 409

``` r
# estimate ancestry coefficients by pop
q_mat <- LEA::Q(data_snmf, K = best_K, run = best_run) 
colnames(q_mat) <- paste0("Q", 1:best_K)
q_matrix <- data.frame(unlist(q_mat))
q_matrix <- q_matrix %>%
  mutate(individual_names, pop_labels)
q_by_pop <- aggregate(q_matrix[, 1:best_K], list(q_matrix$pop_labels), mean)
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

Compute genetic distinctness

``` r
# Average pairwise Fst
pairwise_WCfst <- pairwise.WCfst(data_hierfstat) #this is very slow, took about 16 hours for 15 pops and 90k SNPs in a laptop
ave_pairwise_fst<-rowMeans(pairwise_WCfst, na.rm = TRUE) 

# Population-specific Fst from bayescan

#If you are running RStudio in Windows 10
system2("powershell", arg = c("-file", "bayescan.ps1")) #check out running parameters in ps1 file; bayescan is very slow, took about xxx hours for 15 pops and 90k SNPs using 12 CPUs and standard parameters
source ("bayescan/plot_R.r")
plot_bayescan("bayescan_results/data_bayescan_fst.txt",add_text=T,pos=0.02,size=0.7,FDR=0.05)
```

![](FORGENIUS-pipeline_files/figure-gfm/genetic%20differentiation-1.png)<!-- -->

    ## $outliers
    ## integer(0)
    ## 
    ## $nb_outliers
    ## [1] 0

``` r
output_bayescan<-read.table("bayescan_results/data_bayescan.sel")
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
#    pr_odds = 1000,
#    bayescan.path = "bayescan/BayeScan2.1_linux64bits"
#    )
```

Print table

``` r
# prepare dataframe with basic stats
table_basic_stats<-data.frame(unlist(individuals_per_pop), unlist(nb_loci), unlist(polymoprhic_loci), unlist(Hs), unlist(Fis), unlist(ave_pairwise_fst), unlist(pop_specific_fst), row.names = NULL)
names (table_basic_stats) = c("GCU", "N", "nb_loci", "polymorphic_loci", "Hs", "Fis", "ave_pairwise_fst", "pop_specific_fst")

# merge dataframes in a single table
table_results <- merge(table_basic_stats, q_by_pop, by.x = "GCU", by.y = "Group.1") # merge dataframes by GCU name, it orders the new dataframe alphabetically
write.csv(table_results,"table_results_Aalba_20SNPs.csv") #change species name
```
