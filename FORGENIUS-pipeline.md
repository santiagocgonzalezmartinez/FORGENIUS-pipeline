FORGENIUS-pipeline, Phalepensis
================
SCGM & MW
2024-04-01

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
library(rnaturalearth)
library(dplyr)
library(mapplots)
library(shapefiles)
```

Read data and convert to different formats

``` r
data_vcf <- read.vcfR("data/Phalepensis_SNP_sampleFilt.vcf.gz") #change species name
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7434
    ##   header_line: 7435
    ##   variant count: 14053
    ##   column count: 352
    ## Meta line 1000 read in.Meta line 2000 read in.Meta line 3000 read in.Meta line 4000 read in.Meta line 5000 read in.Meta line 6000 read in.Meta line 7000 read in.Meta line 7434 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 14053
    ##   Character matrix gt cols: 352
    ##   skip: 0
    ##   nrows: 14053
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant: 14053
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
    ##  // 15 populations; 14,053 loci; 28,106 alleles; size: 9.5 Mb
    ## 
    ##  // Basic content
    ##    @tab:  15 x 28106 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 2-2)
    ##    @loc.fac: locus factor for the 28106 columns of @tab
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
write.struct(dat=data_hierfstat,ilab=individual_names,pop=NULL,MARKERNAMES=FALSE,MISSING=-9,fname="data_structure") 

# geno format to external file
struct2geno (input.file="data_structure", ploidy=2, FORMAT=2, extra.row=0, extra.column=1) #writting this file can also take quite long
```

    ## Input file in the STRUCTURE format. The genotypic matrix has 343 individuals and 14053 markers. 
    ## The number of extra rows is 0 and the number of extra columns is 1 .
    ## Missing alleles are encoded as -9 , converted as 9.

    ## Warning in struct2geno(input.file = "data_structure", ploidy = 2, FORMAT = 2, : Monomorphic alleles generated during conversion were removed.

    ## Output files: data_structure.geno  .lfmm.

``` r
# genepop format to external file
#genind_to_genepop(data_genind, output = "data_genepop.txt")
```

Run admixture analysis with LEA

Estimate admixture coefficients

``` r
# processing sNMF results
plot(data_snmf, cex = 1.2, col = "lightblue", pch = 19)
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-1.png)<!-- -->

``` r
best_K = 4  #select best K based on entropy graph
ce <-  cross.entropy(data_snmf, K = best_K)
best_run <- which.min(ce)

# plot barchart for best run of best K
barchart (data_snmf, best_K, best_run, sort.by.Q = T, col = rainbow(best_K),border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients") #use "sort.by.Q = F, names.arg= individual_names, cex.names = 0.4, las=2, " to get ind labels in original order
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-2.png)<!-- -->

    ## $order
    ##   [1] 126 127 128 129 130 132 133 134 136 137 138 139 140 141 142 143 144 146
    ##  [19] 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164
    ##  [37] 165 167 168 169 170 171 172 173 174 175 176   3   4   5   6   7   8   9
    ##  [55]  10  11  13  14  15  16  17  18  19  20  21  22  24  25  26  27  28  29
    ##  [73]  30  31  33  34  35  36  37  38  39  40  42  43  44  45  46  47  48  49
    ##  [91]  50  51  53  54  55  56  57  58  59  60  61  63  64  65  66  67  68  69
    ## [109]  70  71  72  74  75  76  77  78  79  80  85  92 131 166 183 184 185 186
    ## [127] 187 188 189 190 191 192 193 194 195 196 198 199 200 201 202 203 204 205
    ## [145] 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    ## [163] 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241
    ## [181] 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259
    ## [199] 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277
    ## [217] 278 279 280 281 282 283 284 285 286 287 289 290 291 292 293 294 295 296
    ## [235] 326 327 328 330 331 332 333 334 335 336 337 338 339 340 341 342 343  81
    ## [253]  82  84  86  87  88  89  90  91  93  94  95  96  97  98 100 101 104 105
    ## [271] 106 107 108 109 110 111 112 113 115 116 117 118 119 120 121 122 123 124
    ## [289] 177 178 179 180 181 298   1   2  12  23  32  41  52  62  73  83  99 102
    ## [307] 103 114 125 135 145 182 197 288 297 299 300 301 302 303 304 305 306 307
    ## [325] 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325
    ## [343] 329

``` r
# estimate ancestry coefficients by pop
q_mat <- LEA::Q(data_snmf, K = best_K, run = best_run) 
colnames(q_mat) <- paste0("Q", 1:best_K)
q_matrix <- data.frame(unlist(q_mat))
q_matrix <- q_matrix %>%
  mutate(individual_names, pop_labels)
q_by_pop <- aggregate(q_matrix[, 1:best_K], list(q_matrix$pop_labels), mean)

# plot map for best K
coord = read.table("data/Phalepensis_coord.txt", header = T) #change species name; GCU, longitude and latitude in this order
Npop = length(unique(pop_names))
qmatrix2plot <- merge(coord, q_by_pop, by.x = "GCU", by.y = "Group.1")
xlim <- c(-10,30) #adjust if needed
ylim <- c(34,50)  #adjust if needed
basemap(xlim, ylim) 
shape <- read.shapefile("shapefiles/chorological_maps_dataset/Pinus_halepensis_plg_clip") #change species name
draw.shape(shape, col = "grey90")
for (i in 1:Npop){
  z_var=NULL
  for (j in 4:(best_K+3)){z_var <- append(z_var, qmatrix2plot[i,j])}
  add.pie(z = z_var, x = qmatrix2plot[i,2], y = qmatrix2plot[i,3], labels =   "", col = rainbow(best_K))
 }
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-3.png)<!-- -->

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
system2("powershell", arg = c("-file", "bayescan.ps1")) #check out running parameters in ps1 file; bayescan is very slow, took about xfrom 9 pmxx hours for 15 pops and 14k SNPs using 12 CPUs and standard parameters
source ("bayescan/plot_R.r")
plot_bayescan("bayescan_results/data_bayescan_fst.txt",add_text=T,pos=0.02,size=0.7,FDR=0.05)
```

![](FORGENIUS-pipeline_files/figure-gfm/genetic%20differentiation-1.png)<!-- -->

    ## $outliers
    ##   [1]    79   300   514   742   744   758   761   763   769   804   805   873
    ##  [13]   892   893   899   902   938  1027  1108  1279  1280  1281  1282  1283
    ##  [25]  1284  1294  1308  1309  1428  1429  1430  1431  1433  1437  1438  1556
    ##  [37]  1629  1630  1649  1810  1819  1874  2090  2317  2323  2393  2628  2637
    ##  [49]  2682  2683  2709  2738  2773  2781  2821  2823  2826  2903  2905  3061
    ##  [61]  3169  3170  3171  3188  3221  3223  3347  3348  3438  3572  3664  3869
    ##  [73]  3912  4357  4358  4360  4361  4362  4364  4365  4366  4367  4439  4441
    ##  [85]  4442  4443  4444  4495  4655  4656  4695  4740  4867  4868  4869  4870
    ##  [97]  4871  4872  4873  4874  4875  4876  4877  4886  4918  4919  5026  5028
    ## [109]  5145  5158  5161  5162  5171  5387  5627  5700  5751  5767  5768  5821
    ## [121]  5847  5848  5850  5851  5976  5978  5979  5982  6005  6132  6133  6136
    ## [133]  6286  6659  6667  6671  6676  6677  6679  6684  6685  6687  6701  6791
    ## [145]  6927  6928  6937  7084  7087  7102  7170  7300  7317  7321  7322  7323
    ## [157]  7640  7743  7788  7790  7791  7997  8011  8047  8272  8274  8275  8276
    ## [169]  8339  8340  8341  8343  8344  8345  8351  8352  8569  8607  8620  8752
    ## [181]  8773  8840  8890  8899  8904  8979  8980  9018  9023  9038  9072  9096
    ## [193]  9204  9306  9429  9430  9431  9478  9556  9558  9564  9921 10145 10493
    ## [205] 10495 10497 10681 10682 10687 10696 10697 10698 10718 10721 10762 10769
    ## [217] 10770 10771 10776 10814 10815 10940 10942 10996 11104 11233 11243 11244
    ## [229] 11350 11534 11690 11778 11784 11785 11869 11870 11913 11954 11966 11969
    ## [241] 11970 11971 12115 12116 12445 12446 12623 12683 12684 12768 12976 13203
    ## [253] 13215 13216 13505 13506 13627 13677 13695 13740 13916 13919 13920
    ## 
    ## $nb_outliers
    ## [1] 263

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
write.csv(table_results,"table_results_Phalepensis.csv") #change species name
```
