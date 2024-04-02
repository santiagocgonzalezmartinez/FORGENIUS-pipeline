FORGENIUS-pipeline, test 20 SNPs
================
SCGM & MW
2024-04-02

Load packages

``` r
rm(list = ls()) # cleans the working environment

#IMPORTANT: don't change order of libraries as some mask others!
library(vcfR)
library(adegenet)
library(hierfstat)
library(LEA)
library(dplyr)
library(mapplots)
library(shapefiles)
library(miscTools)
```

Read data and convert to different formats

``` r
species_name = "Aalba_20SNPs" #change species name
data_file <- paste("data/", species_name, ".vcf.gz", sep="")
data_vcf <- read.vcfR(data_file) 
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
write.bayescan(dat=data_hierfstat,diploid=TRUE,fn="data_bayescan") #writing this file can take about 6 hours for 90k SNPs; better to store these files once generated!

# structure format to external file
write.struct(dat=data_hierfstat,ilab=individual_names,pop=NULL,MARKERNAMES=FALSE,MISSING=-9,fname="snmf/data_structure") 

# geno format to external file
struct2geno (input.file="snmf/data_structure", ploidy=2, FORMAT=2, extra.row=0, extra.column=1) #writting this file can also take quite long
```

    ## Input file in the STRUCTURE format. The genotypic matrix has 413 individuals and 20 markers. 
    ## The number of extra rows is 0 and the number of extra columns is 1 .
    ## Missing alleles are encoded as -9 , converted as 9.
    ## Output files: snmf/data_structure.geno  .lfmm.

``` r
# genepop format to external file - need to be fixed or try plink
#source("functions/TJ_genind2genepop_function.R")
#genind2genepop(data_genind,file = "data_genepop.txt")
```

Run admixture analysis with LEA

Estimate admixture coefficients

``` r
# processing sNMF results
plot(data_snmf, cex = 1.2, col = "lightblue", pch = 19) # run this line separately to determine best K
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-1.png)<!-- -->

``` r
best_K = 5  #select best K based on entropy graph
ce <-  cross.entropy(data_snmf, K = best_K)
best_run <- which.min(ce)

# plot barchart for best run of best K
barchart (data_snmf, best_K, best_run, sort.by.Q = T, col = rainbow(best_K), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients") #use "sort.by.Q = F, names.arg= individual_names, cex.names = 0.4, las=2, " to get ind labels in original order
```

![](FORGENIUS-pipeline_files/figure-gfm/Estimate%20admixture%20coefficients-2.png)<!-- -->

    ## $order
    ##   [1]  13  28  40  46  47  52  61  73  76  83  85  86  99 101 102 103 105 122
    ##  [19] 132 167 170 174 176 187 192 193 207 209 212 218 230 277 278 280 281 291
    ##  [37] 292 293 296 307 319 331 340 363 371 372 373 382 384 401 404 412   1   2
    ##  [55]   3   6   7   8   9  11  12  15  17  18  19  20  21  23  25  26  29  31
    ##  [73]  35  36  37  38  39  41  42  43  44  45  48  49  50  51  53  55  56  57
    ##  [91]  58  65  66  67  69  70  72  74  75  78  80  81  82  84  87  88  89  90
    ## [109]  91  92  93  94  95  96  97  98 104 106 107 108 109 110 111 112 113 115
    ## [127] 116 117 118 120 121 123 124 125 126 127 128 129 131 133 134 135 137 139
    ## [145] 140 141 144 145 146 147 148 149 150 151 152 153 155 156 157 158 159 160
    ## [163] 161 162 163 166 168 169 173 175 177 178 179 183 184 185 186 190 194 197
    ## [181] 199 200 201 204 205 215 216 217 219 220 221 222 223 224 225 226 227 228
    ## [199] 231 233 234 235 236 237 238 239 242 243 244 245 246 247 248 253 256 259
    ## [217] 261 262 267 268 272 276 283 284 289 290 294 295 297 298 299 303 304 306
    ## [235] 309 312 314 317 321 328 333 334 335 337 339 343 345 346 347 349 352 353
    ## [253] 356 358 361 362 365 367 368 375 378 379 383 386 388 389 390 393 395 398
    ## [271] 402 407 408 409   4   5  10  22  27  33  64  68  71  77 100 172 181 182
    ## [289] 191 203 210 213 229 232 270 271 282 285 287 300 357 366 392  16  24  30
    ## [307]  32  34  54  62  63 114 119 130 136 138 142 154 164 165 180 189 195 198
    ## [325] 202 206 208 211 214 241 252 254 255 257 258 260 263 264 266 273 275 279
    ## [343] 286 301 305 311 313 315 316 318 320 322 323 327 329 330 332 336 344 348
    ## [361] 359 360 364 370 374 387 391 394 396 397 400 403 406 410 411 413  14  59
    ## [379]  60  79 143 171 188 196 240 249 250 251 265 269 274 288 302 308 310 324
    ## [397] 325 326 338 341 342 350 351 354 355 369 376 377 380 381 385 399 405

``` r
# estimate ancestry coefficients by pop
q_mat <- LEA::Q(data_snmf, K = best_K, run = best_run) 
colnames(q_mat) <- paste0("Q", 1:best_K)
q_matrix <- data.frame(unlist(q_mat))
q_matrix <- q_matrix %>%
  mutate(individual_names, pop_labels)
q_by_pop <- aggregate(q_matrix[, 1:best_K], list(q_matrix$pop_labels), mean)

# clean-up
#remove.snmfProject("snmf/data_structure.snmfProject")

# plot map for best K
data_coord <- paste("data/", species_name, "_coord.txt", sep="")
coord = read.table(data_coord, header = T) #GCU, longitude and latitude in this order
Npop = length(unique(pop_names))
qmatrix2plot <- merge(coord, q_by_pop, by.x = "GCU", by.y = "Group.1")
xlim <- c(-10,30) 
ylim <- c(35,55)  #adjust if needed
basemap(xlim, ylim) 
#shape <- read.shapefile("shapefiles/chorological_maps_dataset/Abies_alba_plg") #change species name
shape <- read.shapefile("shapefiles/Europe_merged")
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
system2("powershell", arg = c("-file", "bayescan.ps1")) #check out running parameters in ps1 file; bayescan is very slow, took about 5 hours for 15 pops and 14k SNPs using 12 CPUs and standard parameters
```

    ## Warning in system2("powershell", arg = c("-file", "bayescan.ps1")): error in
    ## running command

``` r
#If you are running RStudio in Linux

source ("bayescan/plot_R.r")
plot_bayescan("bayescan/results/data_bayescan_fst.txt",add_text=T,pos=0.02,size=0.7,FDR=0.05)
```

![](FORGENIUS-pipeline_files/figure-gfm/genetic%20differentiation-1.png)<!-- -->

    ## $outliers
    ## integer(0)
    ## 
    ## $nb_outliers
    ## [1] 0

``` r
output_bayescan<-read.table("bayescan/results/data_bayescan.sel")
output_bayescan$logL<-NULL
pop_specific_fst<-colMeans(output_bayescan, na.rm = TRUE)
```

Print table

``` r
# prepare dataframe with basic stats
table_basic_stats<-data.frame(unlist(individuals_per_pop), unlist(nb_loci), unlist(polymoprhic_loci), unlist(Hs), unlist(Fis), unlist(ave_pairwise_fst), unlist(pop_specific_fst), row.names = NULL)
names (table_basic_stats) = c("GCU", "N", "nb_loci", "polymorphic_loci", "Hs", "Fis", "ave_pairwise_fst", "pop_specific_fst")

# merge dataframes in a single table
table_results <- merge(table_basic_stats, q_by_pop, by.x = "GCU", by.y = "Group.1") # merge dataframes by GCU name, it orders the new dataframe alphabetically
data_print <- paste("results/table_results", species_name, ".csv", sep="")
write.csv(table_results,"results/table_results_Aalba_20SNPs.csv")
```
