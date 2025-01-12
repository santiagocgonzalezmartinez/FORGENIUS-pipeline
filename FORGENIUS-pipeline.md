FORGENIUS-pipeline, test 20 SNPs
================
SCGM & MW
2025-01-12

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
library(dartR)
library(poppr)
library(ggplot2)
library(kableExtra)

# maybe needed, to check:
#library(knitr)
#library(htmltools)
#library(stringr)
#library(radiant.data)


# Function to build tables
kable_mydf <- function(x, boldfirstcolumn,font_size){
  x %>% 
    kable() %>%  
    kable_styling(font_size=font_size,
                  bootstrap_options = c("stripped","hover", "condensed"), full_width = F) %>% 
    {if(boldfirstcolumn == TRUE) column_spec(., 1, bold = T) else .}
}
```

Read data, filtering and convert to different formats

``` r
# keep only random SNPs in VCF file using a bed file to filter in vcftools
# system2("functions/extract_random.sh") # change input and output file names in functions/extract_random.sh

# define file name for input and output
species_name = "Aalba" #change species name if needed
suffix = "_20SNPs" #change suffix if needed

# read data from VCF into genind format
data_file <- paste("data/", species_name, suffix, ".vcf.gz", sep="")
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
data_genind <- vcfR2genind(data_vcf)
pop(data_genind) <- substr(indNames(data_genind), 1,8) #takes pop name from the first 8 digits of sample name, e.g. AUT00215

# read data from genotype file (e.g., from Axiom chips) into genind format
# data_file <- paste("data/", species_name, suffix, ".txt", sep="")
# data <- read.table(data_file, header=T)
# data[1:10,1:10] %>% kable_mydf(boldfirstcolumn = F, font_size = 11) # use for visualization of input data but careful as it may give an error when rendering!
# data_genind <- df2genind(X = data[, -(1:2)], sep="/", pop = data[, 2], ind.names = data[, 1], NA.char ="NA") #change column numbers, separator or missing data string if needed

# Remove loci with too much missing data from genind object
#data_genind <- missingno(data_genind, type = "loci", cutoff = 0.05, quiet = FALSE, freq = FALSE)

# Remove individuals with too much missing data from genind object
#data_genind <- missingno(data_genind, type = "geno", cutoff = 0.05, quiet = FALSE, freq = FALSE)

# Filter genind by locus list
#loci2keep <- scan(file="gea_outliers.txt", what=character()) # read list of loci to keep
#data_genind = data_genind[loc=loci2keep] # filter by locus name in list

# Filter by Minor Allele Frequency (MAF) from genind object
data_genind <- informloci(data_genind, cutoff = 2/nInd(data_genind), MAF = 0.00, quiet = FALSE) # Keep at least 0.00, to remove monomorphic loci in the full dataset
```

    ## cutoff value: 0.484261501210654 % ( 2 samples ).

    ## MAF         : 0

    ## 
    ## All sites polymorphic

``` r
# extract individual names
individual_names<-indNames(data_genind)

# genpop format (don't mistake wiwth genepop format, below)
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
nb_pop<-length(pop_names)

# counting individuals per population
pop_labels <- data_genind@pop
individuals_per_pop <- table(pop_labels)

# hierfstat format
data_hierfstat <- genind2hierfstat(data_genind,pop=NULL)

# bayescan format to external file
write.bayescan(dat=data_hierfstat,diploid=TRUE,fn="bayescan/data_bayescan")

# structure format to external file
write.struct(dat=data_hierfstat,ilab=individual_names,pop=NULL,MARKERNAMES=FALSE,MISSING=-9,fname="snmf/data_structure") 

# geno format to external file
struct2geno (input.file="snmf/data_structure", ploidy=2, FORMAT=2, extra.row=0, extra.column=1) #writing this file can take quite long
```

    ## Input file in the STRUCTURE format. The genotypic matrix has 413 individuals and 20 markers. 
    ## The number of extra rows is 0 and the number of extra columns is 1 .
    ## Missing alleles are encoded as -9 , converted as 9.
    ## Output files: snmf/data_structure.geno  .lfmm.

``` r
# genepop format to external file
data_gl <- gi2gl(data_genind) # convert first to genlight object
```

    ## Starting gi2gl 
    ## Starting gl.compliance.check 
    ##   Processing genlight object with SNP data
    ##   The slot loc.all, which stores allele name for each locus, is empty. 
    ## Creating a dummy variable (A/C) to insert in this slot. 
    ##   Checking coding of SNPs
    ##     SNP data scored NA, 0, 1 or 2 confirmed
    ##   Checking locus metrics and flags
    ##   Recalculating locus metrics
    ##   Checking for monomorphic loci
    ##     No monomorphic loci detected
    ##   Checking for loci with all missing data
    ##     No loci with all missing data detected
    ##   Checking whether individual names are unique.
    ##   Checking for individual metrics
    ##   Warning: Creating a slot for individual metrics
    ##   Checking for population assignments
    ##     Population assignments confirmed
    ##   Spelling of coordinates checked and changed if necessary to 
    ##             lat/lon
    ## Completed: gl.compliance.check 
    ## Completed: gi2gl

``` r
data_genepop <- gl2genepop(
  data_gl,
  outfile = "data_genepop",
  outpath = "./metapop",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = 5
) 
```

    ## Starting gl2genepop 
    ## [dartR vers. 2.9.7 Build = Jody ]
    ##   Processing genlight object with SNP data
    ##   The genepop file is saved as:  ./metapop/data_genepop/
    ## Completed: gl2genepop

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
    ##   [1]   1   2   3   6   7   8   9  12  15  17  18  20  21  24  25  26  29  31
    ##  [19]  35  36  39  41  42  43  44  45  48  49  50  51  53  55  57  58  66  67
    ##  [37]  75  77  78  81  85  87  89  90  91  92  93  94  95  96  97  98 104 106
    ##  [55] 107 108 109 110 111 112 113 115 116 117 118 120 121 123 124 125 126 127
    ##  [73] 128 129 131 133 134 135 137 140 141 144 145 146 147 148 149 150 151 152
    ##  [91] 153 155 156 157 158 159 160 161 162 163 166 168 169 173 175 177 179 183
    ## [109] 184 185 186 197 199 201 203 204 205 215 216 217 219 220 221 222 223 227
    ## [127] 228 231 232 233 236 237 239 244 245 247 253 262 267 268 272 277 278 283
    ## [145] 284 295 298 299 304 306 312 314 321 333 335 337 347 353 356 358 361 362
    ## [163] 368 375 379 383 388 389 390 393 395 402 407  11  19  23  37  38  56  59
    ## [181]  65  69  70  72  74  80  82  84  88 139 171 178 190 194 200 210 224 225
    ## [199] 226 229 234 235 238 242 243 246 248 250 254 256 259 261 270 276 282 285
    ## [217] 289 290 294 297 300 303 308 309 317 325 328 334 339 343 345 346 349 352
    ## [235] 365 367 372 378 386 392 398 408 409  14  16  79  83 143 180 206 208 211
    ## [253] 212 269 273 274 275 279 288 305 315 318 320 338 350 355 370 381   4   5
    ## [271]  10  13  22  27  28  33  40  46  47  52  61  64  68  71  73  76  86  99
    ## [289] 100 101 102 103 105 122 132 167 170 172 174 176 181 182 187 191 192 193
    ## [307] 207 209 213 218 230 271 280 281 287 291 292 296 307 319 331 340 363 366
    ## [325] 371 373 382 384 401 404 412  30  32  34  54  60  62  63 114 119 130 136
    ## [343] 138 142 154 164 165 188 189 195 196 198 202 214 240 241 249 251 252 255
    ## [361] 257 258 260 263 264 265 266 286 293 301 302 310 311 313 316 322 323 324
    ## [379] 326 327 329 330 332 336 341 342 344 348 351 354 357 359 360 364 369 374
    ## [397] 376 377 380 385 387 391 394 396 397 399 400 403 405 406 410 411 413

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
allele_freq <- tab(data_genpop, freq=TRUE)

# Count monomorphic loci and missing loci per population
monomorphic_counts <- rowSums(allele_freq == 1)
nb_loci <- rowSums(allele_freq != "NA")/2
nb_loci_zero <- rowSums(allele_freq == 0)
nb_loci_noNA = nb_loci - (nb_loci_zero-monomorphic_counts)/2
polymoprhic_loci <- 1-(monomorphic_counts/nb_loci_noNA)
```

Compute genetic diversity and inbreeding stats with hierfstat

``` r
# Compute basic statistics
Hs <- Hs(data_hierfstat)
Ho <- Ho(data_hierfstat)
Fis <- 1-(Ho/Hs)
```

Compute population genetic and allelic diversity contributions with
metapop2

``` r
# Genetic and allelic diversity partitions with metapop2
system2("./metapop/metapop", arg = "metapop/data_genepop metapop/config")

# Clean-up
file.copy("res_run_FULL.md", "metapop/results", overwrite = TRUE)
```

    ## [1] TRUE

``` r
file.copy("res_run_p.csv", "metapop/results",  overwrite = TRUE)
```

    ## [1] TRUE

``` r
file.remove("res_run_FULL.md")
```

    ## [1] TRUE

``` r
file.remove("res_run_p.csv")
```

    ## [1] TRUE

``` r
# Read all lines from the file
file_path <- "metapop/results/res_run_FULL.md"
lines <- readLines(file_path)

# GENETIC DIVERSITY CONTRIBUTIONS
# Identify the start of the table
start_line <- grep("Percentage of loss \\(\\+\\) or gain \\(-\\) of gene diversity after removal of each subpopulation", lines) + 4

# Extract the lines containing the table data
table_lines <- lines[start_line:(start_line + (nb_pop-1))]

# Clean up the lines: remove leading/trailing spaces and split by '|' character
table_data <- lapply(table_lines, function(line) {
  # Remove extra spaces and split by '|'
  values <- strsplit(trimws(line), "\\|")[[1]]
  # Further clean the values by trimming and removing any extra characters
  trimws(values)
})

# Convert list of vectors to a data frame
contrib_div <- as.data.frame(do.call(rbind, table_data), stringsAsFactors = FALSE)

# Set column names appropriately
names(contrib_div) <- c("ID1", "cHs", "cDg", "cHt")

# Add pop levels
contrib_div <- contrib_div %>% mutate(GCU=pop_names, .after=ID1)

# Convert to numeric
contrib_div$cHs <- as.numeric(contrib_div$cHs)
contrib_div$cDg <- as.numeric(contrib_div$cDg)
contrib_div$cHt <- as.numeric(contrib_div$cHt)

# ALLELIC DIVERSITY CONTRIBUTIONS
# Identify the start of the table
start_line <- grep("Percentage of loss \\(\\+\\) or gain \\(-\\) of allelic diversity after removal of each subpopulation", lines) + 4

# Extract the lines containing the table data (next 3 lines)
table_lines <- lines[start_line:(start_line + (nb_pop-1))]

# Clean up the lines: remove leading/trailing spaces and split by '|' character
table_data <- lapply(table_lines, function(line) {
  # Remove extra spaces and split by '|'
  values <- strsplit(trimws(line), "\\|")[[1]]
  # Further clean the values by trimming and removing any extra characters
  trimws(values)
})

# Convert list of vectors to a data frame
contrib_a <- as.data.frame(do.call(rbind, table_data), stringsAsFactors = FALSE)

# Set column names appropriately
names(contrib_a) <- c("ID2", "cAs", "cDa", "cAt")

# Add pop levels
contrib_a <- contrib_a %>% mutate(GCU=pop_names, .after=ID2)

# Convert to numeric
contrib_a$cAs <- as.numeric(contrib_a$cAs)
contrib_a$cDa <- as.numeric(contrib_a$cDa)
contrib_a$cAt <- as.numeric(contrib_a$cAt)

# Merge dataframes
contrib <- merge(contrib_div, contrib_a, by = "GCU") # Careful, merge() reorder the populations by alphabetical order (as in the output)
contrib <- subset(contrib, select = -c(ID1, ID2))
str(contrib)
```

    ## 'data.frame':    18 obs. of  7 variables:
    ##  $ GCU: chr  "AUT00179" "AUT00215" "DEU00114" "ESP00339" ...
    ##  $ cHs: num  1.872 1.601 1.619 0.244 0.126 ...
    ##  $ cDg: num  -0.5348 -0.2233 -0.4445 0.6569 -0.0556 ...
    ##  $ cHt: num  1.3376 1.3781 1.1743 0.9005 0.0702 ...
    ##  $ cAs: num  0.854 0.927 0.543 -0.33 0.26 ...
    ##  $ cDa: num  -0.261 -0.186 -0.239 -0.301 0.135 ...
    ##  $ cAt: num  0.593 0.741 0.304 -0.632 0.395 ...

``` r
# Create contribution figures
# Genetic diversity
contrib_temp1 <- subset(contrib, select = c(GCU, cHs))
colnames(contrib_temp1) <- c("GCU", "values")
contrib_temp1 <- contrib_temp1 %>% mutate(Partition="Diversity", .after=values)
contrib_temp2 <- subset(contrib, select = c(GCU, cDg))
colnames(contrib_temp2) <- c("GCU", "values")
contrib_temp2 <- contrib_temp2 %>% mutate(Partition="Differentiation", .after=values)
contrib_fig_div <- rbind(contrib_temp1,contrib_temp2)
ggplot(contrib_fig_div, aes(x = GCU, y = values))+ geom_col(aes(fill = Partition), width = 0.7) + scale_y_continuous() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("GCUs") + ylab("Contribution % (genetic diversity)") + scale_alpha_manual(values = 1) + labs(alpha="")
```

![](FORGENIUS-pipeline_files/figure-gfm/genetic%20and%20allelic%20diversity%20partitions-1.png)<!-- -->

``` r
# Allelic diversity
contrib_temp3 <- subset(contrib, select = c(GCU, cAs))
colnames(contrib_temp3) <- c("GCU", "values")
contrib_temp3 <- contrib_temp3 %>% mutate(Partition="Diversity", .after=values)
contrib_temp4 <- subset(contrib, select = c(GCU, cDa))
colnames(contrib_temp4) <- c("GCU", "values")
contrib_temp4 <- contrib_temp4 %>% mutate(Partition="Differentiation", .after=values)
contrib_fig_allelic <- rbind(contrib_temp3,contrib_temp4)
ggplot(contrib_fig_allelic, aes(x = GCU, y = values))+ geom_col(aes(fill = Partition), width = 0.7) + scale_y_continuous() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("GCUs") + ylab("Contribution % (allelic diversity)") + scale_alpha_manual(values = 1) + labs(alpha="")
```

![](FORGENIUS-pipeline_files/figure-gfm/genetic%20and%20allelic%20diversity%20partitions-2.png)<!-- -->

Compute genetic distinctness

``` r
# Average pairwise Fst
pairwise_WCfst <- pairwise.WCfst(data_hierfstat) #this is very slow, took about 16 hours for 15 pops and 90k SNPs in a laptop
ave_pairwise_fst<-rowMeans(pairwise_WCfst, na.rm = TRUE) 

# Population-specific Fst from bayescan

#If you are running RStudio in Windows 10
#system2("powershell", arg = c("-file", "bayescan.ps1")) #check out running parameters in ps1 file; bayescan is very slow in Windows, took about 5 hours for 15 pops and 14k SNPs using 12 CPUs and standard parameters

#If you are running RStudio in Linux
system2("./bayescan/BayeScan2.1_linux64bits", arg = "bayescan/data_bayescan -n 5000 -thin 10 -nbp 10 -pilot 5000 -burn 10000 -pr_odds 1000 -od bayescan/results")

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
table_results <- merge(table_basic_stats, q_by_pop, by.x = "GCU", by.y = "Group.1")
table_results_contrib <- merge(table_results, contrib, by = "GCU")
# merge dataframes by GCU name, it orders the new dataframe alphabetically
data_print <- paste("results/table_results_contrib_", species_name, suffix, ".csv", sep="")
write.csv(table_results_contrib, data_print)
```
