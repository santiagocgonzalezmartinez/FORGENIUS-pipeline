vcftools --gzvcf data/Aalba_SNP_sampleFilt.vcf.gz --bed data/Aalba_random_regions.txt --recode --stdout | gzip -c > data/Aalba_SNP_sampleFilt_random.vcf.gz

