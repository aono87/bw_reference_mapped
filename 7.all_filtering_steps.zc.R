### some of this has been modified slightly from Emma's code as her's didnt run on my samples and i had to do some de-buggine

##load required libraries
library(readr)
library(scales)
library(reshape2)
library(tidyverse)
library(stringr)


######THIS CODE IS UPDATED WITH NEW EMMA CODE AS OF JULY 2018 ##################################################################
## Step 1: Initial filtering by quality and missing data


# Filter genotypes with genotype quality < 20 and minimum depth of 5
system("vcftools --vcf populations.snps.vcf --out out.1 --minDP 5 --minGQ 20 --recode --recode-INFO-all")


# Filter out sites that were made monomorphic by the previous filter
system("vcftools --vcf out.1.recode.vcf --maf 0.001 --out out.2 --recode --recode-INFO-all")

# Remove sites with more than 50% missing data
system("vcftools --vcf out.2.recode.vcf --out out.3 --max-missing 0.5 --recode --recode-INFO-all")

# Cannot remove loci with extreme allele balance or filter by SNP quality as CH did but this is accounted for by the genotype calling algorithm used by Stacks

# Produce a file with missingness per individual
system("vcftools --vcf out.3.recode.vcf --out out.3 --missing-indv")

################################################################
## Step 2: Filtering and plotting inidvidual samples by missing data

# Plot missingness
# Load the data for the missingness file and sort data by population

out_3_imiss <- read_tsv("out.3.imiss")

pop<-c(rep("med", 35), rep("natl", 80), rep("npac", 25), rep ("satl", 2), rep("spac", 12))
missingness<-cbind(out_3_imiss, pop)
missplot<-melt(missingness[,c(1,5,6)])
Pops<-c("med", "natl", "npac", "satl", "spac")
missplot$INDV<-factor(missplot$INDV, levels = missplot$INDV[order(missplot$pop)])  


# set colours for plotting
myPal4 <- c("green" ,"blue", "yellow", "pink", "orange")

q<-ggplot(missplot, aes(x = INDV, y=value, fill=pop))
q<-q+geom_bar(stat = "identity") + scale_y_continuous(labels=percent_format())+scale_fill_manual(values = myPal4)
q<-q+theme(axis.text.x = element_text(size=6, angle = 90),
             axis.text.y = element_text(size=10),
             legend.text=element_text(size=10),
             legend.title=element_blank())

# save plot as pdf
pdf("missingness_per_sample_by_pop.zcav.final154.pdf")
q
dev.off()

# Select individuals with more than 50% missing data
miss_50 <- filter(out_3_imiss, F_MISS > 0.5) %>% select(INDV)

# Write the individuals to remove to a file
write_delim(miss_50, "remove.3.inds", col_names = FALSE)

# Remove individuals with >50% missing data
system("vcftools --vcf out.3.recode.vcf --out out.4 --remove remove.3.inds --recode --recode-INFO-all")

####################################################################
## Step 3: Filtering loci with high read depth and missing data
# Calculate site depth
system("vcftools --vcf out.4.recode.vcf --site-depth --out out.5")

# Read in the site depth file and calculate mean depth divided by number of samples
site_depth_5 <- read_tsv("out.5.ldepth") %>% mutate(MEAN_DEPTH = SUM_DEPTH / 154)

# Plot a histogram of the mean site depth per locus
qplot(site_depth_5$MEAN_DEPTH, binwidth = 10)

# Filter out loci with a mean site depth > 3x the overall mean
mean_site_depth_5 <- mean(site_depth_5$MEAN_DEPTH)
to_keep_5 <- filter(site_depth_5, MEAN_DEPTH < 3 * mean_site_depth_5)
mean_site_depth_5_filt <- mean(to_keep_5$MEAN_DEPTH)

# Plot the distribution again
qplot(to_keep_5$MEAN_DEPTH)

# Make a list of the sites to filter
to_filter_5 <- filter(site_depth_5, MEAN_DEPTH >= 3 * mean_site_depth_5) %>% select(CHROM, POS)

# Write the sites to remove to a file
write_delim(to_filter_5, "remove.5.sites", col_names = FALSE)

# Remove the sites with VCFtools
system("vcftools --vcf out.4.recode.vcf --out out.5 --exclude-positions remove.5.sites --recode --recode-INFO-all")



#################################################################
## Step 4: Filtering individuals with >25% missing data

# Calculate individual missingness
system("vcftools --vcf out.5.recode.vcf --out out.6 --missing-indv")
       
# Load the data for the out.6.recode.vcf file
out_6_imiss <- read_tsv("out.6.imiss")

# Plot a quick histogram of the data
qplot(out_6_imiss$F_MISS)

# Remove sites with more than 75% missing data
system("vcftools --vcf out.5.recode.vcf --out out.6 --max-missing 0.75 --recode --recode-INFO-all")

# Calculate individual missingness
system("vcftools --vcf out.6.recode.vcf --out out.7 --missing-indv")
       
# Load the data for the out.7.recode.vcf file
out_7_imiss <- read_tsv("out.7.imiss")

# Plot a quick histogram of the data
qplot(out_7_imiss$F_MISS)

# Select individuals with more than 25% missing data
miss_7 <- filter(out_7_imiss, F_MISS > 0.25) %>% select(INDV)

# Write the individuals to remove to a file
write_delim(miss_7, "remove.7.inds", col_names = FALSE)

# Remove the individuals with more than 25% missing genotypes
system("vcftools --vcf out.6.recode.vcf --out out.7 --remove remove.7.inds --recode --recode-INFO-all")

#Remove indels
system("vcftools --vcf out.7.recode.vcf --out out.8 --remove-indels --recode --recode-INFO-all")

######################################################################
## Step 5: Create 'whitelist' of loci to give back to stakcs and then filter using minor allele frequency, examine hwe and estimate error rate (if there are duplicates)

system("cut out.8.recode.vcf -f 3 > whitelist")
system("sed -e 's/_/\t/g' whitelist > whitelistv2")
system("sed '1,15d' whitelistv2 > whitelistv3")
