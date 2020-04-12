---
title: Derived Homo sapiens cis-eQTL regulation - implications for brain evolution
authors: Alejandro Andirkó, Cedric Boeckx
date: 2020
contact: cedric.boeckx@ub.edu
---

This markdown document contains all the code necessary to replicate the manuscript "Differential gene regulation by Homo sapiens-specific alleles in brain tissues: implications for brain evolution", including figures. Materials and other information not included here can be retrieved from the github repository.

#Folder structure
eQTL
  ├─ hsap         #contains the files resulting from combining the GTEX and Homo sapiens variation databases
  ├─ materials    #contains raw and preprocessed files   
  │   └── GTEx_Analysis_v8_eQTL
  │   └── permutations #files for permutation tests (positive selection)
  └─ postprocessing                                             
      ├── clinvar # contains the results of crossing with clinvar data
      ├── clumped # SNP clumping analysis
      ├── gwas    # GWAS data 
           ├── gwascat      # data from the NHGRI-EBI GWAS catalog
      └── plots   # contains files with specific information for plotting
           ├── circos 
           ├── conseqs
           └── TSS



#Materials
Creates some of the parent directories for the project, and download the original list of variants from Kuhlwilm & Boeckx (2019), as well as the list of GTEX significant variants.

```sh
#Create directories
mkdir eQTL
cd eQTL
mkdir materials
mkdir hsap
mkdir postprocessing

cd materials

#Get Kuhlwilm & Boeckx (2019) raw data 
wget https://ndownloader.figshare.com/files/15253052
mv 15253052 NNaall.tsv.gz
gzip -d NNaall.tsv.gz 

#Get GTEX data (v8)
#Listed as "eGene and significant variant-gene associations based on permutations" in https://gtexportal.org/home/datasets
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
tar xf GTEx_Analysis_v8_eQTL.tar && rm GTEx_Analysis_v8_eQTL.tar

#Remove all non-brain tissues (minus Adrenal gland) - see justification in the article of the text
cd GTEx_Analysis_v8_eQTL
rm -v !(Brain*|Adrenal*)
rm -v !(*pairs*)
gzip -d * & cd ..
```


#Filtering/Preprocessing
## Filtering the raw Kuhlwilm & Boeckx (2019) data
Kuhlwilm & Boeckx (2019) (from now on, K&B) provide all their data in one file, including all frequencies of derived alleles. However, the results discussed in their article only concern those variants that are located in the +90% in African, American, East Asian, European and South Asian metapopulations. K&B (2019) define high-frequency mutations as:

> High-frequency (HF) differences: Positions where more than 90% of present-day humans carry a derived allele, while at least the Denisovan and one Neanderthal carry the ancestral allele, accounting for different types of errors and bi-directional gene flow.


We applied the same criteria to the master file provided in the supplementary:

```sh
# From eQTL/materials
awk '{
if ($2 >= 0.90 && $2 != "NA" && $4 != $6 && $6 !~ "*") print $0}' NNaall.tsv | awk '{ 
if ($7 == 1 || $10 == 1 || $13 == 1 && $7 != 0 ) print $0}' |  awk '{
if ($10 == "NA" && $13 == "NA") next;
else if ($7 == "NA" && $13 == "NA") next;
else if ($13 == 0 && $14 == 1) next;
else if ($11 == 0 && $14 == 0) next;
else if ($8 == 0 && $13 == 0) next;
else if ($10 == 0 && $13 == 0) next;
else if ($10 == 0 && $13 == "NA") next;
else if ($10 == "NA" && $13 == "NA") next;
else if ($10 == "NA" && $13 == 0) next;
else if ($7 == "NA" && $10 == "NA") next;
else if ($7 == "NA" && $13 == 0) next;
else if ($5 != $6) next;
else if ($7 == 0 || $10 == 0 && $13 == 0.5) next;
else if ($13 == 0 && $10 == 1 && $7 == 1 ) print $0;
else {print $0}
}' | grep -v "MAF\=.\:0\.[1-9]" > hfsapiens.freqfilt.tsv


#For a reasoning on why these criteria were chosen, please read the methods section in the article
```



##RsIDs for the GTEx files 
This part of the script downloads a lookup table to change the default format of the GTEx data to rsIDs, a less error-prone way to perform analysis accross databases. 

```sh
# From eQTL/materials
cd GTEx_Analysis_v8_eQTL

#Download, decompress, change name for convenience
wget https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
gzip -d GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
mv GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt GTEx_lookup_table.txt

# Pruning out those without rsID, which include an otherwise inexistent "." instead of the pertinent rsid
grep -v "\." GTEx_lookup_table.txt > GTEx_lookup_table_filtd.txt

#Split GTEx_lookup_table_filtd.txt in four more manageable files for the awk array of the next step to work without depleting memory
split --number=l/6 ${fspec} GTEx_lookup_table_filtd.txt
rm GTEx_lookup_table.txt GTEx_lookup_table_filtd.txt
mkdir splitlookup && mv xaa xab xac xad xae xaf splitlookup/

# An awk loop to add the rsID to all files. Divided due to computer memory problems.
for i in *.txt; do
    awk 'FNR==NR{a[$1]=$7;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' splitlookup/xaa $i >> $i.out
    awk 'FNR==NR{a[$1]=$7;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' splitlookup/xab $i >> $i.out
    awk 'FNR==NR{a[$1]=$7;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' splitlookup/xac $i >> $i.out
    awk 'FNR==NR{a[$1]=$7;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' splitlookup/xad $i >> $i.out
    awk 'FNR==NR{a[$1]=$7;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' splitlookup/xae $i >> $i.out
    awk 'FNR==NR{a[$1]=$7;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' splitlookup/xaf $i >> $i.out
done
rm -r splitlookup *.txt
cd ..
```


# 3. Cross the human variation data with the GTEx significant variant database.

```sh
#Here (hsap/) we will store the product of crossing GTEx and the Homo Sapiens high frequency mutations database
# From the main eQTL/ directory:
cp materials/GTEx_Analysis_v8_eQTL/* hsap/
cp materials/hfsapiens.freqfilt.tsv hsap/ && cd hsap

#Print the rsids only and clean 1. missing rsIDS
awk '{print $19}' hfsapiens.freqfilt.tsv | sed 's/-//g' | sed 's/,/\n/g' | sed '/^$/d' > hfsapiens.rsids

for i in *.out; do
    fgrep -w -f hfsapiens.rsids $i > $i.rs
done

rm *.txt.out
rm hf*
```


#Preprocessing 
##Clumping
We then clumped the data on Linkage Disequilibrium blocks based on eQTL mapping p-value. 

```sh
#from main directory of project (eQTL)
cd hsap
for i in *.rs;do
  awk '{print $2, $8}' $i > $i.pval
done

cd ..
mv hsap/*.pval postprocessing
#Add headers - necessary for plink clumping later on
cd postprocessing/
for i in *.rs.pval; do
  echo 'SNP P' | cat - $i > temp && mv temp $i
done

#wget vcfs of chromosomes from the 1000 genomes project ftp server
wget -m -A "*.vcf.gz" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ --reject "index.html*"


#creates a folder for each tissue where we will keep the .vcf for each chromosome
for i in *.rs.pval; do
  dir="${i%%[.]*}"
  mkdir $dir/
  mv $dir.* $dir/
done

#Now for the extraction of the snps from each chromosome (per tissue)
#Remember that the .rs files are already filtered by frequency, as per step 3
for i in *.vcf.gz; do
  vcftools --gzvcf $i --snps Adrenal_Gland/Adrenal_Gland.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Adrenal_Gland/
  vcftools --gzvcf $i --snps Brain_Frontal_Cortex_BA9/ Brain_Frontal_Cortex_BA9.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Frontal_Cortex_BA9/
  vcftools --gzvcf $i --snps Brain_Amygdala/Brain_Amygdala.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Amygdala/
  vcftools --gzvcf $i --snps Brain_Hippocampus/Brain_Hippocampus.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out Brain_Hippocampus/$i
  mv *.recode.vcf Brain_Hippocampus/
  vcftools --gzvcf $i --snps Brain_Anterior_cingulate_cortex_BA24/Brain_Anterior_cingulate_cortex_BA24.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Anterior_cingulate_cortex_BA24/
  vcftools --gzvcf $i --snps Brain_Hypothalamus/Brain_Hypothalamus.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Hypothalamus/
  vcftools --gzvcf $i --snps Brain_Caudate_basal_ganglia/Brain_Caudate_basal_ganglia.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Caudate_basal_ganglia/
  vcftools --gzvcf $i --snps Brain_Nucleus_accumbens_basal_ganglia/Brain_Nucleus_accumbens_basal_ganglia.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Nucleus_accumbens_basal_ganglia/
  vcftools --gzvcf $i --snps Brain_Cerebellar_Hemisphere/Brain_Cerebellar_Hemisphere.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Cerebellar_Hemisphere/
  vcftools --gzvcf $i --snps Brain_Putamen_basal_ganglia/Brain_Putamen_basal_ganglia.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Putamen_basal_ganglia/
  vcftools --gzvcf $i --snps Brain_Cerebellum/Brain_Cerebellum.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Cerebellum/
  vcftools --gzvcf $i --snps Brain_Spinal_cord_cervical_c-1/Brain_Spinal_cord_cervical_c-1.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Spinal_cord_cervical_c-1/
  vcftools --gzvcf $i --snps Brain_Cortex/Brain_Cortex.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Cortex/
  vcftools --gzvcf $i --snps Brain_Substantia_nigra/Brain_Substantia_nigra.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Brain_Substantia_nigra/
  Pituitary/Pituitary.v8.signif_variant_gene_pairs.txt.out.rs.pval --recode --recode-INFO-all --out $i
  mv *.recode.vcf Pituitary/
done
rm *.log *.gz *.gz.tbi

# Name cleaning
for fold in */; do
  rename 's/phase.*.vcf/recode.vcf/' $fold/*.vcf
done

#Download plink 1.9 / unzip / Copy plink file to all tissue folders
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200121.zip
unzip plink_linux_x86_64_20200121.zip
for fold in */; do
  cp plink $fold/
done

# Get .map file
for fold in */; do
  for a in $fold/*.vcf; do
    ./plink --vcf $a --recode --out $a.out
    rm $fold/*.log
  done
done

# Clumping
for fold in */; do
  for a in $fold/*.out.map; do
    input="${a%%[.]map}"
    ./plink --file $input --clump $fold/*.pval --clump-p1 1 --clump-r2 0.99 --out $input
  done
done

#Remove everything unused, leave only clumpled stuff and logs
for fold in */; do
  cd $fold
  rm -r *.vcf *.map *.nosex *.ped *.log
  cd ..
done

mkdir clumped
cp -r * clumped/
rm -r !(clumped)
cd clumped
rm -r plink* toy* LICENSE prettify clumped/

```

Recheck because it really doesn't work

```sh
#Extract clumped rsIDs: Necessary for clinvar and some of the plots
for fold in */; do
  for i in $fold/*.clumped; do
    dir="${i}"
    awk 'NR>2 {print $3}' $i >> $fold/${fold/\//}.rsclump.txt
  done
done

```

#Postprocessing
##Clinvar
```sh
#From eQTL/postprocessing/
mkdir clinvar
cd clinvar

#Downloads clinvar database
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz 
gzip -d clinvar.vcf.gz
#Fixes rsIDs' formating
sed 's/RS\=/rs/g' clinvar.vcf > clinvar_out.vcf

# Move tissue.rsclump.txts
cd ..
for fold in clumped/*/; do
  for i in $fold/*rsclump.txt; do
    mv $i clinvar/
  done
done

for i in *.txt; do
  sed -r '/^\s*$/d' $i > $i.out
  rm $i
done

for i in *.out; do
   grep -w -f $i clinvar_out.vcf > $i.cvar;
   rm $i 
done

rm clinvar_out.vcf clinvar.vcf

```


#If you want to include the unclumped data, here it is
```sh
cd hsap
for i in *.rs;do
  awk '{print $2}' $i > $i.pval
done
cd ..
mkdir postprocessing/clinvar_unclumped
mv hsap/*.pval postprocessing/clinvar_unclumped
cd clinvar_unclumped 

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz 
gzip -d clinvar.vcf.gz
for i in *.pval; do
   grep -w -f $i clinvar_out.vcf > $i.cvar; 
done
rm !(*.cvar)
```


#Visualizations
##Plots 1 and 2.
**Add K-S test, chi-square and the down-up plot vs all
Get dead code out**

```r
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(sanzo)
setwd("~/v8/version4/eQTL/hsap/") #Change to your path

#Import
temp = list.files(pattern="*")
myfiles = lapply(temp,
                 read_table2,
                 col_names = FALSE)

#Functions
divide <- function(X) {
  split(X,  X$X9 > 0 )
}

divided <- lapply(myfiles, divide)
all <- NULL
all$tissue <- c("Adrenal Gland",
                "Amygdala",
                "BA24",
                "Caudate",
                "Cerebellar hemisphere",
                "Cerebellum",
                "Cortex",
                "BA9",
                "Hippocampus",
                "Hypothalamus",
                "Nucleus Accumbens",
                "Putamen",
                "Spinal cord",
                "Substantia Nigra",
                "Pituitary")
all <- as.data.frame(all)
all2 <- NULL
all2$tissue <- all$tissue
all2 <- as.data.frame(all2)

divided <- lapply(divided[], function(x) lapply(x, '[[', 2))

up <- NULL
down <- NULL
unfilt <- NULL
for (i in 1:15){
  up <- append(up, length(divided[[i]][["TRUE"]]))
  down <- append(down, length(divided[[i]][["FALSE"]]))
  unfilt <- append(unfilt, length(divided[[i]][["FALSE"]])+length(divided[[i]][["TRUE"]]))
}

all$Up <- up
all$Down <- down
all2$Unfilt <- unfilt
all.m <- melt(all)
all2.m <- melt(all2)
#all.m$prop <- reorder(all.m$tissue, all.m$value,FUN=function(x) mean(as.numeric(x)))

#Erase
GTExPortal <- read_csv("~/v8/GTExPortal.csv", 
                       col_names = FALSE, skip = 1) #Change to your path

sample <- NULL
sample$tissue <- all$tissue
sample$sampleval <- GTExPortal$X3
sample$eqtl <- unfilt
sample <- as.data.frame(sample)

#Fill plot
ggplot() +
  theme_minimal() +
  geom_bar(aes(y = value, x = tissue, fill = variable), 
           data = all.m, 
           stat= "identity",
           position="fill",
           alpha = 0.75)+
  scale_fill_manual(values = sanzo.duo("c106")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=0.50, linetype="dashed") 


all2.m$sample <- sample$sampleval
  
#Plot
ggplot(all2.m, aes(x = all2.m$sample, y= all2.m$value, group = all2.m$variable, shape = tissue)) +
  theme_minimal() +
  theme(legend.position="top") +
  labs(y= "Number of eQTL", x = "Sample size") +
  geom_point(size = 3.5) +
  geom_smooth(alpha = 0.25) + 
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) +
  stat_cor(method = "pearson", label.x = 175, label.y = 100,color="black")
#geom_text(label=all2.m$tissue) 
#geom_label(label = all.m$tissue) 

#Pearson: sample ~ value correlation 
cor.test(all2.m$sample, all2.m$value, method = "pearson", conf.level = 0.95)

#Binomial test per tissue
binom.test(all$Up[1], c(all$Up[1]+all$Down[1]), 0.5)
binom.test(all$Up[2], c(all$Up[2]+all$Down[2]), 0.5)
binom.test(all$Up[3], c(all$Up[3]+all$Down[3]), 0.5)
binom.test(all$Up[4], c(all$Up[4]+all$Down[4]), 0.5)
binom.test(all$Up[5], c(all$Up[5]+all$Down[5]), 0.5)
binom.test(all$Up[6], c(all$Up[6]+all$Down[6]), 0.5)
binom.test(all$Up[7], c(all$Up[7]+all$Down[7]), 0.5)
binom.test(all$Up[8], c(all$Up[8]+all$Down[8]), 0.5)
binom.test(all$Up[9], c(all$Up[9]+all$Down[9]), 0.5)
binom.test(all$Up[10], c(all$Up[10]+all$Down[10]), 0.5)
binom.test(all$Up[11], c(all$Up[11]+all$Down[11]), 0.5)
binom.test(all$Up[12], c(all$Up[12]+all$Down[12]), 0.5)
binom.test(all$Up[13], c(all$Up[13]+all$Down[13]), 0.5)
binom.test(all$Up[14], c(all$Up[14]+all$Down[14]), 0.5)
binom.test(all$Up[15], c(all$Up[15]+all$Down[15]), 0.5)



```

# Genes ~ SNPs

# Circos
```sh
#extract rsids, move them to folder
for fold in */; do
  for i in $fold/*.clumped; do
    dir="${i}"
    awk 'NR>2 {print $3}' $i >> $fold/${fold/\//}.rsclump.txt
  done
done

cd ..
for fold in clumped/*/; do
  for i in $fold/*rsclump.txt; do
    mv $i clinvar/
  done
done
 
mkdir circos
cp *.txt circos/ && cd circos/
cat *.txt >> allsnps.txt

#from eQTL
for i in materials/GTEx_Analysis_v8_eQTL/*; do
  grep -w -f postprocessing/plots/circos/allsnps.txt $i >> out.txt
done

#Pipe this
grep -w -f postprocessing/plots/circos/allsnps.txt hsap/*.rs | awk '{print $1, $2, $9}' | sed 's/.*chr/hs/g' | sed 's/\_/ /g' | awk '{print $1, $2, $2, $7}' | sed -e 's/\s/\t/g' > postprocessing/plots/circos/circos_data
cd postprocessing/plots/circos/

#Add header
echo -e 'chr\tstart\tend\tvalue'| cat - circos_data > temp && mv temp circos_data
```

Now you can input `circos_data` into `circos`. See the `circos` folder for the rest of configuration files.
**INCLUDE ALSO CONF FILES**

# Permutation tests
```sh
#From eQTL:
cd materials && mkdir permutations && cd permutations
```
This folder (`permutations`) should include the **Peyregne and Racimo** data (see references) - included in the github repository for convenience. You can contrast this data with the originals in the respective articles. Once you have them:

```sh
#From `permutations` folder
sed 's/hs/chr/g' racimo_coords.txt > racimo_perm
sed 's/hs/chr/g' peyregne_coords.txt > pey_perm
cd ../..
cp postprocessing/plots/circos/circos_data materials/permutations/
awk '{print $1, $2, $3}' circos_data | sed 's/hs/chr/g' > snppermut.txt
rm circos_data racimo_coords.txt peyregne_coords.txt
```

Now, in R...

```r
library(readr)
library(regioneR)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
  
setwd("~/eQTL/materials/permutations/")  #Change to your pathr

pey <- read_table2("pey_perm", col_names = FALSE)
pey <- as.data.frame(pey)
pey <- toGRanges(pey)

racimo <- read_table2("racimo_perm", col_names = FALSE)
racimo <- as.data.frame(racimo)
racimo <- toGRanges(racimo)

data <- read_table2("snppermut.txt")
data <- as.data.frame(data)
data <- toGRanges(data)

set.seed(12345) #ensures that the results are always the same
pt1 <- permTest(A=data, B=pey, ntimes=1000,
                randomize.function=circularRandomizeRegions,
                evaluate.function=numOverlaps, mc.set.seed=FALSE)
pt2 <- permTest(A=data, B=racimo, ntimes=1000,
                randomize.function=circularRandomizeRegions,
                evaluate.function=numOverlaps, mc.set.seed=FALSE)

pt1
pt2
plot(pt1, plotType="Tailed")
plot(pt2, plotType="Tailed")
```


# TSS 
```sh
#From eQTL
mkdir postprocessing/plots/TSS 

#Extract TSS information for high frequency variants, no GTEx filtering
awk '{print $2, $4}' materials/GTEx_Analysis_v8_eQTL/* > postprocessing/plots/TSS/TSS_unfilt

#Extract TSS information for high frequency variants, no GTEx filtering
awk '{print $2, $4}' hsap/* > postprocessing/plots/TSS/TSS_gtex
```

```r
library(ggplot2)
library(viridis)
library(readr)

alltss <- read_table2("~/v8/version4/eQTL/postprocessing/plots/TSS/TSS_unfilt") #Change to your path
gtextss <- read_table2("~/v8/version4/eQTL/postprocessing/plots/TSS/TSS_gtex", col_names = FALSE) #Change to your path

ggplot(alltss, aes(x=alltss$tss_distance)) + 
  geom_histogram(fill = viridis(79), bins = 80) +
  theme_minimal()+
  labs(x = "TSS distance", y = "Counts") +
  scale_x_continuous(name = "Counts", labels = scales::comma, breaks =c(-1000000, -500000, 0, 500000, 1000000), limits= c(-1000000, 1000000))

ggplot(gtextss, aes(x=gtextss$X2)) + 
  geom_histogram(fill = viridis(79), bins = 80) +
  theme_minimal()+
  labs(x = "TSS distance", y = "Counts") +
  scale_x_continuous(name = "Counts", labels = scales::comma, breaks =c(-1000000, -500000, 0, 500000, 1000000), limits= c(-1000000, 1000000))

```

# Variant consequences
```r
library(readr)
library(ggplot2)
library(biomaRt)
library(reshape2)
library(plyr)
ensembl = useEnsembl(biomart="snp", dataset = "hsapiens_snp", GRCh = 37)

input <- read_table2("~/eQTL/postprocessing/plots/TSS/TSS_gtex", col_names = FALSE) #Change to your path

#Queries Biomart
result <- getBM(attributes=c('consequence_type_tv', 'refsnp_id'),
                filters = 'snp_filter',
                values = input$X1,
                mart = ensembl)

#Orders types by frequency in a new dataframe
plotinput <- count(result, vars = c("consequence_type_tv"))
#There are a lot of empty cells
plotinput = plotinput[-1,]

#Reorder factor levels for plot
plotinput$consequence_type_tv<-factor(plotinput$consequence_type_tv, levels = plotinput$consequence_type_tv[order(plotinput$freq)])
plotinput$freq <- as.factor(plotinput$freq)

#Renames labels
labelsx <- c("Stop lost",
             "Splice region variant",
             "Synonymous variant",
             "Missense variant",
             "5 prime UTR variant",
             "3 prime UTR variant",
             "Non coding transcript exon variant",
             "NMD transcript variant",
             "Non coding transcript vartiant", 
             "Intron variant")


ggplot(plotinput, aes(x = consequence_type_tv, y = plotinput$freq)) +
  theme_minimal() +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_colour_viridis_d() +
  coord_flip() +
  labs(x = "Variant consequence", y = "Nr. of eQTL loci") +
  scale_x_discrete(labels = labelsx)
```

## Variant dating #Two R codes for these: check which one to use
With unclumped data

```sh
#From eQTL
mkdir postprocessing/dating && cd postprocessing/dating

#run provided datingdownload.sh script in this folder
awk '{print $2}' hsap/* > postprocessing/dating/rsids
cd postprocessing/dating/rsids/
# Removes duplicates
sort -u rsids | uniq > uniqrs.txt

# 
grep -w -f uniqrs.txt *.csv >> dated_uniq.txt
rm uniqrs.txt

#For each tissue
for i in hsap/*.rs;do
  awk '{print $2}' $i > $i.dat
  mv $i.dat postprocessing/dating/
done

cd postprocessing/dating
for file in *.dat; do
  sudo rename 's/v8.*.dat/v8.dat/' $file
done
  
for file in *.dat; do
  sort -u $file | uniq > $file.txt
  grep -w -f $file.txt *.csv >> $file.tiss
  rm $file.txt $file.tiss
done

```

```r
library(readr)
library(ggplot2)
library(viridis)
library(reshape2)
library(kSamples)

setwd("~/v8/version4/eQTL/postprocessing/dating/")

temp = list.files(pattern="*.tiss")
myfiles = lapply(temp,
                 read_table2,
                 col_names = FALSE)

#Prepare the data for plotting
dates <- NULL
dates <- lapply(myfiles, subset, X7 == "Combined,")
dates <- lapply(dates, '[', "X23")
dates <- melt(dates)

#Change the names of variables
dates[dates=="1"]<-"Adrenal gland"
dates[dates=="2"]<-"Amygdala"
dates[dates=="3"]<-"BA24"
dates[dates=="4"]<-"Caudate"
dates[dates=="5"]<-"Cereb. Hemisph."
dates[dates=="6"]<-"Cerebellum"
dates[dates=="7"]<-"Cortex"
dates[dates=="8"]<-"BA9"
dates[dates=="9"]<-"Hippocampus"
dates[dates=="10"]<-"Hypothalamus"
dates[dates=="11"]<-"Nucleus accumbens"
dates[dates=="12"]<-"Putamen"
dates[dates=="13"]<-"Spinal cord"
dates[dates=="14"]<-"Subs. Nigra"
dates[dates=="15"]<-"Pituitary"


ggplot(dates, aes(value*25, color = dates$L1))+
  facet_grid(L1 ~ ., scales = "free") +
  theme_minimal() +
  geom_freqpoly(stat = "bin", bins = 300) +
  theme(strip.text.y = element_blank()) +
  #stat_bin(bin = 40)+
  theme(legend.position="right") +
  #theme(legend.title = element_text("Tissue")) +
  labs(color = "Tissue") +
  labs(x = "Years", y = "Variant count")

```


## Heatmap plot
```r
library(readr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(magrittr)
library(dplyr)

#Change
setwd("~/v8/version4/eQTL/hsap")

temp = list.files(pattern="*.rs")
myfiles = lapply(temp,
                 read_table2,
                 col_names = FALSE)

#Tidy data + give it long format
data <- NULL
data <- lapply(myfiles, '[', "X2")
data <- melt(data)
data <- unique(data)

#Factors
data <- mutate(data, L1 = factor(L1, ordered=T))

#Gets unique elements for each tissue
test <- data %>%  
  full_join(data, by="X2") %>% 
  group_by(X2) %>% 
  filter(n() == 1) 
test <- test %>% 
  group_by(L1.x) %>% 
  summarize(n())

#Gets total number of eqtl, not counting duplicated variants
test2 <- data %>%  
  full_join(data, by="X2") %>% 
  group_by(X2)
test2 <- test2 %>% 
  group_by(L1.x) %>% 
  summarize(n())

#Then, computes a matrix adition 
#A measure of how different two tissues is would be here the sum of of the number of unique elements of every pairwise combination of tissues
#That's column sum
# Then, divide between the total of each two tissues to check the percentage of unique variants
#
ind <- t(combn(nrow(test),2))
out <- apply(ind, 1, function(x) sum(test[x[1], -1] + test[x[2], -1]))
out.total <- apply(ind, 1, function(x) sum(test2[x[1], -1] + test2[x[2], -1]))
what <- cbind(ind, out, out.total)
what <- as.data.frame(what)
what$out3 <- out/out.total
# The previous filter(n() == 1) gets rid also of the comparaison of tissue vs itself
# Since we are measuring number of unique variants per tissue, and all tissues share variants with themselves, this is NA
what <- add_row(what, V1 = 1:15, V2 = 1:15, out = 0)

#Labelling
what[what=="1"]<-"Adrenal gland"
what[what=="2"]<-"Amygdala"
what[what=="3"]<-"BA24"
what[what=="4"]<-"Caudate"
what[what=="5"]<-"Cereb. Hemisph."
what[what=="6"]<-"Cerebellum"
what[what=="7"]<-"Cortex"
what[what=="8"]<-"BA9"
what[what=="9"]<-"Hippocampus"
what[what=="10"]<-"Hypothalamus"
what[what=="11"]<-"Nucleus accumbens"
what[what=="12"]<-"Putamen"
what[what=="13"]<-"Spinal cord"
what[what=="14"]<-"Subs. Nigra"
what[what=="15"]<-"Pituitary"

ggplot(what, aes(x = V2, y = V1, fill = out3)) +
    theme_minimal() +
    geom_tile(colour="black",size=0.25) +
    labs(x="",y="") +
    geom_tile(aes(V1, V2),colour="black",size=0.25) +
    scale_fill_distiller(palette = "YlGnBu", na.value = 'white')  +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text( size=10, angle=90, hjust=0))
  
  

```

**DEADCODE**

```r
#ord <- hclust(dist(what))$order 
#what$V1 <- what[ord,] # re-order matrix accoring to clustering
#what$V1 <- as.factor(what$V1, levels = ord)
#what$V1 <- factor(what$V1, levels = ord)
#what$V1
#ord <- levels(what)
```

# GO analysis
With clumped data to avoid noise.

```sh
for fold in */; do
  for asd in $fold/*.clumpled; do
  nametxt="${fold%%[/]*}"
  awk 'FNR > 1 {print $3}' $fold/*.clumped > $nametxt.txt
  done
done

sed -i '/^$/d' *.txt
```

```r
library(gprofiler2)
library(readr)

setwd("~/v8/version4/eQTL/postprocessing/clumped/")

temp <- list.files(pattern="*.txt")
rsids <- lapply(temp,
                 read_table2,
                 col_names = FALSE)
results <- NULL
for (n in 1:15){
  results[n] <- gost(query = rsids[[n]][[1]],
     organism = "hsapiens", ordered_query = FALSE,
     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
     measure_underrepresentation = FALSE, evcodes = FALSE,
     user_threshold = 0.05, correction_method = "g_SCS",
     domain_scope = "annotated", custom_bg = NULL,
     numeric_ns = "", sources = NULL, as_short_link = FALSE)
}

```

#PCA analysis
```{bash}
#From eQTL
cd hsap
cat *.rs > ../all.txt
```

```r
library(devtools)
library(ggbiplot)
library(reshape2)
library(factoextra)
library(readr)
library(ggbiplot)


all <- read_table2("~/eQTL/all.rs") #Change if you set other path for the project

#Cleaning data
asd <- all[, -c(1:3) ] 

#PCA
asd.pca <- prcomp(asd, center = TRUE,scale. = TRUE)
summary(asd.pca)

#Habillage (grouping)
dup <- ifelse(duplicated(all$rs_id_dbSNP151_GRCh38p7), "Shared", "Unique")

#Plotting
fviz_pca_ind(asd.pca,
             geom = c("point"),
             repel = TRUE,
             select.ind = list(contrib = 10000),
             addEllipses = TRUE, ellipse.level=0.99,
             habillage = dup
              )


fviz_pca_var(asd.pca,
             geom = c("point", "text"),
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             title = "PCA - Variables"
)
```

#GWAS

## GWAS catalog
```sh
#From eQTL
cd postprocessing && mkdir gwas && cd gwas && mkdir gwascat gwasneale
cd gwascat && wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
mv alternative gwas.txt
```

```r
gwas <- read_delim("~/eQTL/postprocessing/gwas/gwascat/gwas.txt", 
+     "\t", escape_double = FALSE, trim_ws = TRUE)

intersect(rsids[[1]][[1]], gwas[[22]])
intersect(rsids[[2]][[1]], gwas[[22]])
intersect(rsids[[3]][[1]], gwas[[22]])
intersect(rsids[[4]][[1]], gwas[[22]])
intersect(rsids[[5]][[1]], gwas[[22]])
intersect(rsids[[6]][[1]], gwas[[22]])
intersect(rsids[[7]][[1]], gwas[[22]])
intersect(rsids[[8]][[1]], gwas[[22]])
intersect(rsids[[9]][[1]], gwas[[22]])
intersect(rsids[[10]][[1]], gwas[[22]])
intersect(rsids[[11]][[1]], gwas[[22]])
intersect(rsids[[12]][[1]], gwas[[22]])
intersect(rsids[[13]][[1]], gwas[[22]])
intersect(rsids[[14]][[1]], gwas[[22]])
intersect(rsids[[15]][[1]], gwas[[22]])
```


##?. Neale
Download all GWAS sumstat results from the Neale webpage: see **attached neale.sh** script.

```sh
cd ~/eQTL/postrocessing/gwasneale 
#Download and run neale.sh (not included here for readability reasons, but provided in the article's github page)

for i in *.bgz; do
  dir="${i%.*}"
  gunzip -c $i > $dir
done

```


```sh
# Positive selection

```
**ADD POSITIVE SELECTION TABLES**


```sh
#From eQTL
# `allrs` is derived from the PCA bash code lines. If you erase it after the analysis, please rerun that block of code.
awk '{print $1, $2, $3}' all.rs | sed 's/_/ /g' | awk 'NR>1 {print $1, $2-1, $2, $6, $7}' | sed 's/\s/\t/g' > coords.bed
# $2-1 due to the reasons explained here: https://bedtools.readthedocs.io/en/latest/content/overview.html#bed-starts-are-zero-based-and-bed-ends-are-one-based  

#Download the positive selections tables provided in the github repository

bedtools intersect -a coords.txt -b pey_coords.txt > ps_pey
bedtools intersect -a coords.bed -b racimo_coords.txt > ps_rac

```

Should I include this section?

```r
library(gprofiler2)
library(readr)

#Peyregne
ps_pey <- read_delim("~/v8/ps_pey", "\t", 
                     escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)

#Residing genes name retrieval
ps_pey$X5 <- gsub("[.]\\d*$", "", ps_pey$X5)
result1 <- gconvert(query = ps_pey$X5,
         organism = "hsapiens",
         target = "ENSG",
         filter_na = TRUE)
ps_pey$eGene <- c(result1$name)
    
result1 <- gconvert(query = ps_pey$X4,
                    organism = "hsapiens",
                    target = "ENSG",
                    filter_na = FALSE)
ps_pey$target <- c(result1$name)


#Racimo
ps_rac <- read_delim("~/v8/ps_rac", "\t", 
                     escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
ps_rac$X5 <- gsub("[.]\\d*$", "", ps_rac$X5)
result1 <- gconvert(query = ps_rac$X5,
                    organism = "hsapiens",
                    target = "ENSG",
                    filter_na = TRUE)
ps_rac$eGene <- c(result1$name)

results1 <- gconvert(query = ps_rac$X4,
         organism = "hsapiens",
         target = "ENSG",
         filter_na = TRUE)


ps_rac$target <- c(result1$name)
asd <- ifelse(ps_rac$X4 %in% results1$input, 
       results1$name, 
       NA)
ps_rac$target <- asd
```
