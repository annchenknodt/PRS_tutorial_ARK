## QC of Base Data
 
# Standard GWAS QC
gunzip -c Height.gwas.txt.gz |\
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > Height.gz

# Duplicate SNPs
gunzip -c Height.gz |\
awk '{ print $3}' |\
sort |\
uniq -d > duplicated.snp

gunzip -c Height.gz  |\
grep -vf duplicated.snp |\
gzip - > Height.nodup.gz

# Ambiguous SNPs
gunzip -c Height.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
       ($4=="T" && $5=="A") || \
       ($4=="G" && $5=="C") || \
       ($4=="C" && $5=="G")) {print}' |\
   gzip > Height.QC.gz

## QC of Target Data

# Standard GWAS QC
plink \
   --bfile EUR \
   --maf 0.01 \
   --hwe 1e-6 \
   --geno 0.01 \
   --mind 0.01 \
   --write-snplist \
   --make-just-fam \
   --out EUR.QC

plink \
   --bfile EUR \
   --keep EUR.QC.fam \
   --extract EUR.QC.snplist \
   --indep-pairwise 200 50 0.25 \
   --out EUR.QC

plink \
   --bfile EUR \
   --extract EUR.QC.prune.in \
   --keep EUR.QC.fam \
   --het \
   --out EUR.QC

Rscript prune_het.R

Rscript mismatchingSNPs.R

# Make a back up
mv EUR.bim EUR.bim.bk
cp EUR.QC.adj.bim EUR.bim

# Sex chromosomes
plink \
   --bfile EUR \
   --extract EUR.QC.prune.in \
   --keep EUR.valid.sample \
   --check-sex \
   --out EUR.QC

Rscript sexCheck.R

# Relatedness
plink \
   --bfile EUR \
   --extract EUR.QC.prune.in \
   --keep EUR.QC.valid \
   --rel-cutoff 0.125 \
   --out EUR.QC

# Generate final QC'd target data file
plink \
   --bfile EUR \
   --make-bed \
   --keep EUR.QC.rel.id \
   --out EUR.QC \
   --extract EUR.QC.snplist \
   --exclude EUR.mismatch

### ARK: also create a file where Step #3 (recoding mismatching SNPs in target file) is NOT performed
### do this by using the original EUR.bim file rather than EUR.QC.adj.bim
cp EUR.bim.bk EUR.bim

# Generate final QC'd target data file
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC.noRecode \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch

