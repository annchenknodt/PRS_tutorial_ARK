# Update Effect Size
Rscript effectSize.R

# Clumping
### ARK, took out clump-field argument bc it was causing plink to not recognize P header
plink \
   --bfile EUR.QC \
   --clump-p1 1 \
   --clump-r2 0.1 \
   --clump-kb 250 \
   --clump Height.QC.Transformed \
   --clump-snp-field SNP \
   --out EUR

awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp

# Generate PRS
awk '{print $3,$8}' Height.QC.Transformed > SNP.pvalue

echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

plink \
   --bfile EUR.QC \
   --score Height.QC.Transformed 3 4 9 header \
   --q-score-range range_list SNP.pvalue \
   --extract EUR.valid.snp \
   --out EUR
### ARK, using column 9 instead of 12 for OR column since that's how it comes out in the sample data

# Accounting for Population Stratification
# First, we need to perform prunning
plink \
   --bfile EUR.QC \
   --indep-pairwise 200 50 0.25 \
   --out EUR
# Then we calculate the first 6 PCs
plink \
   --bfile EUR.QC \
   --extract EUR.prune.in \
   --pca 6 \
   --out EUR

# Finding the "best-fit" PRS
Rscript bestFit.R

#### ARK, now calculate score with non-recoded bim file
plink \
    --bfile EUR.QC.noRecode \
    --score Height.QC.Transformed 3 4 9 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR.noRecode

