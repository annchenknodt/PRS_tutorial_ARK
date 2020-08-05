# following tutorial at choishingwan.github.io/PRS-Tutorial/target using data from the tutorial
# comparing PRS in plink & PRSice with and without the recoding of mismatching SNPs step
# 7/31/20 ARK

workdir <- "H:/Projects/Annchen/DNS/Genetic/PRScheck_Kelly/tutorial/"

## load data
prs.plink <- read.table(paste0(workdir, "EUR.0.5.profile"), header=T)
prs.plink.noRecode <- read.table(paste0(workdir, "EUR.noRecode.0.5.profile"), header=T)
prs.PRSice <- read.table(paste0(workdir, "EUR.all_score"), header=T)
prs.PRSice.noRecode <- read.table(paste0(workdir, "EUR.noRecode.all_score"), header=T)

## check that IDs match and merge lazy way
table(prs.PRSice$IID==prs.plink$IID)
table(prs.PRSice$IID==prs.plink.noRecode$IID)
table(prs.PRSice$IID==prs.PRSice.noRecode$IID)

prs_merged <- data.frame(plink=prs.plink$SCORE, plink.noRecode=prs.plink.noRecode$SCORE, PRSice=prs.PRSice$Pt_0.5, PRSice.noRecode=prs.PRSice.noRecode$Pt_0.5)

cor(prs_merged)

# plot(prs.plink$SCORE, prs.plink.noRecode$SCORE) 
