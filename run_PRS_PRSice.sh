# # Running PRS analysis
# Rscript format_cov.R

# # (ARK added no-regress)
# Rscript /cifs/hariri-long/Scripts/Tools/PRSice_linux/v2.3.2/PRSice.R \
    # --prsice /cifs/hariri-long/Scripts/Tools/PRSice_linux/v2.3.2/PRSice_linux \
    # --no-regress \
    # --base Height.QC.gz \
    # --target EUR.QC \
    # --binary-target F \
    # --pheno EUR.height \
    # --cov EUR.cov \
    # --base-maf MAF:0.01 \
    # --base-info INFO:0.8 \
    # --stat OR \
    # --or \
    # --out EUR

### un-recoded 
# (ARK added no-regress)
Rscript /cifs/hariri-long/Scripts/Tools/PRSice_linux/v2.3.2/PRSice.R \
    --prsice /cifs/hariri-long/Scripts/Tools/PRSice_linux/v2.3.2/PRSice_linux \
    --no-regress \
    --base Height.QC.gz \
    --target EUR.QC.noRecode \
    --binary-target F \
    --pheno EUR.height \
    --cov EUR.cov \
    --base-maf MAF:0.01 \
    --base-info INFO:0.8 \
    --stat OR \
    --or \
    --out EUR.noRecode
