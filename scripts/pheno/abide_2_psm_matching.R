library(MatchIt)

root_p = '/home/surchs/temp/abide_univariate_data/pheno/'
qc_table_p = file.path(root_p, 'ABIDE2_Pheno_QC_pass.tsv')
psm_table_p = file.path(root_p, 'ABIDE2_Pheno_PSM_matched.tsv')
qc_data = read.table(file = qc_table_p, sep = '\t', header = TRUE)
# DX_GROUP==2 -> Control, ==1 -> Autism
# Turn the controls into zeros, leave the autsists as they
qc_data$DX_NUM[qc_data$DX_GROUP == 'Control'] <- 0
qc_data$DX_NUM[qc_data$DX_GROUP == 'Autism'] <- 1
covariates = c('SITE_ID','SUB_ID', 'DX_NUM', 'AGE_AT_SCAN', 'fd_scrubbed')

sites <- unique(qc_data$SITE_ID)
count = 0
for (site in sites) {
  slice_data = subset(qc_data, SITE_ID==site, select=covariates)
  matched = matchit(DX_NUM ~ AGE_AT_SCAN + fd_scrubbed, data=slice_data, method="nearest", caliper=0.8, ratio=1)
  temp_matched = match.data(matched)
  if (count == 0) {
    data_matched = temp_matched
  }
  else {
    data_matched = rbind(data_matched, temp_matched)
  }
  count = count + 1
}
read.table
rownames(data_matched) <- NULL
# Rejoin the PSM data with the remaining phenotype file
psm_data <- merge(x = qc_data, y = data_matched[,c('SUB_ID', 'distance')], by = "SUB_ID")
write.table(psm_data, file = psm_table_p, sep="\t")