# ========================
# === get descriptives ===
# ========================

# set working directory
setwd('/users/philippe/desktop/projects/base2')

# detach 'other packages' if there are any
if (!is.null(names(sessionInfo()$otherPkgs))) {
  invisible(lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE))
}

# activate R environment
if (exists('.rs.restartR', mode = 'function')) { .rs.restartR() }
source('renv/activate.R')
renv::activate(getwd())
renv::restore(prompt = FALSE)

# attach packages to current R session
library(psych)
library(dplyr)

# load data
df = read.delim('code/derivatives/03_brainage_w_phenotypes.txt', sep = '\t', header = TRUE)

# bias-correction of brain age gap variables
df$gmres = resid(lm(brainage_gap_gm_stack ~ sex + age + age2 + TIV, data = df))
df$wmres = resid(lm(brainage_gap_wm_stack ~ sex + age + age2 + TIV, data = df))
df$gwmres = resid(lm(brainage_gap_gwm_stack ~ sex + age + age2 + TIV, data = df))

# age and sex
describe(df$age, type = 2)
table(df$sex)

# define vars of interest
varint = as.data.frame(matrix(c(
  'gmres', 'Grey matter',
  'wmres', 'White matter',
  'gwmres', 'Grey and white matter',
  'Educ_final', 'Years of education',
  'hnetto', 'Household income',
  'MMSE_Summe', 'Mini-mental state examination',
  'GDS_Summe', 'Geriatric depression scale',
  'CESD_Summe', 'CES-Depression',
  'Rauchen_aktuell_inverted', 'Smoking status',
  'Alkohol_haufigkeit', 'Frequency of alcohol intake',
  'Alkohol_Menge', 'Amount of alcohol intake',
  'Alkohol_6Glaser', 'Frequency of 6 glasses of alcohol intake',
  'CH_Diabetes', 'Diabetes diagnosis',
  'HOMAIR', 'HOMA-Insulin resistance',
  'HbA1c', 'HbA1c',
  'BZP1', 'Fasting glucose',
  'BZP2', 'Post-load glucose',
  'BMI', 'Body mass index',
  'RRdi', 'Diastolic blood pressure',
  'RRsy', 'Systolic blood pressure',
  'finalMetLscore', 'Metabolic load factor',
  'GammaGTGGTUL', 'Gamma-glutamyl-transferase',
  'HarnsaeuremgdL', 'Uric acid',
  'TNF1', 'Tumor necrosis factor-alpha',
  'DS2_corr', 'Digit symbol substitution test',
  'EM_final', 'Episodic memory',
  'WM_final', 'Working memory',
  'Gf_final', 'Fluid intelligence',
  'futi_mean', 'Future time perspective',
  'cfc_mean', 'Consideration of future consequences'), ncol = 2, byrow = T))
names(varint) = c('pheno', 'phenoname')
  
# get descriptive statistics
descr = describe(df[,varint$pheno], type = 2)
descr = cbind(descr, as.data.frame(t(sapply(df[,varint$pheno], quantile, na.rm = TRUE))))
descr$pheno = rownames(descr)
descr = descr[,c(ncol(descr),1:(ncol(descr)-1))]
descr = descr[,c('pheno', 'n', 'mean', 'sd', 'min', 'max', '25%', '50%', '75%', 'skew', 'kurtosis')]
names(descr) = c('pheno', 'n', 'mean', 'sd', 'min', 'max', 'Q1', 'Q2', 'Q3', 'skew', 'kurtosis')

# get levels
descr$levels = NA
descr$levels[descr$pheno == 'Rauchen_aktuell_inverted'] = paste0("‘never’: ", sum(df$Rauchen_aktuell_inverted == 0, na.rm = TRUE), ", ‘stopped more than a year ago’: ", sum(df$Rauchen_aktuell_inverted == 1, na.rm = TRUE), ", ‘stopped less than a year ago’: ", sum(df$Rauchen_aktuell_inverted == 2, na.rm = TRUE), ", ‘current smoker’: ", sum(df$Rauchen_aktuell_inverted == 3, na.rm = TRUE))
descr$levels[descr$pheno == 'Alkohol_haufigkeit'] = paste0("‘never’: ", sum(df$Alkohol_haufigkeit == 0, na.rm = TRUE), ", ‘once a month or less’: ", sum(df$Alkohol_haufigkeit == 1, na.rm = TRUE), ", ‘two to four times a month’: ", sum(df$Alkohol_haufigkeit == 2, na.rm = TRUE), ", ‘two to four times a week’: ", sum(df$Alkohol_haufigkeit == 3, na.rm = TRUE), ", ‘four times a week or more’: ", sum(df$Alkohol_haufigkeit == 4, na.rm = TRUE))
descr$levels[descr$pheno == 'Alkohol_Menge'] = paste0("‘one to two glasses’: ", sum(df$Alkohol_Menge == 0, na.rm = TRUE), ", ‘three to four glasses’: ", sum(df$Alkohol_Menge == 1, na.rm = TRUE), ", ‘five to six glasses’: ", sum(df$Alkohol_Menge == 2, na.rm = TRUE), ", ‘seven to nine glasses’: ", sum(df$Alkohol_Menge == 3, na.rm = TRUE), ", ‘ten or more glasses’: ", sum(df$Alkohol_Menge == 4, na.rm = TRUE))
descr$levels[descr$pheno == 'Alkohol_6Glaser'] = paste0("‘never’: ", sum(df$Alkohol_6Glaser == 0, na.rm = TRUE), ", ‘less than once a month’: ", sum(df$Alkohol_6Glaser == 1, na.rm = TRUE), ", ‘once a month’: ", sum(df$Alkohol_6Glaser == 2, na.rm = TRUE), ", ‘once a week’: ", sum(df$Alkohol_6Glaser == 3, na.rm = TRUE), ", ‘daily or almost daily’: ", sum(df$Alkohol_6Glaser == 4, na.rm = TRUE))
descr$levels[descr$pheno == 'CH_Diabetes'] = paste0("controls: ", sum(df$CH_Diabetes == 0, na.rm = TRUE), ", cases: ", sum(df$CH_Diabetes == 1, na.rm = TRUE))

# merge with varint
descr = left_join(varint, descr, 'pheno')

# save
write.table(descr, 'code/tables/descriptives.txt', sep = '\t', quote = FALSE, row.names = FALSE)

