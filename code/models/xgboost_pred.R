#!/usr/bin/env Rscript

# set base2 project folder as working directory
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename = dirname(script.name)
setwd(script.basename)

# activate R environment
setwd('../../')
source("renv/activate.R")
renv::restore(prompt = FALSE)
setwd(script.basename)

# attach packages to current R session
require(xgboost)

# ===================
# === predict age ===
# ===================

# --- Grey matter ---
# load model and data
load('xgb_gm_models.RData')
Test_Samples_pca = read.table('xgb_gm_Test_Samples_pca.txt', sep='\t')

# predict
dtest = xgb.DMatrix(data = as.matrix(Test_Samples_pca), label=as.matrix(c(1:length(Test_Samples_pca[,1]))))
pred = data.frame(cbind(predict(model_tree, dtest), predict(model_linear, dtest)))

# --- White matter ---
# load model and data
load('xgb_wm_models.RData')
Test_Samples_pca = read.table('xgb_wm_Test_Samples_pca.txt', sep='\t')

# predict
dtest = xgb.DMatrix(data = as.matrix(Test_Samples_pca), label=as.matrix(c(1:length(Test_Samples_pca[,1]))))
pred = data.frame(cbind(pred, predict(model_tree, dtest), predict(model_linear, dtest)))

# save models and predictions
write.table(pred, 'xgb_pred.txt', sep = '\t', row.names = F, col.names = F)
message('Age prediction completed.')
