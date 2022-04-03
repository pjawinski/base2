# ================================================================================================================
# === (1) Calculate correlations between brain age gap and outcome variables without age2 serving as covariate ===
# ================================================================================================================

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
library(dplyr)
library(ggplot2)
library(scales)
library(plotly)
library(ppcor)
library(reshape2)
library(tidyverse)
library(htmlwidgets)

# load data
df = read.delim('code/derivatives/03_brainAGE_w_phenotypes.txt', sep = '\t', header = TRUE)

# rename brain age gap variables for reasons of convenience
names(df)[names(df) == "brainage_gap_gm_stack"] = "gm"
names(df)[names(df) == "brainage_gap_wm_stack"] = "wm"
names(df)[names(df) == "brainage_gap_gwm_stack"] = "gwm"

# set independent variables and covariates
indep = c("gm","wm","gwm")
covs = c("sex", "age", "TIV")

# set dependent variables and expected effect directions
dep = as.data.frame(matrix(ncol = 3, byrow = T, data = c(
  'Educ_final', 'Years of education', 0,
  'hnetto', 'Household income', 0,
  'MMSE_Summe', 'Mini-mental state examination', 0,
  'GDS_Summe', 'Geriatric depression scale', 1,
  'CESD_Summe', 'CES-Depression', 1,
  'Rauchen_aktuell_inverted', 'Smoking status', 1,
  'Alkohol_haufigkeit', 'Frequency of alcohol intake', 1,
  'Alkohol_Menge', 'Amount of alcohol intake', 1,
  'Alkohol_6Glaser', 'Frequency of 6 glasses of alcohol intake', 1,
  'CH_Diabetes', 'Diabetes diagnosis', 1,
  'HOMAIR', 'HOMA-Insulin resistance', 1,
  'HbA1c', 'Hemoglobin A1c', 1,
  'BZP1', 'Fasting glucose', 1,
  'BZP2', 'Post-load glucose', 1,
  'BMI', 'Body mass index', 1,
  'RRdi', 'Diastolic blood pressure', 1,
  'RRsy', 'Systolic blood pressure', 1,
  'finalMetLscore', 'Metabolic load factor', 1,
  'GammaGTGGTUL', 'Gamma-glutamyl-transferase', 1,
  'HarnsaeuremgdL', 'Uric acid', 1,
  'TNF1', 'Tumor necrosis factor-alpha', 1,
  'DS2_corr', 'Digit symbol task', 0,
  'EM_final', 'Episodic memory', 0,
  'WM_final', 'Working memory', 0,
  'Gf_final', 'Fluid intelligence', 0,
  'futi_mean', 'Future time perspective', 0,
  'cfc_mean', 'Consideration of future consequences', 0)))

# set column names of dep and variable format
names(dep) = c('pheno', 'phenoname', 'hypothesis')
dep$hypothesis = as.numeric(dep$hypothesis)

##  run correlations
for (i in indep) {
  
  # create temporary data frame for results
  results = data.frame(matrix(NA, nrow = nrow(dep), ncol = 6))
  results = cbind(dep, results)
  names(results) = c(names(dep),c('estimate', 'p.value', 'statistic', 'n', 'gp', 'method'))
    
  # loop across outcome variables
  k = 0
  for (j in results$pheno) {
    k = k+1
    dftemp = df[!is.na(df[,i]) & !is.na(df[,j]),c(i,j,covs)]
    results[k,4:9] = pcor.test(dftemp[i],dftemp[,j],dftemp[,covs], method = 'pearson')
    results$p.value[k] = pt(results$statistic[k], results$n[k]-2-results$gp[k], lower=results$hypothesis[k]==0) # calculate one-tailed p-value
  }
  
  # rename variables and assign individual name for results data.frame
  names(results) = c(names(dep),paste(i, c('estimate', 'p.value', 'statistic', 'n', 'gp', 'method'), sep = "_"))
  assign(paste0("results_",i), results)
}

# merge results and check if n, gp and method is the same across brainAGE association results
results = cbind(results_gm, results_wm[4:9], results_gwm[4:9])
identical(results$gm_n,results$wm_n,results$gwm_n)
identical(results$gm_gp,results$wm_gp,results$gwm_gp)
sum(results$gm_method == results$wm_method) + sum(results$gm_method == results$gwm_method) == nrow(results)*2
   
# make results sparse
results = results[,c('pheno','phenoname','hypothesis', 'gm_n', 'gm_gp',
                     'gm_estimate','gm_statistic','gm_p.value',
                     'wm_estimate','wm_statistic','wm_p.value',
                     'gwm_estimate','gwm_statistic','gwm_p.value')]

names(results) = c('pheno','phenoname','hypothesis', 'n', 'ncovs',
                   'gm_estimate','gm_t.statistic','gm_p.value',
                   'wm_estimate','wm_t.statistic','wm_p.value',
                   'gwm_estimate','gwm_t.statistic','gwm_p.value')

# convert hypothesis variable to factor
results$hypothesis = factor(results$hypothesis, labels = c('-','+'))

# save results
write.table(results, 'code/tables/main_corr_age2dropped.txt', quote = FALSE, sep = '\t', row.names = FALSE)

# =====================================================================
# === Plot correlations between brain age gap and outcome variables ===
# =====================================================================

# load data
results = read.delim('code/tables/main_corr_age2dropped.txt', sep = '\t', header = TRUE)

# convert hypothesis-variable to numeric
results$hypothesis = as.numeric(as.factor(results$hypothesis))
results$hypothesis[results$hypothesis==1] = -1
results$hypothesis[results$hypothesis==2] = 1

# prepare matrices for ggplot data.frame
r = data.frame(matrix(data = c(results$gm_estimate, results$wm_estimate,  results$gwm_estimate), nrow = 3, byrow = TRUE))
p = data.frame(matrix(data = c(results$gm_p.value, results$wm_p.value,  results$gwm_p.value), nrow = 3, byrow = TRUE))
t = data.frame(matrix(data = c(results$gm_t.statistic, results$wm_t.statistic,  results$gwm_t.statistic), nrow = 3, byrow = TRUE))
n = data.frame(matrix(data = rep(results$n, 3), nrow = 3, byrow = TRUE))
h = data.frame(matrix(data = rep(results$hypothesis, 3), nrow = 3, byrow = TRUE))
  
colnames(r) = colnames(p) = colnames(t) = colnames(n) = colnames(h) = results$phenoname
r$brainAGE = p$brainAGE = t$brainAGE = n$brainAGE = h$brainAGE = factor(c('Grey matter brain age gap', 'White matter brain age gap', 'Grey & white matter brain age gap'))
  
# melt matrices
r.m = melt(r, id="brainAGE", variable.name="outcome", value.name="r", na.rm = F)
p.m = melt(p, id="brainAGE", variable.name="outcome", value.name="p", na.rm = F)
t.m = melt(t, id="brainAGE", variable.name="outcome", value.name="t", na.rm = F)
n.m = melt(n, id="brainAGE", variable.name="outcome", value.name="n", na.rm = F)
h.m = melt(h, id="brainAGE", variable.name="outcome", value.name="h", na.rm = F)
  
# combine matrices to ggplot data.frame
dfplot = cbind(r.m, p.m$p, t.m$t, n.m$n, h.m$h)
names(dfplot) = c(names(r.m),"p", "t", "n", "h")
dfplot$brainAGE = factor(dfplot$brainAGE, levels=unique(dfplot$brainAGE))

# color by t-value:
# set t-value to zero where associations are non-significant
# set t-value to NA where effect direction is not consistent with hypothesis
dfplot$tvalue = dfplot$t
for (i in 1:nrow(dfplot)) {
  if (dfplot$p[i] > 0.05) { dfplot$tvalue[i] = 0 }
  if (dfplot$h[i] != sign(dfplot$r[i])) { dfplot$tvalue[i] = NA }
}
dfplot$h = factor(dfplot$h, levels = c(-1, 1), labels = c("-", "+"))
  
# create group variable for facet grid
dfplot$group = c(rep("Replication", 63), rep("Cognition", 12), rep("TH", 6))
dfplot$group = factor(dfplot$group, levels = c("Replication","Cognition","TH"))
  
# set one-tailed level of significance
sign_threshold = qnorm(0.95)

# draw ggplot with beta values and t-value fill gradient
corrplot = ggplot(dfplot, aes(outcome, brainAGE, fill=tvalue)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%0.3f", round(r, 3)+0)), size=1.8, hjust = 1, nudge_x = 0.35) +
  theme_bw(base_size=10) +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 0),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background=element_rect(fill="white", size = 0.5),
        strip.text.x = element_text(size = 7, margin = margin(3, 0, 3, 0)),
        legend.position = "bottom",
        plot.margin = unit(c(7, 30, 1, 5), "mm")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(dfplot$brainAGE))) +
  facet_grid(.~group,scales = "free", space = "free") +
  scale_fill_gradientn(colours = c("#294979", "#D1F8FF", "#EBFFE5", "#EBFFE5", "#EBFFE5", "#FFDAD8", "#9E2013"), values = rescale(x = c(-4, -sign_threshold-1E-6, -sign_threshold, sign_threshold, sign_threshold+1E-6, 4), to = c(0, 1)),  na.value="white", guide = "colourbar", breaks = c(-4, -2, 0, 2, 4), limits = c(-4,4)) +
  guides(fill = guide_colourbar(title = 't-value',  title.position = "top", title.theme = element_text(size = 7, hjust = 0.5, vjust = 1), label.theme = element_text(size = 6, hjust = 1), barwidth = 30, barheight = 0.5, ticks = F, limits = c(-4,4)))
corrplot
  
# save corrplot as png
png('code/figures/main_corr_age2dropped.png', width=10.33, height=3.15, units = "in", res = 300)
corrplot
dev.off()
  
# save corrplot as pdf
pdf('code/figures/main_corr_age2dropped.pdf', width=10.33, height=3.15)
corrplot
dev.off()
 
# =========================================
# === make interactive plot with plotly ===
# =========================================

# prepare variables for ggplotly tooltip
dfplot$x = dfplot$outcome
dfplot$y = dfplot$brainAGE
dfplot$r.orig = dfplot$r
dfplot$r = sprintf("%.6f", dfplot$r)
dfplot$p.orig = dfplot$p
dfplot$p = sprintf("%.6f", dfplot$p.orig)
dfplot[dfplot$p.orig < 0.000001] = sprintf("%.2g", dfplot$p.orig[dfplot$p.orig < 0.000001])
dfplot$t.orig = dfplot$t
dfplot$t = sprintf("%.6f", dfplot$t.orig)
  
# draw ggplotly
corrplotly = dfplot %>%
  split(.$group) %>%
  map(function(x) {
    ggplotly(
      ggplot(x, aes(outcome, brainAGE, fill=tvalue, textx = x, texty = y, texth = h, textn = n, textr = r, textp = p, textt = t)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%0.3f", round(r.orig, 3)+0)), size=1.8, hjust = 1) +
      theme_bw(base_size=10) +
      theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background=element_rect(fill="white", size = 0.5),
            strip.text.x = element_text(size = 7, margin = margin(3, 0, 3, 0)),
            plot.margin = unit(c(7, 30, 1, 5), "mm")) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(limits = rev(levels(dfplot$brainAGE))) +
      facet_grid(.~group,scales = "free", space = "free") +
      scale_fill_gradientn(colours = c("#294979", "#D1F8FF", "#EBFFE5", "#EBFFE5", "#EBFFE5", "#FFDAD8", "#9E2013"), values = rescale(x = c(-4, -sign_threshold-1E-6, -sign_threshold, sign_threshold, sign_threshold+1E-6, 4), to = c(0, 1)),  na.value="white", guide = "colourbar", breaks = c(-4, -2, 0, 2, 4), limits = c(-4,4)) +
      guides(fill = guide_colourbar(title = 't-value',  title.position = "top", title.theme = element_text(size = 7, hjust = 0.5, vjust = 1), label.theme = element_text(size = 6, hjust = 1), barwidth = 30, barheight = 0.5, ticks = F, limits = c(-4,4))),
    tooltip = c("textx", "texty", "texth", "textn", "textr", "textp", "textt"))
  }) %>%
  subplot(margin = 0.005, widths = c(21/27,4/27,2/27)) %>%
  config(displayModeBar = F) %>% 
  layout(xaxis=list(side = "top", fixedrange=TRUE), xaxis2=list(side = "top", fixedrange=TRUE), xaxis3=list(side = "top", fixedrange=TRUE),
         yaxis=list(fixedrange=TRUE), yaxis2=list(ticktext = rep("",3), fixedrange=TRUE), yaxis3=list(ticktext = rep("",3), fixedrange=TRUE))
  
# adjust layout
corrplotly$width = 1000
corrplotly$height = 300
corrplotly$x$layout$margin$t = 180
corrplotly$x$layout$margin$b = 30
corrplotly$x$layout$xaxis$ticks = corrplotly$x$layout$xaxis2$ticks = corrplotly$x$layout$xaxis3$ticks = "outside"
corrplotly$x$layout$xaxis$ticklen = corrplotly$x$layout$xaxis2$ticklen = corrplotly$x$layout$xaxis3$ticklen = 15
corrplotly$x$layout$xaxis$tickcolor = corrplotly$x$layout$xaxis2$tickcolor = corrplotly$x$layout$xaxis3$tickcolor = "white"
corrplotly$x$data[[9]]$marker$colorbar$len = 1.2
corrplotly$x$data[[9]]$marker$colorbar$outlinecolor = "transparent"
corrplotly$x$data[[9]]$marker$colorbar$thickness = 15
  
# fix colorscale of LP subplot
corrplotly$x$data[[7]]$colorscale = data.frame(rbind(c(0,"#EBFFE5"),corrplotly$x$data[[7]]$colorscale, c(1,"#EBFFE5")))
names(corrplotly$x$data[[7]]$colorscale) = NULL
    
# draw corrplotly as html
corrplotly
saveWidget(corrplotly, paste0('code/figures/main_corr_age2dropped.html'), selfcontained = TRUE)
system(paste0('rm -rf code/figures/main_corr_age2dropped_files'))

# ==============================================================================
# === compare association results with vs. without age2 serving as covariate ===
# ==============================================================================

# load data
w.age2 = read.delim('code/tables/main_corr.txt', sep = '\t', header = TRUE)
wout.age2 = read.delim('code/tables/main_corr_age2dropped.txt', sep = '\t', header = TRUE)

# prepare dataframe
dfplot = data.frame(tissue = c(rep('grey matter', nrow(w.age2)), rep('white matter', nrow(w.age2)), rep('grey and white matter', nrow(w.age2))),
                    w.age2 = c(w.age2$gm_t.statistic, w.age2$wm_t.statistic, w.age2$gwm_t.statistic),
                    wout.age2 = c(wout.age2$gm_t.statistic, wout.age2$wm_t.statistic, wout.age2$gwm_t.statistic))
  
dfplot = data.frame(tissue = factor(c(rep(1, nrow(w.age2)), rep(2, nrow(w.age2)), rep(3, nrow(w.age2))), 
                                    levels = c(1:3), labels = c('grey matter', 'white matter', 'grey and white matter')),
                    w.age2 = c(w.age2$gm_t.statistic, w.age2$wm_t.statistic, w.age2$gwm_t.statistic),
                    wout.age2 = c(wout.age2$gm_t.statistic, wout.age2$wm_t.statistic, wout.age2$gwm_t.statistic))

# create plot
tvalComparison = ggplot(dfplot, aes(x = w.age2, y = wout.age2, color=tissue)) +
  geom_point(shape=19, alpha = 0.8, size = 0.8) +
  xlab("t-value\n(covariates: sex, age, age2, TIV)") +
  ylab("t-value\n(covariates: sex, age, and TIV)") +
  scale_color_manual(values=c("#f8766d", "#00ba38", "#619cff")) +
  theme_bw(base_size=10) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(5, 5, 5, 5), "mm"))

tvalComparison = ggplot(dfplot, aes(x = w.age2, y = wout.age2, color=tissue)) +
  geom_point(shape=19, alpha = 0.8, size = 0.8) +
  xlab("t-value\n(covariates: sex, age, age2, TIV)") +
  ylab("t-value\n(covariates: sex, age, and TIV)") +
  scale_x_continuous(expand = c(0,0), limits = c(-4,3.5), breaks = seq(-4,3,2)) + # label = axis.set$CHR % label = c(1:22, 'X', 'Y MT', '') label = c(1:18,'', 20, '', 22, 'X', 'Y MT', '')
  scale_y_continuous(expand = c(0,0), limits = c(-4,3.5), breaks = seq(-4,3,2)) +
  scale_color_manual(values=c("#f8766d", "#00ba38", "#619cff")) +
  geom_segment(aes(x=-4,xend=-4,y=-3,yend=3), colour = "black", size = 0.25) +
  geom_segment(aes(y=-4,yend=-4,x=-3,xend=3), colour = "black", size = 0.25) +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.ticks.length=unit(.15, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 10, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
    axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 5)),
    axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    plot.margin=unit(c(0.25,1,-0.5,0),"cm"))

# save plot
png('code/figures/tvalComparison_age2.png', width=5.6, height=3.7, units = "in", res = 300)
tvalComparison
dev.off()


# ========================================================================================
# === do educational attainment and income independently correlate with brain age gap? ===
# ========================================================================================

# load data
df = read.delim('code/derivatives/03_brainAGE_w_phenotypes.txt', sep = '\t', header = TRUE)

# rename brain age gap variables for reasons of convenience
names(df)[names(df) == "brainage_gap_gm_stack"] = "gm"
names(df)[names(df) == "brainage_gap_wm_stack"] = "wm"
names(df)[names(df) == "brainage_gap_gwm_stack"] = "gwm"

# set independent variables and covariates
indep = c("gm","wm","gwm")
covs = c("sex", "age", "age2", "TIV")

# set dependent variables and expected effect directions
dep = as.data.frame(matrix(ncol = 3, byrow = T, data = c(
  'Educ_final', 'Years of education', 0,
  'hnetto', 'Household income', 0)))

# set column names of dep and variable format
names(dep) = c('pheno', 'phenoname', 'hypothesis')
dep$hypothesis = as.numeric(dep$hypothesis)

##  run correlations
for (i in indep) {
  
  # create temporary data frame for results
  results = data.frame(matrix(NA, nrow = nrow(dep), ncol = 12))
  results = cbind(dep, results)
  names(results) = c(names(dep),c('estimate', 'p.value', 'statistic', 'n', 'gp', 'method', 'uniq.estimate', 'uniq.p.value', 'uniq.statistic', 'uniq.n', 'uniq.gp', 'uniq.method'))
  
  # loop across outcome variables
  k = 0
  for (j in results$pheno) {
    k = k+1
    dftemp = df[!is.na(df[,i]) & !is.na(df[,j]),c(i,j,covs)]
    results[k,4:9] = pcor.test(dftemp[i],dftemp[,j],dftemp[,covs], method = 'pearson')
    results$p.value[k] = pt(results$statistic[k], results$n[k]-2-results$gp[k], lower=results$hypothesis[k]==0) # calculate one-tailed p-value
    
    # now calculate unique association by partialling out effects of all other variables in dependent variable list
    dftemp = df[complete.cases(df[,c(i,dep$pheno,covs)]),c(i,dep$pheno,covs)]
    results[k,10:15] = pcor.test(dftemp[i],dftemp[,j],dftemp[,-which(names(dftemp) %in% c(i,j))], method = 'pearson')
    results$uniq.p.value[k] = pt(results$uniq.statistic[k], results$uniq.n[k]-2-results$uniq.gp[k], lower=results$hypothesis[k]==0) # calculate one-tailed p-value
  }
  
  # rename variables and assign individual name for results data.frame
  names(results) = c(names(dep),paste(i, c('estimate', 'p.value', 'statistic', 'n', 'gp', 'method', 'uniq.estimate', 'uniq.p.value', 'uniq.statistic', 'uniq.n', 'uniq.gp', 'uniq.method'), sep = "_"))
  assign(paste0("results_",i), results)
}

# merge results and check if n, gp and method is the same across brainAGE association results
results = cbind(results_gm, results_wm[4:15], results_gwm[4:15])
identical(results$gm_n,results$wm_n,results$gwm_n)
identical(results$gm_gp,results$wm_gp,results$gwm_gp)
sum(results$gm_method == results$wm_method) + sum(results$gm_method == results$gwm_method) == nrow(results)*2

# make results sparse
results = results[,c('pheno','phenoname','hypothesis', 'gm_n', 'gm_gp','gm_uniq.n', 'gm_uniq.gp',
                     'gm_estimate','gm_statistic','gm_p.value',
                     'wm_estimate','wm_statistic','wm_p.value',
                     'gwm_estimate','gwm_statistic','gwm_p.value',
                     'gm_uniq.estimate','gm_uniq.statistic','gm_uniq.p.value',
                     'wm_uniq.estimate','wm_uniq.statistic','wm_uniq.p.value',
                     'gwm_uniq.estimate','gwm_uniq.statistic','gwm_uniq.p.value')]

names(results) = c('pheno','phenoname','hypothesis', 'n', 'ncovs','uniq.n','uniq.ncovs',
                   'gm_estimate','gm_statistic','gm_p.value',
                   'wm_estimate','wm_statistic','wm_p.value',
                   'gwm_estimate','gwm_statistic','gwm_p.value',
                   'gm_uniq.estimate','gm_uniq.statistic','gm_uniq.p.value',
                   'wm_uniq.estimate','wm_uniq.statistic','wm_uniq.p.value',
                   'gwm_uniq.estimate','gwm_uniq.statistic','gwm_uniq.p.value')

# convert hypothesis variable to factor
output = data.frame(pheno = rep(results$pheno,2),
                    phenoname = rep(results$phenoname,2),
                    hypothesis = rep(results$hypothesis,2),
                    n = c(results$n, results$uniq.n),
                    ncovs = c(results$ncovs, results$uniq.ncovs),
                    gm_estimate = c(results$gm_estimate, results$gm_uniq.estimate),
                    gm_statistic = c(results$gm_statistic, results$gm_uniq.statistic),
                    gm_p.value = c(results$gm_p.value, results$gm_uniq.p.value),
                    wm_estimate = c(results$wm_estimate, results$wm_uniq.estimate),
                    wm_statistic = c(results$wm_statistic, results$wm_uniq.statistic),
                    wm_p.value = c(results$wm_p.value, results$wm_uniq.p.value),                  
                    gwm_estimate = c(results$gwm_estimate, results$gwm_uniq.estimate),
                    gwm_statistic = c(results$gwm_statistic, results$gwm_uniq.statistic),
                    gwm_p.value = c(results$gwm_p.value, results$gwm_uniq.p.value))
                    
                    
output$hypothesis = factor(output$hypothesis, levels = c(0,1), labels = c('-','+'))

# save results
write.table(output, 'code/tables/expl_income_education.txt', quote = FALSE, sep = '\t', row.names = FALSE)






