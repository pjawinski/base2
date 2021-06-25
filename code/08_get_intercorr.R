# ============================================================
# === get bivariate correlations between outcome variables ===
# ============================================================

# set working directory
setwd('/Users/philippe/Desktop/base2')

# activate R environment
source("renv/activate.R")
renv::restore(prompt = FALSE)

# attach packages to current R session
library(htmlwidgets)
library(ggplot2)
library(plotly)
library(poolr)
library(ppcor)
library(reshape2)
library(scales)

# load data
df = read.delim('data/03_brainage_w_phenotypes.txt', sep = '\t', header = TRUE)

# bias-correction of brain age gap variables
df$gmres = resid(lm(brainage_gap_gm_stack ~ sex + age + age2 + TIV, data = df))
df$wmres = resid(lm(brainage_gap_wm_stack ~ sex + age + age2 + TIV, data = df))
df$gwmres = resid(lm(brainage_gap_gwm_stack ~ sex + age + age2 + TIV, data = df))

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
  'GammaGTGGTUL', 'Gamma-glutamyltransferase',
  'HarnsaeuremgdL', 'Urea',
  'TNF1', 'Tumor necrosis factor-alpha',
  'DS2_corr', 'Digit symbol task',
  'EM_final', 'Episodic memory',
  'WM_final', 'Working memory',
  'Gf_final', 'Fluid intelligence',
  'futi_mean', 'Future time perspective',
  'cfc_mean', 'Consideration of future consequences'), ncol = 2, byrow = T))
names(varint) = c('pheno', 'phenoname')

## define function for correlation matrix
makecorr = function(df, varnames, varlabels, covnames = NULL, type) {
  
  # prepare data.frame
  if (!is.null(covnames)) { covs = df[,covnames] }
  df = df[,varnames]
  names(df) = varlabels
  
  # make empty data.frames for rho and pval
  df.rho = df.pval = df.n = data.frame(matrix(data = NA, nrow = length(df), ncol = length(df)))
  names(df.rho) = names(df.pval) = names(df.n) = names(df)
  
  # do calculations
  for (i in 1:(length(df)-1)) {
    for (j in (i+1):length(df)) {
      temp.df = df[!is.na(df[,i]) & !is.na(df[,j]),c(i,j)]
      
      if (is.null(covnames)) { 
        temp.corr = cor.test(temp.df[,1], temp.df[,2], method = type)
      } else {
        temp.corr = pcor.test(temp.df[,1], temp.df[,2], covs[!is.na(df[,i]) & !is.na(df[,j]),], method = type)
      }
      
      df.rho[i,j] = df.rho[j,i] = temp.corr$estimate
      df.pval[i,j] = df.pval[j,i] = temp.corr$p.value
      df.n[i,j] = df.n[j,i] = nrow(temp.df)
      
    }
  }
  
  # fill rho diagonal with ones and p diagonal with zeros
  for (i in 1:(length(df))) {
    df.rho[i,i] = 1 
    df.pval[i,i] = 0
    df.n[i,i] = sum(!is.na(df[,i]))
  }
  
  # add row.names
  df.rho$id = df.pval$id = df.n$id = names(df)
  df.rho = df.rho[,c(length(df.rho),1:(length(df.rho)-1))]
  df.pval = df.pval[,c(length(df.pval),1:(length(df.pval)-1))]
  df.n = df.n[,c(length(df.n),1:(length(df.n)-1))]
  
  # format rho matrix for ggplot
  df.m = melt(df.rho, id="id", variable.name="id_y", value.name="rho", na.rm = F)
  df.m$id = as.character(df.m$id)
  df.m$id = factor(df.m$id, levels=unique(df.m$id))
  df.m$rho = df.m$rho_4color = as.numeric(df.m$rho)
  
  # use rho as fill gradient - do not color if corresponding Bonferroni-corrected p value > 0.05
  df.m$pval = as.numeric(reshape2::melt(df.pval, id="id", variable.name="id_y", value.name="pval", na.rm = F)$pval)
  # df.m$rho_4color[df.m$pval > 0.05/(length(df)*(length(df)-1)/2)] = NA # Bonferroni
  df.m$rho_4color[df.m$pval > 0.05] = NA # Bonferroni
  
  # add variables for ggplotly tooltip
  df.m$x = df.m$id
  df.m$y = df.m$id_y
  df.m$n = as.numeric(reshape2::melt(df.n, id="id", variable.name="id_y", value.name="n", na.rm = F)$n)
  df.m$p = sprintf("%.6f", as.numeric(df.m$pval))
  df.m$p[as.numeric(df.m$pval) < 0.000001] = sprintf("%.2g", as.numeric(df.m$pval[as.numeric(df.m$pval) < 0.000001]))
  df.m$r = sprintf("%.6f", as.numeric(df.m$rho))
  
  # draw plots
  for (type in c('plotly', 'ggplot')) {
    if (type == 'plotly') { nudgex = 0 } else { nudgex = 0.35 }
    
    tempplot = ggplot(df.m, aes(id_y, id, fill=rho_4color, textx = x, texty = y, textn = n, textr = r, textp = p)) + 
      geom_tile() + 
      geom_text(data = df.m, aes(label = sprintf("%0.2f", round(rho, 2))), size=1.8, hjust = 1, nudge_x = nudgex) +
      theme_bw(base_size=10) + 
      theme(legend.position="right",
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 6, angle = 40, hjust = 0),
            axis.text.y = element_text(size = 6),
            panel.border = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),   
            plot.margin = unit(c(5, 25, 5, 7), "mm")) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(limits = rev(levels(df.m$id))) +
      scale_fill_gradientn(colours = c("#9E2013", "white", "#294979"), values = rescale(x = c(-1, 0, 1), to = c(0, 1)),  na.value="white", guide = "colourbar", breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1,1)) +
      guides(fill = guide_colourbar(title = 'rho', direction = 'vertical', title.position = 'top', title.theme = element_text(size = 7, hjust = 0), label.theme = element_text(size = 6, hjust = 0.5), barwidth = 0.5, barheight = 10, ticks = F, limits = c(-1,1)))
    
    if (type == 'plotly') { 
      tempplot = ggplotly(tempplot, tooltip = c("textx", "texty", "textn", "textr", "textp", "textn"), dynamicTicks = FALSE, width = 1100, height = 700)
      tempplot = tempplot %>% config(displayModeBar = F) %>% layout(xaxis=list(side = "top", fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE))
      tempplot$x$data[[3]]$marker$colorbar$outlinecolor = "transparent"
    }
    
    assign(paste0('df.', type), tempplot)
  }
  
  # return results (rho + pval + n + plot)
  results = list("rho" = df.rho, "pval" = df.pval, "n" = df.n, "ggplot" = df.ggplot, "plotly" = df.plotly)
  results
}

# hijack function solve.default and set tolerance = 0 (multicollinearity between covariates age and age2 results in singularity error)
trace("solve.default", tracer = quote(tol <- 1E-20), print = FALSE) # alternatively use trace(base::solve.default,edit=TRUE) and set tol = 0

# calculate correlations
results = makecorr(df = df, varnames = varint$pheno[4:nrow(varint)], varlabels = varint$phenoname[4:nrow(varint)], type = 'pearson')
results_partial = makecorr(df = df, varnames = varint$pheno[4:nrow(varint)], varlabels = varint$phenoname[4:nrow(varint)], covnames = c('sex','age','age2', 'TIV'), type = 'pearson')

# draw ggplots as png
png(width = 10.25, height = 5.97, units = "in", res = 300, filename = 'code/figures/intercorr_outcome.png')
results$ggplot
dev.off()

png(width = 10.25, height = 5.97, units = "in", res = 300, filename = 'code/figures/intercorr_outcome_partial.png')
results_partial$ggplot
dev.off()

# draw ggplots as pdf
pdf(width = 10.25, height = 5.97, file = 'code/figures/intercorr_outcome.pdf')
results$ggplot
dev.off()

pdf(width = 10.25, height = 5.97, file = 'code/figures/intercorr_outcome_partial.pdf')
results_partial$ggplot
dev.off()

# draw plotly as html
saveWidget(results$plotly, 'code/figures/intercorr_outcome.html', selfcontained = TRUE)
system('rm -rf code/figures/intercorr_outcome_files')

saveWidget(results_partial$plotly, 'code/figures/intercorr_outcome_partial.html', selfcontained = TRUE)
system('rm -rf code/figures/intercorr_outcome_partial_files')

# calculate effective number of tests - outcome variables
results_partial = makecorr(df = df, varnames = varint$pheno[4:nrow(varint)], varlabels = varint$phenoname[4:nrow(varint)], covnames = c('sex','age','age2', 'TIV'), type = 'pearson')
meff(results_partial$rho[,2:ncol(results$rho)], method = "liji")

# calculate effective number of tests - brain age gap
results_partial = makecorr(df = df, varnames = c('brainage_gap_gm_stack', 'brainage_gap_wm_stack', 'brainage_gap_gwm_stack'), varlabels = c('grey matter', 'white matter', 'grey + white matter'), covnames = c('sex','age','age2', 'TIV'), type = 'pearson')
meff(results_partial$rho[,2:ncol(results_partial$rho)], method = "liji")

# save table of brain age gap intercorrelations
gap = cbind(results_partial$rho, results_partial$pval, results_partial$n)
gap = gap[,-c(5,9)]
names(gap) = c('var','gm_rho','wm_rho','gwm_rho','gm_pval','wm_pval','gwm_pval','gm_n','wm_n','gwm_n')
write.table(gap, 'code/tables/intercorr_gap.txt', sep = '\t', quote = FALSE, row.names = FALSE)

