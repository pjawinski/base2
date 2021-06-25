# =======================================================================
# === Compare base2 and matched UKB subset regarding model accuracies ===
# =======================================================================

# --------------------
# --- prepare data ---
# --------------------

# set working directory
setwd('/Users/philippe/Desktop/base2')

# activate R environment
source("renv/activate.R")
renv::restore(prompt = FALSE)

# attach packages to current R session
library(boot)
library(cocor)
library(dplyr)
library(ggplot2)
library(MatchIt)
library(psych)

# load base2 dataset
base2 = read.delim('data/03_brainage.txt', header = TRUE, sep  = '\t')
base2 = base2[,c('Row', 'sex', 'age', 'TIV', 'brainage_gm_stack', 'brainage_wm_stack', 'brainage_gwm_stack')]

# load UKB data
for (i in c('gm', 'wm', 'gwm')) {
  message(paste0('Starting with ', i))
  df = read.delim(paste0('data/ukb/gap_', i, '_stack.txt'), header = TRUE, sep = '\t')
  covs = read.delim(paste0('data/ukb/covs.txt'), header = TRUE, sep = '\t')
  df = left_join(df, covs, by = 'IID')
  df = df[, c('IID','sex', 'age', 'TIV', paste0('brainage_gap_',i,'_stk'))]
  df$sex = df$sex*-1+3
  assign(paste0(i),df)
}

wm = wm[,c('IID', 'brainage_gap_wm_stk')]
gwm = gwm[,c('IID', 'brainage_gap_gwm_stk')]

ukb = left_join(gm, wm, by = 'IID')
ukb = left_join(ukb, gwm, by = 'IID')
ukb$brainage_gap_gm_stk = ukb$brainage_gap_gm_stk + ukb$age
ukb$brainage_gap_wm_stk = ukb$brainage_gap_wm_stk + ukb$age
ukb$brainage_gap_gwm_stk = ukb$brainage_gap_gwm_stk + ukb$age

# put into one data.frame
base2$group = 1
ukb$group = 0

names(base2) = c('IID', 'sex', 'age', 'TIV', 'gm', 'wm', 'gwm', 'group')
names(ukb) = c('IID', 'sex', 'age', 'TIV', 'gm', 'wm', 'gwm', 'group')
df = rbind(base2, ukb)
df = df[, c('IID', 'group', 'sex', 'age', 'TIV', 'gm', 'wm', 'gwm')]

# ---------------
# --- MatchIt ---
# ---------------

# match UKB to base2 individuals
set.seed(7339)
m.out = matchit(group ~ age + sex, data = df, distance = "mahalanobis", ratio = 10)
df.matched = match.data(m.out)
df.matched$group = factor(df.matched$group, labels = c('ukb', 'base2'))

# compare min, max, and mean - samples are well-matched
table(df.matched$sex[df.matched$group == 'base2']) 
table(df.matched$sex[df.matched$group == 'ukb']) 

min(df.matched$age[df.matched$group=='base2'])
min(df.matched$age[df.matched$group=='ukb'])
max(df.matched$age[df.matched$group=='base2'])
max(df.matched$age[df.matched$group=='ukb'])
mean(df.matched$age[df.matched$group == 'base2'])
mean(df.matched$age[df.matched$group == 'ukb'])

# show density plot of age - samples are well-matched
ggplot(df.matched, aes(age, fill = group)) + geom_density(alpha = 0.2)

# compare ukb and base2 regarding sex and age
chisq.test(df.matched$sex, df.matched$group)
chisq.test(df.matched$sex, df.matched$group)$statistic[[1]]/nrow(df.matched)

t.test(df.matched$age[df.matched$group=='base2'],df.matched$age[df.matched$group=='ukb'],paired=FALSE, var.equal = TRUE)
cohen.d(df.matched$age,df.matched$group)$cohen.d[[2]]

# plot distributions
ggplot(df.matched, aes(age, gm, color = group)) +
  geom_point(alpha = 0.5) +
  theme_bw()
ggplot(df.matched, aes(age, wm, color = group)) +
  geom_point(alpha = 0.5) +
  theme_bw()
ggplot(df.matched, aes(age, gwm, color = group)) +
  geom_point(alpha = 0.5) +
  theme_bw()

# compare correlations
ukb.matched = df.matched[df.matched$group == 'ukb',]
base2.matched = df.matched[df.matched$group == 'base2',]
c(cor(base2.matched$gm, base2.matched$age), cor(base2.matched$wm, base2.matched$age), cor(base2.matched$gwm, base2.matched$age))
c(cor(ukb.matched$gm, ukb.matched$age), cor(ukb.matched$wm, ukb.matched$age), cor(ukb.matched$gwm, ukb.matched$age))

# compare mean absolute errors
c(mean(abs(base2.matched$gm - base2.matched$age)), mean(abs(base2.matched$wm - base2.matched$age)), mean(abs(base2.matched$gwm - base2.matched$age)))
c(mean(abs(ukb.matched$gm - ukb.matched$age)), mean(abs(ukb.matched$wm - ukb.matched$age)), mean(abs(ukb.matched$gwm - ukb.matched$age)))

# get bootstrapped confidence interval of correlation
bootCorTest = function(data, i){
  d = data[i, ]
  resultsgm = cor.test(d$gm, d$age, method='pearson')
  resultswm = cor.test(d$wm, d$age, method='pearson')
  resultsgwm = cor.test(d$gwm, d$age, method='pearson')
  
  c(gmest = resultsgm$estimate, wmest = resultswm$estimate, gwmest = resultsgwm$estimate)
}

bbase2 = boot(data = base2.matched, bootCorTest,R = 5000)
bukb = boot(data = ukb.matched, bootCorTest,R = 5000)

base2_ci_gm = boot.ci(bbase2, type = "bca", index = 1)
ukb_ci_gm = boot.ci(bukb, type = "bca", index = 1)

base2_ci_wm = boot.ci(bbase2, type = "bca", index = 2)
ukb_ci_wm = boot.ci(bukb, type = "bca", index = 2) 

base2_ci_gwm = boot.ci(bbase2, type = "bca", index = 3)
ukb_ci_gwm = boot.ci(bukb, type = "bca", index = 3)

# compare correlations between brain-predicted and chronological age
cocor.indep.groups(bukb$t0[[1]], bbase2$t0[[1]], nrow(ukb.matched), nrow(base2.matched), alternative = "two.sided")
cocor.indep.groups(bukb$t0[[2]], bbase2$t0[[2]], nrow(ukb.matched), nrow(base2.matched), alternative = "two.sided")
cocor.indep.groups(bukb$t0[[3]], bbase2$t0[[3]], nrow(ukb.matched), nrow(base2.matched), alternative = "two.sided")

# compare brain age estimates for grey matter, white matter, and combined grey and white matter brain age
ttest = t.test(df.matched$gm[df.matched$group=='base2'],df.matched$gm[df.matched$group=='ukb'],paired=FALSE, var.equal = TRUE)
pval = 2*pt(ttest$statistic[[1]], ttest$parameter[[1]], lower=FALSE)
d = psych::cohen.d(df.matched$gm,df.matched$group)$cohen.d[[2]]
ttest$parameter[[1]]; ttest$statistic[[1]]; pval; d; 

ttest = t.test(df.matched$wm[df.matched$group=='base2'],df.matched$wm[df.matched$group=='ukb'],paired=FALSE, var.equal = TRUE)
pval = 2*pt(ttest$statistic[[1]], ttest$parameter[[1]], lower=FALSE)
d = psych::cohen.d(df.matched$wm,df.matched$group)$cohen.d[[2]]
ttest$parameter[[1]]; ttest$statistic[[1]]; pval; d; 

ttest = t.test(df.matched$gwm[df.matched$group=='base2'],df.matched$gwm[df.matched$group=='ukb'],paired=FALSE, var.equal = TRUE)
pval = 2*pt(ttest$statistic[[1]], ttest$parameter[[1]], lower=FALSE)
d = psych::cohen.d(df.matched$gwm,df.matched$group)$cohen.d[[2]]
ttest$parameter[[1]]; ttest$statistic[[1]]; pval; d; 

# write table
df.matched$group = as.numeric(df.matched$group)
write.table(df.matched, file = 'data/04_matched_datasets.txt', sep = '\t', quote = FALSE, row.names = FALSE)
