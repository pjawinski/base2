# ========================
# === plot power curve ===
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
library(pwr)
library(scales)

# load data
df = read.delim('code/tables/descriptives.txt', sep = '\t', header = TRUE)

# create function for power calculations
pwrcurve = function(n, lowerlimit, upperlimit, steps, thresh, alternative) {
  xy = cbind(NULL, NULL)       
  for (i in seq(lowerlimit, upperlimit, length.out = (upperlimit-lowerlimit)*steps+1)){
    calc = pwr.r.test(n = n, r = i, sig.level = thresh, power = NULL, alternative = alternative)
    xy = rbind(xy, cbind(calc$r, calc$power))
  }
  xy = data.frame(xy)
  colnames(xy) = c("r","power")
  xy
}

# settings
lowerlimit = 0
upperlimit = 0.4
steps = 10000
ntests = 50 # corresponds to p = 0.001 after multiple testing correction
thresh = 0.05
threshcorr = 1-(1-thresh)^(1/ntests)
alternative = 'greater'

# calculate pwr for minN and maxN
maxN = pwrcurve(max(df$n), lowerlimit, upperlimit, steps, thresh, alternative)
minN = pwrcurve(min(df$n), lowerlimit, upperlimit, steps, thresh, alternative)

# calculate pwr for minN and maxN after multiple-testing correction
maxNcorr = pwrcurve(max(df$n), lowerlimit, upperlimit, steps, threshcorr, alternative)
minNcorr = pwrcurve(min(df$n), lowerlimit, upperlimit, steps, threshcorr, alternative)

# set xlim and x-axis ticks
xlim = c(lowerlimit,upperlimit)
xticks = seq(lowerlimit,upperlimit,0.1)

# Open pdf
dir.create('code/figures')
for (img in c('pdf', 'png')) {
  
  if (img == 'pdf') { 
    pdf(file = 'code/figures/power.pdf', width=5.98, height=4.48)
  } else {
    png(filename = 'code/figures/power.png', width=5.98, height=4.48, units = "in", res = 600)
  }
  
  par(mgp = c(2, 0.7, 0), lwd=0.5)
  plot(maxN$r,maxN$power,
       type="n",
       main="", 
       ylab=expression(paste("1-", beta)),
       xlab="",
       pch=20,
       cex.lab=0.9,
       col="royalblue3",
       xlim=xlim,
       ylim=c(0,1),
       las=1,
       axes=FALSE)
  
  # Insert lines
  lines(minN$r, minN$power, col="royalblue3", type="l", lty=2)
  lines(maxN$r, maxN$power, col="royalblue3")
  lines(minNcorr$r, minNcorr$power, col="red", type="l", lty=2)
  lines(maxNcorr$r, maxNcorr$power, col="red")
  
  # Adjust axis
  axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), las=1, tck=-0.02, cex.axis=0.75, lwd=0.75)
  par(mgp = c(1, 0.3, 0))
  axis(side=1, at=xticks, labels=xticks, las=1, tck=-0.02, cex.axis=0.75, lwd=0.75)
  
  # X-axis title
  mtext(expression(paste("Pearson ", rho)), side=1, line=1.5, cex=0.9)
  
  # Insert dots at 20/50/80% Power
  doty = c(0.2,0.5,0.8)
  for (i in 1:3){
  dotx = pwr.r.test(n = max(df$n), r = NULL, sig.level = thresh, power = doty[i], alternative = alternative)$r
  text(x = dotx, y = doty[i], labels = bquote(paste("1-", beta," = ", .(doty[i]), ", ", rho, " = ", .(format(round(dotx,digits=3),nsmall=3)))), cex=0.5, pos=4)
  points(x = dotx, y = doty[i], pch=16, col="royalblue3", cex = 0.5)
  }
  
  # Finish plot
  dev.off()
}

