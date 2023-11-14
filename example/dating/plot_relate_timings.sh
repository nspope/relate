set -e
#grep "TIMING" sim_constne/*/benchmarks/benchmarks >sim_constne/timings
#grep "DATING" sim_constne/*/benchmarks/benchmarks >sim_constne/datings
awk '$1 ~ "TIMING" {print FILENAME" "$0}' sim_constne/*/benchmarks/benchmarks >sim_constne/timings
awk '$1 ~ "DATING" {print FILENAME" "$0}' sim_constne/*/benchmarks/benchmarks >sim_constne/datings

echo '
library(ggplot2)
library(cowplot)
library(dplyr)
hi <- read.table("sim_constne/timings")
mcmc_rep <- as.numeric(factor(hi[,1]))
hi <- hi[,2:ncol(hi)]
bye <- read.table("sim_constne/datings")
date_rep <- as.numeric(factor(bye[,1]))
bye <- bye[,2:ncol(bye)]
mcmc1 <- hi[,4]
date1 <- bye[,6]
mcmc3 <- hi[,5]
date3 <- bye[,7]
mcmc10 <- hi[,6]
date10 <- bye[,8]
mcmc32 <- hi[,7]
date32 <- bye[,9]
mcmc100 <- hi[,8]
date100 <- bye[,10]
mcmc316 <- hi[,9]
date316 <- bye[,11]
mcmc1000 <- hi[,10]
date1000 <- bye[,12]
mcmcEP <- hi[,12]
dateEP <- bye[,14]
mcmc <- data.frame()
mcmc_idx <- hi[,3]
date_idx <- bye[,3]
for(i in c(1,3,10,32,100,316,1000)){
  df <- data.frame(idx=mcmc_idx, rep=mcmc_rep, samps=i, seconds=get(paste0("mcmc", i)))
  df2 <- data.frame(idx=date_idx, rep=date_rep, samps=i, age=get(paste0("date", i)), baseline=date1000)
  df2 %>% group_by(idx, rep) %>% filter(baseline > 0) %>% summarise(err = mean(abs(log10(age) - log10(baseline)))) %>% as.data.frame ->df3
  df$err = df3$err
  mcmc <- rbind(mcmc, df)
}
df <- data.frame(idx=mcmc_idx, rep=mcmc_rep, samps=i, seconds=mcmcEP)
df2 <- data.frame(idx=date_idx, rep=date_rep, samps=i, age=dateEP, baseline=date1000)
df2 %>% group_by(idx, rep) %>% filter(baseline > 0) %>% summarise(err = mean(abs(log10(age) - log10(baseline)))) %>% as.data.frame ->df3
df$err = df3$err
ggplot(mcmc) + geom_point(aes(x=seconds, y=err, color=samps), pch=19) + scale_x_log10() + theme_cowplot() + 
 xlab("CPU seconds") + ylab("RMSE (log10 ages)") + 
 theme(legend.position=c(0.6,0.8), panel.background=element_rect(fill="white")) +
 geom_point(data=df, aes(x=seconds, y=err), color="firebrick", pch=19) +
 annotate(geom="text", color="firebrick", label="EP", x=0.5, y=0.02) +
 scale_color_gradient(name="MCMC samples", trans="log10")
  
ggsave("/home/natep/public_html/tsdate-paper/deleteme-timing.png", height=4, width=5., units="in", dpi=300)
' | R --slave
