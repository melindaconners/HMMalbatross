# -------------------------------------------------------------------------------------------------
#  Explore features to determine final set of features used in HMMS for Conners et al 2021:
# "Hidden Markov models identify major movement modes in accelerometer and magnetometer data from four albatross species." Movement Ecology
# Contact M. Conners (connersm@gmail.com) for correspondence
# -------------------------------------------------------------------------------------------------

library(tictoc)
library(tidyverse)
library(Hmisc)
library(corrr)
library(flextable)


### Import data with features extracted: This dataset is the 30s data with the first two hours trimmed for each bird.

m<-read.csv(paste0(~,'/allbirds_featExtr30sec_2hrcut.csv'))

### Select desired columns of candidate features
m<-m %>% dplyr::select(spp,ID,df_heave,hf_heave,m_stheave,sd_stheave,p5_stheave,sd_head,iqr_heave,m_odba) 
colnames(m) <- c('spp', 'ID', 'df', 'hf', 'ms','ss','p5','sh','iqr','mo')
# df_heave   ('df') = dominant frequency in heave acceleration
# hf_heave   ('hf') = highest dominant frequency in heave acceleration
# m_stheave  ('ms') = mean static heave acceleration
# sd_stheave ('ss') = standard deviation of static heave acceleration
# p5_stheave ('p5') = top 5th percentile of static heave acceleration
# sd_head    ('sh') = circular deviation of heading
# iqr_heave  ('iqr') = the inter-quartile range of dynamic heave acceleration
# m_odba     ('mo') = mean ODBA


#### Initial data check
summary(m)
length(which(m$p5<0))/length(m$p5) # check p5 < 0 (all need to be positive numbers)
map(m, ~sum(is.na(.))) # note % of NaNs in full dataset

#### Plot histograms of features for initial exploration
png(filename=paste0(~,'/finalfeat_hist_fulldataset.png'), width = 365, height = 225, units='mm', res = 300)
hist.data.frame(m[,c(3:10)]) 
dev.off

#### Use correlation matrix to identify final set of features
d <- correlate(m[,c(3:10)], quiet = TRUE)
dc<-d %>% 
  fashion()  

dct <- flextable(data = dc) 
dct <- autofit(dct)
dct <- width(dct, j = 3:5, width = 1.5)
dct <- align(dct, align = 'center', part = "body")
dct <- align(dct, align = 'left', part = "all")
dct


#### Select final features in final dataframe to be used in HMM
m<-m %>% dplyr::select(ID,spp,hf,p5,sh,mo)

#### Plot histograms of features across species
mlong<-melt(m, id.vars=c("ID", "spp"))

cairo_ps(filename=paste0(~,'/finalfeat_hist_bySpp.eps'), width=15, height=7)
ggplot(data=mlong, aes(x=value, group=spp, fill=spp))+
  geom_density(adjust=1, alpha=.2)+
  scale_fill_manual(values=c("black", "black", "black", "black"))+
  theme_ipsum() +
  facet_wrap(~variable, nrow=1, scales=c("free")) 
ggtitle("Feature Distributions by Species")
dev.off()


