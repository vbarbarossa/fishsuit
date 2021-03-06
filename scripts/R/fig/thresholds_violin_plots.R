source('config.R')

# read hist niches first
library(foreach);

hist_niches <- foreach(clmod = climate_models,.combine = 'rbind') %do% {
  
  t <- read.csv(paste0('proc/',clmod,'/niches/niches_filtered.csv'))
  t$clmod <- clmod
  return(t[,-1])
  
}

# format table for E2O
vars <- c('Qmi','Qma','Qzf','Tmi','Tma')
th <- c(2.5,100.0,97.5,2.5,97.5)

data2 <- reshape2::melt(hist_niches[hist_niches$Tma_97.5. != 0,], id = c('id_no','clmod'))
data2$clmod <- factor(data2$clmod)

data2$value[data2$variable == 'Qmi_2.5.'] <- (data2$value[data2$variable == 'Qmi_2.5.'])**(1/2)
data2$value[data2$variable == 'Qma_100.0.'] <- (data2$value[data2$variable == 'Qma_100.0.'])**(1/2)

data2$value[data2$variable == 'Tma_97.5.'] <- (data2$value[data2$variable == 'Tma_97.5.']) - 273.15
data2$value[data2$variable == 'Tmi_2.5.'] <- (data2$value[data2$variable == 'Tmi_2.5.']) - 273.15
data2$value[data2$variable == 'Tmi_2.5.'][data2$value[data2$variable == 'Tmi_2.5.'] < 0]  <- 0

levels(data2$clmod) <- c('GFDL','HadGEM','IPSL','MIROC','NorESM')
levels(data2$variable) <- c('Minimum weekly flow (2.5%) [sqrt(x) - m3/s]',
                            'Number of zero flow weeks (97.5%) [-]',
                            'Maximum weekly flow (100%) [sqrt(x) - m3/s]',
                            'Maximum weekly water temperature (97.5%) [oC]',
                            'Minimum weekly water temperature (2.5%) [oC]')


data2 <- data2[!is.na(data2$clmod),]

library(ggplot2)
p <- ggplot(data2,aes(x=clmod,y=value)) + #
  geom_violin(alpha=0.3,width = 0.5,draw_quantiles = c(0.25,0.5,0.75)) +
  # scale_fill_manual(values = brewer.pal(n = 3,'Set2')) +
  # geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(data2$clmod)) +
  ylab(label = ' ') +
  xlab(label = ' ') +  
  # coord_cartesian(ylim = c(-0.1,0.1)) +
  guides(fill = guide_legend(title = ' ')) +
  theme_bw() +
  facet_wrap('variable',ncol=1,scales = 'free_y') + 
  theme(
    # legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    text = element_text(size=20),
    axis.text.x = element_text(color='black'),
    axis.text.y = element_text(color='black')
  ) 
p

ggsave(paste0('figs/thresholds_violins.jpg'),p,width = 220,height = 220,units = 'mm',dpi=600)


