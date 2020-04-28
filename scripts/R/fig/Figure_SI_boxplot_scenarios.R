source('R/4targets/MASTER.R')

library(foreach); library(ggplot2); library(RColorBrewer)

ESH_filename <- paste0(dir_mod,'tabs/ESH_merged_warmingtargets.rds')

if(file.exists(ESH_filename)){
  tab <- readRDS(ESH_filename)
}else{
  
  tab <- foreach(clmod = climate_models,.combine='rbind') %do% {
    
    ESH_tab <- read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab.csv'))
    
    tcl <- reshape2::melt(ESH_tab,id.vars = c(1:2,ncol(ESH_tab)))
    colnames(tcl)[(ncol(tcl)-1):ncol(tcl)] <- c('ESH_type','ESH')
    levels(tcl$ESH_type) <- c('Total','Q','Tw','Q&Tw')
    tcl$GCM <- clmod
    return(tcl)
  }
  tab <- droplevels(tab)
  tab$GCM <- factor(tab$GCM)
  
  tab <- cbind(tab,as.data.frame(do.call('rbind',strsplit(as.character(tab[,3]),'_'))))
  colnames(tab)[7:9] <- c('RCP','warmt','year')
  
  saveRDS(tab,ESH_filename)
}

levels(tab$RCP) <- c('2.6','4.5','6.0','8.5')
levels(tab$GCM) <- c('GFDL','HadGEM','IPSL','MIROC','NorESM')

# based on warming degrees only
p <- ggplot(tab[tab$ESH_type == 'Total',],aes(x=warmt,y=ESH)) + #
  geom_boxplot(varwidth = F,width = 0.5) + #
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(tab$warmt)) +
  scale_y_reverse() +
  xlab(label = ' ') +
  ylab(label = '% of suitable habitat lost') +
  theme_bw() +
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

ggsave(paste0(dir_(paste0(dir_mod,'figs/Figure_SI_boxplots/')),'boxplot_warming_targets.jpg'),p,width = 220, height = 120, units = 'mm', dpi = 600,scale = 1.3)


# based on warming degrees and climate
p <- ggplot(tab[tab$ESH_type == 'Total',],aes(x=warmt,y=ESH)) + #
  geom_boxplot(aes(color = GCM),varwidth = F,width = 0.5) + #
  scale_color_manual(values = brewer.pal(n = 5,'Set1')) +
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(tab$warmt)) +
  scale_y_reverse() +
  xlab(label = ' ') +
  ylab(label = '% of suitable habitat lost') +
  theme_bw() +
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

ggsave(paste0(dir_(paste0(dir_mod,'figs/Figure_SI_boxplots/')),'boxplot_warming_targets_GCM.jpg'),p,width = 220, height = 120, units = 'mm', dpi = 600,scale = 1.3)

# based on warming degrees and scenarios
p <- ggplot(tab[tab$ESH_type == 'Total',],aes(x=warmt,y=ESH)) + #
  geom_boxplot(aes(linetype = RCP),varwidth = F,width = 0.5) + #
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(tab$warmt)) +
  scale_y_reverse() +
  xlab(label = ' ') +
  ylab(label = '% of suitable habitat lost') +
  theme_bw() +
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

ggsave(paste0(dir_(paste0(dir_mod,'figs/Figure_SI_boxplots/')),'boxplot_warming_targets_RCP.jpg'),p,width = 220, height = 120, units = 'mm', dpi = 600,scale = 1.3)

# based on warming degrees RCP and GCM
p <- ggplot(tab[tab$ESH_type == 'Total',],aes(x=warmt,y=ESH)) + #
  geom_boxplot(aes(color = GCM,linetype = RCP),varwidth = F,width = 0.5,alpha=0.5,outlier.size = 0.1) + #
  scale_color_manual(values = brewer.pal(n = 5,'Set1')) +
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(tab$warmt)) +
  scale_y_reverse() +
  # xlab(label = 'RCP Scenario') +
  ylab(label = '% of suitable habitat lost') +
  xlab(label = ' ') +
  # guides(color = guide_legend(title = ' ')) +
  theme_bw() +
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

ggsave(paste0(dir_(paste0(dir_mod,'figs/Figure_SI_boxplots/')),'boxplot_warming_targets_combo.jpg'),p,width = 220, height = 120, units = 'mm', dpi = 600,scale = 1.3)


# one figure for GCM-RCP
d <- droplevels(rbind(cbind(tab[tab$ESH_type == 'Total',],data.frame(g = 'clim')),
           cbind(tab[tab$ESH_type == 'Total',],data.frame(g = 'scen'))
))
levels(d$RCP) <- paste0('RCP ',levels(d$RCP))
colnames(d)[6] <- 'f'
d$f <- as.character(d$f)
d$f[d$g == 'scen'] <- as.character(d$RCP[d$g == 'scen'])
d$f <- as.factor(d$f)

levels(d$g) <- c('GCM','RCP scenario')
# based on warming degrees and scenarios
p <- ggplot(d,aes(x=warmt,y=ESH)) + #
  geom_boxplot(aes(fill = f),varwidth = F,width = 0.5,outlier.size = 0.3,outlier.alpha = 0.3,outlier.colour = 'Grey') + #
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(tab$warmt)) +
  scale_fill_manual(values = c(brewer.pal(n = 5,'Blues'),brewer.pal(4,'Reds'))) +
  # guides(fill = ' ') +
  scale_y_reverse() +
  xlab(label = ' ') +
  ylab(label = '% of suitable habitat lost') +
  facet_wrap('g',ncol = 1) +
  theme_bw() +
  theme(
    # legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    text = element_text(size=20),
    axis.text.x = element_text(color='black'),
    axis.text.y = element_text(color='black'),
    legend.title = element_blank()
  ) 

ggsave(paste0(dir_mod,'figs/Figure_SI_boxplot_warming_targets_combo.jpg'),p,width = 170, height = 150, units = 'mm', dpi = 600,scale = 1.3)


