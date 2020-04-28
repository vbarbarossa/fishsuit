#For guidance, Nature's standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). 
#The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).

source('R/4targets/MASTER.R');

library(foreach); library(dplyr); library(ggplot2); library(RColorBrewer)

#> FUNCTIONS TO BUILD DATATABLES -------------------------------------------------------------------------------------------

build_tab <- function(dispersal = FALSE){
  suffix <- ''
  if(dispersal) suffix <- '_dispersal3'
  return(
    foreach(clmod = climate_models,.combine='rbind') %do% {
      read.csv(paste0('fishsuit_completeRun_warming_4targets/',clmod,'/ESH_tab',suffix,'.csv')) %>%
        reshape2::melt(.,measure.vars = c('ESH_all', 'ESH_Q', 'ESH_T', 'ESH_QnT')) %>%
        rename(ESH_type = variable, ESH = value) %>%
        mutate(ESH_type = forcats::fct_recode(ESH_type, Total='ESH_all', Q = 'ESH_Q', Tw = 'ESH_T', 'Q&Tw' = 'ESH_QnT'),
               GCM = clmod)
      
    } %>%
      droplevels(.) %>%
      mutate(GCM = factor(GCM))%>%
      cbind(.,as.data.frame(do.call('rbind',strsplit(as.character(.[,'comboscen']),'_')))) %>%
      as_tibble() %>%
      rename(RCP = V1,warmt = V2,year = V3) %>%
      mutate(RCP = forcats::fct_recode(RCP,'2.6' = 'rcp2p6','4.5' = 'rcp4p5','6.0' = 'rcp6p0','8.5' ='rcp8p5' ),
             GCM = forcats::fct_recode(GCM,'GFDL'="gfdl",'HadGEM'="hadgem",'IPSL'="ipsl",'MIROC'="miroc",'NorESM'="noresm"))
  )
}

compute_median <- function(tab){
  return(
    foreach(sp = unique(tab$id_no),.combine = 'rbind') %do% {
      tsp <- droplevels(tab[tab$id_no == sp,])
      foreach(wt = paste0(warming_targets,'C'),.combine = 'rbind') %do% {
        
        t <- droplevels(tsp[tsp$warmt == wt,])
        return(
          data.frame(
            id_no = sp,
            ESH_median = median(t$ESH,na.rm=T), #<<<<<<<<<<<<<<<<<<<< NEED TO CHECK WHY MEDIAN HERE
            warmt = wt
          )
        )
      } 
    } %>% as_tibble()
  )
}


#> RC OVERALL VIOLIN PLOTS -----------------------------------------------------------------------------------------

tab <- rbind(
  build_tab() %>% select(id_no,ESH_type,ESH,GCM,RCP,warmt,year) %>% filter(ESH_type == 'Total') %>% 
    droplevels() %>% compute_median() %>% mutate(scenario = 'no dispersal'),
  build_tab(dispersal = TRUE) %>% select(id_no,ESH_type,ESH,GCM,RCP,warmt,year) %>% filter(ESH_type == 'Total') %>% 
    droplevels() %>% compute_median() %>% mutate(scenario = 'maximal dispersal')
) %>%
  mutate(scenario = factor(scenario))


p <- ggplot(tab,aes(x=warmt,y=ESH_median)) + #
  geom_violin(aes(fill = scenario),lwd = .5,color='transparent') + #
  geom_boxplot(aes(color = scenario),fill='white',outlier.color = 'transparent',
               width = 0.08,lwd=0.5,coef=0,notch = 1,position = position_dodge(0.9)) +
  
  scale_fill_manual(values = viridis::viridis(10,option='C')[c(5,8)]) +
  scale_color_manual(values = viridis::viridis(10,option='C')[c(5,8)]) +
  
  # stat_summary(fun.y=mean, geom="point",aes(fill = warmt), 
  #              shape=23, size=2,show.legend = FALSE) +
  
  scale_x_discrete(labels=c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  scale_y_reverse(breaks = c(0,25,50,75,100),limits = c(100,0),labels=paste0(c(0,25,50,75,100),'%')) +
  xlab(label = ' ') +
  ylab(label = 'Range losses') +
  
  theme_bw() +
  coord_cartesian(expand=F) +
  theme(
    legend.position="bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-20,0,0,0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed',color='black'),
    axis.ticks.x = element_blank(),
    text = element_text(size=10),
    axis.text.x = element_text(color='black',vjust = 3),
    axis.text.y = element_text(color='black',angle=90,hjust = 0.5, vjust=1),
    axis.line.y = element_line(color='black'),
    # axis.line.y.right = element_line(),
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    # , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    # , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_text(angle = 0, size = 16)
    
    
  )

p

ggsave(paste0(dir_mod,'figs/violins_overall_RC.jpg'),p,
       width = 89,height = 60,dpi = 600,units = 'mm')
ggsave(paste0(dir_mod,'figs/violins_overall_RC.pdf'),p,
       width = 89,height = 60,units = 'mm')

# save the table for figshare
tab_wide <- tab %>% reshape2::dcast(.,id_no ~ warmt + scenario,value.var = 'ESH_median') %>% as_tibble()
iucn_simplified <- foreach(i = 1:2,.combine = 'rbind') %do% foreign::read.dbf(paste0('data/FW_FISH_PART_',i,'.dbf')) %>%
  as_tibble() %>%
  distinct(binomial,.keep_all = T)

tab_wide <- right_join(iucn_simplified %>% select(id_no,binomial),tab_wide) %>%
  select(-id_no)

write.csv(tab_wide,
          paste0(dir_mod,'figshare/RC_by_species.csv'),row.names = F)


# BOXPLOTS BY RCP AND GCM ------------------------------------------------------------------------------

tab <- rbind(
  build_tab() %>% select(id_no,ESH_type,ESH,GCM,RCP,warmt,year) %>% mutate(scenario = 'no dispersal'),
  build_tab(dispersal = TRUE) %>% select(id_no,ESH_type,ESH,GCM,RCP,warmt,year) %>% mutate(scenario = 'maximal dispersal')
) %>%
  mutate(scenario = factor(scenario)) %>%
  droplevels()

# one figure for GCM-RCP
d <- rbind(
  tab %>% filter(ESH_type == 'Total') %>% mutate(g = 'GCM') %>% mutate(f = factor(GCM)),
  tab %>% filter(ESH_type == 'Total') %>% mutate(g = 'RCP scenario') %>% 
    mutate( RCP = forcats::fct_recode(RCP, 'RCP 2.6' = '2.6', 'RCP 4.5' = "4.5", 'RCP 6.0' = "6.0", 'RCP 8.5' = "8.5" )) %>%
    mutate(f = factor(RCP))
  ) %>% 
  mutate(g = factor(g)) %>%
  droplevels()

# based on warming degrees and scenarios
p <- ggplot(d,aes(x=warmt,y=ESH)) + #
  geom_boxplot(aes(fill = f),varwidth = F,width = 0.5,outlier.size = 0.3,outlier.alpha = 0.3,outlier.colour = 'Grey') + #
  geom_hline(yintercept=0, linetype="dotted") +
  scale_x_discrete(labels=levels(tab$warmt)) +
  scale_fill_manual(values = c(brewer.pal(n = 5,'Blues'),brewer.pal(4,'Reds'))) +
  # guides(fill = ' ') +
  scale_y_reverse() +
  xlab(label = ' ') +
  ylab(label = 'Range losses [%]') +
  facet_grid(g ~ scenario) +
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
p

ggsave(paste0(dir_mod,'figs/boxplot_RCP_and_GCM.jpg'),p,width = 170, height = 170, units = 'mm', dpi = 600,scale = 1.3)



