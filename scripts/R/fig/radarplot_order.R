
valerioUtils::libinv('dplyr'); library(fmsb); library(foreach); library(grid); library(RColorBrewer)
source('config.R')

df <- read.csv('figshare/RC_by_species.csv') %>%
  as_tibble()
colnames(df) <- c('binomial',paste0('RC',rep(c('1.5','2.0','3.2','4.5'),each=2),c('_dsp','')))

taxa <- rfishbase::load_taxa() %>% select(binomial = Species,Order)

df_raw <- df %>%
  left_join(taxa)

df <- df_raw %>% select(-ends_with('_dsp'))
df_dsp <- df_raw %>% select(binomial,ends_with('_dsp'),Order)
colnames(df_dsp) <- colnames(df)

# abbreviate names with first 7 letters (cannot do less or Cyprini-Cyprino get confused)
nchar_names <- 7
# set orders with less than 20 species to Other
min_sample_size <- 20

df$order_ <- lapply(df$Order,function(x) strsplit(x,'')[[1]][1:7] %>% paste0(collapse='') %>% paste0('.')) %>% do.call('c',.)
df$order_[df$order_ %in% (table(df$order_)[table(df$order_) < min_sample_size] %>% names)] <- 'Other'
df$order_ <- as.factor(df$order_)

df_dsp$order_ <- lapply(df_dsp$Order,function(x) strsplit(x,'')[[1]][1:7] %>% paste0(collapse='') %>% paste0('.')) %>% do.call('c',.)
df_dsp$order_[df_dsp$order_ %in% (table(df_dsp$order_)[table(df_dsp$order_) <= min_sample_size] %>% names)] <- 'Other'
df_dsp$order_ <- as.factor(df_dsp$order_)

# abbreviation table for SI
write.csv(df %>% select(ORDER_NAME = Order, ABBREVIATION_USED = order_) %>% distinct(),
           'tabs/order_names_abbreviation.csv',row.names = F)

df <- df %>% select(-Order)
df_dsp <- df_dsp %>% select(-Order)



cat <- 'order_'
n_min <- 1 #n_minv[i]
nchar_names <- NA #nchar_namesv[i]
size_names <- 0.5
tit = NULL

library(foreach)
data <- foreach(wt = warming_targets,.combine = 'rbind') %do% {
  #split based on categories
  ts <- split(df,df[,cat])
  # calculate the mean proportion of species losing >50% of habitat at each scenario
  res <- do.call('cbind',
                 lapply(ts,function(x) {
                   x = x[,paste0('RC',wt)]
                   x = x[complete.cases(x),]
                   # calculate the proportion of species losing >50% of habitat for each scenario
                   # then average the values across the scenarios
                   c(nrow(x),mean(apply(x,2,function(c) {sum(c > 50)/length(c)*100})) )
                 }))
  row.names(res) <- c('n',wt)
  if(wt == '1.5'){
    return(as.data.frame(res))
  }else{
    res2 <- t(as.data.frame(res[2,]))
    row.names(res2) <- wt
    return(res2)
  }
  
}

data <- data[,data[1,] >= n_min]

new_names <- colnames(data)
if(!is.na(nchar_names)) new_names <- paste0(as.character(sapply(colnames(data),function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')

colnames(data) <- paste0(new_names,'\n(',as.integer(data['n',]),')')

data <- rbind(rep(100,ncol(data)),rep(0,ncol(data)),data[2:nrow(data),][rev(row.names(data[2:nrow(data),])),])

library(RColorBrewer)
custom_pal <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)])
colors_border<- custom_pal
colors_in <- do.call('c',lapply(custom_pal,function(x) rgb(t(col2rgb(x)) ,alpha = 180,maxColorValue = 255)))


# par(mar = c(0,0,0,0),oma = c(0,0,0,0))
library(ggplot2); library(ggplotify)
p <- ggplotGrob(
  as.ggplot(
    ~radarchart( data  , axistype=1 , pty=32, seg=4,
                 #custom polygon
                 pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
                 #custom the grid
                 cglty=1, cglwd=1, cglcol="black", axislabcol="black", calcex = 0.7, caxislabels=c('0%','','50%','',''), #caxislabels=paste0(seq(0,100,25),'%'),
                 #custom labels
                 vlcex=size_names,title = tit)
  ) + theme(plot.margin = unit(c(-0.5, -0.5, -2, -2), "cm"))
)

# dispersal
data <- foreach(wt = warming_targets,.combine = 'rbind') %do% {
  #split based on categories
  ts <- split(df_dsp,df_dsp[,cat])
  # calculate the mean proportion of species losing >50% of habitat at each scenario
  res <- do.call('cbind',
                 lapply(ts,function(x) {
                   x = x[,paste0('RC',wt)]
                   x = x[complete.cases(x),]
                   # calculate the proportion of species losing >50% of habitat for each scenario
                   # then average the values across the scenarios
                   c(nrow(x),mean(apply(x,2,function(c) {sum(c > 50)/length(c)*100})) )
                 }))
  row.names(res) <- c('n',wt)
  if(wt == '1.5'){
    return(as.data.frame(res))
  }else{
    res2 <- t(as.data.frame(res[2,]))
    row.names(res2) <- wt
    return(res2)
  }
  
}

data <- data[,data[1,] >= n_min]

new_names <- colnames(data)
if(!is.na(nchar_names)) new_names <- paste0(as.character(sapply(colnames(data),function(x) strsplit(x, paste0("(?<=.{",nchar_names,"})"), perl = TRUE)[[1]][1])),'.')

colnames(data) <- paste0(new_names,'\n(',as.integer(data['n',]),')')

data <- rbind(rep(100,ncol(data)),rep(0,ncol(data)),data[2:nrow(data),][rev(row.names(data[2:nrow(data),])),])

library(RColorBrewer)
custom_pal <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)])
colors_border<- custom_pal
colors_in <- do.call('c',lapply(custom_pal,function(x) rgb(t(col2rgb(x)) ,alpha = 180,maxColorValue = 255)))


# par(mar = c(0,0,0,0),oma = c(0,0,0,0))
library(ggplot2); library(ggplotify)
p_dsp <- ggplotGrob(
  as.ggplot(
    ~radarchart( data  , axistype=1 , pty=32, seg=4,
                 #custom polygon
                 pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
                 #custom the grid
                 cglty=1, cglwd=1, cglcol="black", axislabcol="black", calcex = 0.7, caxislabels=c('0%','','50%','',''), #caxislabels=paste0(seq(0,100,25),'%'),
                 #custom labels
                 vlcex=size_names,title = tit)
  ) + theme(plot.margin = unit(c(-0.5, -0.5, -2, -2), "cm"))
)
grid.draw(p_dsp)


#------------------------------------------------------------
#> LEGEND
p_leg <- ggplot(df %>% reshape2::melt(id.vars = c('binomial','order_')),aes(x = binomial,y = value,color = variable)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(40)[c(10,20,30,40)],
                     labels = c(expression('1.5'^o*C),expression('2.0'^o*C),expression('3.2'^o*C),expression('4.5'^o*C))) +
  guides(colour = guide_legend(title=NULL,override.aes = list(size = 3))) +
  theme_bw() +
  theme(legend.direction = 'horizontal',
        text = element_text(size = 15))

tmp <- ggplot_gtable(ggplot_build(p_leg)) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 


library(ggpubr)
fig <- ggarrange(
  ggarrange(p,p_dsp,ncol= 2, nrow = 1, labels = c('No dispersal','Maximal dispersal'),
            font.label = list(size = 12, color = "black", face = "plain", family = NULL),
            label.x = 0,label.y=1),
  legend,
  ncol = 1,nrow = 2,heights = c(1,0.05))
ggsave(paste0('figs/radarplot_order.jpg'),fig,width = 180,height = 90,units='mm',dpi = 600)


