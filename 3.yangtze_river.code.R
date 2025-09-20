###load packages
library(DECIPHER)
library(reshape2)
library(ggpubr)
library(vegan)
library(ggmap)
library(leaflet)
library(ggridges)
require(ggplotify)
library(dplyr)
library(ggrepel)
library(forcats)
library(RColorBrewer)
library(ggplot2)
library(ggprism)
##################step1: sequence composition
##################step2: Calculate the DNA extraction efficiency of the sample
##################step3: Quantitative calibration
##################step4: compare the adjustment effect
###set wd
setwd("../yangtze-river/")
###################################################
##################step1: sequence composition
###################################################
###read NGS otu table
otu <- read.csv("ngs.otu.table.csv" , row.names = 1 ) ; head( otu )
fishreads <- sum( otu[ grepl("Actinopteri" , otu$taxonomy) , !grepl("taxonomy" , colnames( otu ) ) ])
nonfishreads <- sum( otu[ (!grepl("Actinopteri" , otu$taxonomy) ) & (!grepl("IPS" , otu$taxonomy) ) , !grepl("taxonomy" , colnames( otu ) ) ] )
ipsreads <- sum( otu[ grepl("IPS" , otu$taxonomy) , !grepl("taxonomy" , colnames( otu ) ) ] )
###make sequence summary table
dfs <- data.frame( categroy = c("IPS1" , "IPS2.1" , "IPS2.2" , "IPS2.3" , "Fish" , "Non-fish"), reads = c(rowSums( otu[ 1:4 , !grepl("taxonomy" , colnames( otu ) )] , na.rm = T)  , fishreads , nonfishreads) ) ; head( dfs )
dfs$prop <- round( dfs$reads/sum( dfs$reads) , 4) ; dfs
dfs$type <- "Total"
dfs2 <- dfs[ grepl("IPS" , dfs$categroy), ]
dfs2$prop <- round( dfs2$reads/sum( dfs2$reads) , 4) ; dfs2
dfs2$type <- "IPS"
seqdf <- rbind( dfs , dfs2)
seqdf$type <- factor( seqdf$type , levels = c("Total" , "IPS"))
seqdf$categroy <- factor( seqdf$categroy , levels = c("IPS1" , "IPS2.1" , "IPS2.2" , "IPS2.3" , "Fish" , "Non-fish"))
###create label information
lb <- data.frame(type = c("Total" , "IPS") , count = c( sum(dfs$reads) , sum( dfs2$reads)))
###plot sequence composition
library(ggalluvial)
library(ggprism)
psc <- ggplot( data = seqdf , aes(x = type , y = prop , fill = categroy) ) +
  geom_bar( stat = "identity" , position = "fill") +
  geom_flow( aes(alluvium = categroy ), alpha = 0.5) + 
  scale_fill_manual( values = c("IPS1" = "#db6968",
                                "IPS2.1" ="#4d97cd",
                                "IPS2.2" = "#f8984e",
                                "IPS2.3" ="#459943",
                                "Fish" = brewer.pal(9, "PuBu")[7],
                                "Non-fish" = brewer.pal(9, "PuBu")[5] ) 
                     )+
  geom_text( data = lb , aes(x = type , y=1, label = paste0( count ) ) , inherit.aes = FALSE , vjust = -0.2 , size =3)+
  theme_light()+
  labs( x = "Reads" , y = "Proportion")+
  guides(fill = guide_legend(title = ""))+
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
psc

######count the ips2.1 ips2.2 and ips2.3 reads
otu <- read.csv("ngs.otu.table.csv" , row.names = 1 ) ; head( otu )
ipscount <- as.data.frame( t( otu[ grepl("IPS" , rownames( otu ) ), !grepl("taxonomy" , colnames( otu ) )] ) ) ; head( ipscount )
df.ips <- melt( ipscount ) ; head( df.ips )
ips2_compare <- list(c("IPS2.1","IPS2.2"), c("IPS2.2","ISP2.3") )

p.ips2 <- ggplot( df.ips[ !df.ips$variable %in% "IPS1", ] , aes(x = variable , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = variable), pch = 21,size = 2, position = position_jitter(0.2))+
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#98d09d")) + 
  labs(x = NULL , y = "Reads") +
  theme_prism() +
  guides(fill = guide_legend(title = ""))+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10 , angle = 45 , hjust = 0.5),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p.ips2

###################################################
##################step2: Calculate the DNA extraction efficiency of the sample
###################################################
###read IPS1 qPCR data
ips1.qpcr <- read.csv("ips1.qpcr.csv") ; head( ips1.qpcr )
ips1.qpcr <- ips1.qpcr[ !ips1.qpcr$ct_value %in% "Undetermined" , ]
ips1.qpcr$ct_value <- as.numeric( ips1.qpcr$ct_value )
##plate1
plate1.qpcr <- ips1.qpcr[ ips1.qpcr$plate == "plate1", ]
plate2.qpcr <- ips1.qpcr[ ips1.qpcr$plate == "plate2", ]
plate3.qpcr <- ips1.qpcr[ ips1.qpcr$plate == "plate3", ]
plate4.qpcr <- ips1.qpcr[ ips1.qpcr$plate == "plate4", ]

library(dplyr)
library(Rmisc)
## calculate CT  mean
plate1.qpcr.ctmean <- plate1.qpcr %>% group_by(sample.1) %>% dplyr::summarise(mean_value = mean(ct_value))
plate1.qpcr.ctmean <- as.data.frame( plate1.qpcr.ctmean )
plate1.qpcr.ctmean <- plate1.qpcr.ctmean[ !plate1.qpcr.ctmean$sample.1 %in% "Standard curve", ]

plate2.qpcr.ctmean <- plate2.qpcr %>% group_by(sample.1) %>% dplyr::summarise(mean_value = mean(ct_value))
plate2.qpcr.ctmean <- as.data.frame( plate2.qpcr.ctmean )
plate2.qpcr.ctmean <- plate2.qpcr.ctmean[ !plate2.qpcr.ctmean$sample.1 %in% "Standard curve", ]

plate3.qpcr.ctmean <- plate3.qpcr %>% group_by(sample.1) %>% dplyr::summarise(mean_value = mean(ct_value))
plate3.qpcr.ctmean <- as.data.frame( plate3.qpcr.ctmean )
plate3.qpcr.ctmean <- plate3.qpcr.ctmean[ !plate3.qpcr.ctmean$sample.1 %in% "Standard curve", ]

plate4.qpcr.ctmean <- plate4.qpcr %>% group_by(sample.1) %>% dplyr::summarise(mean_value = mean(ct_value))
plate4.qpcr.ctmean <- as.data.frame( plate4.qpcr.ctmean )
plate4.qpcr.ctmean <- plate4.qpcr.ctmean[ !plate4.qpcr.ctmean$sample.1 %in% "Standard curve", ]

ips1.conc <- 30
## linear fitting model for plate1
stand1 <- plate1.qpcr[plate1.qpcr$sample.1 == "Standard curve" , ]
stand1$copies <- log10( stand1$copies)
model1 <- lm( stand1$copies ~ stand1$ct_value )
plate1.qpcr.ctmean$copies <- 10^( coef(model1)[1] + (plate1.qpcr.ctmean$mean_value * coef(model1)[2]) )
plate1.qpcr.ctmean$logcopy <- log10( plate1.qpcr.ctmean$copies )
plate1.qpcr.ctmean$extract_rate <- plate1.qpcr.ctmean$copies/ips1.conc
plate1.qpcr.ctmean
## linear fitting model for plate2
stand2 <- plate2.qpcr[plate2.qpcr$sample.1 == "Standard curve" , ]
stand2$copies <- log10( stand2$copies )
model2 <- lm( stand2$copies ~ stand2$ct_value )
plate2.qpcr.ctmean$copies <- 10^(coef(model2)[1] + (plate2.qpcr.ctmean$mean_value * coef(model2)[2]))
plate2.qpcr.ctmean$logcopy <- log10( plate2.qpcr.ctmean$copies )
plate2.qpcr.ctmean$extract_rate <- plate2.qpcr.ctmean$copies/ips1.conc
plate2.qpcr.ctmean
## linear fitting model for plate3
stand3 <- plate3.qpcr[plate3.qpcr$sample.1 == "Standard curve" , ]
stand3$copies <- log10( stand3$copies )
model3 <- lm( stand3$copies ~ stand3$ct_value )
plate3.qpcr.ctmean$copies <- 10^(coef(model3)[1] + (plate3.qpcr.ctmean$mean_value * coef(model3)[2]))
plate3.qpcr.ctmean$logcopy <- log10( plate3.qpcr.ctmean$copies )
plate3.qpcr.ctmean$extract_rate <- plate3.qpcr.ctmean$copies/ips1.conc
plate3.qpcr.ctmean
## linear fitting model for plate4
stand4 <- plate4.qpcr[plate4.qpcr$sample.1 == "Standard curve" , ]
stand4$copies <- log10( stand4$copies )
model4 <- lm( stand4$copies ~ stand4$ct_value )
plate4.qpcr.ctmean$copies <- 10^(coef(model2)[1] + (plate4.qpcr.ctmean$mean_value * coef(model2)[2]))
plate4.qpcr.ctmean$logcopy <- log10( plate4.qpcr.ctmean$copies )
plate4.qpcr.ctmean$extract_rate <- plate4.qpcr.ctmean$copies/ips1.conc
plate4.qpcr.ctmean
##combine all data and save table
df.ips1.copy <- rbind( plate1.qpcr.ctmean , plate2.qpcr.ctmean , plate3.qpcr.ctmean , plate4.qpcr.ctmean)
df.ips1.copy$sample.1 <- gsub("-",".",df.ips1.copy$sample.1) ; head( df.ips1.copy )
dir.create("results")
write.csv( df.ips1.copy , file = "results/ips1.qpcr.copies.csv" , row.names = F)
## plot the linear fitting model
library(RColorBrewer)
p.fitmodel1 <- ggplot(stand1 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = plate1.qpcr.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(plate1.qpcr.ctmean$logcopy), xend =plate1.qpcr.ctmean$mean_value[plate1.qpcr.ctmean$logcopy == min(plate1.qpcr.ctmean$logcopy)], yend = min(plate1.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(plate1.qpcr.ctmean$logcopy), xend =plate1.qpcr.ctmean$mean_value[plate1.qpcr.ctmean$logcopy == max(plate1.qpcr.ctmean$logcopy)], yend = max(plate1.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(plate1.qpcr.ctmean$mean_value), y = -Inf, xend = min(plate1.qpcr.ctmean$mean_value), yend = plate1.qpcr.ctmean$logcopy[plate1.qpcr.ctmean$mean_value == min(plate1.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(plate1.qpcr.ctmean$mean_value), y = -Inf, xend = max(plate1.qpcr.ctmean$mean_value), yend = plate1.qpcr.ctmean$logcopy[plate1.qpcr.ctmean$mean_value == max(plate1.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(plate1.qpcr$ct_value) , y = min(plate1.qpcr.ctmean$logcopy)+0.2 , label = round(min(plate1.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(plate1.qpcr$ct_value) , y = max(plate1.qpcr.ctmean$logcopy)+0.2 , label = round(max(plate1.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "Plate1 copies (log10)") +
  theme_classic2() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel1

p.fitmodel2 <- ggplot(stand2 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = plate2.qpcr.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(plate2.qpcr.ctmean$logcopy), xend =plate2.qpcr.ctmean$mean_value[plate2.qpcr.ctmean$logcopy == min(plate2.qpcr.ctmean$logcopy)], yend = min(plate2.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(plate2.qpcr.ctmean$logcopy), xend =plate2.qpcr.ctmean$mean_value[plate2.qpcr.ctmean$logcopy == max(plate2.qpcr.ctmean$logcopy)], yend = max(plate2.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(plate2.qpcr.ctmean$mean_value), y = -Inf, xend = min(plate2.qpcr.ctmean$mean_value), yend = plate2.qpcr.ctmean$logcopy[plate2.qpcr.ctmean$mean_value == min(plate2.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(plate2.qpcr.ctmean$mean_value), y = -Inf, xend = max(plate2.qpcr.ctmean$mean_value), yend = plate2.qpcr.ctmean$logcopy[plate2.qpcr.ctmean$mean_value == max(plate2.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(plate2.qpcr$ct_value) , y = min(plate2.qpcr.ctmean$logcopy)+0.2 , label = round(min(plate2.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(plate2.qpcr$ct_value) , y = max(plate2.qpcr.ctmean$logcopy)+0.2 , label = round(max(plate2.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "Plate2 copies (log10)") +
  theme_classic2() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel2

## plot the linear fitting model
p.fitmodel3 <- ggplot(stand3 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = plate3.qpcr.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(plate3.qpcr.ctmean$logcopy), xend =plate3.qpcr.ctmean$mean_value[plate3.qpcr.ctmean$logcopy == min(plate3.qpcr.ctmean$logcopy)], yend = min(plate3.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(plate3.qpcr.ctmean$logcopy), xend =plate3.qpcr.ctmean$mean_value[plate3.qpcr.ctmean$logcopy == max(plate3.qpcr.ctmean$logcopy)], yend = max(plate3.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(plate3.qpcr.ctmean$mean_value), y = -Inf, xend = min(plate3.qpcr.ctmean$mean_value), yend = plate3.qpcr.ctmean$logcopy[plate3.qpcr.ctmean$mean_value == min(plate3.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(plate3.qpcr.ctmean$mean_value), y = -Inf, xend = max(plate3.qpcr.ctmean$mean_value), yend = plate3.qpcr.ctmean$logcopy[plate3.qpcr.ctmean$mean_value == max(plate3.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(plate3.qpcr$ct_value) , y = min(plate3.qpcr.ctmean$logcopy)+0.2 , label = round(min(plate3.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(plate3.qpcr$ct_value) , y = max(plate3.qpcr.ctmean$logcopy)+0.2 , label = round(max(plate3.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "Plate3 copies (log10)") +
  theme_classic2() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel3

## plot the linear fitting model
p.fitmodel4 <- ggplot(stand4 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = plate4.qpcr.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(plate4.qpcr.ctmean$logcopy), xend =plate4.qpcr.ctmean$mean_value[plate4.qpcr.ctmean$logcopy == min(plate4.qpcr.ctmean$logcopy)], yend = min(plate4.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(plate4.qpcr.ctmean$logcopy), xend =plate4.qpcr.ctmean$mean_value[plate4.qpcr.ctmean$logcopy == max(plate4.qpcr.ctmean$logcopy)], yend = max(plate4.qpcr.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(plate4.qpcr.ctmean$mean_value), y = -Inf, xend = min(plate4.qpcr.ctmean$mean_value), yend = plate4.qpcr.ctmean$logcopy[plate4.qpcr.ctmean$mean_value == min(plate4.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(plate4.qpcr.ctmean$mean_value), y = -Inf, xend = max(plate4.qpcr.ctmean$mean_value), yend = plate4.qpcr.ctmean$logcopy[plate4.qpcr.ctmean$mean_value == max(plate4.qpcr.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(plate4.qpcr$ct_value) , y = min(plate4.qpcr.ctmean$logcopy)+0.2 , label = round(min(plate4.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(plate4.qpcr$ct_value) , y = max(plate4.qpcr.ctmean$logcopy)+0.2 , label = round(max(plate4.qpcr.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "Plate4 copies (log10)") +
  theme_classic2() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel4

p.fit <- ggarrange(p.fitmodel1 , p.fitmodel2 , p.fitmodel3 , p.fitmodel4 ) ; p.fit

dir.create("figures")
###save the linear fitting model
pdf("figures/ips1.fit.pdf")
p.fit
dev.off()

###correlation between ips1 copies and NGS read
otu <- read.csv("ngs.otu.table.csv" , row.names = 1 ) ; head( otu )
ipscount <- as.data.frame( t( otu[ grepl("IPS" , rownames( otu ) ), !grepl("taxonomy" , colnames( otu ) )] ) ) ; head( ipscount )
ipscount$sample <- rownames( ipscount ) ; head( ipscount )
ips1copy <- read.csv("results/ips1.qpcr.copies.csv") ; head( ips1copy ) ; dim( ips1copy )

library(reshape2)
df.ips1.com <- merge( ips1copy , ipscount , by.x = "sample.1" , by.y = "sample") ; head( df.ips1.com )

p.copy.reads <- ggplot(df.ips1.com , aes(x = IPS1 , y = extract_rate )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  labs(x = "IPS1 NGS reads", y = "DNA extraction efficiency") +
  theme_classic2() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.., sep = "~~~")), formula = y ~ x, parse = TRUE, color = "red") +  # 显示拟合公式、R² 和 p值
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson") #+
  #stat_regline_equation(size = 3, color = "red",label.y.npc = 1)
p.copy.reads

p1 <- ggarrange(psc , p.ips2 , ncol = 2 , widths = c(1.2, 1) , labels = c("a" , "b")) ; p1
p2 <- ggarrange( p1 ,  p.copy.reads , ncol = 1 , labels = c("a" , "c")) ; p2

##save the figure
pdf("figures/Figure4.pdf")
p2
dev.off()
###################################################
##################step2: Calculate the PCR amplification efficiency of the sample
###################################################
## read NGS reads table
otu <- read.csv("ngs.otu.table.csv" , row.names = 1 ) ; head( otu ) ; dim(otu)
ipscount <- as.data.frame( t( otu[ grepl("IPS" , rownames( otu ) ), !colnames(otu) %in% c("YCS11.RT" ,"YCS15.LT","YCS13.RM" , "taxonomy")] ) ) ; head( ipscount )

yangtze <- otu[ !grepl( "IPS" , otu$taxonomy ) , !colnames(otu) %in% c("YCS11.RT" ,"YCS15.LT","YCS13.RM" )] ; head(yangtze.clean) ; dim( yangtze )
yangtze.clean <- otu[ !grepl( "IPS" , otu$taxonomy ) , !colnames(otu) %in% c("YCS11.RT" ,"YCS15.LT","YCS13.RM" , "taxonomy")] ; head(yangtze.clean) ; dim( yangtze.clean )

yangtze.pcr <- data.frame( matrix( ncol = ncol( yangtze.clean ) , nrow = nrow( yangtze.clean ) ) )
rho <- c() ; pvalue <- c() ; formula <- c() ; ips1adj <- c() ; slope <- c()
for( i in 1: nrow( ipscount ) ){
  x <- c(0 , 5.5 , 32.75 , 66.83)
  y <- c( 0 , as.numeric( ipscount[i, colnames( ipscount ) %in% c("IPS2.1" , "IPS2.2" , "IPS2.3")] ) )
  model <- lm( y ~ x )
  cts <- summary( model )
  rho[i] <- cts$r.squared
  slope[i] <- coef(model)[2]
  pvalue[i] <- round(cts$coefficients[2,4] , 4)
  formula[i] <- paste( "y =", round(coef(model)[1] , 2) , "+", paste0(round(coef(model)[2] , 2)  , "x") )
  #ips1adj[i] <- round( ( ips1 - coef(model)[1][[1]])/coef(model)[2][[1]] , 0)
  #yangtze.pcr[ ,i] <- round( ( yangtze.clean[,i] - coef(model)[1][[1]])/coef(model)[2][[1]] , 0)
  ips1adj[i] <- round( ipscount$IPS1 /coef(model)[2][[1]] , 0)
  yangtze.pcr[ ,i] <- round( 100*yangtze.clean[,i] /coef(model)[2][[1]] , 0)
}
yangtze.pcr[ yangtze.pcr < 0 ] <- 0
rownames( yangtze.pcr ) <- rownames( yangtze.clean )
colnames( yangtze.pcr ) <- colnames( yangtze.clean )
yangtze.pcr$taxonomy <- yangtze$taxonomy ; head( yangtze.pcr )

df1 <- ipscount ; head( df1 )
df1$R = rho
df1$pvalue = pvalue
df1$slope  = slope
df1$formula = formula
df1
#save the PCR amplification efficiency
write.csv( df1 , file = "results/IPS2.fitting.model.csv" , row.names = T)
#write.csv( yangtze.pcr , file = "results/ngs.reads.PCRadj.sample.csv" )
###################################################
##################step3: Quantitative calibration
###################################################
library(dplyr)
##read DNA extraction rate
ips1copy <- read.csv("results/ips1.qpcr.copies.csv" , row.names = 1) ; head( ips1copy ) ; dim( ips1copy )
##read PCR amplification efficiency
pcramp <- read.csv("results/IPS2.fitting.model.csv" , row.names = 1) ; head( pcramp ) ; dim( pcramp )
##read NGS reads table
otu <- read.csv("ngs.otu.table.csv" , row.names = 1) ; head( otu ) ; dim( otu )
yangtze <- otu[ !grepl("IPS" , otu$taxonomy ) , ] ; head( yangtze )
yangtze.clean <- yangtze[ ,!colnames(yangtze) %in% c("YCS11.RT" ,"YCS15.LT","YCS13.RM")] ; head(yangtze.clean) ; dim( yangtze.clean )
yangtze.clean <- yangtze.clean[ , colnames( yangtze.clean ) %in% rownames(ips1copy) ] ; dim( yangtze.clean )
yangtze.clean$taxonomy <- yangtze$taxonomy ; head( yangtze.clean )
write.csv(yangtze.clean , file = "results/ngs.reads.rawdata.csv")
##PCR amplification efficiency correction
pcr.adj <- yangtze.clean[ , !grepl("taxonomy" , colnames( yangtze.clean ) )] ; head( pcr.adj )
for(i in 1: ncol( pcr.adj ) ){
  pcr.adj[,i] <- round( 100* pcr.adj[ , i]/pcramp$slope[rownames( pcramp ) %in% colnames( pcr.adj )[i] ] , 0)
}
pcr.adj$taxonomy <- yangtze.clean$taxonomy; head( pcr.adj ) ; dim( pcr.adj )
write.csv( pcr.adj , file = "results/ngs.reads.PCRadj.sample.csv" , row.names = T) 

##DNA extraction correction
dna.adj <- yangtze.clean[ , !grepl("taxonomy" , colnames( yangtze.clean ) )] ; head( dna.adj )
for(i in 1: ncol( dna.adj ) ){
  dna.adj[,i] <- round( dna.adj[ , i]/ips1copy$extract_rate[rownames(ips1copy) %in% colnames( dna.adj )[i]] , 0 )
}
dna.adj$taxonomy <- yangtze.clean$taxonomy ; head( dna.adj ) ; dim( dna.adj )
write.csv( dna.adj , file = "results/ngs.reads.DNAadj.sample.csv" , row.names = T)

##both DNA and PCR correction
dual.adj <- yangtze.clean[ , !grepl("taxonomy" , colnames( yangtze.clean ) )] ; head( dual.adj )
for(i in 1: ncol( dual.adj ) ){
  r1 <- ips1copy$extract_rate[rownames(ips1copy) %in% colnames( dual.adj )[i]]
  r2 <- pcramp$slope[rownames( pcramp ) %in% colnames( dual.adj )[i]]
  dual.adj[,i] <- round( 100* dual.adj[ , i] / (r1*r2) , 0 )
}
dual.adj$taxonomy <- yangtze.clean$taxonomy ; head( dual.adj ) ; dim( dual.adj )
write.csv( dual.adj , file = "results/ngs.reads.PCRadj.DNAadj.sample.csv" , row.names = T)

###################################################
##################step4: compare the adjustment effect
###################################################
topnspecies <- function( x , n){
  r <- 100*rowSums( x[ , !grepl("taxonomy" , colnames( x ))] )/sum( x[ , !grepl("taxonomy" , colnames( x ))] )
  names( r ) <- x$taxonomy
  r2 <- sort( r , decreasing = T)
  r2sp <- data.frame(species = c(gsub(".*; ","",names(r2)[1:n]) , "Other") , prop = c(as.numeric(r2[1:n]) , 100-sum(as.numeric(r2[1:n]))))
  return( r2sp )
}

topnmeansd <- function( x , n){
  x.clean <- t( x[ , !grepl( "taxonomy" , colnames(x))] )
  x.prop <- 100*t( x.clean/rowSums( x.clean ) )
  r <- rowMeans( x.prop )
  names( r ) <- x$taxonomy
  r2 <- sort( r , decreasing = T)
  r2sp <- data.frame(species = c(gsub(".*; ","",names(r2)[1:n]) , "Other") , prop = c(as.numeric(r2[1:n]) , 100-sum(as.numeric(r2[1:n]))))
  return( r2sp )
}

library(metabarcoding)
raw <- read.csv("results/ngs.reads.rawdata.csv" , row.names = 1) ; head( raw ) ; dim( raw )
raw <- raw[ grepl("Actinopteri" , raw$taxonomy) & ! grepl("; NA" , raw$taxonomy), ] ; dim( raw )
raw$taxonomy <- gsub("\\|" , "; " , raw$taxonomy ) ; head( raw )
raw <- raw[ , !grepl("B" , colnames( raw) )] ; head( raw ) ; dim( raw )
raw <- raw[ , grepl("T" , colnames( raw) )| grepl("taxonomy" , colnames( raw) )] ; head( raw ) ; dim( raw )

dna <- read.csv("results/ngs.reads.DNAadj.sample.csv" , row.names = 1 ) ; head( dna )
dna <- dna[ grepl("Actinopteri" , dna$taxonomy) & ! grepl("; NA" , dna$taxonomy), ] ; dim( dna )
dna$taxonomy <- gsub("\\|" , "; " , dna$taxonomy ) ; head( dna )
dna <- dna[ , !grepl("B" , colnames( dna) )] ; head( dna ) ; dim( dna )
dna <- dna[ , grepl("T" , colnames( dna) ) | grepl("taxonomy" , colnames( dna) )] ; head( dna ) ; dim( dna )

both <- read.csv( "results/ngs.reads.PCRadj.DNAadj.sample.csv" , row.names = 1 ) ; head( both ) ; dim( both )
both <- both[ grepl("Actinopteri" , both$taxonomy) & ! grepl("; NA" , both$taxonomy), ] ; dim( both )
both$taxonomy <- gsub("\\|" , "; " , both$taxonomy ) ; head( both )
both <- both[ , !grepl("B" , colnames( both) )] ; head( both ) ; dim( both )
both <- both[ , grepl("T" , colnames( both) ) | grepl("taxonomy" , colnames( both) )] ; head( both ) ; dim( both )

df.raw <- data.frame(  species = gsub(".*; ","",raw$taxonomy) , reads = rowSums(raw[, -ncol(raw)]) )
df.raw <- df.raw[ order( df.raw$reads , decreasing = T),] ; head( df.raw )
df.raw$prop <- round( 100*df.raw$reads/sum( df.raw$reads ) , 4) ; df.raw
dfraw <- df.raw[ 1:5, ] ; dfraw
dfrawtop5 <- rbind( dfraw , data.frame( species=  "other" , reads =sum( df.raw[ 6:nrow( df.raw),]$reads ) , prop = sum( df.raw[ 6:nrow( df.raw),]$prop ) ) )
dfrawtop5$group <- "raw" ; dfrawtop5

df.both <- data.frame(  species = gsub(".*; ","",both$taxonomy) , reads = rowSums(both[, -ncol(both)]) )
df.both <- df.both[ order( df.both$reads , decreasing = T),] ; head( df.both )
df.both$prop <- round( 100*df.both$reads/sum( df.both$reads ) , 4) ; head( df.both )
dfboth <- df.both[ 1:5, ] ; head(dfboth)
dfbothtop5 <- rbind( dfboth , data.frame( species=  "other" , reads =sum( df.both[ 6:nrow( df.both),]$reads ) , prop = sum( df.both[ 6:nrow( df.both),]$prop ) ) )
dfbothtop5$group <- "both" ; dfbothtop5

rawtab <- otutabs(otutable =  raw , sep = "; ")$order.summary
rawtab <- rawtab[ rawtab$Reads > 1000 | rawtab$OTU > 2, ] ; rawtab

dnatab <- otutabs(otutable =  dna , sep = "; ")$order.summary
bothtab <- otutabs(otutable =  both , sep = "; ")$order.summary

df <- merge( rawtab[ , colnames(rawtab) %in% c("order" , "Reads")] , bothtab[,colnames(bothtab) %in% c("order" , "Reads")] , by = "order")
colnames( df ) <- c("order" , "raw" , "calibration") ; df

df.prop <- data.frame( order = df$order , 
                       raw = df$raw/sum( df$raw),
                       calibration = df$calibration/sum( df$calibration))

df.prop <- melt( df.prop ) ; df.prop
df.prop$order <- factor(df.prop$order , levels = c("Perciformes","Anabantiformes" , "Clupeiformes" ,"Gobiiformes" ,"Cypriniformes", "Siluriformes"))

library(ggalluvial)
library(ggprism)
pvs <- ggplot( data = df.prop , aes(x = variable , y = value , fill = order) ) +
  geom_bar( stat = "identity" , position = "fill") +
  geom_flow( aes(alluvium = order ), alpha = 0.5) + 
  scale_fill_manual( values = c("Anabantiformes" = "#db6968",
                                "Clupeiformes" ="#4d97cd",
                                "Cypriniformes" = "#f8984e",
                                "Gobiiformes" ="#459943",
                                "Perciformes" = brewer.pal(9, "PuBu")[7],
                                "Siluriformes" = brewer.pal(9, "PuBu")[5] ) 
                     )+
  theme_light()+
  labs( x = "Reads" , y = "Proportion")+
  guides(fill = guide_legend(title = ""))+
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
pvs
###########top10 species
raw10 <- topnspecies(raw , 10)
both10 <- topnspecies(both , 10)
top10 <- merge( raw10 , both10 , by = "species" )
colnames( top10 ) <- c("species" , "Raw" , "Calibrated")
df.top10 <- melt( top10 ) ; head( df.top10 )
# 绘制分组柱形图
p.top10sp.reads <- ggplot(df.top10, aes(x = species, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs( x = NULL, y = "Relative abundance(%)", fill = NULL)+
  scale_fill_manual( values = c("Raw" = "#db6968",
                                "Calibrated" ="#4d97cd" ) )+
  theme_bw()+
  theme( strip.text = element_text(size = 12),
         axis.line = element_line(color = "black",size = 0.4),
         axis.text.y = element_text(color = "black",size = 8),
         axis.text.x = element_text(color = "black",size = 10 , angle = 45 , hjust = 1),
         axis.title = element_text(color = "black",size = 12),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
         legend.position = c(0.25, 0.6),
         legend.text = element_text(size = 8),
         legend.title = element_text(size = 8),
         legend.key.size = unit(1, "lines"))
p.top10sp.reads

# mean reads for each site
#define a function for site combine
site.combine <- function( otutable ){
  site <- unique( gsub("\\..*","",colnames( otutable )[grepl("YCS" , colnames( otutable ) )]) ) ; site
  otu <- otutable[ , !grepl( "taxonomy" , colnames( otutable ) ) ]
  mt <- data.frame( matrix( ncol = length( site ) , nrow = nrow( otutable ) ) )
  for( i in 1: length( site ) ){
    if( sum(grepl(site[i] , colnames( otu ) ) ) > 1 ){
      mt[,i] <- round( apply( otu[ , grepl(site[i] , colnames( otu ) )] , 1 , mean) , 0)
    } else {
      mt[,i] <- otu[ , grepl(site[i] , colnames( otu ) )]
    }
  }
  colnames( mt ) <- site ; rownames( mt ) <- rownames( otu )
  mt$taxonomy <- otutable$taxonomy
  return( mt )
}

tax.split <- function( otutable , level){
  tax <- sapply( strsplit( otutable$taxonomy , "; ") , "[" , level)
}

raw <- read.csv("results/ngs.reads.rawdata.csv" , row.names = 1) ; head( raw ) ; dim( raw )
raw <- raw[ grepl("Actinopteri" , raw$taxonomy) & ! grepl("; NA" , raw$taxonomy), ] ; dim( raw )
raw$taxonomy <- gsub("\\|" , "; " , raw$taxonomy ) ; head( raw )
raw <- raw[ , !grepl("B" , colnames( raw) )] ; head( raw ) ; dim( raw )
raw <- raw[ , grepl("T" , colnames( raw) )| grepl("taxonomy" , colnames( raw) )] ; head( raw ) ; dim( raw )

both <- read.csv( "results/ngs.reads.PCRadj.DNAadj.sample.csv" , row.names = 1 ) ; head( both ) ; dim( both )
both <- both[ grepl("Actinopteri" , both$taxonomy) & ! grepl("; NA" , both$taxonomy), ] ; dim( both )
both$taxonomy <- gsub("\\|" , "; " , both$taxonomy ) ; head( both )
both <- both[ , !grepl("B" , colnames( both) )] ; head( both ) ; dim( both )
both <- both[ , grepl("T" , colnames( both) ) | grepl("taxonomy" , colnames( both) )] ; head( both ) ; dim( both )

raw.site <- site.combine( raw ) ; head( raw.site )
both.site <- site.combine( both ) ; head( both.site  )

df.site.reads <- data.frame(
  site = colnames( raw.site )[1:19],
  raw = as.numeric( colSums(raw.site[ , !grepl("taxonomy" , colnames( raw.site))]) ),
  calibarted = as.numeric( colSums(both.site[ , !grepl("taxonomy" , colnames( both.site))]) )
)

df.site.reads.com <- data.frame(
  site = colnames( raw.site )[1:19],
  raw = as.numeric( colSums(raw.site[ !grepl("Leiocassis longirostris" ,raw.site$taxonomy), !grepl("taxonomy" , colnames( raw.site))]) ),
  calibarted = as.numeric( colSums(both.site[ !grepl("Leiocassis longirostris" ,both.site$taxonomy), !grepl("taxonomy" , colnames( both.site))]) )
)

write.csv(df.site.reads , file = "results/site.reads.sum.csv" , row.names = F)
write.csv(df.site.reads.com , file = "results/site.reads.sum.no.ll.csv" , row.names = F)

##################top 10 fish 
site <- unique( gsub( "\\..*", "" , gsub( "YCS" , "" , colnames( raw ) ) ) )
site <- site[ ! grepl("taxonomy" , site ) ] ; site
dfs <- c()
for(i in 1: nrow( raw ) ){
  tem <- data.frame( matrix( ncol = 7 , nrow =  length(site)))
  colnames( tem ) <- c("Site" , "Raw" , "RawSE" , "DNA" , "DNASE" , "PCR" , "PCRSE")
  tem$Site <- as.numeric( site )
  for(j in 1: length( site ) ){
    tem1 <- as.numeric( raw[i , grepl( site[j] , colnames( raw ) ) ] )
    tem2 <- as.numeric( dna[i , grepl( site[j] , colnames( dna ) ) ] )
    tem3 <- as.numeric( both[i , grepl( site[j] , colnames( both ) ) ] )
    tem[ j , 2] <- mean( tem1 )
    tem[ j , 3] <- sd( tem1 )/sqrt( length( tem1 ) )
    tem[ j , 4] <- mean( tem2 )
    tem[ j , 5] <- sd( tem2 )/sqrt( length( tem2 ) )
    tem[ j , 6] <- mean( tem3 )
    tem[ j , 7] <- sd( tem3 )/sqrt( length( tem3 ) )
  }
  tem$species = gsub(".*; ","", raw$taxonomy[i] )
  tem$taxonomy <- raw$taxonomy[i]
  tem[ is.na(tem) ] <- 0
  dfs <- rbind( dfs , tem )
}
##
temfunc <- function( temdf ){
  tem1 <- data.frame( reads = rep( temdf$Site , sqrt(temdf$Raw) ) , type = "Raw")
  tem2 <- data.frame( reads = rep( temdf$Site , sqrt(temdf$DNA) ) , type = "DNA")
  tem3 <- data.frame( reads = rep( temdf$Site , sqrt(temdf$PCR) ) , type = "PCR")
  return( rbind( tem1 , tem2 , tem3 ) )
}

splist <- raw10$species[1:10]
dfspsite <- c()
for( i in 1:length( splist )){
  temdf <- dfs[ dfs$species == splist[i], ]
  dfspsite <- rbind( dfspsite , temdf )
}

p3 <- ggplot(  dfspsite , aes( x = as.numeric(Site) ) ) + 
  geom_line(aes(y = Raw , colour="Site") , colour=brewer.pal(7,"Set1")[2] , size=1 , lty = 20 )+
  geom_line(aes(y =  PCR * max(dfspsite$Raw) / max( dfspsite$PCR) , colour="Site") , colour=brewer.pal(7,"Set1")[1],size=1 ) +
  geom_errorbar(aes(ymin=Raw-RawSE, ymax=Raw+RawSE), width=0.1,color=brewer.pal(7,"Set1")[2],position=position_dodge(0.2))+
  geom_errorbar(aes(ymin=(PCR-PCRSE) * max(dfspsite$Raw) / max( dfspsite$PCR) ,
                    ymax=(PCR+PCRSE) * max(dfspsite$Raw) / max( dfspsite$PCR) ), width=0.1,color=brewer.pal(7,"Set1")[1],position=position_dodge(0.2))+
  scale_x_continuous( name= NULL,breaks = c(1:19),labels = 1:19 )+
  scale_y_continuous( name='Raw reads' ,
                      sec.axis = sec_axis(~./(max(dfspsite$Raw) / max( dfspsite$PCR)),
                                          name = "Calibrated reads") 
                      )+
  theme_light() +
  facet_wrap(vars(species) , scales = "free" , ncol = 2)+
  theme(
    legend.position="none",
    axis.text.x = element_text(size=10),
    axis.text=element_text(size=10),
    axis.title.x = element_text(size=0),
    plot.title = element_text(hjust = 0.5, size= 10)
  )
p3

################PCA
library(FactoMineR)
library(factoextra)
raw.clean <- t( raw[!grepl("Leiocassis longirostris" , raw$taxonomy) , !grepl("taxonomy" , colnames( raw ))] )
colnames( raw.clean ) <- gsub( ".*; ","",raw$taxonomy[!grepl("Leiocassis longirostris" , raw$taxonomy)])
raw.pc <- 100*raw.clean/rowSums( raw.clean )
both.clean <- t( both[ !grepl("Leiocassis longirostris" , both$taxonomy), !grepl("taxonomy" , colnames(both))] )
colnames( both.clean ) <- gsub( ".*; ","",both$taxonomy[!grepl("Leiocassis longirostris" , both$taxonomy)])
#sitegroup <- as.numeric( gsub("YCS","", gsub( "\\..*","",rownames( raw.clean ) ) ) )
sitegroup <- read.csv("site.group.csv")
raw.site.clean <- t(raw.site[ !grepl("Leiocassis longirostris" , raw.site$taxonomy), !grepl("taxonomy" , colnames( raw.site ) )] )
both.site.clean <- t(both.site[ !grepl("Leiocassis longirostris" , both.site$taxonomy), !grepl("taxonomy" , colnames( both.site ) )] )
colnames( raw.site.clean ) <- gsub( ".*; ","",raw.site$taxonomy[!grepl("Leiocassis longirostris" , raw.site$taxonomy)])
colnames( both.site.clean ) <- gsub( ".*; ","",both.site$taxonomy[!grepl("Leiocassis longirostris" , both.site$taxonomy)])

raw.pca <- PCA( raw.site.clean, graph = T , ncp = 5, scale.unit = F)
raw.eig <- get_eigenvalue(raw.pca) ; head( raw.eig)
raw.site <- as.data.frame( get_pca_ind( raw.pca )$coord ) ; head( raw.site )
raw.site$site <- rownames( raw.site ) ; head( raw.site )
raw_site <- merge( raw.site , sitegroup , by = "site") ; head( raw_site )
raw.spe <- raw.pca$var$coord ; head( raw.spe )
raw.spe.abund1 <- raw.spe[ order(abs(raw.spe[,1]) , decreasing = T) , c(1:2) ][1:3, ] ; raw.spe.abund1
raw.spe.abund2 <- raw.spe[ order(abs(raw.spe[,2]) , decreasing = T) , c(1:2) ][1:3, ] ; raw.spe.abund2
raw.spe.abund <- round( as.data.frame( raw.spe[ rownames(raw.spe) %in% c( rownames( raw.spe.abund1) , rownames( raw.spe.abund2)), 1:2]) , 2) ; raw.spe.abund
raw.spe.abund$x = 0
raw.spe.abund$y = 0
raw.spe.abund <- raw.spe.abund/(max(raw.spe.abund)/(max( raw.site[, colnames( raw.site ) %in% c("Dim.1" , "Dim.2" ) ] )*0.7))
raw.spe.abund
# 执行 PERMANOVA 分析
raw_distance_bray <- vegdist( raw.site.clean, method = "bray")
raw_permanova_bray <- adonis2(raw_distance_bray ~ sitegroup$group ) ; raw_permanova_bray

p.pca.raw <- ggplot(data = raw_site , aes( x = Dim.1 , y = Dim.2  ) ) +
  geom_point( data = raw_site , aes(color = group ), size = 4) +
  stat_ellipse( aes(color = group , fill = group ) , level = 0.6 , geom = "polygon", alpha = 0.2 , linetype = "dashed") +
  geom_text(data = raw.spe.abund, aes(x = Dim.1, y = Dim.2, label = rownames(raw.spe.abund) ), vjust = -1, hjust = 1, size=2) +
  scale_color_manual(values = c("#db6968","#4d97cd","#98d09d" ,"#f8984e" )) +
  scale_fill_manual(values = c("#db6968","#4d97cd","#98d09d" ,"#f8984e" )) +
  theme_bw()+
  geom_segment(data = raw.spe.abund, aes(x = x, y = y, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm") ), color = "blue")+
  labs(x = paste0("PC1 (" , round( raw.eig[1,2] , 2) , "%)") ,
       y = paste0("PC2 (" , round( raw.eig[2,2] , 2) , "%)") ,fill = NULL , color = NULL)+
  annotate("text", x = Inf, y = Inf, label = paste("PERMANOVA R²:", round(raw_permanova_bray$R2[1],2), " p-value:", round(raw_permanova_bray$`Pr(>F)`[1], 4)),
           hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic")+
  theme(
    strip.text = element_text(size = 12),
    axis.line = element_line(color = "black",size = 0.4),
    axis.text.y = element_text(color = "black",size = 10),
    axis.title = element_text(color = "black",size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
    legend.position = c(0.8 ,0.7),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"))
p.pca.raw

both.pca <- PCA(  both.site.clean, graph = T , ncp = 5, scale.unit = F)
both.eig <- get_eigenvalue(both.pca) ; head( both.eig)
both.site <- as.data.frame( get_pca_ind( both.pca )$coord ) ; head( both.site )
both.site$site <- rownames( both.site ) ; head( both.site )
both_site <- merge( both.site , sitegroup , by = "site") ; head( both_site )
both.spe <- both.pca$var$coord ; head( both.spe )
both.spe.abund1 <- both.spe[ order(abs(both.spe[,1]) , decreasing = T) , c(1:2) ][1:3, ] ; both.spe.abund1
both.spe.abund2 <- both.spe[ order(abs(both.spe[,2]) , decreasing = T) , c(1:2) ][1:3, ] ; both.spe.abund2
both.spe.abund <- round( as.data.frame( both.spe[ rownames(both.spe) %in% c( rownames( both.spe.abund1) , rownames( both.spe.abund2)), 1:2]) , 2) ; both.spe.abund
both.spe.abund$x = 0
both.spe.abund$y = 0
both.spe.abund <- both.spe.abund/(max(both.spe.abund)/(max( both.site[, colnames( both.site ) %in% c("Dim.1" , "Dim.2" ) ] )*0.5))
both.spe.abund
# 执行 PERMANOVA 分析
both_distance_bray <- vegdist( both.site.clean, method = "bray")
both_permanova_bray <- adonis2(both_distance_bray ~ sitegroup$group ) ; both_permanova_bray

p.pca.both <- ggplot(data = both_site , aes( x = Dim.1 , y = Dim.2  ) ) +
  geom_point( data = both_site , aes(color = group ), size = 4) +
  stat_ellipse( aes(color = group , fill = group ) , level = 0.6 , geom = "polygon", alpha = 0.2 , linetype = "dashed") +
  geom_text(data = both.spe.abund, aes(x = Dim.1, y = Dim.2, label = rownames(both.spe.abund) ), vjust = -1, hjust = 1, size=2) +
  scale_color_manual(values = c("#db6968","#4d97cd","#98d09d" ,"#f8984e" )) +
  scale_fill_manual(values = c("#db6968","#4d97cd","#98d09d" ,"#f8984e" )) +
  theme_bw()+
  geom_segment(data = both.spe.abund, aes(x = x, y = y, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm") ), color = "blue")+
  labs(x = paste0("PC1 (" , round( both.eig[1,2] , 2) , "%)") ,
       y = paste0("PC2 (" , round( both.eig[2,2] , 2) , "%)") ,fill = NULL , color = NULL)+
  annotate("text", x = Inf, y = Inf, label = paste("PERMANOVA R²:", round(both_permanova_bray$R2[1],2), " p-value:", round(both_permanova_bray$`Pr(>F)`[1], 4)),
           hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic")+
  theme(
    strip.text = element_text(size = 12),
    axis.line = element_line(color = "black",size = 0.4),
    axis.text.y = element_text(color = "black",size = 10),
    axis.title = element_text(color = "black",size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
    legend.position = c(0.8 ,0.7),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"))
p.pca.both

p.blank <- ggplot() + theme_void()

ptem <- ggarrange( p.blank, p.top10sp.reads , widths = c(1,4)) ; ptem
p4 <- ggarrange( ptem , p.pca.raw, p.pca.both , ncol = 1 , labels = c("d" ,"e" , "f") , heights = c(1, 1 ,1)) ; p4
p5 <- ggarrange(p4 ,p3 ,labels = c("d" ,"g")  , widths = c(1, 2.3))

########make the sampling map
load("map/chinadf.rda")
china_shp <- "map/china.json"
nine <- "map/nahai_nine.json"
river <- "map/yangtze_river.json"
river2 <- "map/yangtze_river2.json"
china <- sf::read_sf(china_shp)
nine_line <- sf::read_sf(nine)
river <- sf::read_sf(river)
river2 <- sf::read_sf(river2)
library(ggspatial)
library(cowplot)

map <- ggplot() + 
  geom_sf(data = china_shp, fill=NA ) + 
  geom_sf(data = nine_line, color='gray50', size=.8)+
  geom_sf(data = river, color= "blue3") + 
  geom_sf(data = river2, color= "blue3") +
  coord_sf(ylim = c(-2387082 , 1654989), crs="+proj=laea +lat_0=40 +lon_0=104")+
  annotation_scale(location = "bl", width_hint = 0.1) +
  annotation_north_arrow(location = "tr", 
                         which_north = "false",
                         style = north_arrow_fancy_orienteering(text_size = 6),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"))+
  theme_minimal()+
  theme(aspect.ratio = 1,
        text = element_text(size = 10,face = "bold"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey80",size=.2),
        legend.key = element_rect(fill = "white"),
        legend.position = "bottom",
        plot.margin=unit(c(0,0,0,0),"mm")
  )
map

nine_map <- ggplot() +
  geom_sf(data = china, fill='NA') + 
  geom_sf(data = nine_line, color='gray70', size=1.)+
  coord_sf(ylim = c(-4028017,-1777844), xlim = c(117131.4, 2115095),crs="+proj=laea +lat_0=40 +lon_0=104")+
  theme_minimal()+
  theme(
    aspect.ratio = 1, #调节长宽比
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill=NA),
    plot.margin=unit(c(0,0,0,0),"mm"))

p.blank <- ggplot() + theme_void()

china_map <- ggarrange(map , ggarrange(p.blank , nine_map , p.blank, ncol = 1 , heights = c(1,1,0.02) ) , p.blank, ncol = 3, widths = c(1.7,0.5,0.4) )

# 包导入
library(ggspatial)
library(ggplot2)
library(sf)
library(pacman)
p_load(sf, dplyr, ggplot2)

# 读取shp文件，导入底图
river_shp <- st_read("map/Export_Output.shp")
river_shp <- st_transform(river_shp, 4326)
# 读取点位信息和对应数据
points_df <- read.csv(file="map/site.csv",fileEncoding = "GBK"); points_df
points_df$lon_copy <- points_df$longitude
points_df$lat_copy <- points_df$latitude
# 将数据框转换为 sf 对象，并指定坐标参考系统
points_sf <- st_as_sf(points_df, coords = c("longitude", "latitude"), crs = 4326)
df.points_sf <- melt( points_sf , measure.vars = c("Rawdata" ,"Calibrated" ,"raw_noll" ,"calibrat_noll"))
df.points_sf.all <- df.points_sf[ df.points_sf$variable %in% c("Rawdata" ,"Calibrated"), ]
df.points_sf.no <- df.points_sf[ !df.points_sf$variable %in% c("Rawdata" ,"Calibrated"), ]
df.points_sf.no$variable <- gsub("raw_noll","Rawdata",df.points_sf.no$variable)
df.points_sf.no$variable <- gsub("calibrat_noll","Calibrated",df.points_sf.no$variable)
# 给河流数据添加一个唯一标识符
river_shp$river_id <- seq_along(river_shp$geometry)
# 合并点位数据到河流数据
river_with_species <- st_join(river_shp, points_sf, join = st_nearest_feature)
library(reshape2 )
colnames(river_with_species)[2] <- "type"
river_with_species$type <- gsub("长江","river",river_with_species$type)
river_with_species$type <- gsub("海洋","marine",river_with_species$type)

site_type <- data.frame(start_x = c(118.537 , 119.370 , 119.852 ,121.044),
                        end_x = c(119.074 , 119.797 , 120.819, 121.776) , 
                        start_y = 32.425 , end_y = 32.425 , variable = "Calibrated")

site_type2 <- data.frame(start_x = c(118.537 , 119.370 , 119.852 ,121.044, 119.074 , 119.797 , 120.819, 121.776),
                        end_x = c(118.537 , 119.370 , 119.852 ,121.044, 119.074 , 119.797 , 120.819, 121.776) , 
                        start_y = 32.425 , end_y = 32.405 , variable = "Calibrated")

site_type_label <- data.frame(x = (site_type$start_x+site_type$end_x)/2 ,
                              y = 32.525 , 
                              label = c("Nanjing" , "Yangzhou" ,"Taizhou" ,"Nantong"), 
                              variable = "Calibrated")
# 绘制地图并根据“类型”列填充不同颜色，同时添加点并根据size = raw的值调整点大小
pmapa <- ggplot() +
  geom_sf(data = river_with_species , aes(fill = type)) +
  geom_point(data = df.points_sf.all, aes(x = lon_copy, y = lat_copy, size = value ), color = "red", alpha = 0.6) +
  scale_fill_manual(values = c("river" = "#cccccc", "marine" = "#97dbf2")) +
  geom_text(data = df.points_sf.all, aes(x = lon_copy, y = lat_copy, label = site), vjust = -0.2 ,
            nudge_y = -0.1, size = 2, color = "black") +  # 添加采样点名称
  guides(fill = FALSE) +
  geom_segment(data = site_type , aes(x = start_x , y = start_y , xend = end_x, yend = end_y ), linetype = 1 , lwd = 0.4 ) +
  geom_segment(data = site_type2 , aes(x = start_x , y = start_y , xend = end_x, yend = end_y ), linetype = 1 , lwd = 0.4 ) +
  geom_text(data = site_type_label, aes(x = x, y = y, label = label), vjust = -0.5 ,
            nudge_y = -0.1, size = 3, color = "black") +  # 添加点位分组
  scale_size_continuous(range = c(1, 10) ) +  # 调整点大小范围
  facet_wrap( vars(variable) , ncol = 1 , strip.position = "right" ) +
  theme_minimal() +
  labs(title = "", fill = "type", size = "Reads") +
  theme(strip.text = element_text(size = 12 ),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.9,0.2),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.2, "lines") ) +
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         style = north_arrow_fancy_orienteering(text_size = 6),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm")) + # 添加指北针
  annotation_scale(location = "bl", width_hint = 0.1) + # 添加比例尺
  coord_sf(xlim = c(118.5, 122.5), ylim = c(31.2, 32.5))  # 设置经度和纬度显示范围

pmapa

pmapb <- ggplot() +
  geom_sf(data = river_with_species , aes(fill = type)) +
  geom_point(data = df.points_sf.no, aes(x = lon_copy, y = lat_copy, size = value ), color = "red", alpha = 0.6) +
  scale_fill_manual(values = c("river" = "#cccccc", "marine" = "#97dbf2")) +
  geom_text(data = df.points_sf.no, aes(x = lon_copy, y = lat_copy, label = site), vjust = -0.2 ,
            nudge_y = -0.1, size = 2, color = "black") +  # 添加采样点名称
  guides(fill = FALSE) +
  geom_segment(data = site_type , aes(x = start_x , y = start_y , xend = end_x, yend = end_y ), linetype = 1 , lwd = 0.4 ) +
  geom_segment(data = site_type2 , aes(x = start_x , y = start_y , xend = end_x, yend = end_y ), linetype = 1 , lwd = 0.4 ) +
  geom_text(data = site_type_label, aes(x = x, y = y, label = label), vjust = -0.5 ,
            nudge_y = -0.1, size = 3, color = "black") +  # 添加点位分组
  scale_size_continuous(range = c(1, 10) ) +  # 调整点大小范围
  facet_wrap( vars(variable) , ncol = 1 , strip.position = "right" ) +
  theme_minimal() +
  labs(title = "", fill = "type", size = "Reads") +
  theme(strip.text = element_text(size = 12 ),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.9,0.2),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.2, "lines") ) +
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         style = north_arrow_fancy_orienteering(text_size = 6),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm")) + # 添加指北针
  annotation_scale(location = "bl", width_hint = 0.1) + # 添加比例尺
  coord_sf(xlim = c(118.5, 122.5), ylim = c(31.2, 32.5))  # 设置经度和纬度显示范围

pmapb


p6 <- ggarrange( china_map , pmapa , pmapb , ncol = 1 , heights = c(0.7, 1, 1) , labels = c("a" , "b" , "c") )

pdf("figures/Figure5.pdf" , width = 20, height = 11 )
ggarrange( p6 , p5 , widths = c(1.2, 3) )
dev.off()

          

