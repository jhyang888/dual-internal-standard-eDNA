####load necessary R packkages
library(xlsx)
library(plyr)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggprism)
library(RColorBrewer)
############################step1: check the basic NGS reads
############################step2: PCR amplification efficiency calibration
############################step3: DNA extraction efficiency correction
############################step4: only correct DNA extract efficiency
############################step5: plot to display correction effect
############################step6: species qPCR copies calibration
############################step7: plot summary
## set wd
setwd("lab-microcosm/")
############################
############################step1: check the basic NGS reads
## read NGS OTUs reads table
lab <- read.csv("ngs.reads.csv" , row.names = 1) ; head( lab )
lab.clean <- lab[ !grepl("Ctrl" , rownames( lab ) ), ] ; head(lab.clean)##remove control sample
## determine the filter volume
lab$volume <- gsub( "200","", gsub( "150","", gsub( "100" ,"", gsub("-.*" , "" , rownames( lab ) ) ) ) )
## species biomass group
lab$volume <- factor( lab$volume , levels = c("L" ,"M" , "H" , "Ctrl")) ; head( lab )
df.lab <- melt( lab ) ; head( df.lab )
df.lab$vector <- gsub("IPS2..*" ,"IPS2" , df.lab$variable)
my_comp <- list(c("L", "M"), c("L", "H"), c("L", "Ctrl"))
## IPS sequences composition
seqcom <- data.frame( colSums( lab.clean ) ) ; seqcom
colnames( seqcom ) <- "reads";seqcom
seqcom$categroy <- c( "IPS1","IPS2.1","IPS2.2","IPS2.3" , "MP" ,"HN","PD","CA","Other" );seqcom
seqcom$prop <- round(seqcom$reads/sum( seqcom$reads ) , 4)
seqcom$type <- "Total"; seqcom

seqcomp <- seqcom[1:4, ]
seqcomp$prop <- round(seqcomp$reads/sum(seqcomp$reads),4)
seqcomp$type <- "IPS";seqcomp

seqdf <- rbind( seqcom , seqcomp) ; seqdf
seqdf$categroy <- factor(seqdf$categroy , levels = unique(seqdf$categroy) )
seqdf$type <- factor( seqdf$type , levels = unique(seqdf$type) )

lb <- data.frame(type = c("Total" , "IPS") , count = c( sum(seqcom$reads) , sum( seqcomp$reads)))

library(ggalluvial)
pf <- ggplot( data = seqdf , aes(x = type , y = prop , fill = categroy) ) +
  geom_bar( stat = "identity" , position = "fill") +
  geom_flow( aes(alluvium = categroy ), alpha = 0.5) + 
  scale_fill_manual( values = c("IPS1" = "#db6968",
                                "IPS2.1" ="#4d97cd",
                                "IPS2.2" = "#f8984e",
                                "IPS2.3" ="#459943",
                                "MP" = "#98d09d" , 
                                "HN" = "#d7e698" , 
                                "PD" = "#fbf398" , 
                                "CA" = brewer.pal(9, "PuBu")[7],
                                "Other" = brewer.pal(9, "PuBu")[5] )
                     )+
  geom_text( data = lb , aes(x = type , y=1, label = paste0( count ) ) , inherit.aes = FALSE , vjust = -0.2 , size =3)+
  theme_light()+
  labs( x = "Reads" , y = "Proportion")+
  theme_prism() +
  guides(fill = guide_legend(title = ""))+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 10),        # 调整图例文本大小
        legend.title = element_text(size = 10),       # 调整图例标题大小
        legend.key.size = unit(1, "lines") )
pf

## fish sequence composition
library(dplyr)
dfseq <- data.frame( biomass = substr(rownames( lab.clean ) , 1, 1), reads = rowSums( lab.clean[ , 5:9]) )
summary <- dfseq %>%
  group_by( biomass ) %>%
  summarise(
    mean = mean(reads),
    se = sd(reads) / sqrt(n() )
  )
summary$biomass <- factor( summary$biomass , levels = c("L" , "M" , "H"))

pfish <- ggplot(summary, aes(x = biomass, y = mean, fill = biomass)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_text(aes(label = round(mean,0) ), vjust = -0.5, size = 3) +
  labs(x = "Biomass", y = "Fish reads") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e")) + 
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")
pfish

## IPS1 reads in different biomass groups
ips1reads <- df.lab[ (!df.lab$volume == "Ctrl") & df.lab$vector %in% "IPS1",]
p.ips1 <- ggplot(ips1reads , aes(x = volume , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = volume), pch = 21,size = 2, position = position_jitter(0.2))+
  #facet_wrap( .~variable , scales="free") +
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e")) + 
  labs(x = NULL,y = "IPS1 reads") +
  stat_compare_means(comparisons=list(c("L", "M"), c("M", "H"), c("L", "H")), method = "t.test",label="p.signif")+
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")
p.ips1 

## IPS2 Reads in different biomass groups
ips2reads <- df.lab[ (!df.lab$volume == "Ctrl") & df.lab$vector %in% "IPS2",]
p.ips2.1 <- ggplot( ips2reads , aes(x = variable , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = volume), pch = 21,size = 2, position = position_jitter(0.2))+
  #facet_wrap( .~vector , scales="free") +
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e")) + 
  labs(x = NULL,y = "IPS2 reads") +
  theme_prism() +
  guides(fill = guide_legend(title = "Biomass"))+
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 10),        # 调整图例文本大小
        legend.title = element_text(size = 10),       # 调整图例标题大小
        legend.key.size = unit(1, "lines") )
p.ips2.1

p.ips2.2 <- ggplot(df.lab[(!df.lab$volume == "Ctrl") & df.lab$vector %in% c("IPS2" ),] , aes(x = volume , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = volume), pch = 21,size = 2, position = position_jitter(0.2))+
  facet_wrap( .~variable , scales="free") +
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e")) + 
  labs(x = "Biomass group",y = "Reads") +
  stat_compare_means(comparisons=list(c("L", "M"), c("M", "H"), c("L", "H")), method="t.test", label="p.signif")+
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")
p.ips2.2

###########################
############################step2: PCR amplification efficiency calibration
## The IPS2 linear fitting equation corrected the PCR amplification efficiency of the samples
## correction by each sample IPS2 fitting model
## read NGS reads table
lab <- read.csv("ngs.reads.csv" , row.names = 1) ; head( lab )
##remove control sample
lab.clean <- lab[ !grepl("Ctrl" , rownames( lab ) ), ] ; head(lab.clean)
##find IPS2 reads
lab.ips2 <- lab[ !grepl("Ctrl" , rownames( lab )) , grepl("IPS2" , colnames( lab ) ) ] ; head( lab.ips2 )
lab.ips2 <- log(lab.ips2) ; lab.ips2
lab.ips2$biomass <- substr( rownames( lab.ips2 ) , 1, 1) ; head( lab.ips2 )
lab.ips2$volume <- substr( rownames( lab.ips2 ) , 2, 4) ; head( lab.ips2 )
##calculate the liner model
rho <- c() ; pvalue <- c() ; formula <- c()
ips1adj <- c() ; mpadj <- c() ; hnadj <- c() ; pdadj <- c() ; caadj <- c() 
for( i in 1: nrow( lab.ips2 ) ){
  x <- c(22, 65.5, 134)###ips2 concentration
  y <- c(as.numeric( lab.ips2[i,1:3] ) )
  temd <- data.frame( conc = x , reads = y)
  model <- lm( y ~ x )
  cts <- summary( model )
  rho[i] <- cts$r.squared
  pvalue[i] <- round(cts$coefficients[2,4] , 4)
  formula[i] <- paste( "y =",round(coef(model)[1] , 2) , "+", paste0(round(coef(model)[2] , 2)  , "x") )
  if(cts$r.squared >= 0.9 ){
    ips1adj[i] <- round( log(lab$IPS1[i]-coef(model)[1][[1]])/coef(model)[2][[1]] ,0)
    mpadj[i] <- round( (log(lab$Mylopharyngodon_piceus[i])-coef(model)[1][[1]])/coef(model)[2][[1]],0)
    hnadj[i] <- round( (log(lab$Hypophthalmichthys_nobilis[i])-coef(model)[1][[1]])/coef(model)[2][[1]] , 0)
    pdadj[i] <- round( (log(lab$Paramisgurnus_dabryanus[i])-coef(model)[1][[1]])/coef(model)[2][[1]] ,0)
    caadj[i] <- round( (log(lab$Channa_argus[i])-coef(model)[1][[1]])/coef(model)[2][[1]],0)
  } else {
    ips1adj[i] <- log(lab$IPS1[i])
    mpadj[i] <- log(lab$Mylopharyngodon_piceus[i])
    hnadj[i] <- log(lab$Hypophthalmichthys_nobilis[i])
    pdadj[i] <- log(lab$Paramisgurnus_dabryanus[i])
    caadj[i] <- log(lab$Channa_argus[i])
  }
}

df1 <- data.frame(R = rho , pvalue = pvalue , formula = formula , IPS1 = ips1adj , MP = mpadj , HN = hnadj , PD = pdadj , CA = caadj)
dfadj <- cbind( lab.ips2 , df1 ) ; dfadj
dfadj <- dfadj[ order(rownames(dfadj)), ] ; dfadj

dir.create("results")
write.csv( dfadj , file = "results/ngs.reads.PCRadj.sample.csv" )

## show the PCR efficiency correct effects
## read NGS raw  data
lab <- read.csv("ngs.reads.csv" , row.names = 1) ; head( lab )
rawdata <- lab[ !grepl("Ctrl" , rownames( lab ) ), colnames(lab) %in% c("Mylopharyngodon_piceus" , "Hypophthalmichthys_nobilis" , "Paramisgurnus_dabryanus" , "Channa_argus" )]
colnames(rawdata) <- c( "MP" , "HN" , "PD" , "CA") ; rawdata
rawdata <- rawdata[order(rownames(rawdata)) ,] ; rawdata
grouptype <- data.frame(
  group = gsub( "200","",gsub( "150","",gsub( "100" ,"", gsub("-.*" , "" , rownames( rawdata  ) ) ) ) ),
  volume = substr( rownames( rawdata) , 2 , 4) )
pcr.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1)[,10:13] ; head( pcr.adj.sample )
pcr.adj.sample <- reshape2::melt(cbind( pcr.adj.sample , grouptype ) , variable.name = "species" , value.name = "reads")
pcr.adj.sample$biogroup <- paste( pcr.adj.sample$group , pcr.adj.sample$species) ; pcr.adj.sample 
pcr.adj.sample <- merge(pcr.adj.sample , spbiomas , by = "biogroup")
pcr.adj.sample$type <- "PCR" ; head( pcr.adj.sample ) ; head( pcr.adj.sample )

## correction by biomass group IPS2 fitting model
lab.ips2 <- lab[ !grepl("Ctrl" , rownames( lab )) , grepl("IPS2" , colnames( lab ) ) ] ; head( lab.ips2 )
lab.ips2$biomass <- substr( rownames( lab.ips2 ) , 1, 1) ; head( lab.ips2 )
lab.ips2$volume <- substr( rownames( lab.ips2 ) , 2, 4) ; head( lab.ips2 )

rho <- c() ; pvalue <- c() ; formula <- c()
ips1adj <- c() ; mpadj <- c() ; hnadj <- c() ; pdadj <- c() ; caadj <- c() 
for( i in 1: nrow( lab.ips2 ) ){
  ips2df <- data.frame(conc = c(22, 65.5, 134) , variable = c("IPS2.1" , "IPS2.2" , "IPS2.3"))
  tem <- lab.ips2[ lab.ips2$biomass== lab.ips2$biomass[i] &lab.ips2$volume== lab.ips2$volume[i], 1:3]
  temd <- merge( melt( tem ) , ips2df , by = "variable" )
  model <- lm( temd$value ~ temd$conc )
  cts <- summary( model )
  rho[i] <- cts$adj.r.squared
  pvalue[i] <- round(cts$coefficients[2,4] , 4)
  formula[i] <- paste( "y =",round(coef(model)[1] , 2) , "+", paste0(round(coef(model)[2] , 2)  , "x") )
  if( cts$coefficients[2,4] <= 0.1 ){
    ips1adj[i] <- round( (lab$IPS1[i]-coef(model)[1][[1]])/coef(model)[2][[1]] ,0)
    mpadj[i] <- round( (lab$Mylopharyngodon_piceus[i]-coef(model)[1][[1]])/coef(model)[2][[1]],0)
    hnadj[i] <- round( (lab$Hypophthalmichthys_nobilis[i]-coef(model)[1][[1]])/coef(model)[2][[1]] , 0)
    pdadj[i] <- round( (lab$Paramisgurnus_dabryanus[i]-coef(model)[1][[1]])/coef(model)[2][[1]] ,0)
    caadj[i] <- round( (lab$Channa_argus[i]-coef(model)[1][[1]])/coef(model)[2][[1]],0)
  } else {
    ips1adj[i] <- lab$IPS1[i]
    mpadj[i] <- lab$Mylopharyngodon_piceus[i]
    hnadj[i] <- lab$Hypophthalmichthys_nobilis[i]
    pdadj[i] <- lab$Paramisgurnus_dabryanus[i]
    caadj[i] <- lab$Channa_argus[i]
  }
}

df1 <- data.frame(R = rho , pvalue = pvalue , formula = formula , IPS1 = ips1adj , MP = mpadj , HN = hnadj , PD = pdadj , CA = caadj)
dfadj <- cbind( lab.ips2 , df1 ) ; dfadj
write.csv( dfadj , file = "results/ngs.reads.PCRadj.biomass.csv" )

## plot the ips2 linear fitting for biomass group
library( reshape2 )
df.ips2 <- data.frame(conc = c(22, 65.5, 134) , variable = c("IPS2.1" , "IPS2.2" , "IPS2.3"))
df.ips2 <- merge( melt(lab.ips2,id.vars = )  , df.ips2 , by = "variable") ; head( df.ips2 )
df.ips2$biomass <- factor( df.ips2$biomass , levels = c("L" , "M" , "H"))

ips2 <- ggplot(df.ips2 , aes(x = conc , y = value )) +
  geom_point(aes(fill = conc), pch = 21, size = 3)+
  labs(x = "IPS2 concentration", y = "IPS2 NGS Reads") +
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  facet_grid(volume~biomass , scales = "free")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
ips2

pdf("figures/IPS2.linear.fitting.biomass.pdf")
ips2
dev.off()
############################
############################step3: DNA extraction efficiency correction
library(dplyr)
## Correct the DNA extraction efficiency by qPCR
## read qPCR data
lab.qpcr1 <- read.csv("ips1.qpcr.valume100-150.csv") ; head( lab.qpcr1 )
lab.qpcr1 <- lab.qpcr1[ !(lab.qpcr1$sample == "L150-2.3" ), ] ; head( lab.qpcr1 ) # remove cv > 3.5%
lab.qpcr100 <- lab.qpcr1[!lab.qpcr1$type == "Standard curve", colnames(lab.qpcr1 ) %in% c("ct_value" , "type")]
lab.qpcr2 <- read.csv("ips1.qpcr.valume200.csv") ; head( lab.qpcr2 )
lab.qpcr2 <- lab.qpcr2[ !(lab.qpcr2$sample %in% c("ctrl_4.1" ,"ctrl_5.2" , "L200-3.3") ), ] ; head( lab.qpcr2 ) # remove cv > 3.5%
lab.qpcr200 <- lab.qpcr2[!lab.qpcr2$type == "Standard curve", colnames(lab.qpcr2 ) %in% c("ct_value" , "type")]
## calculate CT  mean
lab.qpcr100.ctmean <- lab.qpcr100 %>% group_by(type) %>% summarise(mean_value = mean(ct_value))
lab.qpcr100.ctmean <- as.data.frame( lab.qpcr100.ctmean )

lab.qpcr200.ctmean <- lab.qpcr200 %>% group_by(type) %>% summarise(mean_value = mean(ct_value))
lab.qpcr200.ctmean <- as.data.frame( lab.qpcr200.ctmean )
## linear fitting model 100 and 150 ml
stand1 <- lab.qpcr1[lab.qpcr1$type == "Standard curve" , ]
stand1$copies <- log10( stand1$copies)

## calculate IPS1 copy number 100 and 150 ml
model1 <- lm( stand1$copies ~ stand1$ct_value )
cts1 <- summary( model1 )
lab.qpcr100.ctmean$copies <- 10^(coef(model1)[1] + (lab.qpcr100.ctmean$mean_value * coef(model1)[2]))
lab.qpcr100.ctmean$logcopy <- log10( lab.qpcr100.ctmean$copies )
lab.qpcr100.ctmean$extract_rate <- lab.qpcr100.ctmean$copies/160.2
## plot the linear fitting model 100 and 150 ml
p.fitmodel100 <- ggplot(stand1 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = lab.qpcr100.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(lab.qpcr100.ctmean$logcopy), xend =lab.qpcr100.ctmean$mean_value[lab.qpcr100.ctmean$logcopy == min(lab.qpcr100.ctmean$logcopy)], yend = min(lab.qpcr100.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(lab.qpcr100.ctmean$logcopy), xend =lab.qpcr100.ctmean$mean_value[lab.qpcr100.ctmean$logcopy == max(lab.qpcr100.ctmean$logcopy)], yend = max(lab.qpcr100.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(lab.qpcr100.ctmean$mean_value), y = -Inf, xend = min(lab.qpcr100.ctmean$mean_value), yend = lab.qpcr100.ctmean$logcopy[lab.qpcr100.ctmean$mean_value == min(lab.qpcr100.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(lab.qpcr100.ctmean$mean_value), y = -Inf, xend = max(lab.qpcr100.ctmean$mean_value), yend = lab.qpcr100.ctmean$logcopy[lab.qpcr100.ctmean$mean_value == max(lab.qpcr100.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(lab.qpcr1$ct_value) , y = min(lab.qpcr100.ctmean$logcopy)+0.2 , label = round(min(lab.qpcr100.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(lab.qpcr1$ct_value) , y = max(lab.qpcr100.ctmean$logcopy)+0.2 , label = round(max(lab.qpcr100.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "IPS1 copies (log10)") +
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel100

## linear fitting model 200 ml
stand2 <- lab.qpcr2[lab.qpcr2$type == "Standard curve" , ]
stand2$copies <- log10( stand2$copies)

## calculate IPS1 copy number 200 ml
model2 <- lm( stand2$copies ~ stand2$ct_value )
cts2 <- summary( model2 )
lab.qpcr200.ctmean$copies <- 10^(coef(model2)[1] + (lab.qpcr200.ctmean$mean_value * coef(model2)[2]))
lab.qpcr200.ctmean$logcopy <- log10( lab.qpcr200.ctmean$copies )
lab.qpcr200.ctmean$extract_rate <- lab.qpcr200.ctmean$copies/160.2
## plot the linear fitting model 100 and 150 ml
p.fitmodel200 <- ggplot(stand2 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = lab.qpcr200.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(lab.qpcr200.ctmean$logcopy), xend =lab.qpcr200.ctmean$mean_value[lab.qpcr200.ctmean$logcopy == min(lab.qpcr200.ctmean$logcopy)], yend = min(lab.qpcr200.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(lab.qpcr200.ctmean$logcopy), xend =lab.qpcr200.ctmean$mean_value[lab.qpcr200.ctmean$logcopy == max(lab.qpcr200.ctmean$logcopy)], yend = max(lab.qpcr200.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(lab.qpcr200.ctmean$mean_value), y = -Inf, xend = min(lab.qpcr200.ctmean$mean_value), yend = lab.qpcr200.ctmean$logcopy[lab.qpcr200.ctmean$mean_value == min(lab.qpcr200.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(lab.qpcr200.ctmean$mean_value), y = -Inf, xend = max(lab.qpcr200.ctmean$mean_value), yend = lab.qpcr200.ctmean$logcopy[lab.qpcr200.ctmean$mean_value == max(lab.qpcr200.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(lab.qpcr2$ct_value) , y = min(lab.qpcr200.ctmean$logcopy)+0.2 , label = round(min(lab.qpcr200.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(lab.qpcr2$ct_value) , y = max(lab.qpcr200.ctmean$logcopy)+0.2 , label = round(max(lab.qpcr200.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "IPS1 copies (log10)") +
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel200

ps4 <- ggarrange( p.fitmodel100 , p.fitmodel200 , labels = c("a" , "b")) ; ps4

pdf("figures/lab.DNA.extract.rate.pdf" , width = 8, height = 4)
ps4
dev.off()

lab.dna.extract.rate <- rbind(lab.qpcr100.ctmean , lab.qpcr200.ctmean) ; lab.dna.extract.rate
lab.dna.extract.rate <- lab.dna.extract.rate[order(lab.dna.extract.rate$type) , ] ; lab.dna.extract.rate
write.csv( lab.dna.extract.rate , file = "results/lab.dna.extract.rate.qcpr.csv" , row.names = F)

#read DNA extraction calculated by qPCR
dnarate <- read.csv("results/lab.dna.extract.rate.qcpr.csv") ; head( dnarate )
reads.adj.biomass <- read.csv("results/ngs.reads.PCRadj.biomass.csv") ; head( reads.adj.biomass )
temdf <- merge(reads.adj.biomass[, colnames(reads.adj.biomass) %in% c("X" , "IPS1") ] , dnarate, by.x = "X" , by.y = "type" )
temdf$biomass <- substr( temdf$X , 1,1) ; temdf
temdf <- rbind( temdf , data.frame(X = "blank" , IPS1 = 0 , mean_value = 0, copies = 0 , logcopy = 0 , extract_rate = 0, biomass = "L") )

#plot the correlation between qPCR and NGS reads
p.dna <- ggplot(temdf , aes(y = IPS1 , x = extract_rate )) +
  geom_point(aes(fill = biomass), pch = 21, size = 3)+
  labs(y = "IPS1 reads", x = "DNA extraction efficiency" , fill = "Biomass") +
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.4, "cm") )+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value..,sep = "~~~")), formula = y ~ x, parse = TRUE, color = "black") 
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation( size = 4,label.y.npc = 1)

p.dna

#correct the DNA extraction efficiency by qPCR
reads.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1) ; head( reads.adj.sample )
reads.adj.sample <- reads.adj.sample[ order(rownames(reads.adj.sample)), ] ; rownames(reads.adj.sample)
reads.adj.sample.dna <- reads.adj.sample[ , colnames(reads.adj.sample) %in% c("IPS1" , "MP" , "HN" , "PD" , "CA")]/dnarate$extract_rate
write.csv(reads.adj.sample.dna , file = "results/ngs.reads.PCRadj.DNAadj.sample.qPCR.csv")

#Correct the DNA extraction efficiency by NGS reads
reads.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1) ; head( reads.adj.sample )
reads.adj.sample <- reads.adj.sample[ order(rownames(reads.adj.sample)), ] ; reads.adj.sample
reads.adj.sample.ngs <- reads.adj.sample[ , colnames(reads.adj.sample) %in% c( "MP" , "HN" , "PD" , "CA")]/(reads.adj.sample$IPS1/max(reads.adj.sample$IPS1))
write.csv(reads.adj.sample.ngs , file = "results/ngs.reads.PCRadj.DNAadj.sample.ngs.csv")

############################
############################step4: only correct DNA extract efficiency
##correct DNA efficiency by qPCR
## read raw NGS reads table
lab <- read.csv("ngs.reads.csv" , row.names = 1) ; head( lab )
lab.clean <- lab[ !grepl("Ctrl" , rownames( lab ) ), ] ; head( lab.clean )
lab.clean <- lab.clean[order(rownames(lab.clean)) , ] ; head( lab.clean )

#read qPCR DNA extraction data
dnarate <- read.csv("results/lab.dna.extract.rate.qcpr.csv") ; dnarate

lab.clean.adj.dna.qpcr <- lab.clean[ , colnames(lab.clean) %in% c("Mylopharyngodon_piceus" , "Hypophthalmichthys_nobilis" , "Paramisgurnus_dabryanus" , "Channa_argus" )]/dnarate$extract_rate
colnames(lab.clean.adj.dna.qpcr) <- c( "MP" , "HN" , "PD" , "CA")
write.csv(lab.clean.adj.dna.qpcr , file = "results/ngs.reads.DNAadj.sample.qpcr.csv")

##correct DNA efficiency by IPS1 NGS reads
## read raw NGS reads table
lab <- read.csv("ngs.reads.csv" , row.names = 1) ; head( lab )
lab.clean <- lab[ !grepl("Ctrl" , rownames( lab ) ), ]
lab.clean <- lab.clean[order( rownames(lab.clean)) , ] ; lab.clean

reads.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1) ; head( reads.adj.sample )
reads.adj.sample <- reads.adj.sample[order(rownames(reads.adj.sample)) ,] ; reads.adj.sample
#calculate DNA rate
reads.adj.sample$IPS1/max(reads.adj.sample$IPS1)

lab.clean.adj.dna.ngs <- lab.clean[ , colnames(lab.clean) %in% c("Mylopharyngodon_piceus" , "Hypophthalmichthys_nobilis" , "Paramisgurnus_dabryanus" , "Channa_argus" )]/(reads.adj.sample$IPS1/max(reads.adj.sample$IPS1))
colnames(lab.clean.adj.dna.ngs) <- c( "MP" , "HN" , "PD" , "CA") ; head(lab.clean.adj.dna.ngs)
lab.clean.adj.dna.ngs <- lab.clean.adj.dna.ngs[order(rownames(lab.clean.adj.dna.ngs)) , ] ; lab.clean.adj.dna.ngs
write.csv(lab.clean.adj.dna.ngs , file = "results/ngs.reads.DNAadj.sample.ngs.csv")
############################
############################step5: plot to display correction effect
spbiomas <- read.csv("lab.species.biomass.csv") ; spbiomas
## read raw NGS data
lab <- read.csv("ngs.reads.csv" , row.names = 1) ; head( lab )
rawdata <- lab[ !grepl("Ctrl" , rownames( lab ) ), colnames(lab) %in% c("Mylopharyngodon_piceus" , "Hypophthalmichthys_nobilis" , "Paramisgurnus_dabryanus" , "Channa_argus" )]
colnames(rawdata) <- c( "MP" , "HN" , "PD" , "CA") ; rawdata
rawdata <- rawdata[order(rownames(rawdata)) ,] ; rawdata
###group type
grouptype <- data.frame(
  group = gsub( "200","",gsub( "150","",gsub( "100" ,"", gsub("-.*" , "" , rownames( rawdata  ) ) ) ) ),
  volume = substr( rownames( rawdata) , 2 , 4) )
##combine
rawdata <- reshape2::melt(cbind( rawdata , grouptype ) , variable.name = "species" , value.name = "reads")
rawdata$biogroup <- paste( rawdata$group , rawdata$species) ; rawdata
rawdata <- merge(rawdata , spbiomas , by = "biogroup")
rawdata$type <- "Rawdata" ; head( rawdata )
## DNA correct only
dna.adj.qpcr <- read.csv("results/ngs.reads.DNAadj.sample.qpcr.csv" , row.names = 1) ; head( dna.adj.qpcr )
dna.adj.qpcr <- reshape2::melt(cbind( dna.adj.qpcr , grouptype ) , variable.name = "species" , value.name = "reads") ; head( dna.adj.qpcr )
dna.adj.qpcr$biogroup <- paste( dna.adj.qpcr$group , dna.adj.qpcr$species) ; dna.adj.qpcr
dna.adj.qpcr <- merge(dna.adj.qpcr , spbiomas , by = "biogroup")
dna.adj.qpcr$type <- "DNA calibrated" ; head( dna.adj.qpcr )
## PCR correct only
pcr.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1)[,10:13] ; head( pcr.adj.sample )
pcr.adj.sample <- reshape2::melt(cbind( pcr.adj.sample , grouptype ) , variable.name = "species" , value.name = "reads")
pcr.adj.sample$biogroup <- paste( pcr.adj.sample$group , pcr.adj.sample$species) ; pcr.adj.sample 
pcr.adj.sample <- merge(pcr.adj.sample , spbiomas , by = "biogroup")
pcr.adj.sample$type <- "PCR calibrated" ; head( pcr.adj.sample )
## both correct
both.adj.sample.qpcr <- read.csv("results/ngs.reads.PCRadj.DNAadj.sample.qPCR.csv" , row.names = 1)[,2:5] ; head( both.adj.sample.qpcr )
both.adj.sample.qpcr <- reshape2::melt(cbind( both.adj.sample.qpcr , grouptype ) , variable.name = "species" , value.name = "reads")
both.adj.sample.qpcr$biogroup <- paste( both.adj.sample.qpcr$group , both.adj.sample.qpcr$species) ; both.adj.sample.qpcr
both.adj.sample.qpcr <- merge(both.adj.sample.qpcr , spbiomas , by = "biogroup")
both.adj.sample.qpcr$type <- "PCR&DNA calibrated" ; head(both.adj.sample.qpcr)
####combine all table
dfall <- rbind(rawdata , dna.adj.qpcr , pcr.adj.sample , both.adj.sample.qpcr) ; head( dfall )
dfall$type <- factor(dfall$type , levels = unique(dfall$type) )
dfall$species <- factor(dfall$species , levels = c("MP" ,"PD" ,"HN" , "CA"))
p200 <- ggplot(dfall[dfall$volume == 200,], aes(x = biomass , y = reads )) +
  geom_point(aes(fill = species), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "NGS reads") +
  facet_grid( type~species , scales = "free" )+
  #theme_prism() +
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943")) + 
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), 
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3,label.y.npc = 1)
p200

p150 <- ggplot(dfall[dfall$volume == 150 ,], aes(x = biomass , y = reads )) +
  geom_point(aes(fill = species), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "NGS reads") +
  facet_grid( type~species , scales = "free" )+
  #theme_prism() +
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943")) + 
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), 
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3,label.y.npc = 1)
p150

p100 <- ggplot(dfall[dfall$volume == 100 ,], aes(x = biomass , y = reads )) +
  geom_point(aes(fill = species), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "NGS reads") +
  facet_grid( type~species , scales = "free" )+
  #theme_prism() +
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943")) + 
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), 
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3,label.y.npc = 1)
p100

pdf("figures/both.adj.100ml.pdf" , height = 9, width = 7)
p100
dev.off()

pdf("figures/both.adj.150ml.pdf" , height = 9, width = 7)
p150
dev.off()

pdf("figures/both.adj.200ml.pdf" , height = 9, width = 7)
p200
dev.off()

#####calculate r square between different filter volume
##new function for Rho collect
rhocal <- function(x , newlabel){
  vo <- unique(x$volume)
  dfr <- matrix(ncol = 3 , nrow = 4)
  for(i in 1:3){
    tem <- x[x$volume == vo[i],]
    sp <- unique( tem$species )
    for(j in 1: 4){
      tem2 <- tem[ tem$species == sp[j], ]
      ct <- cor.test( tem2$reads , tem2$biomass)
      dfr[j,i] <- ct$estimate
    }
  }
  colnames( dfr ) <- vo
  rownames( dfr ) <- unique( tem$species )
  dfrboth <- melt( dfr ) ; head( dfrboth )
  colnames( dfrboth ) <- c( "species" , "volume" , "Rho")
  dfrboth <- as.data.frame( dfrboth )
  dfrboth$label <- newlabel
  return( dfrboth )
}

dfr <- rbind( rhocal(rawdata , "Rawdata"), rhocal(both.adj.sample.qpcr , "DNA&PCR"), rhocal(dna.adj.qpcr , "DNA") , rhocal(pcr.adj.sample , "PCR") )
dfr$volume <- factor( dfr$volume , levels = c("100" , "150" , "200"))
dfr$label <- factor( dfr$label , levels = c("Rawdata" , "DNA" , "PCR" , "DNA&PCR"))
dfr$species <- factor( dfr$species , levels = c("MP" , "PD" , "HN" , "CA"))
head( dfr )
pr <- ggplot( dfr[ !dfr$label == "PCR", ] , aes(x = volume , y = Rho , group = volume) ) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = species), pch = 21,size = 2, position = position_jitter(0.2))+
  facet_wrap( .~label , nrow = 1) +
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943")) + 
  labs(x = "Filter volume (mL)",y = "Linear fitting Rho") +
  stat_compare_means(comparisons=list(c("100", "150"), c("150", "200"), c("100", "200")),  label="p.signif" )+
  theme_prism() +
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 8),        # 调整图例文本大小
        legend.title = element_text(size = 8),       # 调整图例标题大小
        legend.key.size = unit(0.8, "lines"))
pr
############################
############################step6: species qPCR copies calibration
## read species qPCR data and calculate species copy number
sp.qpcr <- read.csv("species.qpcr.valume200.csv") ; head( sp.qpcr )

## Paramisgurnus dabryanus
sp.qpcr.pd <- sp.qpcr[ sp.qpcr$target == "PD" & !sp.qpcr$type == "Standard curve" , colnames(sp.qpcr) %in% c("type" ,"ct_value")] ; sp.qpcr.pd
sp.qpcr.pd.ctmean <- sp.qpcr.pd %>% group_by(type) %>% summarise(mean_value = mean(ct_value))

stand.pd <- sp.qpcr[ sp.qpcr$target == "PD" & sp.qpcr$type == "Standard curve" , ] ; stand.pd
stand.pd$copies <- log10( stand.pd$copies) ; stand.pd
model.pd <- lm( stand.pd$copies ~ stand.pd$ct_value )

sp.qpcr.pd.ctmean$copies <- 10^(coef(model.pd)[1] + (sp.qpcr.pd.ctmean$mean_value * coef(model.pd)[2]))
sp.qpcr.pd.ctmean$logcopy <- log10( sp.qpcr.pd.ctmean$copies )

p.qpcr.pd <- ggplot(stand.pd , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = sp.qpcr.pd.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(sp.qpcr.pd.ctmean$logcopy), xend =sp.qpcr.pd.ctmean$mean_value[sp.qpcr.pd.ctmean$logcopy == min(sp.qpcr.pd.ctmean$logcopy)], yend = min(sp.qpcr.pd.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(sp.qpcr.pd.ctmean$logcopy), xend =sp.qpcr.pd.ctmean$mean_value[sp.qpcr.pd.ctmean$logcopy == max(sp.qpcr.pd.ctmean$logcopy)], yend = max(sp.qpcr.pd.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(sp.qpcr.pd.ctmean$mean_value), y = -Inf, xend = min(sp.qpcr.pd.ctmean$mean_value), yend = sp.qpcr.pd.ctmean$logcopy[sp.qpcr.pd.ctmean$mean_value == min(sp.qpcr.pd.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(sp.qpcr.pd.ctmean$mean_value), y = -Inf, xend = max(sp.qpcr.pd.ctmean$mean_value), yend = sp.qpcr.pd.ctmean$logcopy[sp.qpcr.pd.ctmean$mean_value == max(sp.qpcr.pd.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(stand.pd$ct_value) , y = min(sp.qpcr.pd.ctmean$logcopy)+0.1 , label = round(min(sp.qpcr.pd.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(stand.pd$ct_value) , y = max(sp.qpcr.pd.ctmean$logcopy)+0.1 , label = round(max(sp.qpcr.pd.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
    labs(x = "qPCR CT value", y = "PD copies (log10)") +
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
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3, color = "red",label.y.npc = 1)
p.qpcr.pd

## Mylopharyngodon piceus
sp.qpcr.mp <- sp.qpcr[ sp.qpcr$target == "MP" & !sp.qpcr$type == "Standard curve" , colnames(sp.qpcr) %in% c("type" ,"ct_value")] ; sp.qpcr.mp
sp.qpcr.mp.ctmean <- sp.qpcr.mp %>% group_by(type) %>% summarise(mean_value = mean(ct_value))

stand.mp <- sp.qpcr[ sp.qpcr$target == "MP" & sp.qpcr$type == "Standard curve" , ] ; stand.mp
stand.mp$copies <- log10( stand.mp$copies) ; stand.mp
model.mp <- lm( stand.mp$copies ~ stand.mp$ct_value )

sp.qpcr.mp.ctmean$copies <- 10^(coef(model.mp)[1] + (sp.qpcr.mp.ctmean$mean_value * coef(model.mp)[2]))
sp.qpcr.mp.ctmean$logcopy <- log10( sp.qpcr.mp.ctmean$copies )

p.qpcr.mp <- ggplot(stand.mp , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = sp.qpcr.mp.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(sp.qpcr.mp.ctmean$logcopy), xend =sp.qpcr.mp.ctmean$mean_value[sp.qpcr.mp.ctmean$logcopy == min(sp.qpcr.mp.ctmean$logcopy)], yend = min(sp.qpcr.mp.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(sp.qpcr.mp.ctmean$logcopy), xend =sp.qpcr.mp.ctmean$mean_value[sp.qpcr.mp.ctmean$logcopy == max(sp.qpcr.mp.ctmean$logcopy)], yend = max(sp.qpcr.mp.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(sp.qpcr.mp.ctmean$mean_value), y = -Inf, xend = min(sp.qpcr.mp.ctmean$mean_value), yend = sp.qpcr.mp.ctmean$logcopy[sp.qpcr.mp.ctmean$mean_value == min(sp.qpcr.mp.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(sp.qpcr.mp.ctmean$mean_value), y = -Inf, xend = max(sp.qpcr.mp.ctmean$mean_value), yend = sp.qpcr.mp.ctmean$logcopy[sp.qpcr.mp.ctmean$mean_value == max(sp.qpcr.mp.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(stand.mp$ct_value) , y = min(sp.qpcr.mp.ctmean$logcopy)+0.1 , label = round(min(sp.qpcr.mp.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(stand.mp$ct_value) , y = max(sp.qpcr.mp.ctmean$logcopy)+0.1 , label = round(max(sp.qpcr.mp.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "MP copies (log10)") +
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
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3, color = "red",label.y.npc = 1)
p.qpcr.mp

## Channa argus
sp.qpcr.ca <- sp.qpcr[ sp.qpcr$target == "CA" & !sp.qpcr$type == "Standard curve" , colnames(sp.qpcr) %in% c("type" ,"ct_value")] ; sp.qpcr.ca
sp.qpcr.ca.ctmean <- sp.qpcr.ca %>% group_by(type) %>% summarise(mean_value = mean(ct_value))

stand.ca <- sp.qpcr[ sp.qpcr$target == "CA" & sp.qpcr$type == "Standard curve" , ] ; stand.ca
stand.ca$copies <- log10( stand.ca$copies) ; stand.ca
model.ca <- lm( stand.ca$copies ~ stand.ca$ct_value )

sp.qpcr.ca.ctmean$copies <- 10^(coef(model.ca)[1] + (sp.qpcr.ca.ctmean$mean_value * coef(model.ca)[2]))
sp.qpcr.ca.ctmean$logcopy <- log10( sp.qpcr.ca.ctmean$copies )

p.qpcr.ca <- ggplot(stand.ca , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = sp.qpcr.ca.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(sp.qpcr.ca.ctmean$logcopy), xend =sp.qpcr.ca.ctmean$mean_value[sp.qpcr.ca.ctmean$logcopy == min(sp.qpcr.ca.ctmean$logcopy)], yend = min(sp.qpcr.ca.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(sp.qpcr.ca.ctmean$logcopy), xend =sp.qpcr.ca.ctmean$mean_value[sp.qpcr.ca.ctmean$logcopy == max(sp.qpcr.ca.ctmean$logcopy)], yend = max(sp.qpcr.ca.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(sp.qpcr.ca.ctmean$mean_value), y = -Inf, xend = min(sp.qpcr.ca.ctmean$mean_value), yend = sp.qpcr.ca.ctmean$logcopy[sp.qpcr.ca.ctmean$mean_value == min(sp.qpcr.ca.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(sp.qpcr.ca.ctmean$mean_value), y = -Inf, xend = max(sp.qpcr.ca.ctmean$mean_value), yend = sp.qpcr.ca.ctmean$logcopy[sp.qpcr.ca.ctmean$mean_value == max(sp.qpcr.ca.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(stand.ca$ct_value) , y = min(sp.qpcr.ca.ctmean$logcopy)+0.1 , label = round(min(sp.qpcr.ca.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(stand.ca$ct_value) , y = max(sp.qpcr.ca.ctmean$logcopy)+0.1 , label = round(max(sp.qpcr.ca.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "CA copies (log10)") +
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
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3, color = "red",label.y.npc = 1)
p.qpcr.ca

## Hypophthalmichthys nobilis
sp.qpcr.hn <- sp.qpcr[ sp.qpcr$target == "HN" & !sp.qpcr$type == "Standard curve" , colnames(sp.qpcr) %in% c("type" ,"ct_value")] ; sp.qpcr.hn
sp.qpcr.hn.ctmean <- sp.qpcr.hn %>% group_by(type) %>% summarise(mean_value = mean(ct_value))

stand.hn <- sp.qpcr[ sp.qpcr$target == "HN" & sp.qpcr$type == "Standard curve" , ] ; stand.hn
stand.hn$copies <- log10( stand.hn$copies) ; stand.hn
model.hn <- lm( stand.hn$copies ~ stand.hn$ct_value )

sp.qpcr.hn.ctmean$copies <- 10^(coef(model.hn)[1] + (sp.qpcr.hn.ctmean$mean_value * coef(model.hn)[2]))
sp.qpcr.hn.ctmean$logcopy <- log10( sp.qpcr.hn.ctmean$copies )

p.qpcr.hn <- ggplot(stand.hn , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = sp.qpcr.hn.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(sp.qpcr.hn.ctmean$logcopy), xend =sp.qpcr.hn.ctmean$mean_value[sp.qpcr.hn.ctmean$logcopy == min(sp.qpcr.hn.ctmean$logcopy)], yend = min(sp.qpcr.hn.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(sp.qpcr.hn.ctmean$logcopy), xend =sp.qpcr.hn.ctmean$mean_value[sp.qpcr.hn.ctmean$logcopy == max(sp.qpcr.hn.ctmean$logcopy)], yend = max(sp.qpcr.hn.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(sp.qpcr.hn.ctmean$mean_value), y = -Inf, xend = min(sp.qpcr.hn.ctmean$mean_value), yend = sp.qpcr.hn.ctmean$logcopy[sp.qpcr.hn.ctmean$mean_value == min(sp.qpcr.hn.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(sp.qpcr.hn.ctmean$mean_value), y = -Inf, xend = max(sp.qpcr.hn.ctmean$mean_value), yend = sp.qpcr.hn.ctmean$logcopy[sp.qpcr.hn.ctmean$mean_value == max(sp.qpcr.hn.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(stand.hn$ct_value) , y = min(sp.qpcr.hn.ctmean$logcopy)+0.1 , label = round(min(sp.qpcr.hn.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(stand.hn$ct_value) , y = max(sp.qpcr.hn.ctmean$logcopy)+0.1 , label = round(max(sp.qpcr.hn.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "HN copies (log10)") +
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
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3, color = "red",label.y.npc = 1)
p.qpcr.hn

## plot standard curve and species copies
ps5 <- ggarrange(p.qpcr.mp , p.qpcr.pd , p.qpcr.hn , p.qpcr.ca , nrow = 2, ncol = 2  ) ; ps5
## DNA extraction calibration
lab.dna.rate <- read.csv("results/lab.dna.extract.rate.qcpr.csv") ; lab.dna.rate
lab.dna.rate <- lab.dna.rate[ grepl("200" , lab.dna.rate$type), ] ; lab.dna.rate

sp.qpcr.mp.ctmean$species <- "MP"
sp.qpcr.pd.ctmean$species <- "PD"
sp.qpcr.hn.ctmean$species <- "HN"
sp.qpcr.ca.ctmean$species <- "CA"

sp.qpcr.mp.ctmean$DNA <- sp.qpcr.mp.ctmean$copies/lab.dna.rate$extract_rate ; sp.qpcr.mp.ctmean
sp.qpcr.pd.ctmean$DNA <- sp.qpcr.pd.ctmean$copies/lab.dna.rate$extract_rate ; sp.qpcr.pd.ctmean
sp.qpcr.hn.ctmean$DNA <- sp.qpcr.hn.ctmean$copies/lab.dna.rate$extract_rate ; sp.qpcr.hn.ctmean
sp.qpcr.ca.ctmean$DNA <- sp.qpcr.ca.ctmean$copies/lab.dna.rate$extract_rate ; sp.qpcr.ca.ctmean

sp.qpcr.all <- as.data.frame( rbind( sp.qpcr.mp.ctmean , sp.qpcr.pd.ctmean ,sp.qpcr.hn.ctmean , sp.qpcr.ca.ctmean) )
sp.qpcr.all$group <- substr( sp.qpcr.all$type , 1, 1)
sp.qpcr.all$biogroup <- paste( sp.qpcr.all$group , sp.qpcr.all$species) ; sp.qpcr.all

spcopy <- melt(sp.qpcr.all , id.vars = c("species" ,"group" , "biogroup") , measure.vars = c("copies" , "DNA"))
spbiomas <- read.csv("lab.species.biomass.csv") ; spbiomas
spcopy <- merge( spcopy , spbiomas , by = "biogroup") ; spcopy 
spcopy$species <- factor(spcopy$species , levels = c("MP" , "PD" ,"HN" ,"CA"))
spcopy$variable <- gsub("copies" ,"Rawdata" , spcopy$variable)
spcopy$variable <- gsub("DNA" ,"DNA calibrated" , spcopy$variable)
spcopy$variable <- factor( spcopy$variable , levels = c("Rawdata" , "DNA calibrated"))

p.qpcr <- ggplot(spcopy , aes(x = biomass , y = value , group = species)) +
  geom_point(aes(fill = group), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "DNA copies") +
  facet_grid( variable~species , scales = "free" )+
  #scale_y_log10()+
  theme_light()+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, label.y = 0.8, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3)
p.qpcr

p.qpcr.mp2 <- ggplot(spcopy[spcopy$species == "MP",] , aes(x = biomass , y = value , group = species)) +
  geom_point(aes(fill = group), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "DNA copies") +
  facet_wrap( vars(variable) , scales = "free", ncol = 2)+
  #scale_y_log10()+
  theme_light()+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, label.y = 0.8, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3)
p.qpcr.mp2

p.qpcr.pd2 <- ggplot(spcopy[spcopy$species == "PD",] , aes(x = biomass , y = value , group = species)) +
  geom_point(aes(fill = group), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "DNA copies") +
  facet_wrap( vars(variable) , scales = "free", ncol = 2)+
  #scale_y_log10()+
  theme_light()+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, label.y = 0.8, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3)
p.qpcr.pd2

p.qpcr.hn2 <- ggplot(spcopy[spcopy$species == "HN",] , aes(x = biomass , y = value , group = species)) +
  geom_point(aes(fill = group), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "DNA copies") +
  facet_wrap( vars(variable) , scales = "free", ncol = 2)+
  #scale_y_log10()+
  theme_light()+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, label.y = 0.8, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3)
p.qpcr.hn2

p.qpcr.ca2 <- ggplot(spcopy[spcopy$species == "CA",] , aes(x = biomass , y = value , group = species)) +
  geom_point(aes(fill = group), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "DNA copies") +
  facet_wrap( vars(variable) , scales = "free", ncol = 2)+
  #scale_y_log10()+
  theme_light()+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, label.y = 0.8, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3)
p.qpcr.ca2

############################
############################step7: plot summary
# 创建一个空白的ggplot对象，不包含图层
p.blank <- ggplot() + theme_void()
ps1 <- ggarrange( pfish , pf ,  ncol = 2 , labels = c("a" , "b") , widths = c(1,1.4)) ; ps1
ps2 <- ggarrange(p.ips1 , p.ips2.1 , widths = c(1,1.4) , labels = c("c" , "d")) ; ps2
ps3 <- ggarrange(ps1 , p.blank, ps2 , heights = c(1,0.1,1) , labels = c("a" ,"", "c") , nrow = 3) ; ps3
pk2 <- ggarrange(p.dna , ggarrange( p.blank , pr , ncol = 1 , heights = c(0.1,1)) , nrow = 1, labels = c("e" , "f") , widths = c(1,1.4)) ; pk2
pk3 <- ggarrange(ps3 , pk2 , nrow = 2 , labels = c("a" , "e") , heights = c(1.5,1) ) ; pk3

pdf("figures/Figure1.pdf" , width = 17 , height = 10)
ggarrange(pk3 , p150 , nrow = 1 , labels = c("a" , "g") , widths = c(1,1.2) )
dev.off()

ps5 <- ggarrange(p.blank , p.qpcr.mp , p.qpcr.pd , p.qpcr.hn , p.qpcr.ca , p.blank , nrow = 1 , widths = c(0.1,1,1,1,1,0.1) ) ; ps5
ps6 <- ggarrange(p.blank , p.qpcr.mp , p.qpcr.pd , p.qpcr.hn , p.qpcr.ca , ncol = 1 , heights = c(0.1,1,1,1,1) ) ; ps6
ps7 <- ggarrange( ps5 , p.qpcr , nrow = 2 , heights = c(1,2) , labels = c("a" ,"b")) ; ps7
ps8 <- ggarrange( p.qpcr.mp2 , p.qpcr.pd2 ,p.qpcr.hn2, p.qpcr.ca2, nrow = 4 , labels = c("b" ,"c" , "d" , "e")) ; ps8
ps9 <- ggarrange( ps6 , ps8 , nrow = 1 , widths = c(1,2) , labels = c("a" ,"b")) ; ps9

pdf("figures/Figure2.pdf" , width = 9 , height = 9)
ps9
dev.off()


