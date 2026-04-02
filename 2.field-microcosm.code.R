#### Load necessary R packages
library(xlsx)
library(plyr)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggprism)
library(ggalluvial)  # Additional package for alluvial plots

############################
############################ STEP 1: Initial assessment of NGS reads
############################
## Purpose: Examine the composition and distribution of raw NGS reads across pond samples
## Input: ngs.reads.csv - Raw NGS reads table with samples as rows and species/IPS as columns
## Output: Figures showing read composition and distribution across different ponds

### Set working directory

setwd("field-microcosm/")

############################
### Read NGS reads table
field <- read.csv("ngs.reads.csv" , row.names = 1) ; head( field )
field$type <- gsub("_.*" , "" , rownames( field ) )
field$type <- factor( field$type , levels = c("Pond1" ,"Pond2" , "Pond3" , "Pond4", "Pond5", "Ctrl"))
head( field )
df.field <- melt( field , id.vars = "type") ; head( df.field )
my_comp <- list(c("Pond1", "Pond2"), c("Pond1", "Pond3"), c("Pond2", "Pond3"))
field.clean <- field[ !grepl("Ctrl" , rownames( field ) ), 1:6 ] ; head(field.clean)

###### IPS sequences composition analysis
seqcom <- data.frame( colSums( field.clean ) ) ; seqcom
colnames( seqcom ) <- "reads";seqcom
seqcom$categroy <- c( "IPS1","IPS2.1","IPS2.2","IPS2.3" , "Cyprinus carpio","Other" );seqcom
seqcom$prop <- round(seqcom$reads/sum( seqcom$reads ) , 4)
seqcom$type <- "Total";seqcom

seqcomp <- seqcom[1:4, ]
seqcomp$prop <- round(seqcomp$reads/sum(seqcomp$reads),4)
seqcomp$type <- "IPS";seqcomp

seqdf <- rbind( seqcom , seqcomp) ; seqdf
seqdf$categroy <- factor(seqdf$categroy , levels = unique(seqdf$categroy) )
seqdf$type <- factor( seqdf$type , levels = unique(seqdf$type) )

lb <- data.frame(type = c("Total" , "IPS") , count = c( sum(seqcom$reads) , sum( seqcomp$reads)))

## Figure: Composition of sequencing reads (Total and IPS only)
pf <- ggplot( data = seqdf , aes(x = type , y = prop , fill = categroy) ) +
  geom_bar( stat = "identity" , position = "fill") +
  geom_flow( aes(alluvium = categroy ), alpha = 0.5) + 
  scale_fill_manual( values = c("IPS1" = "#db6968",
                                "IPS2.1" ="#4d97cd",
                                "IPS2.2" = "#f8984e",
                                "IPS2.3" ="#459943",
                                "Cyprinus carpio" = brewer.pal(9, "PuBu")[7],
                                "Other" = brewer.pal(9, "PuBu")[5] )
  )+
  geom_text( data = lb , aes(x = type , y=1, label = paste0( count ) ) , inherit.aes = FALSE , vjust = -0.2 , size =3)+
  theme_light()+
  labs( x = "Reads" , y = "Proportion")+
  guides(fill = guide_legend(title = ""))+
  theme_prism() +
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.8, "lines"))
pf

### Fish sequence composition across ponds
library(dplyr)
dfseq <- data.frame( biomass = substr(rownames( field.clean ) , 1, 5), reads = rowSums( field.clean[ , !grepl("IPS",colnames( field.clean) ) ] ) )
summary <- dfseq %>%
  group_by(biomass) %>%
  summarise(
    mean = mean(reads),
    se = sd(reads) / sqrt(n())
  )

## Figure: Fish reads distribution across ponds
pfish <- ggplot(summary, aes(x = biomass, y = mean, fill = biomass)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_text(aes(label = round(mean,0) ), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Fish reads") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#98d09d")) + 
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10 , angle = 45, hjust = 0.5),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")
pfish

###### Figure: IPS1 reads distribution across different ponds
ips1reads <- df.field[ (!df.field$type == "Ctrl") & df.field$variable %in% "IPS1",] ; ips1reads

p.ips1 <- ggplot(ips1reads , aes(x = type , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = type), pch = 21,size = 2, position = position_jitter(0.2))+
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#98d09d")) + 
  labs(x = NULL,y = "IPS1 reads") +
  stat_compare_means(comparisons=list(c("Pond1", "Pond2") ,c("Pond2", "Pond4"),c("Pond1", "Pond5")), method = "t.test" , label="p.signif")+
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4) ,
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10 , angle = 45 , hjust = 0.5),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")
p.ips1

###### Figure: IPS2 reads distribution across different IPS2 variants
ips2reads <- df.field[ (!df.field$type == "Ctrl") & grepl("IPS2" , df.field$variable) , ] ; ips2reads
p.ips2 <- ggplot( ips2reads , aes(x = variable , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = type), pch = 21,size = 2, position = position_jitter(0.2))+
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#98d09d")) + 
  labs(x = NULL, y = "Reads") +
  stat_compare_means(comparisons=list(c("IPS2.1","IPS2.2"), c("IPS2.2","IPS2.3"), c("IPS2.1","IPS2.3")), method="t.test", label="p.signif")+
  theme_prism() +
  guides(fill = guide_legend(title = ""))+
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10 , angle = 45 , hjust = 0.5),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p.ips2

## Figure: IPS2 reads distribution by pond and IPS2 variant
p.ips2.2 <- ggplot(df.field[(!df.field$type == "Ctrl") & grepl("IPS2",df.field$variable) ,] , aes(x = type , y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_point(aes(fill = type), pch = 21,size = 2, position = position_jitter(0.2))+
  facet_wrap( .~variable , scales="free") +
  scale_x_discrete(guide = "prism_bracket") +
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#98d09d")) + 
  labs(x = NULL , y = "Reads") +
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10 , angle = 45 , hjust = 0.5),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")
p.ips2.2 

###########################
############################ STEP 2: PCR amplification efficiency calibration
############################
## Purpose: Correct for PCR amplification bias using IPS2 spike-in standards
## Method: Linear regression of log-transformed IPS2 reads vs. known concentrations
## Equation: y = (x - b)/a, where y is calibrated reads, x is raw reads,
##           a is slope, b is intercept from IPS2 linear model
## Input: ngs.reads.csv - Raw NGS reads
## Output: PCR-calibrated reads table

## Read NGS reads table
field <- read.csv("ngs.reads.csv" , row.names = 1) ; head( field )
field.clean <- field[ !grepl("Ctrl" , rownames( field ) ), ] ; head(field.clean)

## Extract IPS2 reads for PCR efficiency calibration
field.ips2 <- field[ !grepl("Ctrl" , rownames( field )) , grepl("IPS2" , colnames( field ) ) ] ; head( field.ips2 )
field.ips2 <- log(field.ips2)  # Log-transform for linear modeling
field.ips2$biomass <- substr( rownames( field.ips2 ) , 1, 5) ; head( field.ips2 )

## Sample-specific PCR efficiency correction
## For each sample, build linear model using the three IPS2 concentrations
rho <- c() ; pvalue <- c() ; formula <- c()
ips1adj <- c() ; mpadj <- c() ; oadj <- c()
for( i in 1: nrow( field.ips2 ) ){
  x <- c(5.5, 32.75, 66.83)  # IPS2 known concentrations (copies/μL) for field samples
  y <- c(as.numeric( field.ips2[i,1:3] ) )  # Log-transformed IPS2 reads
  model <- lm( y ~ x )
  cts <- summary( model )
  rho[i] <- cts$r.squared
  pvalue[i] <- round(cts$coefficients[2,4] , 4)
  formula[i] <- paste( "y =",round(coef(model)[1] , 2) , "+", paste0(round(coef(model)[2] , 2)  , "x") )
  
  # Apply correction only for samples with R² ≥ 0.9
  if( cts$r.squared >= 0.9 ){
    ips1adj[i] <- round( ( log(field$IPS1[i])-coef(model)[1][[1]])/coef(model)[2][[1]] ,0)
    mpadj[i] <- round( ( log(field$Cyprinus.carpio[i])-coef(model)[1][[1]])/coef(model)[2][[1]],0)
    oadj[i] <- round( ( log(field$other[i])-coef(model)[1][[1]])/coef(model)[2][[1]],0)
  } else {
    # If model fit is poor, retain log-transformed raw reads
    ips1adj[i] <- log(field$IPS1[i])
    mpadj[i] <- log(field$Cyprinus.carpio[i])
    oadj[i] <- log(field$other[i])
  }
}

df1 <- data.frame(R = rho , pvalue = pvalue , formula = formula , IPS1 = ips1adj , Cyprinus.carpio = mpadj, other = oadj )
dfadj <- cbind( field.ips2 , df1 ) ; dfadj

dir.create("results", showWarnings = F)
write.csv( dfadj , file = "results/ngs.reads.PCRadj.sample.csv" )

#### Figure: IPS2 linear fitting for each sample
library( reshape2 )
field.ips2$sample <- rownames( field.ips2 ) ; head( field.ips2 )
df.ips2 <- data.frame(conc = c(5.5, 32.75, 66.83) , variable = c("IPS2.1" , "IPS2.2" , "IPS2.3"))
df.ips2 <- merge( melt(field.ips2,id.vars = c("sample" , "biomass") )  , df.ips2 , by = "variable") ; head( df.ips2 )
df.ips2$biomass <- factor( df.ips2$biomass , levels = paste0("Pond" , 1:5) )

ips2.100 <- ggplot(df.ips2 , aes(x = conc , y = value )) +
  geom_point(aes(fill = conc), pch = 21, size = 3)+
  labs(x = "IPS2 concentration (copies/μL)", y = "IPS2 NGS Reads (log-transformed)") +
  theme_prism() +
  theme(strip.text = element_text(size = 12),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 10),
        axis.text.x = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  facet_wrap( .~sample , scales="free" , ncol = 3) +
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
ips2.100

pdf("figures/IPS2.linear.fitting.sample.pdf" , height = 10 , width = 8)
ips2.100
dev.off()

############################
############################ STEP 3: DNA extraction efficiency correction
############################
## Purpose: Correct for DNA extraction efficiency variation using IPS1 qPCR data
## Method: Calculate extraction efficiency as Dext/D0, where Dext is IPS1 copies by qPCR,
##         and D0 = 30 is the theoretical spiked IPS1 copy number for field samples
## Equation: Rq = R0 / (Dext/D0) or Rq = R0 × (D0/Dext)
## Input: ips1.qpcr.csv - qPCR data for IPS1
## Output: DNA extraction efficiency table and DNA-corrected reads

## Read qPCR data for IPS1
field.qpcr <- read.csv("ips1.qpcr.csv") ; head( field.qpcr )
field.qpcr100 <- field.qpcr[!field.qpcr$type == "Standard curve", colnames(field.qpcr ) %in% c("ct_value" , "type")]

## Calculate mean CT values for each sample
field.qpcr100.ctmean <- field.qpcr100 %>% group_by(type) %>% summarise(mean_value = mean(ct_value))
field.qpcr100.ctmean <- as.data.frame( field.qpcr100.ctmean )

## Linear fitting for qPCR standard curve
stand1 <- field.qpcr[field.qpcr$type == "Standard curve" , ]
stand1$copies <- log10( stand1$copies)

## Calculate IPS1 copy number
ips1.conc <- 30  # Theoretical spiked IPS1 copy number for field samples
model1 <- lm( stand1$copies ~ stand1$ct_value )
field.qpcr100.ctmean$copies <- 10^(coef(model1)[1] + (field.qpcr100.ctmean$mean_value * coef(model1)[2]))
field.qpcr100.ctmean$logcopy <- log10( field.qpcr100.ctmean$copies )
field.qpcr100.ctmean$extract_rate <- field.qpcr100.ctmean$copies/ips1.conc  # Extraction efficiency = Dext/D0

## Figure: Linear fitting model for qPCR standard curve
p.fitmodel100 <- ggplot(stand1 , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = field.qpcr100.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(field.qpcr100.ctmean$logcopy), xend =field.qpcr100.ctmean$mean_value[field.qpcr100.ctmean$logcopy == min(field.qpcr100.ctmean$logcopy)], yend = min(field.qpcr100.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(field.qpcr100.ctmean$logcopy), xend =field.qpcr100.ctmean$mean_value[field.qpcr100.ctmean$logcopy == max(field.qpcr100.ctmean$logcopy)], yend = max(field.qpcr100.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(field.qpcr100.ctmean$mean_value), y = -Inf, xend = min(field.qpcr100.ctmean$mean_value), yend = field.qpcr100.ctmean$logcopy[field.qpcr100.ctmean$mean_value == min(field.qpcr100.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(field.qpcr100.ctmean$mean_value), y = -Inf, xend = max(field.qpcr100.ctmean$mean_value), yend = field.qpcr100.ctmean$logcopy[field.qpcr100.ctmean$mean_value == max(field.qpcr100.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(field.qpcr$ct_value) , y = min(field.qpcr100.ctmean$logcopy)+0.2 , label = round(min(field.qpcr100.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(field.qpcr$ct_value) , y = max(field.qpcr100.ctmean$logcopy)+0.2 , label = round(max(field.qpcr100.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "IPS1 copies (log10)") +
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
  stat_cor(size = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson") +
  stat_regline_equation(size = 4, color = "red",label.y.npc = 1)
p.fitmodel100

## Save DNA extraction efficiency data
field.dna.extract.rate <- field.qpcr100.ctmean ; field.dna.extract.rate
field.dna.extract.rate <- field.dna.extract.rate[order(field.dna.extract.rate$type) , ] ; field.dna.extract.rate
write.csv( field.dna.extract.rate , file = "results/field.dna.extract.rate.qcpr.csv" , row.names = F)

## Correct DNA extraction efficiency by qPCR (applied to PCR-corrected data)
dnarate <- read.csv("results/field.dna.extract.rate.qcpr.csv") ; head( dnarate )
reads.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1) ; head( reads.adj.sample )
reads.adj.sample.dna <- reads.adj.sample[ , colnames(reads.adj.sample) %in% c("IPS1" , "Cyprinus.carpio" , "other")]/dnarate$extract_rate
write.csv(reads.adj.sample.dna , file = "results/ngs.reads.PCRadj.DNAadj.sample.qPCR.csv")

############################
############################ STEP 4: DNA extraction efficiency correction only (without PCR correction)
############################
## Purpose: Apply DNA extraction efficiency correction directly to raw NGS reads
## Methods: 1) Using qPCR-derived extraction efficiency
##          2) Using IPS1 NGS reads as proxy for extraction efficiency
## Input: ngs.reads.csv - Raw NGS reads, field.dna.extract.rate.qcpr.csv - qPCR data
## Output: DNA-corrected reads tables (without PCR correction)

## Correct DNA efficiency by qPCR (applied to raw reads)
## Read raw NGS reads table
field <- read.csv("ngs.reads.csv" , row.names = 1) ; head( field )
field.clean <- field[ !grepl("Ctrl" , rownames( field ) ), ]
field.clean <- field.clean[order(rownames(field.clean)) , ] ; field.clean

# Read qPCR DNA extraction data
dnarate <- read.csv("results/field.dna.extract.rate.qcpr.csv") ; dnarate

# Apply qPCR-based DNA correction to raw reads
field.clean.adj.dna.qpcr <- field.clean[ , colnames(field.clean) %in% c("Cyprinus.carpio" , "other" )]/dnarate$extract_rate
write.csv(field.clean.adj.dna.qpcr , file = "results/ngs.reads.DNAadj.sample.qpcr.csv")

## Correct DNA efficiency by IPS1 NGS reads (using IPS1 as proxy for extraction efficiency)
## Read raw NGS reads table
field <- read.csv("ngs.reads.csv" , row.names = 1) ; head( field )
field.clean <- field[ !grepl("Ctrl" , rownames( field ) ), ]
field.clean <- field.clean[order(rownames(field.clean)) , ] ; field.clean

reads.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1) ; head( reads.adj.sample )

# Calculate DNA rate as relative IPS1 reads
reads.adj.sample$IPS1/max(reads.adj.sample$IPS1)

# Apply IPS1-based DNA correction to raw reads
field.clean.adj.dna.ngs <- field.clean[ , colnames(field.clean) %in% c("Cyprinus.carpio" , "other" )]/(reads.adj.sample$IPS1/max(reads.adj.sample$IPS1))
write.csv(field.clean.adj.dna.ngs , file = "results/ngs.reads.DNAadj.sample.ngs.csv")

############################
############################ STEP 5: Visualization of correction effects
############################
## Purpose: Compare different correction methods and evaluate their effectiveness
## Methods: Compare raw data, DNA-only correction, PCR-only correction, combined correction,
##          and a novel IPS1 read-based calibration
## Input: All corrected data files from previous steps
## Output: Figures showing linear relationships between reads and biomass for each method

spbiomas <- read.csv("field.species.biomass.csv") ; spbiomas

## Read raw NGS data
field <- read.csv("ngs.reads.csv" , row.names = 1) ; head( field )
rawdata <- field[ !grepl("Ctrl" , rownames( field ) ), colnames(field) %in% c("Cyprinus.carpio" , "other" )]
rawdata$biogroup <- substr(rownames(rawdata) , 1,5) ; rawdata

rawdata <- merge(rawdata , spbiomas , by = "biogroup");rawdata
rawdata$type <- "Rawdata" ; head( rawdata )

## DNA correction only (qPCR-based)
dna.adj.qpcr <- read.csv("results/ngs.reads.DNAadj.sample.qpcr.csv" , row.names = 1) ; head( dna.adj.qpcr )
dna.adj.qpcr$biogroup <- substr(rownames(dna.adj.qpcr) , 1,5) ; dna.adj.qpcr
dna.adj.qpcr <- merge(dna.adj.qpcr , spbiomas , by = "biogroup") ; dna.adj.qpcr
dna.adj.qpcr$type <- "DNA calibrated" ; head( dna.adj.qpcr )

## PCR correction only
pcr.adj.sample <- read.csv("results/ngs.reads.PCRadj.sample.csv" , row.names = 1)[, 9:10] ; head( pcr.adj.sample )
pcr.adj.sample$biogroup <- substr(rownames(pcr.adj.sample) , 1,5) ; pcr.adj.sample
pcr.adj.sample <- merge(pcr.adj.sample , spbiomas , by = "biogroup") ; pcr.adj.sample
pcr.adj.sample$type <- "PCR calibrated" ; head( pcr.adj.sample )

## Both PCR and DNA correction (qPCR-based)
both.adj.sample.qpcr <- read.csv("results/ngs.reads.PCRadj.DNAadj.sample.qPCR.csv" , row.names = 1)[,2:3] ; head( both.adj.sample.qpcr )
both.adj.sample.qpcr$biogroup <- substr(rownames(both.adj.sample.qpcr) , 1,5) ; both.adj.sample.qpcr
both.adj.sample.qpcr <- merge(both.adj.sample.qpcr , spbiomas , by = "biogroup")
both.adj.sample.qpcr$type <- "PCR+DNA calibrated" ; head(both.adj.sample.qpcr)

# ====================== NEW METHOD: IPS1 read-based calibration ======================
## IPS1 read-based calibrated: Corrected Cyprinus.carpio reads = Raw Cyprinus.carpio reads × 
##                             (IPS1 qPCR spike-in copy number / IPS1 NGS reads)
## This method uses only NGS data and the known IPS1 spike-in concentration

# Read raw data
field_raw <- read.csv("ngs.reads.csv" , row.names = 1)
ips1_qpcr_add <- 30  # IPS1 qPCR spike-in copy number (constant, same as ips1.conc in Step 3)

# Extract raw reads for Cyprinus.carpio and IPS1 from non-Ctrl samples
ips1_calib <- field_raw[!grepl("Ctrl", rownames(field_raw)), c("Cyprinus.carpio", "IPS1")]

# Apply IPS1-based correction: Corrected = Raw × (IPS1_spike / IPS1_NGS)
ips1_calib$Cyprinus.carpio <- ips1_calib$Cyprinus.carpio * ips1_qpcr_add / ips1_calib$IPS1

# Keep other reads unchanged and remove IPS1 column
ips1_calib$other <- field_raw[rownames(ips1_calib), "other"]
ips1_calib$IPS1 <- NULL

# Add biomass group information
ips1_calib$biogroup <- substr(rownames(ips1_calib), 1, 5)
ips1_calib <- merge(ips1_calib, spbiomas, by = "biogroup")

# Label the correction method
ips1_calib$type <- "IPS1 read-based calibrated"
# ==============================================================================

#### Combine all datasets (original 4 methods + new IPS1-based method)
dfall <- rbind(rawdata , dna.adj.qpcr , pcr.adj.sample , both.adj.sample.qpcr, ips1_calib) ; head( dfall )
dfall$type <- factor(dfall$type , levels = c("Rawdata", "DNA calibrated", "PCR calibrated", "PCR+DNA calibrated", "IPS1 read-based calibrated") )

## Figure: Comparison of all correction methods for Cyprinus.carpio reads
p200 <- ggplot(dfall, aes(x = biomass , y = Cyprinus.carpio )) +
  geom_point(aes(fill = biogroup), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "NGS reads") +
  facet_wrap( vars(type) , scales = "free" , nrow = 2)+  # 2 rows, 5 plots automatically arranged
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#74A9CF")) + 
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), 
           color = "black", method = "pearson", label.y.npc = "top") +
  stat_regline_equation(size = 3,label.y.npc = 1)
p200

############################
############################ STEP 6: Species qPCR copies calibration
############################
## Purpose: Quantify Cyprinus.carpio DNA copies using species-specific qPCR
## Method: Convert qPCR CT values to copy numbers using standard curve
## Input: species.qpcr.csv - qPCR data for Cyprinus.carpio
## Output: Species copy numbers and DNA extraction efficiency correction

## Read species qPCR data
sp.qpcr <- read.csv("species.qpcr.csv") ; head( sp.qpcr )

## Cyprinus carpio qPCR analysis
sp.qpcr.pd <- sp.qpcr[ sp.qpcr$target == "Cyprinus carpio" & !sp.qpcr$type == "Standard curve" , colnames(sp.qpcr) %in% c("type" ,"ct_value")] ; sp.qpcr.pd
sp.qpcr.pd.ctmean <- sp.qpcr.pd %>% group_by(type) %>% summarise(mean_value = mean(ct_value))

stand.pd <- sp.qpcr[ sp.qpcr$target == "Cyprinus carpio" & sp.qpcr$type == "Standard curve" , ] ; stand.pd
stand.pd$copies <- log10( stand.pd$copies) ; stand.pd
model.pd <- lm( stand.pd$copies ~ stand.pd$ct_value )

sp.qpcr.pd.ctmean$copies <- 10^(coef(model.pd)[1] + (sp.qpcr.pd.ctmean$mean_value * coef(model.pd)[2]))
sp.qpcr.pd.ctmean$logcopy <- log10( sp.qpcr.pd.ctmean$copies )

## Figure: Cyprinus.carpio qPCR standard curve and sample copy numbers
p.qpcr.pd <- ggplot(stand.pd , aes(x = ct_value , y = copies )) +
  geom_point(aes(fill = copies), pch = 21, size = 3)+
  geom_point(data = sp.qpcr.pd.ctmean , aes(x= mean_value , y = logcopy ), pch = 21, size = 4 , bg = "red3")+
  geom_segment(aes(x = -Inf, y = min(sp.qpcr.pd.ctmean$logcopy), xend =sp.qpcr.pd.ctmean$mean_value[sp.qpcr.pd.ctmean$logcopy == min(sp.qpcr.pd.ctmean$logcopy)], yend = min(sp.qpcr.pd.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = -Inf, y = max(sp.qpcr.pd.ctmean$logcopy), xend =sp.qpcr.pd.ctmean$mean_value[sp.qpcr.pd.ctmean$logcopy == max(sp.qpcr.pd.ctmean$logcopy)], yend = max(sp.qpcr.pd.ctmean$logcopy) ), color = "#e77381", linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = min(sp.qpcr.pd.ctmean$mean_value), y = -Inf, xend = min(sp.qpcr.pd.ctmean$mean_value), yend = sp.qpcr.pd.ctmean$logcopy[sp.qpcr.pd.ctmean$mean_value == min(sp.qpcr.pd.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) + 
  geom_segment(aes(x = max(sp.qpcr.pd.ctmean$mean_value), y = -Inf, xend = max(sp.qpcr.pd.ctmean$mean_value), yend = sp.qpcr.pd.ctmean$logcopy[sp.qpcr.pd.ctmean$mean_value == max(sp.qpcr.pd.ctmean$mean_value)]), color = brewer.pal(8, "Set1")[2], linetype = "dashed" , lwd = 0.4 ) +
  geom_text( x = min(stand.pd$ct_value) , y = min(sp.qpcr.pd.ctmean$logcopy)+0.1 , label = round(min(sp.qpcr.pd.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  geom_text( x = min(stand.pd$ct_value) , y = max(sp.qpcr.pd.ctmean$logcopy)+0.1 , label = round(max(sp.qpcr.pd.ctmean$copies),1) , color = "#e77381" , hjust = 0 , size = 3)+
  labs(x = "qPCR CT value", y = "DNA copies (log10)") +
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
  stat_cor(size = 3, label.y = 5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1),
           color = "black", method = "pearson") +
  stat_regline_equation(size = 3, color = "red",label.y.npc = 1)
p.qpcr.pd

## Apply DNA extraction efficiency correction to qPCR data
lab.dna.rate <- read.csv("results/field.dna.extract.rate.qcpr.csv") ; lab.dna.rate
sp.qpcr.pd.ctmean$DNA <- sp.qpcr.pd.ctmean$copies/lab.dna.rate$extract_rate ; sp.qpcr.pd.ctmean

sp.qpcr.all <- as.data.frame( sp.qpcr.pd.ctmean )
sp.qpcr.all$biogroup <- substr( sp.qpcr.all$type , 1, 5) ; sp.qpcr.all

spcopy <- melt(sp.qpcr.all , id.vars = "biogroup" , measure.vars = c("copies" , "DNA")) ; head(spcopy)
spbiomas <- read.csv("field.species.biomass.csv") ; spbiomas
spcopy <- merge( spcopy , spbiomas , by = "biogroup") ; spcopy 

spcopy$variable <- gsub("copies" ,"Rawdata" , spcopy$variable)
spcopy$variable <- gsub("DNA" ,"DNA calibrated" , spcopy$variable)
spcopy$variable <- factor( spcopy$variable , levels = c("Rawdata" , "DNA calibrated"))

## Figure: Comparison of raw and DNA-corrected qPCR data
p.qpcr <- ggplot(spcopy , aes(x = biomass , y = value )) +
  geom_point(aes(fill = biogroup), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "DNA copies") +
  facet_wrap( vars(variable) , scales = "free" )+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#74A9CF")) +
  scale_y_log10()+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3,  aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "black", method = "pearson", label.x.npc = "left") +
  stat_regline_equation(size = 3,label.y.npc = 1)
p.qpcr

############################
############################ STEP 7: Summary plots for publication
############################
## Purpose: Generate publication-ready figures combining all analyses
## Note: Figures are saved as individual PDF files for flexible use in manuscripts

library(ggpubr)
p.blank <- ggplot() + theme_void()

# Figures a-d: Initial assessment plots
ps1 <- ggarrange( pfish , pf ,  ncol = 2 , labels = c("a" , "b") , widths = c(1 , 1.5)) ; ps1
ps2 <- ggarrange(p.ips1 , p.ips2 , widths = c(1,1.5) , labels = c("c" , "d")) ; ps2
ps3 <- ggarrange( ps1 , p.blank , ps2 , heights = c(1 , 0.05 , 1) , nrow = 3 ) ; ps3

# Figures e-f: DNA extraction efficiency and qPCR standard curve
pk2 <- ggarrange(p.fitmodel100 , p.qpcr.pd , nrow = 1, labels = c("e" , "f") , widths = c(1,1)) ; pk2
pk3 <- ggarrange(ps3 , pk2 , nrow = 2 , heights = c(1.8,1) ) ; pk3

# ==================== Create new figure g (first 4 methods only) ====================

# Create dataset with only the first 4 correction methods (excluding IPS1 read-based)
dfall_g <- rbind(rawdata , dna.adj.qpcr , pcr.adj.sample , both.adj.sample.qpcr) ; head(dfall_g)
dfall_g$type <- factor(dfall_g$type , levels = c("Rawdata", "DNA calibrated", "PCR calibrated", "PCR+DNA calibrated"))
dfall_g$biogroup <- factor(dfall_g$biogroup , levels = paste0("Pond",1:5))

# Create new figure g (only the first 4 methods)
p200_g <- ggplot(dfall_g, aes(x = biomass , y = Cyprinus.carpio )) +
  geom_point(aes(fill = biogroup), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "NGS reads") +
  facet_wrap( vars(type) , scales = "free" , nrow = 2)+  # 2 rows, 4 subplots
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#74A9CF")) + 
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), 
           color = "black", method = "pearson", label.y.npc = "top") +
  stat_regline_equation(size = 3,label.y.npc = 1)

# ==================== Save the 5th method separately (IPS1 read-based calibrated) ====================

# Create dataset with only the 5th method
dfall_ips <- ips1_calib
dfall_ips$type <- factor(dfall_ips$type , levels = c("IPS1 read-based calibrated"))
dfall_ips$biogroup <- factor(dfall_ips$biogroup , levels = paste0("Pond",1:5))

# Create figure for the 5th method
p200_ips <- ggplot(dfall_ips, aes(x = biomass , y = Cyprinus.carpio )) +
  geom_point(aes(fill = biogroup), pch = 21, size = 3)+
  labs(x = "Biomass (g)", y = "NGS reads") +
  facet_wrap( vars(type) , scales = "free" , nrow = 1)+  # Single row
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  scale_fill_manual(values = c("#db6968","#4d97cd","#f8984e","#459943","#74A9CF")) + 
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), 
           color = "black", method = "pearson", label.y.npc = "top") +
  stat_regline_equation(size = 3,label.y.npc = 1)

# Save the 5th method figure separately
pdf("figures/Figure_IPS1_based_calibrated.pdf" , width = 8, height = 4)
print(p200_ips)
dev.off()

# ==================== Create new Figure 3====================

# Adjust internal proportions for g and h for better clarity
pk4_new <- ggarrange(p.blank , p200_g , p.blank , p.qpcr , nrow = 4 , labels = c("","" ,"", "h") , heights = c(0.05,2.2,0.1,1.2) ) ; pk4_new

# Key: Set widths=c(2,3) to achieve left panel 2/5, right panel (g+h) 3/5 proportion
pdf("figures/Figure3.pdf" , width = 16 , height = 10)  # Width adjusted to 16 for better display of g+h
ggarrange(pk3 , pk4_new , nrow = 1 , labels = c("a" , "g") , widths = c(2,3) )
dev.off()

