library(plyr)
library(ggplot2)

## Read in the latest titer data
dat <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/HI_titers_paired.csv")
dat$index <- seq(1,nrow(dat),by=1)

dat <- unique(dat[,c("Visit","Virus","Participant_ID")])
####################################
## Do some data cleaning here
dat[dat$Visit == "v1","Visit"] <- "V1"
dat <- dat[dat$Visit != "20",]
dat[dat$HI_Titer == 1,"HI_Titer"] <- 10
######################################

## Convert to log scale
dat[dat$HI_Titer == 0,"HI_Titer"] <- 5
dat$logTiter <- log2(dat$HI_Titer/5)
##

## Get counts for each visit
HI_counts <- ddply(dat[,c("Visit","Virus","logTiter")], c("Visit","Virus"), count)
ddply(HI_counts, "Virus", function(x) reshape2::dcast(x,Visit~logTiter))

## Plot histograms by visit and virus
p1 <- ggplot(dat) + 
  geom_histogram(aes(x=logTiter,y=..density..),binwidth=1) + 
  facet_grid(Virus~Visit) + 
  theme_bw() + 
  scale_x_continuous(limits = c(-0.5,8),breaks=seq(0,8,by=1),labels=seq(0,8,by=1)) +
  scale_y_continuous(limits=c(0,0.75))


######################################
## Useful stats
######################################
## Number by visit
visit_ns <- ddply(dat[,c("Visit","Participant_ID")], "Visit",function(x) length(unique(x$Participant_ID)))
colnames(visit_ns) <- c("Visit","Unique_Participants")

## How many visits each individual participated in
visit_key <- c("V1","V2","V3","V4")
visits <- ddply(dat[,c("Visit","Participant_ID")], "Participant_ID", function(x){
  as.numeric(visit_key %in% unique(x$Visit))})
visits$total <- rowSums(visits[,2:5])
visit_summary <- count(visits$total)

