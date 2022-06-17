#========================CHICTABA Analysis Script==========================
# This script plots figures for the CHICTABA samples' nitrate isotope analysis
# Written by Pete D Akers, May 2019. Updated and streamlined July 2021.
# Prepared for GitHub June 2022.

library(ggplot2)
library(Rmisc)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(lubridate)
library(tidyverse)

#### Additional Functions ####
#Geom Stepribbon https://github.com/adibender/pammtools/blob/master/R/ggplot-extensions.R
geom_stepribbon <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
  na.rm       = FALSE,
  show.legend = NA,
  inherit.aes = TRUE, ...) {
  
  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepribbon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, ... )
  )
  
}

#' @rdname geom_stepribbon
#' @format NULL
#' @usage NULL
#' @export

GeomStepribbon <- ggproto(
  "GeomStepribbon", GeomRibbon,
  
  extra_params = c("na.rm"),
  
  draw_group = function(data, panel_scales, coord, na.rm = FALSE) {
    
    if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
    data   <- rbind(data, data)
    data   <- data[order(data$x), ]
    data$x <- c(data$x[2:nrow(data)], NA)
    data   <- data[complete.cases(data["x"]), ]
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)
    
  }
  
)

lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

darken <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col - (col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#### END Additional Functions ####

#Reading in CHICTABA nitrate data
chic <- read.csv("chictaba_NO3data.csv", header=TRUE, sep=",")
chic$id <- as.character(chic$id)
chic$site <- as.character(chic$site)
chic$group <- as.character(chic$group)
chic$group <- factor(chic$group, levels=c("skl", "dl", "p1", "p2", "p3", "p4", "p5")) #Making so the order is always correct

#Splitting database by sampling method into pits and skin layer (skl) + 1 m depth layer (dl)
chic.pit <- chic[chic$group != "skl" & chic$group != "dl",]
chic.pit$site[chic.pit$group == "p1"] <- "CHIC-01.p1" #Renaming sites to better distinguish pit data
chic.pit$site[chic.pit$group == "p2"] <- "CHIC-05.p2"
chic.pit$site[chic.pit$group == "p3"] <- "CHIC-11.p3"
chic.pit$site[chic.pit$group == "p4"] <- "ABN.p4"
chic.pit$site[chic.pit$group == "p5"] <- "ABN.p5"
pit.site.index <- c("CHIC-01.p1", "CHIC-05.p2", "CHIC-11.p3", "ABN.p4", "ABN.p5")
chic.pit$site <- factor(chic.pit$site, levels=pit.site.index)
chic.bypit <- split(chic.pit, f=chic.pit$site) #Making a list where each pit is separated out as a dataframe

chic.skd <- chic[chic$group == "skl" | chic$group == "dl",] #skd = skin layer + depth layer
chic.skd$group <- factor(chic.skd$group, levels=c("dl", "skl")) #Making it so the order is always correctly displayed
chic.byskd <- split(chic.skd, f=chic.skd$site) #Making a list where each unique site is separated out a dataframe
chic.byskd <- chic.byskd[lapply(chic.byskd, nrow) > 1] #Keeping only sites that have both a skin layer and depth layer sample

#Making an alternate dataframe that excludes the SKL pseudo data (ie, the duplicate and extrapolated data)
  #The skin layer sample for P2-P5 were duplicated for simpler analysis from the skin layer sample taken at the same site.
  #There was no skin layer sample taken at CHIC-01 site for P1, so it was extrapolated based on the relationship between
    #other skin layers and local SMB.
  #1 m depth layer values were also created from the pit samples by aggregrating one complete yearly cycle. These
    #pseudo depth layers were not included in comparisons between pits vs skl vs dl, but were otherwise used
    #as actual 1 m depth layer samples.
  #Statistical analyses proceed on data with pseudo values removed.
chic.no.pseudo <- chic[!(chic$pseudo == 1),] #For analyses with all pseudo removed
chic.no.pseudo.dl <- chic[!(chic$pseudo == 1 & chic$group != "dl"),] #For analyses when pseudo dl samples are included

#Setting plotting and analysis variables
vbl.clm <- c("NO3", "d18O", "D17O", "d15N") #columns with plotting data
vbl.clm.index <- NA
for (i in 1:length(vbl.clm)) {
  vbl.clm.index[i] <- which(colnames(chic) == vbl.clm[i])
}

#Creating color-iso variable table
color.index <- c("firebrick", "darkorange", "turquoise3", "mediumorchid4")
iso.color.index <- data.frame(colnames(chic[vbl.clm.index]), color.index[seq(1,length(vbl.clm.index))])
colnames(iso.color.index) <- c("iso", "color")

#Calculating means of different group types (excluding pseudo values)
chic.bygroup.means <- chic.no.pseudo %>%
  group_by(group) %>%
  summarize(NO3.mean = mean(NO3, na.rm=T),
            NO3.ci = CI(na.omit(NO3))['upper']-mean(NO3, na.rm=T),
            d18O.mean = mean(d18O, na.rm=T),
            d18O.ci = CI(na.omit(d18O))['upper']-mean(d18O, na.rm=T),
            D17O.mean = mean(D17O, na.rm=T),
            D17O.ci = CI(na.omit(D17O))['upper']-mean(D17O, na.rm=T),
            d15N.mean = mean(d15N, na.rm=T),
            d15N.ci = CI(na.omit(d15N))['upper']-mean(d15N, na.rm=T))

chic.allpits.means <- chic.no.pseudo %>% #Means and CI for all five pits combined
  filter(!group %in% c("dl", "skl")) %>%
  summarize(NO3.mean = mean(NO3, na.rm=T),
            NO3.ci = CI(na.omit(NO3))['upper']-mean(NO3, na.rm=T),
            d18O.mean = mean(d18O, na.rm=T),
            d18O.ci = CI(na.omit(d18O))['upper']-mean(d18O, na.rm=T),
            D17O.mean = mean(D17O, na.rm=T),
            D17O.ci = CI(na.omit(D17O))['upper']-mean(D17O, na.rm=T),
            d15N.mean = mean(d15N, na.rm=T),
            d15N.ci = CI(na.omit(d15N))['upper']-mean(d15N, na.rm=T))

chic.allsamples.means <- chic.no.pseudo %>% #Means and CI for all samples combined
  summarize(NO3.mean = mean(NO3, na.rm=T),
            NO3.ci = CI(na.omit(NO3))['upper']-mean(NO3, na.rm=T),
            d18O.mean = mean(d18O, na.rm=T),
            d18O.ci = CI(na.omit(d18O))['upper']-mean(d18O, na.rm=T),
            D17O.mean = mean(D17O, na.rm=T),
            D17O.ci = CI(na.omit(D17O))['upper']-mean(D17O, na.rm=T),
            d15N.mean = mean(d15N, na.rm=T),
            d15N.ci = CI(na.omit(d15N))['upper']-mean(d15N, na.rm=T))


#================Correlations between variables, pit by pit==============================
pit.cor <- list()
pit.cor.pvalue <- list()
pit.cor.rsd <- list() #Variables detrended for linear change with depth
pit.cor.pvalue.rsd <- list() #Variables detrended for linear change with depth
pit.rsd <- list()
for(i in 1:length(chic.bypit)) {
  rsd.iter <- data.frame(matrix(nrow=nrow(chic.bypit[[i]]), ncol=length(vbl.clm.index)))
  for(j in 1:length(vbl.clm.index)) {
    coef.iter <- coef(lm(chic.bypit[[i]][,vbl.clm.index[j]]~chic.bypit[[i]]$depth.mean))
    rsd.iter[,j] <- chic.bypit[[i]][,vbl.clm.index[j]] - (chic.bypit[[i]]$depth.mean*coef.iter[2] + coef.iter[1])
  }
  pit.rsd[[i]] <- rsd.iter
  pit.cor[[i]] <- round(cor(chic.bypit[[i]][vbl.clm.index], use="complete.obs"), 2)
  pit.cor.pvalue[[i]] <- round(cor.mtest(as.matrix(chic.bypit[[i]][vbl.clm.index]), use="complete.obs")$p,2)
  pit.cor.rsd[[i]] <- round(cor(pit.rsd[[i]], use="complete.obs"), 2)
  pit.cor.pvalue.rsd[[i]] <- round(cor.mtest(as.matrix(pit.rsd[[i]]), use="complete.obs")$p,2)
}

#Collapsing the detrended pit values into a single tibble
pit.all.rsd <- pit.rsd %>%
  map_df(as_tibble)
colnames(pit.all.rsd) <- colnames(chic.pit[vbl.clm.index])

#Correlations for all pits combined
pit.cor.all <- round(cor(chic.pit[vbl.clm.index], use="complete.obs"), 2)
pit.cor.all.pvalue <- round(cor.mtest(as.matrix(chic.pit[vbl.clm.index]), use="complete.obs")$p,3)
pit.cor.rsd.all <- round(cor(pit.all.rsd, use="complete.obs"), 2)
pit.cor.rsd.all.pvalue <- round(cor.mtest(as.matrix(pit.all.rsd), use="complete.obs")$p,3)
quartz()
corrplot(corr = as.matrix(pit.cor.all), p.mat=as.matrix(pit.cor.all.pvalue), method="ellipse",
         type="upper", tl.col="black", tl.srt=0, sig.level = 0.05,
         tl.cex = 1.5, cl.cex=1.3, diag=FALSE, addCoef.col = "black", number.cex = 1.5)
quartz() #residual
corrplot(corr = as.matrix(pit.cor.rsd.all), p.mat=as.matrix(pit.cor.rsd.all.pvalue), method="ellipse",
         type="upper", tl.col="black", tl.srt=0, sig.level = 0.05,
         tl.cex = 1.5, cl.cex=1.3, diag=FALSE, addCoef.col = "black", number.cex = 1.5)

#=====T tests, LM summaries, and Wilcoxon tests for whether the skin vs depth layers differ in mean values
skl <- chic.no.pseudo.dl[chic.no.pseudo.dl$group == "skl",]
dl <- chic.no.pseudo.dl[chic.no.pseudo.dl$group == "dl",]

skl.ttest <- list()
dl.ttest <- list()
pit.ttest <- list()
full.ttest <- list()
skl.lmstat.inv <- list() #More accurate regression (wNO3 or log(dX +1)  vs. 1/SMB)
dl.lmstat.inv <- list() #More accurate regression (wNO3 or log(dX +1)  vs. 1/SMB)
skd.wilcox <- list()
skl.lm.stats <- data.frame(matrix(nrow=length(vbl.clm.index), ncol=7))
dl.lm.stats <- data.frame(matrix(nrow=length(vbl.clm.index), ncol=7))
for (i in 1:length(vbl.clm.index)) {
  skl.ttest[[i]] <- t.test(skl[vbl.clm.index[i]])
  dl.ttest[[i]] <- t.test(dl[vbl.clm.index[i]])
  pit.ttest[[i]] <- t.test(chic.pit[vbl.clm.index[i]])
  full.ttest[[i]] <- t.test(chic.no.pseudo.dl[vbl.clm.index[i]])
  if (colnames(skl)[vbl.clm.index[i]] == "NO3") {
    skl.lmstat.inv[[i]] <- summary(lm(skl[,vbl.clm.index[i]] ~ I(1/skl$smb)))
  } else {
    skl.lmstat.inv[[i]] <- summary(lm(log(skl[,vbl.clm.index[i]]+1) ~ I(1/skl$smb)))
  }
  skl.lm.stats[i,1] <- names(skl)[vbl.clm.index[i]]
  skl.lm.stats[i,2] <- skl.lmstat.inv[[i]]$coefficients[2]
  skl.lm.stats[i,3] <- skl.lmstat.inv[[i]]$coefficients[4]
  skl.lm.stats[i,4] <- skl.lmstat.inv[[i]]$coefficients[8]
  skl.lm.stats[i,5] <- skl.lmstat.inv[[i]]$coefficients[1]
  skl.lm.stats[i,6] <- skl.lmstat.inv[[i]]$coefficients[3]
  skl.lm.stats[i,7] <- skl.lmstat.inv[[i]]$r.squared
  if (colnames(skl)[vbl.clm.index[i]] == "NO3") {
    dl.lmstat.inv[[i]] <- summary(lm(dl[,vbl.clm.index[i]] ~ I(1/dl$smb)))
  } else {
    dl.lmstat.inv[[i]] <- summary(lm(log(dl[,vbl.clm.index[i]]+1) ~ I(1/dl$smb)))
  }
  dl.lm.stats[i,1] <- names(dl)[vbl.clm.index[i]]
  dl.lm.stats[i,2] <- dl.lmstat.inv[[i]]$coefficients[2]
  dl.lm.stats[i,3] <- dl.lmstat.inv[[i]]$coefficients[4]
  dl.lm.stats[i,4] <- dl.lmstat.inv[[i]]$coefficients[8]
  dl.lm.stats[i,5] <- dl.lmstat.inv[[i]]$coefficients[1]
  dl.lm.stats[i,6] <- dl.lmstat.inv[[i]]$coefficients[3]
  dl.lm.stats[i,7] <- dl.lmstat.inv[[i]]$r.squared
  skd.wilcox[[i]] <- wilcox.test(skl[,vbl.clm.index[i]], dl[,vbl.clm.index[i]])
}
names(skl.ttest) <- names(skl[vbl.clm.index])
names(dl.ttest) <- names(dl[vbl.clm.index])
names(pit.ttest) <- names(chic[vbl.clm.index])
names(full.ttest) <- names(skl[vbl.clm.index])
names(skl.lmstat.inv) <- names(skl[vbl.clm.index])
names(dl.lmstat.inv) <- names(dl[vbl.clm.index])
names(skd.wilcox) <- names(skl[vbl.clm.index])
colnames(skl.lm.stats) <- c("var", "slope", "slope.se", "p.value", "intercept", "intercept.se", "r2")
colnames(dl.lm.stats) <- c("var", "slope", "slope.se", "p.value", "intercept", "intercept.se", "r2")

#=======================Epsilon Calculations=================
#==Calculating apparent fractionation coefficiants based on log(d.snow + 1) = eps * log(NO3)
#=Calculating just for sk-d samples
eps.skd.d18O.lm <- list()
eps.skd.d18O <- NA
eps.skd.d18O.se <- NA
eps.skd.D17O.lm <- list()
eps.skd.D17O <- NA
eps.skd.D17O.se <- NA
eps.skd.d15N.lm <- list()
eps.skd.d15N <- NA
eps.skd.d15N.se <- NA
site <- NA
smb <- NA
#Calculating for sk-d
for (i in 1:length(chic.byskd)) {
  if(nrow(chic.byskd[[i]]) > 1) {
    eps.skd.d18O.lm[[i]] <- lm(log(chic.byskd[[i]]$d18O+1) ~ log(chic.byskd[[i]]$NO3/10^9))
    eps.skd.d18O[i] <- eps.skd.d18O.lm[[i]]$coef[2]*10^3
    eps.skd.d18O.se[i] <- summary(eps.skd.d18O.lm[[i]])$coef[4]*10^3
    eps.skd.D17O.lm[[i]] <- lm(log(chic.byskd[[i]]$D17O+1) ~ log(chic.byskd[[i]]$NO3/10^9))
    eps.skd.D17O[i] <- eps.skd.D17O.lm[[i]]$coef[2]*10^3
    eps.skd.D17O.se[i] <- summary(eps.skd.D17O.lm[[i]])$coef[4]*10^3
    eps.skd.d15N.lm[[i]] <- lm(log(chic.byskd[[i]]$d15N+1) ~ log(chic.byskd[[i]]$NO3/10^9))
    eps.skd.d15N[i] <- eps.skd.d15N.lm[[i]]$coef[2]*10^3
    eps.skd.d15N.se[i] <- summary(eps.skd.d15N.lm[[i]])$coef[4]*10^3
    site[i] <- as.character(chic.byskd[[i]]$site[1])
    smb[i] <- chic.byskd[[i]]$smb[1]
  }
}

#Making a list of all the epsilon linear models
eps.skd.lm <- list(eps.skd.d18O.lm, eps.skd.D17O.lm, eps.skd.d15N.lm)

#Table of epsilon results for skin layer-depth layer samples (NOTE: no stnd error or p values
# because only two points in regression)
eps.skd.chic.pit <- data.frame(site, eps.skd.d18O, eps.skd.D17O, eps.skd.d15N, smb)
eps.skd.chic.pit <- eps.skd.chic.pit[order(-eps.skd.chic.pit$smb),]
eps.skd.chic.pit[,-1] <- round(eps.skd.chic.pit[,-1], 2)
print(eps.skd.chic.pit)

#Calculating the linear regressions of eps versus SMBs
eps.clm <- c("eps.skd.d18O", "eps.skd.D17O", "eps.skd.d15N") #columns with needed data
eps.clm.index <- NA
for (i in 1:length(eps.clm)) {
  eps.clm.index[i] <- which(colnames(eps.skd.chic.pit) == eps.clm[i])
}

eps.smb.lmstat <- list()
for (i in 1:length(eps.clm.index)) {
  eps.smb.lmstat[[i]] <- summary(lm(eps.skd.chic.pit[,eps.clm.index[i]] ~ eps.skd.chic.pit$smb))
}
names(eps.smb.lmstat) <- names(eps.skd.chic.pit[eps.clm.index])


#=Calculating for pits using full pit data (NOTE: Results rejected because of poor regressions)
eps.pit.d18O.lm <- list()
eps.pit.d18O <- NA
eps.pit.d18O.se <- NA
eps.pit.d18O.p <- NA
eps.pit.D17O.lm <- list()
eps.pit.D17O <- NA
eps.pit.D17O.se <- NA
eps.pit.D17O.p <- NA
eps.pit.d15N.lm <- list()
eps.pit.d15N <- NA
eps.pit.d15N.se <- NA
eps.pit.d15N.p <- NA
site <- NA
smb <- NA

#Calculating for pits
for (i in 1:length(chic.bypit)) {
  if(nrow(chic.bypit[[i]]) > 1) {
    eps.pit.d18O.lm[[i]] <- lm(log(chic.bypit[[i]]$d18O+1) ~ log(chic.bypit[[i]]$NO3/10^9))
    eps.pit.d18O[i] <- eps.pit.d18O.lm[[i]]$coef[2]*10^3
    eps.pit.d18O.se[i] <- summary(eps.pit.d18O.lm[[i]])$coef[4]*10^3
    eps.pit.d18O.p[i] <- summary(eps.pit.d18O.lm[[i]])$coef[8]
    eps.pit.D17O.lm[[i]] <- lm(log(chic.bypit[[i]]$D17O+1) ~ log(chic.bypit[[i]]$NO3/10^9))
    eps.pit.D17O[i] <- eps.pit.D17O.lm[[i]]$coef[2]*10^3
    eps.pit.D17O.se[i] <- summary(eps.pit.D17O.lm[[i]])$coef[4]*10^3
    eps.pit.D17O.p[i] <- summary(eps.pit.D17O.lm[[i]])$coef[8]
    eps.pit.d15N.lm[[i]] <- lm(log(chic.bypit[[i]]$d15N+1) ~ log(chic.bypit[[i]]$NO3/10^9))
    eps.pit.d15N[i] <- eps.pit.d15N.lm[[i]]$coef[2]*10^3
    eps.pit.d15N.se[i] <- summary(eps.pit.d15N.lm[[i]])$coef[4]*10^3
    eps.pit.d15N.p[i] <- summary(eps.pit.d15N.lm[[i]])$coef[8]
    site[i] <- as.character(chic.bypit[[i]]$site[1])
    smb[i] <- chic.bypit[[i]]$smb[1]
  }
}

#Making a list and predictions of all the epsilon linear models
eps.pit.lm <- list(eps.pit.d18O.lm, eps.pit.D17O.lm, eps.pit.d15N.lm)

#Table of epsilon results for pits
eps.chic.pit <- data.frame(site, eps.pit.d18O, eps.pit.d18O.se, eps.pit.d18O.p, eps.pit.D17O, eps.pit.D17O.se,
                           eps.pit.D17O.p, eps.pit.d15N, eps.pit.d15N.se, eps.pit.d15N.p, smb)
eps.chic.pit <- eps.chic.pit[order(-eps.chic.pit$smb),]
eps.chic.pit[,-1] <- round(eps.chic.pit[,-1], 2)

#### END Epsilons of Pits and SKD ####


#=======================Plots===============================
#### Violin plots ####
#This makes a plot of the general structure of data sets by the different groups. It includes
# the pit pseudo depth layer samples in the depth layer section and the skl samples in the
# pits, but does not have duplicate samples in any one group.

#Removing duplicate skl samples used to match with pit pseudo depth layers
chic.violin <- chic[!(chic$pseudo ==1 & chic$group == "skl"), ]

p.NO3 <- ggplot(data=chic.violin) +
  theme_classic() +
  geom_violin(aes(x=group, y=log(NO3), color=group, fill=group), linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
  geom_violin(aes(x=group, y=log(NO3), color=group), fill='transparent', draw_quantiles = c(0.5)) +
  scale_color_manual(values=brewer.pal(length(unique(chic.violin$group)), "Set2"), name = "Sample Type", labels=toupper(levels(chic.violin$group))) +
  scale_fill_manual(values=alpha(brewer.pal(length(unique(chic.violin$group)), "Set2"), 0.4), name = "Sample Type", labels=toupper(levels(chic.violin$group))) +
  scale_y_continuous(name="NO3 (ng g-1)", breaks = log(c(50, seq(100,600,100))), labels = c(50,seq(100,600,100)), position = "left") +
  scale_x_discrete(labels=toupper(levels(chic.violin$group)), position = "bottom", name=NULL) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color=brewer.pal(length(unique(chic.violin$group)), "Set2")),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

p.d18O <- ggplot(data=chic.violin) +
  theme_classic() +
  geom_violin(aes(x=group, y=d18O, color=group, fill=group), linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
  geom_violin(aes(x=group, y=d18O, color=group), fill='transparent', draw_quantiles = c(0.5)) +
  scale_color_manual(values=brewer.pal(length(unique(chic.violin$group)), "Set2"), name = "Sample Type", labels=toupper(levels(chic.violin$group))) +
  scale_fill_manual(values=alpha(brewer.pal(length(unique(chic.violin$group)), "Set2"), 0.4), name = "Sample Type", labels=levels(chic.violin$group)) +
  scale_y_continuous(name="d18O (‰)", position = "left") +
  scale_x_discrete(labels=toupper(levels(chic.violin$group)), position = "bottom", name=NULL) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color=brewer.pal(length(unique(chic.violin$group)), "Set2")),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

p.D17O <- ggplot(data=chic.violin) +
  theme_classic() +
  geom_violin(aes(x=group, y=D17O, color=group, fill=group), linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
  geom_violin(aes(x=group, y=D17O, color=group), fill='transparent', draw_quantiles = c(0.5)) +
  scale_color_manual(values=brewer.pal(length(unique(chic.violin$group)), "Set2"), name = "Sample Type", labels=levels(chic.violin$group)) +
  scale_fill_manual(values=alpha(brewer.pal(length(unique(chic.violin$group)), "Set2"), 0.4), name = "Sample Type", labels=levels(chic.violin$group)) +
  scale_y_continuous(name="D17O (‰)", position = "left") +
  scale_x_discrete(labels=toupper(levels(chic.violin$group)), position = "bottom", name=NULL) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color=brewer.pal(length(unique(chic.violin$group)), "Set2")),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

p.d15N <- ggplot(data=chic.violin) +
  theme_classic() +
  geom_violin(aes(x=group, y=d15N, color=group, fill=group), linetype='dotted', draw_quantiles = c(0.25, 0.75)) +
  geom_violin(aes(x=group, y=d15N, color=group), fill='transparent', draw_quantiles = c(0.5)) +
  scale_color_manual(values=brewer.pal(length(unique(chic.violin$group)), "Set2"), name = "Sample Type", labels=levels(chic.violin$group)) +
  scale_fill_manual(values=alpha(brewer.pal(length(unique(chic.violin$group)), "Set2"), 0.4), name = "Sample Type", labels=levels(chic.violin$group)) +
  scale_y_continuous(name="d15N (‰)", position = "left") +
  scale_x_discrete(labels=toupper(levels(chic.violin$group)), position = "bottom", name=NULL) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color=brewer.pal(length(unique(chic.violin$group)), "Set2")),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))

quartz(height=8, width=15)
plot_grid(p.NO3, p.d18O, p.D17O, p.d15N, ncol=2, align = "hv")

#Output as PDF
#pdf("Figures/CHIC Violin Plots.pdf", height=8, width=15)
#plot_grid(p.NO3, p.d18O, p.D17O, p.d15N, ncol=2, align = "hv")
#dev.off()

#### END Violin Plots ####

#### Depth Series Plots ####
chic.bypit.ts.plot <- list()
chic.ts.plot <- list()
ylims <- data.frame(c(30,45,20,-18),c(265,100,42,100))
colnames(ylims) <- c("ymin", "ymax")
pit.index <- c("p1", "p2", "p3", "p4", "p5")

#Plot creation
for (j in 1:length(chic.bypit)) {
  for (i in 1:length(vbl.clm.index)) 
    local({ #Have to do local or else plot later only does each plot in most recent i
      i <- i
      color.select <- as.character(iso.color.index[iso.color.index$iso == colnames(chic.bypit[[j]][vbl.clm.index[i]]),2])
      chic.ts.plot[[i]] <<- ggplot() +
        theme_classic() +
        {if(colnames(chic.skd)[vbl.clm.index[i]] == "NO3") {
          geom_smooth(aes(x=chic.bypit[[j]]$depth.top, y=chic.bypit[[j]][,vbl.clm.index[i]]), method="glm", formula=y~x, method.args=list(family=gaussian(link='inverse')),
                      linetype='longdash', color="gray60", fill=NA) 
        } else {
          geom_smooth(aes(x=chic.bypit[[j]]$depth.top, y=chic.bypit[[j]][,vbl.clm.index[i]]*1000), method="glm", formula=y~x, method.args=list(family=gaussian(link='inverse')),
                      linetype='longdash', color="gray60", fill=NA) 
        }} +
        {if(colnames(chic.skd)[vbl.clm.index[i]] == "NO3") {
          geom_stepribbon(aes(x=chic.bypit[[j]]$depth.top, ymax=chic.bypit[[j]][,vbl.clm.index[i]]+chic.bypit[[j]][,vbl.clm.index[i]+5],
                              ymin=chic.bypit[[j]][,vbl.clm.index[i]]-chic.bypit[[j]][,vbl.clm.index[i]+5]), fill=color.index[i], alpha=0.2) 
        } else {
          geom_stepribbon(aes(x=chic.bypit[[j]]$depth.top, ymax=chic.bypit[[j]][,vbl.clm.index[i]]*1000+chic.bypit[[j]][,vbl.clm.index[i]+5]*1000,
                              ymin=chic.bypit[[j]][,vbl.clm.index[i]]*1000-chic.bypit[[j]][,vbl.clm.index[i]+5]*1000), fill=color.index[i], alpha=0.2) 
        }} +
        {if(colnames(chic.skd)[vbl.clm.index[i]] == "NO3") {
          geom_step(aes(x=chic.bypit[[j]]$depth.top, y=chic.bypit[[j]][,vbl.clm.index[i]]), color=color.select)
        } else {
          geom_step(aes(x=chic.bypit[[j]]$depth.top, y=chic.bypit[[j]][,vbl.clm.index[i]]*1000), color=color.select)
        }} +
        scale_x_continuous(name="Depth (cm)") +
        scale_y_continuous(name=colnames(chic.bypit[[j]][vbl.clm.index[i]])) +
        coord_cartesian(ylim=c(ylims[i,1], ylims[i,2])) +
        theme(axis.title.y = element_text(size=15, color=color.select),
              axis.text.y = element_text(size=12, color=color.select),
              axis.line.y = element_line(color=color.select),
              axis.ticks.y = element_line(color=color.select),
              axis.text.x = element_text(size=12, color="gray35"),
              axis.title.x = element_text(size=15, color="gray35"),
              axis.line.x = element_line(color="gray35"),
              plot.title = element_text(hjust=0.5, size=24, color=color.select))
    })
  chic.bypit.ts.plot[[j]] <- chic.ts.plot
}
#NOTE: The inverse regression trends don't match well for d15N, but that's because the data doesn't support it.
# I tested with fake data that looks like the assumed "rise, then plateau" trend, and the script matched that perfectly.

plotgrid.list <- list()
for (j in 1:length(pit.index)) {
  plotgrid.list[[j]] <- plot_grid(plotlist=chic.bypit.ts.plot[[j]], ncol=1, align = "v")
}

quartz(height=8, width=15)
plot_grid(plotlist=plotgrid.list, ncol=length(pit.index))

#Output as PDF
#pdf("Figures/CHIC Pit Depth Series.pdf", height=8, width=15)
#plot_grid(plotlist=plotgrid.list, ncol=length(pit.index))
#dev.off()

#### END Depth Series Plots ####

##### Skin & Depth by SMB Plots ####
skd.ts.plot <- list()
breaks.yaxis <- list(seq(0,500,100), seq(50,90,10), seq(25,35,5), seq(-25,75,25))
limits.yaxis <- list(c(0,510), log(c(50,92)/1000+1), log(c(25,38)/1000+1),
                     log(c(-25,80)/1000+1))

for (i in 1:length(vbl.clm.index)) 
  local({ #Have to do local or else plot later only does each plot in most recent i
    i <- i
    skd.ts.plot[[i]] <<- ggplot(data = chic.skd, aes(x=smb, y=chic.skd[,vbl.clm.index[i]], color=group)) +
      theme_classic() +
      geom_point(shape=16) +
      {if (colnames(chic.skd)[vbl.clm.index[i]] == "NO3") {
        geom_smooth(method="lm", aes(fill=group), formula = y~I(1/x))
      } else {
        geom_smooth(method="lm", aes(fill=group), formula = log(y+1)~I(1/x))
      }} +
      scale_color_manual(values=c(lighten(as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2]), factor=0.3),
                                  darken(as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2]))), limits = c("skl", "dl")) +
      scale_fill_manual(values=c(lighten(as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2]), factor=0.3),
                                 darken(as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2]))), limits = c("skl", "dl")) +
      scale_x_continuous(name="SMB (kg m-2 a-1)") +
      {if (colnames(chic.skd)[vbl.clm.index[i]] == "NO3") {
        scale_y_continuous(name=colnames(chic.skd[vbl.clm.index[i]]), breaks=breaks.yaxis[[i]], limits = limits.yaxis[[i]])
      } else {
        scale_y_continuous(name=colnames(chic.skd[vbl.clm.index[i]]), breaks=log(breaks.yaxis[[i]]/1000+1),
                           labels=function(y) round(1000*(exp(y)-1)), limits = limits.yaxis[[i]])
      }} +
      theme(axis.title.y = element_text(size=15, color=as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2])),
            axis.text.y = element_text(size=12, color=as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2])),
            axis.line.y = element_line(color=as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2])),
            axis.ticks.y = element_line(color=as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2])),
            axis.text.x = element_text(size=12, color="gray35"),
            axis.title.x = element_text(size=13, color="gray35"),
            axis.line.x = element_line(color="gray35"),
            plot.title = element_text(hjust=0.5, size=24, color=as.character(iso.color.index[iso.color.index$iso == colnames(chic.skd[vbl.clm.index[i]]),2])))
  })

quartz(height=8, width=5)
plot_grid(skd.ts.plot[[1]], skd.ts.plot[[4]], skd.ts.plot[[2]],
          skd.ts.plot[[3]], ncol=1, align = "v")

#Output as PDF
#pdf("Figures/SKD SMB Trends.pdf", height=8, width=5)
#plot_grid(skd.ts.plot[[1]], skd.ts.plot[[4]], skd.ts.plot[[2]],
#          skd.ts.plot[[3]], ncol=1, align = "v")
#dev.off()

#### END SKL & DL by SMB ####


####========Plotting the CHICTABA elevation and SMB profiles====####
#Data import
chic.rema <- as_tibble(read.csv("chictaba_rema_profile.csv", header=TRUE, sep=","))
chic.MAR <- as_tibble(read.csv("chictaba_SMBMAR_profile.csv", header=TRUE, sep=","))

#quartz(height=6, width=12)
#pdf("Figures/chic_rema_profile.pdf", height=6, width=12)
p.chic.rema <- ggplot(data=chic.rema) +
  theme_classic() +
  geom_step(aes(x=plot.dist.mean, y=elevation.rema), direction="vh") +
  scale_y_continuous(name="REMA Elevation (m a.s.l.)") +
  scale_x_continuous(name="Distance Along Transect (km)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))
#dev.off()

#quartz(height=6, width=12)
#pdf("Figures/chic_MAR_profile.pdf", height=6, width=12)
p.chic.MAR <- ggplot(data=chic.MAR) +
  theme_classic() +
  geom_step(aes(x=plot.dist.mean, y=smb.mar), direction="vh") +
  scale_y_continuous(name="SMB (kg m-2 a-1)") +
  scale_x_continuous(name="Distance Along Transect (km)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=16, color='gray40'),
        axis.text.x =  element_text(size=14, color='gray40'),
        axis.ticks.x = element_line(color='gray40'),
        axis.line.x = element_line(color='gray40'),
        axis.title.y = element_text(size=16, color='gray40'),
        axis.text.y = element_text(size=14, color='gray40'),
        axis.ticks.y = element_line(color='gray40'),
        axis.line.y = element_line(color='gray40'))
#dev.off()

quartz(height=8, width=12)
#pdf("Figures/chic_rema_MAR_profile.pdf", height=8, width=12)
plot_grid(plotlist=list(p.chic.rema,p.chic.MAR), ncol=1, align="hv")
#dev.off()

####====END PLOTTING REMA SMB PROFILES====