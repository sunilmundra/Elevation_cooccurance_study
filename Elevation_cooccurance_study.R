# I joined succesfully! Miki

library(mgcv) ## (Wood, 2011)
library(vegan) ## (Oksanen et al., 2013)
library(mvabund)
library(boral)

getwd()
setwd("C:/Users/admin/Desktop/elevation")
getwd()

#load("elevation.image.RData")
#names(OTUs.fit.nb)
#ls()

## 1. Import and organize community data
raw.abundance = read.csv(file="Final_OTU_matrix.csv",
                         sep=",", header=T, row.names=1)

raw.abundance = raw.abundance[order(row.names(raw.abundance)),]
raw.abundance

## read the experimental variables
metadata = read.csv(file="Metadata.csv",
                    header=T, row.names=1)
metadata
attach(metadata)

save.image("elevation.image.RData")
## set the order for the experimental factors
#metadata$Location
#experiment = factor(metadata$Location,
#                    levels=c("ava", "boy", "car", "den", "fbk", "ftm", "gpr", "hay", "lov", "mel", "por"))
#experiment


#Then import the stand column separately
#stand.data<-read.csv(file="stand.csv",
#                     sep=",", header=T, row.names=1)
#attach(stand.data)
#names(stand.data)
#Converting variable from numeric to factor:
#mode(Stand)
#stand<-as.factor(Stand)
#stand
#mode(stand) 
#is.factor(stand)

### first: northern garden, natural conditions
experiment.time = factor(metadata$season,
                         levels=c("May", "September"))

experiment.loc = factor(metadata$reg,
                        levels=c("D", "S"))

experiment.veg = factor(metadata$mh,
                        levels=c("B", "V"))


#history1 = factor(metadata$history1)
#history2 = metadata$history2
#FD = metadata$FD
#SD = metadata$SD
#H = metadata$H
#MaxL = metadata$MaxL


## Uneven sequencing depth may have an impact
readNumber = apply(raw.abundance,1,sum) # summing read per samples
readNumber
otuCount =  apply(raw.abundance,2,sum) # summing read per OTU
otuCount
richness = specnumber(raw.abundance, MARGIN = 1)
richness

##################################################
# 1. Relationship between read number and OTU richness
#################################################
LM1 = lm(richness~readNumber)
summary(LM1)
LM1
plot (readNumber,richness)
abline(LM1)

## OTU abundances
png(file="reads_vs_richness.png", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)
par(mar = c(4,4,1,1))
plot (readNumber,richness,xaxt="n", xlab="Reads numbers", ylab="OTU richness", main=NA, col="black",pch=16)
abline(LM1)
mylegend = legend("topright", c("t value = 8.47", "p value < 0.001"), bty="n")
dev.off()

## OTU abundances
png(file="histogram.png", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)
par(mar = c(4,4,1,1))
hist(log10(otuCount), xaxt="n", xlab="Abundance of OTUs (DNA reads)", main=NA, col="grey")
ticks.3 = seq(1,6, by=1)
labels.3 <- sapply(ticks.3, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(1,2,3,4,5,6), labels=labels.3)
dev.off()



###########################################
# 2. Env variable corlation test
##########################################
library(corrplot)

#Normalizing env variables
open.scale = scale(open)
needle.scale = scale(needle)
grass.scale = scale(grass)
erica.scale = scale(erica)
junip.scale = scale(junip)
shrubs.scale = scale(shrubs)
soil_mean.log = soil_mean
pH.log = pH
T_mean.log = T_mean
CO.scale=scale(CO)

model.season = model.matrix(~season)
model.season
season.scale=scale(model.season)
season.scale


microhabitat_data<-data.frame(elevation_dGPS,open.scale, needle.scale, grass.scale, erica.scale, junip.scale, shrubs.scale, 
                              soil_mean.log, pH.log, T_mean.log,CO.scale)
microhabitat_data

X<-data.frame(season,reg,elevation_dGPS, needle.scale,grass.scale, shrubs.scale, 
              soil_mean.log, pH.log, T_mean.log)

#mtcars <- metadata[,4:10]
# Corleation test
M <- cor(microhabitat_data, use = "na.or.complete")

png(file="env_cor_circle.png", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)
corrplot(M, method="circle")
dev.off()
#corrplot(M, type="upper")
png(file="env_cor_mix_plot.png", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)
corrplot.mixed(M)
dev.off()
#corrplot(M, order="hclust", addrect=2,  bg="gold2")
#corrplot(M, order="FPC", addrect=2,  bg="gold2")
#cor.test(needle,season, method="kendall")


#shapiro.test(needle) # p value < 0.05 data differ significantly from normal distribution
## Plot using a qqplot
#qqnorm(needle);qqline(needle, col = 2)
#hist(open, col='grey')
#hist(needle, col='grey', freq = FALSE)
#xbar=mean(needle)
#S=sd(needle)
#curve(dnorm(needle, xbar, S), col = 2, add=TRUE)



#principal component analysis
library(stats) 
prin_comp <- prcomp(microhabitat_data, scale = TRUE)
prin_comp <- prcomp(~ elevation_dGPS+ open.scale + needle.scale + grass.scale + shrubs.scale + 
                      soil_mean.log + pH.log +T_mean.log+CO.scale, scale = TRUE)

names(prin_comp)
#outputs the mean of variables
prin_comp$center

#outputs the standard deviation of variables
prin_comp$scale
prin_comp$rotation
prin_comp$rotation[1:9,1:4]
biplot(prin_comp, scale = 0)
#compute standard deviation of each principal component
std_dev <- prin_comp$sdev
#compute variance
pr_var <- std_dev^2
#check variance of first 10 components
pr_var[1:10]
#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]
#scree plot
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")


# PCA analyisis
library("factoextra")
#pr_comp <- prcomp(~ open.scale + needle.scale + grass.scale + shrubs.scale + 
#               soil_mean.log + pH.log +T_mean.log, scale = TRUE)

#prin_comp <- princomp(~ open.scale + needle.scale + grass.scale + shrubs.scale + 
#         soil_mean.log + pH.log +T_mean.log, scale = TRUE, scores = TRUE)

res.pca <- prcomp(~ elevation_dGPS+open.scale + needle.scale + grass.scale + erica.scale+junip.scale+ shrubs.scale + 
                    soil_mean.log + pH.log +T_mean.log+CO.scale, scale = TRUE)

names(res.pca)
head(res.pca$sdev)
head(unclass(res.pca$rotation)[, 1:4])
# Eigenvalues
eig <- (res.pca$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)
summary(res.pca)

eig.val <- get_eigenvalue(res.pca)
head(eig.val)

barplot(eig.decathlon2.active[, 2], names.arg=1:nrow(eig.decathlon2.active), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.decathlon2.active), 
      eig.decathlon2.active[, 2], 
      type="b", pch=19, col = "red")

fviz_screeplot(res.pca, ncp=10)
fviz_screeplot(res.pca, ncp=10, choice="eigenvalue")
var <- get_pca_var(res.pca)
var
# Coordinates of variables
var$coord[, 1:4]
# Helper function : 
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- res.pca$rotation
sdev <- res.pca$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord[, 1:4])

# Plot the correlation circle
a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = "PC1",  ylab = "PC2")
abline(h = 0, v = 0, lty = 2)
# Add active variables
arrows(0, 0, var.coord[, 1], var.coord[, 2], 
       length = 0.1, angle = 15, code = 2)
# Add labels
text(var.coord, labels=rownames(var.coord), cex = 1, adj=1)

fviz_pca_var(res.pca)

var.cos2 <- var.coord^2
head(var.cos2[, 1:4])

fviz_pca_var(res.pca, col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=55) + theme_minimal()

comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])

fviz_pca_var(res.pca, col.var="contrib") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=50) + theme_minimal()

ind.coord <- res.pca$x
head(ind.coord[, 1:4])


################################################
# 3. Define core OTUs
################################################
## Summarize reads
TotCount = apply(raw.abundance,2,sum)
TotCount

## The average read number of OTUs
MeanCount=apply(raw.abundance,2,function(vec) mean(vec[vec>0]))
MeanCount
## In how many samples is an OTU present?
TotPresent = apply(raw.abundance,2,function(vec) sum(vec>0))
TotPresent
## The highest read number of an OTU in a sample
MaxCount=apply(raw.abundance,2,max)
MaxCount
## Plotting incidence against abundance
plot(TotPresent, MaxCount, xlab="Incidence",
     ylab="Maximum Abundance", pch=20)

plot(TotPresent, log(MaxCount), xlab="Incidence",
     ylab="log(Maximum Abundance)", pch=20)

## Create a smoothed trendline
gam.control(nthreads=1,irls.reg=0.0,epsilon = 1e-07, maxit = 500,
            mgcv.tol=1e-7,mgcv.half=15, trace = FALSE,
            rank.tol=.Machine$double.eps^0.5,
            nlm=list(),optim=list(),newton=list(),
            outerPIsteps=0,idLinksBases=TRUE,scalePenalty=TRUE,
            keepData=FALSE,scale.est="pearson") 

gam1 = gam(log(MaxCount)~s(TotPresent))
summary(gam1)



png(file="GAM_model.png", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)
plot(gam1, residuals=T, shade=T, rug=F, cex=2.6,
     xlab="Incidence", ylab="logMean Abundance") # , xaxp=c(0,150,15)
dev.off()

## keep core OTUs
IsFreq = TotPresent > 143
fun.some = raw.abundance[,IsFreq]
# write.table(fun.some, file = "core_OTUs.csv", sep = ",")

## How many of these in differnt season?
some.May = fun.some[experiment.time == "May",]
length(colnames(some.May)[apply(some.May,2,sum) > 0])

some.Sept = fun.some[experiment.time == "September",]
length(colnames(some.Sept)[apply(some.Sept,2,sum) > 0])

## How many of these in differnt Location?
some.D = fun.some[experiment.loc == "D",]
length(colnames(some.D)[apply(some.D,2,sum) > 0])

some.S = fun.some[experiment.loc == "S",]
length(colnames(some.S)[apply(some.S,2,sum) > 0])

## How many of these in differnt vegetation type?
some.B = fun.some[experiment.veg == "B",]
length(colnames(some.B)[apply(some.B,2,sum) > 0])

some.V = fun.some[experiment.veg == "V",]
length(colnames(some.V)[apply(some.V,2,sum) > 0])

## Model selection for path core OTUs####
fun.mvabund_core_path = mvabund(fun.some)

#null model
fun.m = manyglm(fun.mvabund_core_path ~ sqrt(readNumber), 
                family="negative.binomial", show.residuals=T)
fun.m

fun.m1 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + reg,
                 family="negative.binomial", show.residuals=T)
fun.m1
#ANOVA to test difference
m.m1.anova = anova(fun.m, fun.m1, nBoot=100, test="LR")
m.m1.anova


fun.m2 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + season,
                 family="negative.binomial", show.residuals=T)
fun.m2
#ANOVA to test difference
m.m2.anova = anova(fun.m, fun.m2, nBoot=100, test="LR")
m.m2.anova

fun.m3 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + elevation_dGPS,
                 family="negative.binomial", show.residuals=T)
fun.m3
#ANOVA to test difference
m.m3.anova = anova(fun.m, fun.m3, nBoot=100, test="LR")
m.m3.anova

fun.m4 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + needle.scale,
                 family="negative.binomial", show.residuals=T)
fun.m4
#ANOVA to test difference
m.m4.anova = anova(fun.m, fun.m4, nBoot=100, test="LR")
m.m4.anova

fun.m5 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + shrubs.scale,
                 family="negative.binomial", show.residuals=T)
fun.m5
#ANOVA to test difference
m.m5.anova = anova(fun.m, fun.m5, nBoot=100, test="LR")
m.m5.anova

fun.m6 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + grass.scale,
                 family="negative.binomial", show.residuals=T)
fun.m6
#ANOVA to test difference
m.m6.anova = anova(fun.m, fun.m6, nBoot=100, test="LR")
m.m6.anova

fun.m7 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + erica.scale,
                 family="negative.binomial", show.residuals=T)
fun.m7
#ANOVA to test difference
m.m7.anova = anova(fun.m, fun.m7, nBoot=100, test="LR")
m.m7.anova

fun.m8 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + junip.scale,
                 family="negative.binomial", show.residuals=T)
fun.m8
#ANOVA to test difference
m.m8.anova = anova(fun.m, fun.m8, nBoot=100, test="LR")
m.m8.anova


fun.m9 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + Forest_type,
                 family="negative.binomial", show.residuals=T)
fun.m9
#ANOVA to test difference
m.m9.anova = anova(fun.m, fun.m9, nBoot=100, test="LR")
m.m9.anova

fun.m10 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + soil_mean.log,
                  family="negative.binomial", show.residuals=T)
fun.m10
#ANOVA to test difference
m.m10.anova = anova(fun.m, fun.m10, nBoot=100, test="LR")
m.m10.anova


fun.m11 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + T_mean.log,
                  family="negative.binomial", show.residuals=T)
fun.m11
#ANOVA to test difference
m.m11.anova = anova(fun.m, fun.m11, nBoot=100, test="LR")
m.m11.anova

fun.m12 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + pH.log,
                  family="negative.binomial", show.residuals=T)
fun.m12
#ANOVA to test difference
m.m12.anova = anova(fun.m, fun.m12, nBoot=100, test="LR")
m.m12.anova

fun.m13 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + open.scale,
                  family="negative.binomial", show.residuals=T)
fun.m13
#ANOVA to test difference
m.m13.anova = anova(fun.m, fun.m13, nBoot=100, test="LR")
m.m13.anova

fun.m2.1 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + season + shrubs.scale,
                   family="negative.binomial", show.residuals=T)
fun.m2.1
#ANOVA to test difference
m.m2.1.anova = anova(fun.m, fun.m2, fun.m2.1, nBoot=10, test="LR")
m.m2.1.anova

fun.m2.2 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + season + grass.scale,
                   family="negative.binomial", show.residuals=T)
fun.m2.2
#ANOVA to test difference
m.m2.2.anova = anova(fun.m, fun.m2, fun.m2.2, nBoot=100, test="LR")
m.m2.2.anova

fun.m2.3 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + season + Forest_type,
                   family="negative.binomial", show.residuals=T)
fun.m2.3
#ANOVA to test difference
m.m2.3.anova = anova(fun.m, fun.m2, fun.m2.3, nBoot=100, test="LR")
m.m2.3.anova

fun.m3.1 = manyglm(fun.mvabund_core_path ~ sqrt(readNumber) + season + shrubs.scale + grass.scale,
                   family="negative.binomial", show.residuals=T)
fun.m3.1
#ANOVA to test difference
m.m3.1.anova = anova(fun.m, fun.m2, fun.m2.1, fun.m3.1, nBoot=100, test="LR")
m.m3.1.anova



#If full model fits better:
#see individual p-Values for each OTU 
## Analysis of variance explained by the predictors
fun.anova.mc_path = anova(fun.m2.1, nBoot=100, test="LR", p.uni="adjusted", 
                          resamp="montecarlo")
fun.anova.mc_path

# ## OTUs significantly explained by the affected by the experiment (at p<=0.01)
p.ind.anova.mc <- as.data.frame(fun.anova.mc_path$uni.p)
fun.exp.mc = colnames(p.ind.anova.mc)[p.ind.anova.mc[3,]<=0.01]

## Visualization of experimental and genetic factor interactions
fun.m6.coef = as.data.frame(fun.m2.1$coefficients)
fun.m6.coef

## Extreme outlier
filter.temp = fun.m6.coef == min(fun.m6.coef)
filter.temp
which(filter.temp, arr.ind=T)
names(fun.m6.coef)[458]

## The value of the extreme outlier
fun.m6.coef[4,458]
data.frame(rownames(fun.m6.coef), fun.m6.coef$OTU_1417)

## coefficients without outlier
fun.m6.coef.noout = fun.m6.coef[colnames(fun.m6.coef) != "OTU_1417"]

## mean-centering the contrasts
coef.mean.contrast = fun.m6.coef.noout - apply(fun.m6.coef.noout,2,mean)

## plot mean-centered contrasts
par(mar=c(4,4,1,1))
plot(c(0,5), c(-50,100), type="n", xaxt="n", 
     ylab="Model coefficients", xlab="")
axis(1, at=c(0:5), labels=c("South", "North\nheated", "North\nnatural", 
                            "South", "North\nheated", "North\nnatural"),
     hadj=1, las=2)
for (i in colnames(coef.mean.contrast)) {
  lines(c(0:5), 
        c(0,
          coef.mean.contrast[colnames(coef.mean.contrast) == i][3,],
          coef.mean.contrast[colnames(coef.mean.contrast) == i][4,],
          coef.mean.contrast[colnames(coef.mean.contrast) == i][5,],
          coef.mean.contrast[colnames(coef.mean.contrast) == i][7,],
          coef.mean.contrast[colnames(coef.mean.contrast) == i][8,]), lwd=0.3)
}

## plot mean-centered contrasts
par(mar=c(4,4,1,1))
plot(c(0,4), c(-50,100), type="n", xaxt="n", 
     ylab="Model coefficients", xlab="")
axis(1, at=c(0:3), labels=c("May", "September", "May", "September"),
     hadj=1, las=2)
for (i in colnames(coef.mean.contrast)) {
  lines(c(0:3), 
        c(0,
          coef.mean.contrast[colnames(coef.mean.contrast) == i][3,],
          coef.mean.contrast[colnames(coef.mean.contrast) == i][4,],
          coef.mean.contrast[colnames(coef.mean.contrast) == i][5,],lwd=0.3)
}


# Boral ordination
# LVM of horizons
# simple model, similar to the LVM ordination
OTU.m0 = manyglm(OTUSummed.mva ~
                   SummedExpReord$reads,
                 family = "negative.binomial")


# diagnostics
plot(fun.m3.1, which = c(1:3))

# distribution of the overdispersion parameters
hist(fun.m3.1$theta)
hist(log(fun.m3.1$theta)) # a few OTUs have weird dispersion parameters

# OTUs with weird overdispersion
names(fun.some[,fun.m3.1$theta > 30])

# LVM model
# Set overdispersion prior
set.prior = list(type = c("normal","normal","normal","uniform"),
                 hypparams = c(100, 20, 100, 20))
model.var = data.frame(sqrt(readNumber), season, grass.scale)
model.var
# LVM ordination of horizons
horizon.comm.ord = boral(fun.some[,fun.m3.1$theta < 30],
                         X = model.var,
                         family = "negative.binomial",
                         prior.control = set.prior, num.lv = 2, 
                         n.burnin = 10,
                         n.iteration = 100, n.thin = 1)



plot(c(1:21), c(1:21), lwd=3,
     col=(SummedExpReord$depth_nominal+1)*10, pch = 
       SummedExpReord$depth_nominal+1)

# depth_names = unique(factor(ExpPredictor$depth))
# pdf(file = "LVM_technical_replicates.pdf", width = 10)
par(mfrow = c(1,1), mar = c(4,4,1,1))
hor.ordicomm = ordiplot(horizon.comm.ord$lv.median, choices = c(1,2), 
                        type = "none", cex =0.5,
                        display = "sites")
points(hor.ordicomm,"sites", lwd=4,
       col=(SummedExpReord$depth_nominal+1)*10,
       pch = SummedExpReord$depth_nominal+1,
       cex = 1.5)
legend(-0.05, 0.08, SummedExpReord$depth, border="white",
       col=(SummedExpReord$depth_nominal+1)*10,
       pch = SummedExpReord$depth_nominal+1,
       cex = 0.9, bty="n", lwd = 2, ncol = 2)
ordisurf(hor.ordicomm, POP_short$DDE4.4,
         add=T, col = "red", lwd = 1.5, cex = 1.3)
ordisurf(hor.ordicomm, POP_short$depth,
         add=T, col = "grey", lwd = 1.5, cex = 1.3)


### 4. Visualize differences in community composition (core communities)
## run NMDS
MDS.all <- metaMDS(fun.some)
MDS.all <- metaMDS(fun.some, previous = MDS.all)
gnmds1 <- MDS.all$points[,1]
gnmds2 <- MDS.all$points[,2]
mds.df<-data.frame(gnmds1,gnmds2)
mds.df
fit<-envfit(mds.df,metadata, na.rm = TRUE,999)
fit

###############################################################
###############################################################

## Example 3 - Extend example 2 to demonstrate grouped and individual
## covariate selection
#dat <- data.frame(sqrt(readNumber),season, CO.scale, shrubs.scale, grass.scale)
#dat
#X <- model.matrix(~ sqrt(readNumber)+ season + CO.scale + shrubs.scale + grass.scale, data=dat)
#X
#X <- X[,-1] ## Remove the intercept as boral includes it automatically

X <- data.frame(sqrt(readNumber),season.scale, CO.scale, shrubs.scale, grass.scale)
X
X <- X[,-2]
X
# LVM model
# Set overdispersion prior
set.prior = list(type = c("normal","normal","normal","uniform"),
                 hypparams = c(100, 20, 100, 20))


OTUs.fit.nb <- boral(fun.some, X = X, family = "negative.binomial", 
                     num.lv = 2, n.burnin = 100, n.iteration = 100, n.thin = 1,
                     calc.ics = FALSE, ssvs.index = c(1,2,3,4,5),save.model = TRUE)

summary(OTUs.fit.nb) 

lvsplot(OTUs.fit.nb,biplot = TRUE,ind.spp = NULL) ## Biplot of the latent variables, 
## which can be interpreted in the same manner as an ordination plot.

# pdf(file = "LVM_technical_replicates.pdf", width = 10)
par(mfrow = c(1,1), mar = c(4,4,1,1))

hor.ordicomm = ordiplot(OTUs.fit.nb$lv.median, choices = c(1,2), 
                        type = "none", cex =0.5,
                        display = "sites")

points(hor.ordicomm,"sites", lwd=4,pch=16, cex = 1.5)

#legend(-0.05, 0.08, SummedExpReord$depth, border="white",
#       col=(SummedExpReord$depth_nominal+1)*10,
#       pch = SummedExpReord$depth_nominal+1,
#       cex = 0.9, bty="n", lwd = 2, ncol = 2)
ordisurf(hor.ordicomm, CO.scale,
         add=T, col = "red", lwd = 1.5, cex = 1.5)
ordisurf(hor.ordicomm, grass.scale,
         add=T, col = "grey", lwd = 1.5, cex = 1.3)
ordisurf(hor.ordicomm, shrubs.scale,
         add=T, col = "blue", lwd = 1.5, cex = 1.3)


# extraction of residuals from the model
# residual correlation due to other sources of covariation
res.cors <- get.residual.cor(OTUs.fit.nb, est = "median", prob = 0.95)
res.cors
corrplot(res.cors$sig.cor, title = "Residual correlations", 
         type = "lower", diag = FALSE)
res.cors$sig.cor
res.cors$trace

# extration of correlations between species due to environmental response
env.cors <- get.enviro.cor(OTUs.fit.nb, est = "median", prob = 0.95)
corrplot(env.cors$sig.cor, type = "lower", diag = FALSE,
         title = "Correlations due to covariates",  mar = c(3,0.5,2,1), tl.srt = 45)









data(spider)
y <- spider$abun
y
n <- nrow(y); p <- ncol(y); 
n
p
## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
spider.fit.nb <- boral(y, family = "negative.binomial", num.lv = 2, 
                       row.eff = "fixed", n.burnin = 10, n.iteration = 100, 
                       n.thin = 1, calc.ics = FALSE)

summary(spider.fit.nb)

plot(spider.fit.nb, ask = FALSE, mfrow = c(2,2)) ## Plots used in residual analysis, 
## Used to check if assumptions such an mean-variance relationship 
## are adequately satisfied.

lvsplot(spider.fit.nb) ## Biplot of the latent variables, 
## which can be interpreted in the same manner as an ordination plot.


## Example 2 - model with no latent variables, no site effects, 
## 	and environmental covariates
X <- scale(spider$x)
X
spider.fit.nb <- boral(y, X = X, family = "negative.binomial", 
                       num.lv = 0, n.burnin = 10, n.iteration = 100, n.thin = 1)

summary(spider.fit.nb) 

## The results can be compared with the default example from 
## the manyglm() function in mvabund. Hopefully they are similar =D


## Example 3 - Extend example 2 to demonstrate grouped and individual
## covariate selection
spider.fit.nb2 <- boral(y, X = X, family = "negative.binomial", 
                        num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1,
                        calc.ics = FALSE, ssvs.index = c(-1,-1,-1,0,1,2))
X
summary(spider.fit.nb2) 
lvsplot(spider.fit.nb2)

## Example 3 - model fitted to presence-absence data, no site effects, and
## two latent variables
data(tikus)
y <- tikus$abun
y[y > 0] <- 1
y <- y[1:20,] ## Consider only years 1981 and 1983
y <- y[,apply(y > 0,2,sum) > 2] ## Consider only spp with more than 2 presences

tikus.fit <- boral(y, family = "binomial", num.lv = 2, 
                   n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE)

lvsplot(tikus.fit, biplot = FALSE) 
## A strong location between the two sampling years 


## Example 4 - model fitted to count data, no site effects, and
## two latent variables, plus traits included to explain environmental responses
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
## Include only traits 1, 2, and 5
traits <- as.matrix(cbind(1,antTraits$traits[,c(1,2,5)]))
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)
## Just for fun, the regression coefficients for the second column of X
## will be estimated separately and not regressed against traits.
which.traits[[3]] <- 0

fit.traits <- boral(y, X = X, traits = traits, which.traits = which.traits, 
                    family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, 
                    n.thin = 1, calc.ics = FALSE)

summary(fit.traits)


# LVM ordination of horizons
horizon.comm.ord = boral(fun.some,
                         family = "negative.binomial",
                         prior.control = set.prior, num.lv = 2, 
                         n.burnin = 10,
                         n.iteration = 100, n.thin = 1)

summary(horizon.comm.ord)
plot(horizon.comm.ord, ask = FALSE, mfrow = c(2,2))
lvsplot(horizon.comm.ord)




