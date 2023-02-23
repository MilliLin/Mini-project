# Prepare	work space
rm(list=ls())
setwd("D:/study/R week/mini project")
library(car)
library(ggplot2)
library(plyr)
library(moments)

# Load the data into R
data<-read.table("data.txt", header = TRUE)
str(data)

# Calculate the mean value of 3 repeats
data$Mean_concentration<-rowMeans(data[5:7], na.rm = TRUE)

# See the distribution of data and check the normal distribution
hist(data$Mean_concentration, main="", xlab="Mean corticosterone concentration (pg/mg)", col="grey")

hist(data$Mean_concentration, main="", xlab="Mean corticosterone concentration (pg/mg)", col="grey",
     prob=TRUE, ylim = c(0, 0.09)) 
lines(density(data$Mean_concentration,na.rm=TRUE), # density plot
      lwd = 2)
abline(v = mean(data$Mean_concentration, na.rm = TRUE), col = "red",lwd = 2)
abline(v = mean(data$Mean_concentration, na.rm = TRUE)-sd(data$Mean_concentration, na.rm = TRUE), col = "blue",lwd = 2, lty=5)
abline(v = mean(data$Mean_concentration, na.rm = TRUE)+sd(data$Mean_concentration, na.rm = TRUE), col = "blue",lwd = 2, lty=5)

# Create QQ plot
qqnorm(data$Mean_concentration)
qqline(data$Mean_concentration)

# See the distribution of log data and check the normal distribution
data$log_concentration<- log(data$Mean_concentration)

hist(data$log_concentration, main="", xlab="Mean corticosterone concentration (pg/mg)", col="grey")

hist(data$log_concentration, main="", xlab="Mean corticosterone concentration (pg/mg)", col="grey")

hist(data$log_concentration, main="", xlab="Mean corticosterone concentration (pg/mg)", col="grey",
     prob=TRUE, ylim = c(0, 1.5)) 
lines(density(data$log_concentration,na.rm=TRUE), # density plot
      lwd = 2)
abline(v = mean(data$log_concentration, na.rm = TRUE), col = "red",lwd = 2)
abline(v = mean(data$log_concentration, na.rm = TRUE)-sd(data$log_concentration, na.rm = TRUE), col = "blue",lwd = 2, lty=5)
abline(v = mean(data$log_concentration, na.rm = TRUE)+sd(data$log_concentration, na.rm = TRUE), col = "blue",lwd = 2, lty=5)

# Create histogram
hist(data$log_concentration)

# Create QQ plot
qqnorm(data$log_concentration)
qqline(data$log_concentration)

# Run Shapiro-Wilk test
shapiro.test(data$log_concentration)

# General linear model fitting
model1<-glm(log_concentration~Lateralization+Age+Sex+Lateralization*Age+Lateralization*Sex, data=data)
summary(model1)
par(mfrow=c(2,2))
plot(model1)
hist(model1$residuals)

skewness(model1$residuals)

dev.off()

# Outlying points test
outlierTest(model1)

# Strong impact points
cutoff <- 4/(nrow(data)-length(model1$coefficients)-2) 
plot(model1, which=4, cook.levels=cutoff) 
abline(h=cutoff, lty=2, col="red")

# Violin plot 
ggplot(data, aes(x=Lateralization, y=log_concentration)) + 
  geom_violin(fill="#56B4E9", 
               color="black", notch=TRUE)+ 
  geom_point(position="jitter", color="blue", alpha=.5)+ 
  geom_rug(side="l", color="black")

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

data_summary <- summarySE(data, measurevar="log_concentration", groupvars=c("Lateralization","Sex"))
ggplot(data=data, aes(x=Lateralization, y=log_concentration, fill=Sex)) + 
  geom_split_violin(trim=FALSE,color="white") + #split violin
  geom_point(data = data_summary,aes(x=Lateralization, y=log_concentration),pch=19,position=position_dodge(0.9),size=1.5)+ #Plot the mean as a dot plot
  geom_errorbar(data = data_summary,aes(ymin = log_concentration-ci, ymax=log_concentration+ci), #Error bars indicate 95% confidence intervals
                width=0.1,
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Log_concentration")+xlab("Lateralization") 

data_summary1 <- summarySE(data, measurevar="log_concentration", groupvars=c("Lateralization","Age"))
ggplot(data=data, aes(x=Lateralization, y=log_concentration, fill=Age)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = data_summary1,aes(x=Lateralization, y=log_concentration),pch=19,position=position_dodge(0.9),size=1.5)+ 
  geom_errorbar(data = data_summary1,aes(ymin = log_concentration-ci, ymax=log_concentration+ci),
                width=0.1, 
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Log_concentration")+xlab("Lateralization") 

# ¼ÇµÃÉ¾µô

# Strong impact points
cutoff <- 4/(nrow(data)-length(model1$coefficients)-2) 
plot(model1, which=4, cook.levels=cutoff) 
abline(h=cutoff, lty=2, col="red")

influencePlot(model1, id.method="identify", main="Influence Plot", 
              sub="Circle size is proportional to Cook's distance")

avPlots(model1, ask=FALSE, id.method="identify")

data1<-data[-c(66,19,17,43,89,65,95,91,35,37,20,22,93,34,74,32,84,101,79,88,105,103,15,107,64,40,41,44,23,102,42,1,28,61,62,51),]
model2<-glm(log_concentration~Lateralization+Age+Sex+Lateralization*Age+Lateralization*Sex, data=data1)
summary(model2)

cutoff <- 4/(nrow(data)-length(model2$coefficients)-2) 
plot(model2, which=4, cook.levels=cutoff) 
abline(h=cutoff, lty=2, col="red")

influencePlot(model2, id.method="identify", main="Influence Plot", 
              sub="Circle size is proportional to Cook's distance")

avPlots(model2, ask=FALSE, id.method="identify")