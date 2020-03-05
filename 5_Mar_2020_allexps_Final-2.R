####Load in necessary packages####
library(MuMIn)#for r squared calculations and seeing contribution of random effects to the model
library(plotly)#generating plots
library(orca)#exporting high definition images
library(emmeans)#calculating estimated marginal means
library(lmerTest)#calculating linear mixed effect models
library(lmtest)#calculating linear models

####Experiment 1####
#Looking at the differences in the log10 of clam hiding time between the No Cover Control, Transparent, and Black 
#covers after the 6 hour experimental timescale 

#Set working directory to data location 
setwd("C:/Users/svydr/FBQ Moorea 2020/Clam Photosynthesis/data_Final")

#Read in the data sheet 
exp1 <- read.csv("exp_1_Rcode_data_29FEB2020.csv")

#Reorder the treatments so that "no cover" is the reference level for analysis
exp1$Treatment = factor(exp1$Treatment,levels=c("No Cover","Transparent","Black"))

#Linear mixed model with size, treatment, and treatment order as factors with individual clam as a random effect and the response variable of the log10 of hiding time
m1 <- lmer(log10_after~ Observer + Size_cm + Treatment_Number + Treatment + (1|Clam), data=exp1)
summary(m1) 


#Check the residuals to assess proper model fit
r1 <- residuals(m1)
f1 <- fitted(m1)

#The histogram should approximate normal distribution
hist(r1)

#The qq plot should be approximate the diagnol line
qqnorm(r1)

#Check for heteroskedasticity with plot of residuals versus fitted points
plot(f1,r1)

#Calculate emmeans
e1 <- emmeans(m1, "Treatment")
e1

#Calculate effect size
eff1 <- eff_size(e1, sigma = sigma(m1), edf = 71)
eff1

sigma(m1)

#Make a bar graph with the emmeans data with the desired settings and values from the emmeans data
fon <- list(family = "helvetica", size = 23, color = 'black')
xform <- list(categoryorder = "array",
              categoryarray = c("No Cover",
                                "Transparent",
                                "Black"),
              title = "Treatment")
marker_style <- list(line = list(width = 1,
                                 color = 'rgb(0, 0, 0)'));
yaxis <- list(
  title = "Estimated Marginal Means +/- SE",
  tickvals = list(0, 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8),
  tickformat = (".1f"),
  linecolor = toRGB("black"),
  showgrid = FALSE)
Treatmentexp1SE = c("0.0508","0.0529","0.0528")
x <- c("No Cover", "Transparent", "Black")
y <- c("1.29", "1.33", "1.27")
#data <- data.frame(Treatment, Estimated Marginal Means)
pexp1 <- plot_ly(x=x, y=y, type="bar", color = I("white"),marker = marker_style,
             error_y = list(array = ~Treatmentexp1SE, color='black'))%>%
  layout(autosize=F, yaxis=yaxis)%>%
  layout(xaxis = xform)%>%
  layout(font = fon)
pexp1

#Export graph as high def image
orca(pexp1, "pexp1_.svg")

#Calculate the r^2 of linear mixed model to determine the contribution of random effect to explanation of variance
r.squaredGLMM(m1)


####Experiment 2####
#Looking at the differences in the log10 of clam hiding time between the No Cover Control, Transparent, and Black covers post experiment
#First, days 1 and 2 are analyzed (the 42-46 hr timescale) 
#Then, days 1, 2, and 3 are analyzed together (day 3 was a 104 hr experiment)

#Read the data sheet for day 1 and 2 treatments only)
exp2 <- read.csv("exp_2_day1and2_Rcode_data_29FEB2020.csv")

#Reorder the treatments so they show up no cover -> transparent -> black in models and plots
exp2$Treatment = factor(exp2$Treatment,levels=c("No Cover","Transparent","Black"))

#Linear mixed model using size, treatment order and treatment as variables 
#with a random effect accounting for each individual clam (1|Clam) to see individual effect on the log10 of your hiding time after
m2 <- lmer(log10_after~ Size_cm + Treatment_Number + Treatment + (1|Clam), data=exp2)
summary(m2)

#Calculate the r^2 of linear mixed model to determine the contribution of random effect to explanation of variance
r.squaredGLMM(m2)

#Calculate emmeans
e2 <- emmeans(m2, "Treatment")
e2

#Plot the emmeans and the standard errors
fon <- list(family = "helvetica", size = 24, color = 'black')
xform <- list(categoryorder = "array",
              categoryarray = c("No Cover",
                                "Transparent",
                                "Black"),
              title = "Treatment")
marker_style <- list(line = list(width = 1,
                                 color = 'rgb(0, 0, 0)'));
yaxis <- list(
  title = "Estimated Marginal Means +/- SE",
  tickvals = list(0, 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8),
  tickformat = (".1f"),
  linecolor = toRGB("black"),
  showgrid = FALSE)
Treatment1and2SE = c("0.0715","0.0735","0.0715")
x <- c("No Cover", "Transparent", "Black")
y <- c("1.09", "1.28", "1.58")
#data <- data.frame(Treatment, Estimated Marginal Means)
p1and2 <- plot_ly(x=x, y=y, type="bar", color = I("white"),marker = marker_style,
             error_y = list(array = ~Treatment1and2SE, color='black'))%>%
  layout(autosize=F, yaxis=yaxis)%>%
  layout(xaxis = xform)%>%
  layout(font = fon)
p1and2

#Export image
orca(p1and2, "p1and2.svg")

#Calculate effect size
eff2 <- eff_size(e2, sigma = sigma(m2), edf = 54)
eff2

#Check the residuals to assess proper model fit
r2 <- residuals(m2)
f2 <- fitted(m2)

#The histogram should approximate normal distribution
hist(r2)

#The qq plot should be approximate the diagnol line
qqnorm(r2)

#Check for heteroskedasticity with plot of residuals versus fitted points
plot(f2,r2)

#To test for significance between black and transparent, change ref level to black
exp2$Treatment <- relevel(exp2$Treatment, ref = "Black")

m3 <- lmer(log10_after~ Size_cm + Treatment_Number + Treatment + (1|Clam), data=exp2)
summary(m3)

#Experiment 2 models testing for all treatments, all days
#Read in data of all treatments for all days of experiment 2
exp2_3day <- read.csv("exp_2_alldays_Rcode_data_29FEB2020.csv")

#Reorder the treatments so they show up no cover -> transparent -> black in models and plots
exp2_3day$Treatment = factor(exp2_3day$Treatment,levels=c("No Cover","Transparent","Black"))

m3 <- lmer(log10_after~ Size_cm + Treatment_Number + Treatment + (1|Clam), data=exp2_3day)
summary(m3)


#Calculate the r^2 of linear mixed model to determine the contribution of random effect to explanation of variance
r.squaredGLMM(m3)

#Check the residuals to assess proper model fit
r3 <- residuals(m3)
f3 <- fitted(m3)

#The histogram should approximate normal distribution
hist(r3)

#The qq plot should be approximate the diagnol line
qqnorm(r3)

#Check for heteroskedasticity with plot of residuals versus fitted points
plot(f3,r3)

#test for significance between black and transparent
exp2_3day$Treatment <- relevel(exp2_3day$Treatment, ref = "Black")

m4 <- lmer(log10_after~ Size_cm + Treatment_Number + Treatment + (1|Clam), data=exp2)
summary(m4)

#Calculate emmeans
e3 <- emmeans(m3, "Treatment")
e3

#Calculate effect size
eff3 <- eff_size(e3, sigma = sigma(m3), edf = 75.7)
eff3

#Plot the emmeans and standard errors
fon <- list(family = "helvetica", size = 24, color = 'black')
xform <- list(categoryorder = "array",
              categoryarray = c("No Cover",
                                "Transparent",
                                "Black"),
              title = "Treatment")
marker_style <- list(line = list(width = 1,
                                 color = 'rgb(0, 0, 0)'));
yaxis <- list(
  title = "Estimated Marginal Means +/- SE",
  tickvals = list(0, 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8),
  tickformat = (".1f"),
  linecolor = toRGB("black"),
  showgrid = FALSE)
TreatmentSE = c("0.0585","0.0622","0.0612")
x <- c("No Cover", "Transparent", "Black")
y <- c("1.16", "1.33", "1.59")
#data <- data.frame(Treatment, Estimated Marginal Means)
p <- plot_ly(x=x, y=y, type="bar", color = I("white"),marker = marker_style,
             error_y = list(array = ~TreatmentSE, color='black'))%>%
  layout(autosize=F, yaxis=yaxis)%>%
  layout(xaxis = xform)%>%
  layout(font = fon)
p

#Export image
orca(p, "p1and2and3.svg")



####Experiment 3####
#Looking at the differences in the log10 of clam hiding time between the No Cover Control, Transparent, and Black 
#covers after the 134 hour deprivation experiment

#Read in data for experiment 3
exp3 <- read.csv("exp_3_Rcode_data_29FEB2020.csv")

#Set the order for the ref lvl to be no cover
exp3$Treatment = factor(exp3$Treatment, levels = c("No Cover", "Transparent", "Black"))

#Linear model with log10hiding time after as a function of size, treatment
m5 <- lm(log10_after~ Size_cm + Treatment, data = exp3)
summary(m5)

#Calculate emmeans
e5 <- emmeans(m5, "Treatment")
summary(e5)

#Calculate effect sizes
eff5 <- eff_size(e5, sigma = sigma(m5), edf = 7)
summary(eff5)

#Make emmeans bar graph
fon <- list(family = "helvetica", size = 24, color = 'black')
xform <- list(categoryorder = "array",
              categoryarray = c("No Cover",
                                "Transparent",
                                "Black"),
              title = "Treatment")
marker_style <- list(line = list(width = 1,
                                 color = 'rgb(0, 0, 0)'));
yaxis <- list(
  title = "Estimated Marginal Means +/- SE",
  tickvals = list(0, 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8),
  tickformat = (".1f"),
  linecolor = toRGB("black"),
  showgrid = FALSE)
exp3SE = c("0.156","0.181","0.157")
x <- c("No Cover", "Transparent", "Black")
y <- c("1.11", "1.42", "1.55")

#data <- data.frame(Treatment, Estimated Marginal Means)

pexp3 <- plot_ly(x=x, y=y, type="bar", color = I("white"),marker = marker_style,
                  error_y = list(array = ~exp3SE, color='black'))%>%
  layout(autosize=F, yaxis=yaxis)%>%
  layout(xaxis = xform)%>%
  layout(font = fon)

pexp3

#export image
orca(pexp3, "pexp3.svg")




