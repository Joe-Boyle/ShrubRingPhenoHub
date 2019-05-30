########################################
# Dissertation code                  ###
# Joe Boyle                          ###
# April 2017                         ###
# Adapted from Sandra Angers-Blondin ### 
########################################

# Set working directory
setwd("/Users/joeboyle/Desktop/Uni/Blissertation")

# Libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(lme4)
library(nlme)
library(effects) 
library(MuMIn)
library(tidyr)
library(dplR)
library(broom)
library(wesanderson)
library(clusterSim)
library(reshape2)
library(stringr)
library(ggmap)
library(rworldmap)
library(praise)

# Functions ----
# Create function which ignores NA
cumsum2 <- function(x){
  x[is.na(x)] <- 0
  cumsum(x)
}

# Create a theme for plotting
theme_JB <- function(){
  theme_classic()+
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16),             
          axis.title.y=element_text(size=16),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=16, vjust=1, hjust=0.5),
          legend.text = element_text(size=10),          
          legend.title = element_text(size=10),                              
          legend.position=c(0.9, 0.9))
}

theme_JBangled <- function(){
  theme_classic()+
    theme(axis.text.x=element_text(size=12, angle = 45, hjust = 1),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16),             
          axis.title.y=element_text(size=16),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=16, vjust=1, hjust=0.5),
          legend.text = element_text(size=10),          
          legend.title = element_text(size=10),                              
          legend.position=c(0.9, 0.9))
}

# Import data -------------------------------------------------------------

# Ring width data
dendro <- read.csv("RingWidthData.csv") #put your own path and change the separator to a comma

# Pith data
pithdata <- read.csv("PithData.csv")

# Microscope calibration data
cal_rad <- read.csv("CalData.csv") # this is the magnification of each radius
cal_table <- read.csv("CalTable.csv") # this is the conversion factor, in pixels/mm

# Climate data
load("climQHI.Rdata")
QikTempStart <- read.csv("Qiktemp.csv")

load("cru_pre_shrub_site_df.RData")
precipdata <- cru.pre.shrub.site.df %>% filter(site_id == "s3")

# Load the pheno data
load("qiki_phen.Rda")

# Age Data
AgeData <- read.csv("AgeData.csv")
AgeData$Plot <- as.factor(AgeData$Plot)

# NDVI data
load("/Users/joeboyle/Downloads/MODIS_qik(1).RData")
qiki_modis <- greenup.qik %>% filter(Plot == "HV") %>%
  dplyr::select(Year, ndvi.date.max)

# Sea ice data
sea_ice_data <- read.csv("seaice.csv")
sea_ice_data <- sea_ice_data %>% na.omit() %>% filter (year < 2016 & year > 1990)

# determine max and 85% of sea ice extent
max_extent <- max(sea_ice_data$sea_ice_extent, na.rm = T)
min.extent <- min(sea_ice_data$sea_ice_extent, na.rm = T)
extent_85 <- max_extent * 0.85

# determine earliest date when annual minimum is reached
seaice <- sea_ice_data %>% group_by(year) %>% summarise(min.extent = min(sea_ice_extent, na.rm = T))
seaice$min.doy <- sapply(seaice$year, function(year){
  min.doy <- min(sea_ice_data[sea_ice_data$year == year &
                                sea_ice_data$sea_ice_extent == seaice[
                                  seaice$year == year,]$min.extent
                              ,]$doy, na.rm = T)
})

# find last day at which the 85 extend is esceeded prior maximum is reached and add 1
seaice$onset.melt <- sapply(seaice$year, function(year){
  onset.melt <- max(sea_ice_data[sea_ice_data$year == year &
                                   sea_ice_data$doy < seaice[seaice$year == year,]$min.doy &
                                   sea_ice_data$sea_ice_extent > extent_85,]$doy, na.rm = T) + 1
})

colnames(seaice)[1] <- "Year"
seaice <- seaice %>% dplyr::select(-min.doy)


# Data preparation --------------------------------------------------------

# If you did not use the 0.5X magnification lens on the microscope, need to change the calibration:
cal_table <- cal_table %>% mutate(Conversion = Conversion*2)

# Combine the RW and calibration data to get actual measurements - for that the Site, Plot, Indiv, and Radii headers in the datasets need to be spelled exactly the same and all info must be consistent
dendro <- merge(dendro, cal_rad) # this adds the magnification for each radius
dendro <- merge(dendro, cal_table) # this adds the conversion factor
# potentially also merge a dataframe with establishment date for each individual

# Calculate the ring width in um
dendro <- mutate(dendro, rw = rwpix/Conversion * 1000) %>%
  dplyr::select(Plot, Individual, Radius, Year, rw) # ditch the columns we don't need

# Average the ring width at the individual level (mean of all radii)
dendro_av <- dendro %>% group_by(Plot, Individual, Year) %>% summarise(rw = mean(rw))

# Ring Count Data
dendro_av <- dendro_av %>%
  group_by(Plot,Individual) %>%
  mutate(count = seq(by = 1, length.out = n()))

# Sample Depth Plot
dendro_av$Plot <- as.factor(dendro_av$Plot)
ggplot(dendro_av, aes(Year, fill = Plot)) +
  geom_histogram(binwidth = 1) +  
  scale_x_reverse() +
  #geom_vline(xintercept=2000.5, alpha = 0.5, linetype = 3) +
  theme_classic() +
  scale_fill_manual(values = wes_palette("Moonrise3")) +
  labs(y= "Number of individuals") +
  guides(fill=guide_legend(title="Transect"))

# Create a wide rw df  
dendro_wide <- dendro_av %>% ungroup() %>%
  mutate(Individual = paste(Plot, "I", Individual, sep = "")) %>%
  dplyr::select(Individual, Year, rw) %>%
  spread(Individual, rw)
dendro_wide <- as.data.frame(dendro_wide)
row.names(dendro_wide) <- dendro_wide$Year
dendro_wide <- dplyr::select(dendro_wide, -Year)
dendro_rwl <- dendro_wide %>% as.data.frame()

# Create ring area data
pithdata <-  mutate(pithdata, pw = Pithpix/Conversion * 1000) %>%
  mutate(pr = pw/2) %>% # Diameter to radius
  dplyr::select(Plot, Individual, pr)

area_av <- merge(pithdata,dendro_av)

area_wide <- area_av %>%
  ungroup() %>%
  mutate(Individual = paste(Plot,"I",Individual)) %>%
  dplyr::select(Individual, Year, rw) %>%
  spread(Individual, rw)
row.names(area_wide) <- area_wide$Year
area_wide <- dplyr::select(area_wide, -Year)
cs_wide <- cumsum2(area_wide)
cs_wide[cs_wide == 0] <- NA
cs_long <- cs_wide %>%
  gather(na.rm = TRUE) %>%
  separate(key, c("Plot", "Individual"), sep = "I", remove = TRUE, convert = TRUE)

area_cs <- cbind(area_av, cs_long[3]) %>%
  mutate(area = (pi * (pr + value) ^ 2) - (pi * (pr + value - rw) ^ 2)) %>%
  dplyr::select(-pr, -rw, -value)

dendro_av <- merge(dendro_av, area_cs) %>%
  mutate(Individual = paste(Plot, "I", Individual, sep = ""))

# Create a wide area df  
area_wide <- dendro_av %>% ungroup() %>%
  dplyr::select(Individual, Year, area) %>%
  spread(Individual, area)

# Check that all the variables are the class they need to be
dendro_av$Year <- as.factor(dendro_av$Year) # convert Year into a factor
dendro_av$Individual <- as.factor(dendro_av$Individual) # convert Individual into a factor
dendro_av$Plot <- as.factor(dendro_av$Plot) # convert Plot into a factor

# Filter out individuals with fewer than 6 years' total data
##INTERIM SOLUTION ONLY
dendro_av <- filter(dendro_av, !(Year == 2016 & count < 7 | Year == 2015 & count < 6 | 
                                   Year == 2014 & count < 5 | Year == 2013 & count < 4))

# Detrending rw
# recreate wide df
dendro_wide <- dendro_av %>% ungroup() %>%
  dplyr::select(Individual, Year, rw) %>%
  spread(Individual, rw)
dendro_wide <- as.data.frame(dendro_wide)
row.names(dendro_wide) <- dendro_wide$Year
dendro_wide <- dplyr::select(dendro_wide, -Year)
dendro_wide <- dendro_wide %>% as.data.frame()

# Test relationship shape
rwcountplot <- ggplot(dendro_av, aes(count, rw)) + geom_smooth(method = loess)+ theme_classic()
print(rwcountplot)

detrendedrw <- detrend(dendro_rwl, make.plot = TRUE, method = "ModNegExp")
detrendedrwlong <- detrendedrw %>% tibble::rownames_to_column() %>%
  gather(Individual, drw,-rowname) %>% filter(drw < 4)
names(detrendedrwlong)[1] <- "Year"
dendro_av <- merge(detrendedrwlong, dendro_av)

# Detrending area
area_wide <- as.data.frame(area_wide)
row.names(area_wide) <- area_wide$Year
area_wide <- dplyr::select(area_wide, -Year)
area_rwl <- area_wide %>% as.data.frame()

# Test relationship shape
areacountplot <- ggplot(dendro_av, aes(count, area)) + geom_smooth(method = loess) + theme_classic()
print(areacountplot)

detrendedarea <- detrend(area_rwl, make.plot = TRUE, method = "Spline")
detrendedarealong <- detrendedarea %>% tibble::rownames_to_column() %>%
  gather(Individual, darea, -rowname) %>%
  na.omit()
names(detrendedarealong)[1] <- "Year"
dendro_av <- merge(detrendedarealong, dendro_av)

# Filter out first two years' data for each individual
dendro_av <- filter(dendro_av, count > 2)

# Prepare the pheno data
pheno <- qiki_phen %>% filter(SPP == "SALARC") %>%
  group_by(Year) %>% summarise(P1 = mean(P1, na.rm=TRUE), P2 = mean(P2, na.rm=TRUE),
                               P3 = mean(P3, na.rm=TRUE), P4 = mean(P4, na.rm=TRUE),
                               P5 = mean(P5, na.rm=TRUE), P6 = mean(P6, na.rm=TRUE),
                               P7 = mean(P7, na.rm=TRUE)) %>%
  mutate(Px = P5 - P2, pPx = lag(Px))
pheno[is.na(pheno)] <- NA

#Sorting the daily data into monthly
QikTemp <- dplyr::select(QikTempStart, -X, -doy)
QikTempMonthly <- aggregate(temp ~ Month + Year , QikTemp , mean)
climQHI <- spread(QikTempMonthly, Month, temp)
colnames(climQHI) <- c("Year", "tjan", "tfeb", "tmar", "tapr", "tmay", "tjun", "tjul", "taug", "tsep", "toct", "tnov", "tdec")
climQHI <- climQHI %>% mutate(ptjun = lag(tjun, 1), ptjul = lag(tjul, 1),
                              ptaug = lag(taug, 1), ptsep = lag(tsep, 1),
                              ptoct = lag(toct, 1), ptnov = lag(tnov, 1),
                              ptdec = lag(tdec, 1))

#Climate data collated into seasons
climQHI$tsummer <- (climQHI$tjun + climQHI$tjul)/2
climQHI$tautumn <- (climQHI$taug + climQHI$tsep)/2
climQHI$ptsummer <- (climQHI$ptjun + climQHI$ptjul)/2
climQHI$ptautumn <- (climQHI$ptaug + climQHI$ptsep)/2
climQHI$tspring <- (climQHI$tapr + climQHI$tmay)/2
climQHI$twinter <- (climQHI$ptoct + climQHI$ptnov + climQHI$ptdec + climQHI$tjan +
                      climQHI$tfeb + climQHI$tmar)/6
climyears <- filter(climQHI, Year > 1990)

# Merge the clim data and the ring width data in one dataframe
dendroclim <- merge(dendro_av, climyears, by = "Year", all = TRUE) # this adds the corresponding climate variables to each ring of each year; the by.x and by.y arguments are the column name by which you wish to merge in both datasets. all = TRUE preserves the rows even when there are NAs in one of the datasets. ** You need to make sure that the joining variable (Year in this case) is the same class in both dataframes: here it is a factor.

# Add in the phenology data
dendroclim <- merge(dendroclim, pheno, by = "Year", all = TRUE)

# Add precipitation data
precipdata <- spread(precipdata, month, value) %>% dplyr::select(-site_id, -day, -variable)
names(precipdata) <- c("Year", "pjan", "pfeb", "pmar", "papr", "pmay", "pjun", "pjul", "paug",
                       "psep", "poct", "pnov", "pdec")
precipdata <- precipdata %>% mutate(ppjun = lag(pjun, 1), ppjul = lag(pjul, 1),
                                    ppaug = lag(paug, 1), ppsep = lag(psep, 1),
                                    ppoct = lag(poct, 1), ppnov = lag(pnov, 1),
                                    ppdec = lag(pdec, 1))
precipdata$psummer <- precipdata$pjun + precipdata$pjul
precipdata$pautumn <- precipdata$paug + precipdata$psep
precipdata$ppsummer <- precipdata$ppjun + precipdata$ppjul
precipdata$ppautumn <- precipdata$ppaug + precipdata$ppsep
precipdata$pspring <- precipdata$papr + precipdata$pmay
precipdata$pwinter <- precipdata$ppoct + precipdata$ppnov + precipdata$ppdec + precipdata$pjan +
  precipdata$pfeb + precipdata$pmar
precipyears <- filter(precipdata, Year > 1990)
dendroclim2 <- merge(precipyears, dendroclim, by = "Year", all = TRUE)

# Add modis data
dendroclim2 <- merge(qiki_modis, dendroclim2, by = "Year", all = TRUE)

# Add sea ice data
dendroclim2 <- merge(seaice, dendroclim2, by="Year", all = TRUE) %>%
  filter(Year > 1990)

# Select key variables
dendroclim2 <- dplyr::select(dendroclim2, Year, Plot, Individual, drw, rw, area, darea, count,
                             P2, P5, Px, pPx, ptsummer, ptautumn, twinter, tspring, tsummer, tautumn, ppsummer, ppautumn, 
                             pwinter, pspring, psummer, pautumn, ndvi.date.max, min.extent, onset.melt)
names(dendroclim2)[25] <- "NDVImodis"

## Remove 2016 data (growth not complete) and prior to 2001 (insufficient rw)
dendroclim2 <- filter(dendroclim2, Year > 2000 & Year < 2016)

## Filter out years without full data for final AICs (Do you need all these?!)
dendroclimAIC <- filter(dendroclim2, Year > 2001 & Year < 2016)

# Filter out individuals with fewer than 4 years' data
##INTERIM SOLUTION ONLY
dendroclimAIC <- filter(dendroclimAIC, !(Year == 2012 & count < 4 | Year == 2011 & count < 3 | 
                                           Year == 2010 & count < 2))

# Scaled DF
dendroclimAICvariables1 <- dendroclimAIC %>% dplyr::select(-Year, -Plot, -Individual,
                                                           -drw, -rw, -count,
                                                           -area, -darea)
dendroclimAICvariables1scaled <- data.Normalization(dendroclimAICvariables1, type="n5")
dendroclimAICvariables2 <- dendroclimAIC %>% dplyr::select(Year, Plot, Individual, count, drw,
                                                           rw, area, darea)
dendroclimAICscaled <- cbind(dendroclimAICvariables1scaled, dendroclimAICvariables2)
dendroclimAICscaledlong <- gather(dendroclimAICscaled, Variable, Value, c(1:19))


# Data exploration --------------------------------------------------------
# Plot all individual curves
ggplot(dendroclimAIC, aes(x = count, y = drw, group = Individual)) +
  geom_smooth(stat = "identity") +
  #facet_wrap(~Site) + # use if you want to plot sites or plots separately
  labs(x = "Year", y = "Ring width") +
  theme_JB()

ggplot(dendroclimAIC, aes(x = count, y = darea, group = Individual)) +
  geom_smooth(stat = "identity") +
  #facet_wrap(~Site) + # use if you want to plot sites or plots separately
  labs(x = "Year", y = "Ring area") +
  theme_JB()

dendroclimAIC <- dplyr::select(dendroclimAIC, -count)

# Age Distribution across sites
ggplot(AgeData, aes(RingCount, fill = Plot)) +
  geom_histogram(binwidth = 4, center = 2) +  
  theme_classic() +
  scale_fill_manual(values = wes_palette("Moonrise3"))

# Ring width distribution (raw, log transformed) DOES THIS MATTER?!?

#WRITE A BEAUTIFUL LOOP FOR ALL PLOTS IN A PANEL PLUS A SEXY MATRIX 
#MAKE SURE IT'S TESTING UNIQUE VALUES, NOT SKEWED BY THE LATER YEARS HAVING MORE RINGS
dendroclimscaledlongnorm <- dendroclimAICscaled %>% dplyr::select(-Plot ,-Individual, -rw, -area, 
                                                                  -count) %>%
  gather(Variable, Value, c(1:19, 21:22)) %>%
  filter(Value != "NA") %>% unique(.)

data.normality <- data.frame("Shapiro_p.value"=NA)
var_norm <- unique(dendroclimscaledlongnorm$Variable)
for (i in 1:length(var_norm)) {
  normdata <- dendroclimscaledlongnorm[dendroclimscaledlongnorm$Variable == var_norm[i],]
  r <- diff(range(normdata$Value))
  #  normplot <- ggplot(normdata, aes(Value)) + 
  #    geom_histogram(binwidth = r/8) + theme_classic() + labs(x = var_norm[i])
  #  print(normplot)
  qqnorm <- qplot(sample = Value, data = normdata) + geom_abline(intercept = 0, slope = 1)
  print(qqnorm)
  data.normality[i,1] <- shapiro.test(normdata$Value)[2]
}
data.normality$Shapiro_p.value <- as.numeric(data.normality$Shapiro_p.value)
rownames(data.normality) <- var_norm

qplot(sample = Value, data = normdata) + geom_abline(intercept = 0, slope = 1)

# Crossdating
plot_list <- unique(dendro_av$Plot)
for (i in 1:length(plot_list)) {                                   
  plotdata <- dendro_av[dendro_av$Plot == plot_list[i],]
  plotplot <- ggplot(plotdata, aes(Year, drw, colour = factor(Individual))) +
    geom_line(aes(group = Individual)) + theme_classic() 
  print(plotplot)
  plotplot <- ggplot(plotdata, aes(count, drw, colour = factor(Individual))) +
    geom_line() + theme_classic() 
  print(plotplot)
  
  plotplot <- ggplot(plotdata, aes(Year, darea, colour = factor(Individual))) +
    geom_line(aes(group = Individual)) + theme_classic() 
  print(plotplot)
  plotplot <- ggplot(plotdata, aes(count, darea, colour = factor(Individual))) +
    geom_line() + theme_classic() 
  print(plotplot)
}

### MIXED-MODEL ANALYSES ----
# Decide on random effects
# ggplot(dendroclimAIC, aes(Plot, darea, group = Plot)) + geom_boxplot() + theme_classic()
# boxplot(darea~Plot, dendroclimAIC, range=5)
anova(aov(darea~Plot, dendroclimAIC))
plots <- c(1:5)
randoms <- data.frame("p-value"=NA, "dof"=NA, "F-value"=NA)[numeric(0), ]
for (i in 1:length(plots)) {
  individualdata <- dendroclimAIC[dendroclimAIC$Plot == plots[i],]
  #   individualsplot <- ggplot(individualdata, aes(Individual, darea, group = Individual)) + geom_boxplot() + theme_classic()
  #  individualsplot <- boxplot(darea~Individual, individualdata, range=5)
  #  print(individualsplot)
  randoms[i,1] <- anova(aov(darea~Individual, individualdata))[1,5]
  randoms[i,2] <- anova(aov(darea~Individual, individualdata))[2,1]
  randoms[i,3] <- anova(aov(darea~Individual, individualdata))[1,4]
}
#  ggplot(dendroclimAIC, aes(Year, darea, group = Year)) + geom_boxplot() + theme_classic()
# boxplot(darea~Year, dendroclimAIC, range=5)
anova(aov(darea~Year, dendroclimAIC))

## Overall ----

# Overall linear models
models1scaled <- dendroclimAICscaledlong %>%
  group_by(Variable) %>%
  do(fit = lmer(darea ~ (Value) + (1|Plot) + (1|Individual) + (1|Year), data = ., REML = TRUE)) %>%
  tidy(fit) %>%
  filter(term == "Value")
rownames(models1scaled) <- (c("Minimum sea ice extent", "MODIS NDVI", "Melt Onset Date", "Emergence", "Senescence", "Autumn precipitation", "Previous autumn precipitation",
                              "Previous summer precipitation", "Previous growing season",  "Spring precipitation", 
                              "Summer precipitation", "Previous autumn temperature", "Previous summer temperature",
                              "Winter precipitation", "Growing season",
                              "Autumn temperature", "Spring temperature", "Summer temperature",  "Winter temperature"))
models1scaled <- models1scaled %>% dplyr::select(estimate, std.error, statistic, Variable)

# AIC and Distributions of resids
results <- data.frame("AIC"=NA, "pseudo.R2m"=NA, "pseudo.R2c"=NA, "LR p"=NA, "LR X"=NA)[numeric(0), ]
AIC_var <- unique(dendroclimAICscaledlong$Variable)
var_names <- c("Emergence", "Senescence", "Growing season", "Previous growing season", "Previous summer temperature", 
               "Previous autumn temperature", "Winter temperature", "Spring temperature", "Summer temperature", 
               "Autumn temperature", "Previous summer precipitation", "Previous autumn precipitation", 
               "Winter precipitation", "Spring precipitation", "Summer precipitation", "Autumn precipitation",
               "MODIS NDVI", "Sea ice minimum", "Onset melt doy")
for (i in 1:length(AIC_var)) {
  modeldata <- dendroclimAICscaledlong[dendroclimAICscaledlong$Variable == AIC_var[i],]
  AIC_nullREML <- lmer(darea ~ 1 + (1|Plot) + (1|Individual) + (1|Year), data = modeldata, REML = TRUE)
  AIC_null <- lmer(darea ~ 1 + (1|Plot) + (1|Individual) + (1|Year), data = modeldata, REML = FALSE)
  AIC_modelREML <- lmer(darea ~ Value + (1|Plot) + (1|Individual) + (1|Year), data = modeldata, REML = TRUE)
  AIC_model <- lmer(darea ~ Value + (1|Plot) + (1|Individual) + (1|Year), data = modeldata, REML = FALSE)
  results[i+1,1] <- summary(AIC_model)$AIC[1]
  results[i+1,2] <- r.squaredGLMM(AIC_modelREML)[1]
  results[i+1,3] <- r.squaredGLMM(AIC_modelREML)[2]
  results[i+1,4] <- anova(AIC_model, AIC_null) [2,8]
  results[i+1,5] <- anova(AIC_model, AIC_null) [2,6]
  modelplot <- ggplot(AIC_model, aes(Value, darea)) + geom_point() +  geom_smooth(method=lm) +
    theme_JB() + labs(x = var_names[i])
  diagplot <- ggplot(AIC_model, aes(.fitted, .resid))  + geom_point() + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = 3) +
    theme_JB()
  autocor <- acf(residuals(AIC_modelREML), plot = FALSE)
  acfdf <- with(autocor, data.frame(lag, acf))
  corplot <- ggplot(data = acfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB() # NEED TO ADD THE BLUE DASHED LINES ~0.145 intercept
  grid.arrange(main=textGrob(var_names[i], gp=gpar(fontsize=20)), modelplot, diagplot, corplot)
}
results[1,1] <- summary(AIC_null)$AIC[1]
results[1,2] <- r.squaredGLMM(AIC_nullREML)[1]
results[1,3] <- r.squaredGLMM(AIC_nullREML)[2]
results[1,4] <- "-"
results[1,5] <- "-"
rownames(results) <- c("null", AIC_var)

#Panel figure
P2model <- lmer(darea ~ P2 + (1|Plot) + (1|Individual) + (1|Year), data = dendroclimAIC, REML = FALSE)
  P2plot <- ggplot(P2model, aes(P2, darea)) + geom_point(colour="#7200a3", alpha=0.5) +  
    geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
    labs(x = "Emergence", y = "Relative growth")
P5model <- lmer(darea ~ P5 + (1|Plot) + (1|Individual) + (1|Year), data = dendroclimAIC, REML = FALSE)
  P5plot <- ggplot(P5model, aes(P5, darea)) + geom_point(colour="#7200a3", alpha=0.5) +  
    geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
    labs(x = "Senescence", y = "Relative growth")
Pxmodel <- lmer(darea ~ Px + (1|Plot) + (1|Individual) + (1|Year), data = dendroclimAIC, REML = FALSE)
  Pxplot <- ggplot(Pxmodel, aes(Px, darea)) + geom_point(colour="#7200a3", alpha=0.5) +  
    geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
    labs(x = "GSL", y = "Relative growth")
pPxmodel <- lmer(darea ~ pPx + (1|Plot) + (1|Individual) + (1|Year), data = dendroclimAIC, REML = FALSE)
  pPxplot <- ggplot(pPxmodel, aes(pPx, darea)) + geom_point(colour="#7200a3", alpha=0.5) +  
    geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
    labs(x = "Previous GSL", y = "Relative growth")
grid.arrange(P2plot, P5plot, Pxplot, pPxplot)
  
# Plot overall figure
groupings_order <- c("ice", "ndvi", "ice", rep("pheno", 2), rep("precip", 3), "pheno", 
                     rep("precip", 2), rep("temp", 2), "precip", "pheno",
                     rep("temp", 4))

models1scaled$type <- as.factor(groupings_order)
models1scaled$Variable <- as.factor(models1scaled$Variable)
ggplot(models1scaled, aes(colour = type, fill = type, x = Variable, y = estimate)) + 
  geom_hline(yintercept=0, alpha=0.5, linetype = 3) +
  geom_crossbar(aes(ymin = estimate - std.error, ymax = estimate + std.error, alpha = 0.8)) +
  theme_JBangled() + ylab("Standardised Effect Size") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_x_discrete(limits=c("P2","P5","Px","pPx","ptsummer","ptautumn","twinter","tspring","tsummer",
                            "tautumn","ppsummer","ppautumn","pwinter","pspring","psummer","pautumn",
                            "NDVImodis","min.extent","onset.melt"), 
                   labels =  c("Emergence", "Senescence", "GSL", "Previous GSL", expression(paste("Previous T"[summer])), 
                               expression(paste("Previous T"[autumn])), expression(paste("T"[winter])), expression(paste("T"[spring])), 
                               expression(paste("T"[summer])), expression(paste("T"[autumn])), expression(paste("Previous P"[summer])), 
                               expression(paste("Previous P"[autumn])), expression(paste("P"[winter])), expression(paste("P"[spring])), 
                               expression(paste("P"[summer])), expression(paste("P"[autumn])),
                              "MODIS NDVI", "Minimum sea ice extent", "Sea ice melt onset date")) +
  scale_fill_manual(values = c("#65c0ed","#F2AD00","#7200a3","#00A08A","#ce0000")) +
  scale_colour_manual(values = c("#65c0ed","#F2AD00","#7200a3","#00A08A","#ce0000")) 

#senescence tsummer correlation
ggplot(dendroclimAICscaled, aes(P5, tsummer)) + geom_point()
cor.test(dendroclimAICscaled$P5, dendroclimAICscaled$tsummer)
cor.test(dendroclimAICscaled$ptautumn, dendroclimAICscaled$tsummer)
cor.test(dendroclimAICscaled$ptautumn, dendroclimAICscaled$P5)

### Map ----
map <- ggmap(get_map(location = c(-170, 45, -100, 80), maptype = "terrain")) + 
  geom_rect(aes(xmin = -139.4, xmax = -138.7, ymin = 69.45, ymax = 69.7),
            fill = "transparent", color = "red", size = 0.3) +
  labs(x="Longitude", y="Latitude")

minimap <- ggmap(get_map(location = c(lon = -139.1, lat = 69.55), maptype = "terrain")) +
  scale_x_continuous(limits = c(-139.42, -138.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(69.3, 69.8), expand = c(0, 0)) + 
  geom_rect(aes(xmin = -139.4, xmax = -138.7, ymin = 69.45, ymax = 69.7),
            fill = "transparent", color = "red", size = 2.5) +
  geom_rect(aes(xmin = -139.42, xmax = -138.7, ymin = 69.701, ymax = 69.8), # Positioning
            fill = "white", color = "white", size = 0.5) +
  geom_rect(aes(xmin = -139.42, xmax = -138.7, ymin = 69.3, ymax = 69.45),
            fill = "white", color = "white", size = 0.5) +
  geom_rect(aes(xmin = -139.42, xmax = -139.4, ymin = 69.3, ymax = 69.7),
            fill = "white", color = "white", size = 0.5) +
  geom_rect(aes(xmin = -138.8, xmax = -138.82, ymin = 69.48, ymax = 69.48), # Scale bar
            fill = "transparent", color = "black", size = 0.1) +
  geom_rect(aes(xmin = -138.82, xmax = -138.84, ymin = 69.48, ymax = 69.48),
            fill = "transparent", color = "white", size = 0.5) +
  geom_rect(aes(xmin = -138.84, xmax = -138.86, ymin = 69.48, ymax = 69.48),
            fill = "transparent", color = "black", size = 0.5) +
  geom_rect(aes(xmin = -138.86, xmax = -138.88, ymin = 69.48, ymax = 69.48),
            fill = "transparent", color = "white", size = 0.5) +
  geom_rect(aes(xmin = -138.88, xmax = -138.9, ymin = 69.48, ymax = 69.48),
            fill = "transparent", color = "black", size = 0.5) +
  geom_rect(aes(xmin = -138.9, xmax = -139, ymin = 69.48, ymax = 69.48),
            fill = "transparent", color = "white", size = 0.5) +
  annotate("text", x = -138.8, y = 69.485, label = "10", size = 2) +
  annotate("text", x = -138.9, y = 69.485, label = "0", size = 2) +
  annotate("text", x = -139, y = 69.485, label = "10", size = 2) +
  annotate("text", x = -138.9, y = 69.475, label = "km", size = 2) +
  geom_rect(aes(xmin = -138.9, xmax = -138.9, ymin = 69.57, ymax = 69.57), # Site position
            fill = "transparent", color = "red", size = 1.5) +
  annotate("text", x = -138.81, y = 69.56, label = "Study site", size = 3) +
  annotate("text", x = -139.05, y = 69.43, label = "Expanded section", size = 4) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

grid.arrange(map, minimap, ncol = 2, padding = unit(2.5))
