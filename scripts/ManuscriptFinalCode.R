########################################
# Phenology/growth Manuscript code   ###
# Joe Boyle                          ###
# 04/11/2020                         ###
# Adapted from Sandra Angers-Blondin ### 
########################################

# Libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gridGraphics)
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
library(ggeffects)

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

setwd("Macintosh HD/Users/joeboyle/Documents/GitHub/ShrubRingPhenoHub")

# Ring width data
dendro <- read.csv("data/RingWidthData.csv")

# Pith data
pithdata <- read.csv("data/PithData.csv")

# Microscope calibration data
cal_rad <- read.csv("data/CalData.csv") # this is the magnification of each radius
cal_table <- read.csv("data/CalTable.csv") # this is the conversion factor, in pixels/mm

# Climate data
load("data/climQHI.Rdata")
QikTempStart <- read.csv("data/QikTemp.csv")
ERA <- read.csv("data/ERA_Qik.csv")

# Load the pheno data
load("data/qiki_phen.Rda")

# Age Data
AgeData <- read.csv("data/AgeData.csv")
AgeData$Plot <- as.factor(AgeData$Plot)

# NDVI data
load("data/MODIS6_ShrubHub_ITEX.RData")
qiki_modis <- MODIS %>% filter(site_name == "QHI") %>% na.omit() %>% group_by(year) %>% summarise(NDVI = max(NDVI))
colnames(qiki_modis)[1] <- "Year"

# Sea ice data
sea_ice_data <- read.csv("data/SeaIce.csv")
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

# find last day at which the 85 extent is reached and add 1
seaice$onset.melt <- sapply(seaice$year, function(year){
  onset.melt <- max(sea_ice_data[sea_ice_data$year == year &
                                   sea_ice_data$doy < seaice[seaice$year == year,]$min.doy &
                                   sea_ice_data$sea_ice_extent > extent_85,]$doy, na.rm = T) + 1
})

colnames(seaice)[1] <- "Year"
seaice <- seaice %>% dplyr::select(-min.doy)



cal_table <- cal_table %>% mutate(Conversion = Conversion*2)

# Combine the RW and calibration data to get actual measurements
dendro <- merge(dendro, cal_rad)
dendro <- merge(dendro, cal_table)

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
  theme_classic() +
  scale_fill_manual(values = wes_palette("Moonrise3")) +
  labs(y= "Number of individuals") +
  guides(fill=guide_legend(title="Transect"))
ggsave("SampleDepth.pdf", width = 20, height = 20, units = "cm")

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
  mutate(pr = pw/2) %>%
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

dendro_av$Year <- as.factor(dendro_av$Year) 
dendro_av$Individual <- as.factor(dendro_av$Individual) 
dendro_av$Plot <- as.factor(dendro_av$Plot) 

# Filter out individuals with fewer than 6 years' total data
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
DendroClim <- merge(dendro_av, climyears, by = "Year", all = TRUE) # this adds the corresponding climate variables to each ring of each year; the by.x and by.y arguments are the column name by which you wish to merge in both datasets. all = TRUE preserves the rows even when there are NAs in one of the datasets. ** You need to make sure that the joining variable (Year in this case) is the same class in both dataframes: here it is a factor.

# Add in the phenology data
DendroClim <- merge(DendroClim, pheno, by = "Year", all = TRUE)

# Add precipitation data
ERAdata <- spread(ERA,month,value)
names(ERAdata) <- c("Year", "pjan", "pfeb", "pmar", "papr", "pmay", "pjun", "pjul", "paug",
                    "psep", "poct", "pnov", "pdec")
precip <- ERAdata %>% mutate(ppjun = lag(pjun, 1), ppjul = lag(pjul, 1),
                                 ppaug = lag(paug, 1), ppsep = lag(psep, 1),
                                 ppoct = lag(poct, 1), ppnov = lag(pnov, 1),
                                 ppdec = lag(pdec, 1))
precip$psummer <- precip$pjun + precip$pjul
precip$pautumn <- precip$paug + precip$psep
precip$ppsummer <- precip$ppjun + precip$ppjul
precip$ppautumn <- precip$ppaug + precip$ppsep
precip$pspring <- precip$papr + precip$pmay
precip$pwinter <- precip$ppoct + precip$ppnov + precip$ppdec + precip$pjan +
  precip$pfeb + precip$pmar
precipyears <- filter(precip, Year > 1990)
DendroClim2 <- merge(precipyears, DendroClim, by = "Year", all = TRUE)

# Add modis data
DendroClim2 <- merge(qiki_modis, DendroClim2, by = "Year", all = TRUE)

# Add sea ice data
DendroClim2 <- merge(seaice, DendroClim2, by="Year", all = TRUE) %>%
  filter(Year > 1990)

# Select key variables
DendroClim2 <- dplyr::select(DendroClim2, Year, Plot, Individual, drw, rw, area, darea, count,
                             P1, P2, P5, Px, pPx, ptsummer, ptautumn, twinter, tspring, tsummer, tautumn, ppsummer, ppautumn, 
                             pwinter, pspring, psummer, pautumn, NDVI, min.extent, onset.melt)
names(DendroClim2)[26] <- "NDVImodis"

## Remove 2016 data (growth not complete) and prior to 2001 (insufficient rw)
DendroClim2 <- filter(DendroClim2, Year > 2000 & Year < 2016)

## Filter out years without full data for final AICs
DendroClimAIC <- filter(DendroClim2, Year > 2001 & Year < 2016)

# Filter out individuals with fewer than 4 years' data
DendroClimAIC <- filter(DendroClimAIC, !(Year == 2012 & count < 4 | Year == 2011 & count < 3 | 
                                           Year == 2010 & count < 2))

# Scaled DF
DendroClimAICvariables1 <- DendroClimAIC %>% dplyr::select(-Year, -Plot, -Individual,
                                                           -drw, -rw, -count,
                                                           -area, -darea)
DendroClimAICvariables1Sc <- data.Normalization(DendroClimAICvariables1, type="n5")
DendroClimAICvariables2 <- DendroClimAIC %>% dplyr::select(Year, Plot, Individual, count, drw,
                                                           rw, area, darea)
DendroClimAICSc <- cbind(DendroClimAICvariables1Sc, DendroClimAICvariables2)
DendroClimAICScLong <- gather(DendroClimAICSc, Variable, Value, c(1:20))


# Data exploration --------------------------------------------------------
# Plot all individual curves
ggplot(DendroClimAIC, aes(x = count, y = drw, group = Individual)) +
  geom_smooth(stat = "identity") +
  labs(x = "Year", y = "Ring width") +
  theme_JB()

ggplot(DendroClimAIC, aes(x = count, y = darea, group = Individual)) +
  geom_smooth(stat = "identity") +
  labs(x = "Year", y = "Ring area") +
  theme_JB()

DendroClimAIC <- dplyr::select(DendroClimAIC, -count)

# Age Distribution across sites
ggplot(AgeData, aes(RingCount, fill = Plot)) +
  geom_histogram(binwidth = 5, center = 2) +  
  theme_classic() +
  scale_fill_manual(values = wes_palette("Moonrise3")) +
  labs(x = "Years of data", y = "Count") +
  guides(fill=guide_legend(title="Transect"))

# Ring width distribution (raw, log transformed)
DendroClimSclongnorm <- DendroClimAICSc %>% 
  dplyr::select(-Plot ,-Individual, -rw, -area, -count) %>%
  gather(Variable, Value, c(1:20, 22:23)) %>%
  filter(Value != "NA") %>% unique(.)

data.normality <- data.frame("Shapiro_p.value"=NA)
var_norm <- unique(DendroClimSclongnorm$Variable)
for (i in 1:length(var_norm)) {
  normdata <- DendroClimSclongnorm[DendroClimSclongnorm$Variable == var_norm[i],]
  r <- diff(range(normdata$Value))
  normplot <- ggplot(normdata, aes(Value)) + 
    geom_histogram(binwidth = r/8) + theme_classic() + labs(x = var_norm[i])
  print(normplot)
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
  
  plotplot <- ggplot(plotdata, aes(Year, drw, colour = factor(Individual))) +
    geom_line(aes(group = Individual)) + theme_classic() 
  print(plotplot)
  plotplot <- ggplot(plotdata, aes(count, drw, colour = factor(Individual))) +
    geom_line() + theme_classic() 
  print(plotplot)
}

### MIXED-MODEL ANALYSES ----
anova(aov(drw~Plot, DendroClimAIC))
plots <- c(1:5)
randoms <- data.frame("p-value"=NA, "dof"=NA, "F-value"=NA)[numeric(0), ]
for (i in 1:length(plots)) {
  individualdata <- DendroClimAIC[DendroClimAIC$Plot == plots[i],]
  individualsplot <- ggplot(individualdata, aes(Individual, drw, group = Individual)) + geom_boxplot() + theme_classic()
  individualsplot <- boxplot(drw~Individual, individualdata, range=5)
  print(individualsplot)
  randoms[i,1] <- anova(aov(drw~Individual, individualdata))[1,5]
  randoms[i,2] <- anova(aov(drw~Individual, individualdata))[2,1]
  randoms[i,3] <- anova(aov(drw~Individual, individualdata))[1,4]
}
ggplot(DendroClimAIC, aes(Year, drw, group = Year)) + geom_boxplot() + theme_classic()
boxplot(drw~Year, DendroClimAIC, range=5)
anova(aov(drw~Year, DendroClimAIC))

## Overall ----

# Overall linear models
models1Sc <- DendroClimAICScLong %>%
  group_by(Variable) %>%
  do(estimate = unlist(coef(summary(lmer(drw ~ (Value) + (1|Year), data = ., REML = TRUE)))[2][[1]]), std.error = unlist(coef(summary(lmer(drw ~ (Value) + (1|Year), data = ., REML = TRUE)))[4][[1]]))
  rownames(models1Sc) <- (c("Minimum sea ice extent", "MODIS maximum NDVI", "Sea ice melt onset date", "Date snow free", 
                                "Emergence", "Senescence", "Autumn precipitation", "Previous autumn precipitation",
                                "Previous summer precipitation", "Previous growing season",  "Spring precipitation", 
                                "Summer precipitation", "Previous autumn temperature", "Previous summer temperature",
                                "Winter precipitation", "Growing season",
                                "Autumn temperature", "Spring temperature", "Summer temperature",  "Winter temperature"))

#models1Sc <- models1Sc %>% dplyr::select(estimate, std.error)

# AIC and Distributions of resids ----
results <- data.frame("AIC"=NA, "pseudo.R2m"=NA, "pseudo.R2c"=NA, "LR p"=NA, "LR X"=NA)[numeric(0), ]
AIC_var <- unique(DendroClimAICScLong$Variable)
var_names <- c("Date snow free", "Emergence", "Senescence", "Growing season", "Previous growing season", "Previous summer temperature", 
               "Previous autumn temperature", "Winter temperature", "Spring temperature", "Summer temperature", 
               "Autumn temperature", "Previous summer precipitation", "Previous autumn precipitation", 
               "Winter precipitation", "Spring precipitation", "Summer precipitation", "Autumn precipitation",
               "MODIS maximum NDVI", "Sea ice minimum", "Onset melt doy")
for (i in 1:length(AIC_var)) {
  modeldata <- DendroClimAICScLong[DendroClimAICScLong$Variable == AIC_var[i],]
  AIC_nullREML <- lmer(drw ~ 1 + (1|Year), data = modeldata, REML = TRUE)
  AIC_null <- lmer(drw ~ 1 + (1|Year), data = modeldata, REML = FALSE)
  AIC_modelREML <- lmer(drw ~ Value + (1|Year), data = modeldata, REML = TRUE)
  AIC_model <- lmer(drw ~ Value + (1|Year), data = modeldata, REML = FALSE)
  results[i+1,1] <- summary(AIC_model)$AIC[1]
  results[i+1,2] <- r.squaredGLMM(AIC_modelREML)[1]
  results[i+1,3] <- r.squaredGLMM(AIC_modelREML)[2]
  results[i+1,4] <- anova(AIC_model, AIC_null) [2,8]
  results[i+1,5] <- anova(AIC_model, AIC_null) [2,6]
  modelplot <- ggplot(AIC_model, aes(Value, drw)) + geom_point() +  geom_smooth(method=lm) +
    theme_JB() + labs(x = var_names[i])
  diagplot <- ggplot(AIC_model, aes(.fitted, .resid))  + geom_point() + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = 3) +
    theme_JB()
  autocor <- acf(residuals(AIC_modelREML), plot = FALSE)
  acfdf <- with(autocor, data.frame(lag, acf))
  corplot <- ggplot(data = acfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_hline(aes(yintercept = 0.145)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
  grid.arrange(main=textGrob(var_names[i], gp=gpar(fontsize=20)), modelplot, diagplot, corplot)
}
results[1,1] <- summary(AIC_null)$AIC[1]
results[1,2] <- r.squaredGLMM(AIC_nullREML)[1]
results[1,3] <- r.squaredGLMM(AIC_nullREML)[2]
results[1,4] <- "-"
results[1,5] <- "-"
rownames(results) <- c("null", AIC_var)

#autocorrelation figure ----
P1AIC_modelREML <- lmer(drw ~ P1 + (1|Year), data = DendroClimAICSc, REML = TRUE)
P1autocor <- acf(residuals(P1AIC_modelREML), plot = FALSE)
P1acfdf <- with(P1autocor, data.frame(lag, acf))
P1corplot <- ggplot(data = P1acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="Snow melt date", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
P2AIC_modelREML <- lmer(drw ~ P2 + (1|Year), data = DendroClimAICSc, REML = TRUE)
P2autocor <- acf(residuals(P2AIC_modelREML), plot = FALSE)
P2acfdf <- with(P2autocor, data.frame(lag, acf))
P2corplot <- ggplot(data = P2acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="Emergence", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
P5AIC_modelREML <- lmer(drw ~ P5 + (1|Year), data = DendroClimAICSc, REML = TRUE)
P5autocor <- acf(residuals(P5AIC_modelREML), plot = FALSE)
P5acfdf <- with(P5autocor, data.frame(lag, acf))
P5corplot <- ggplot(data = P5acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="Senescence", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
PxAIC_modelREML <- lmer(drw ~ Px + (1|Year), data = DendroClimAICSc, REML = TRUE)
Pxautocor <- acf(residuals(PxAIC_modelREML), plot = FALSE)
Pxacfdf <- with(Pxautocor, data.frame(lag, acf))
Pxcorplot <- ggplot(data = Pxacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="GSL", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pPxAIC_modelREML <- lmer(drw ~ pPx + (1|Year), data = DendroClimAICSc, REML = TRUE)
pPxautocor <- acf(residuals(pPxAIC_modelREML), plot = FALSE)
pPxacfdf <- with(pPxautocor, data.frame(lag, acf))
pPxcorplot <- ggplot(data = pPxacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="Previous GSL", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ptsummerAIC_modelREML <- lmer(drw ~ ptsummer + (1|Year), data = DendroClimAICSc, REML = TRUE)
ptsummerautocor <- acf(residuals(ptsummerAIC_modelREML), plot = FALSE)
ptsummeracfdf <- with(ptsummerautocor, data.frame(lag, acf))
ptsummercorplot <- ggplot(data = ptsummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("Previous T"[summer])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ptautumnAIC_modelREML <- lmer(drw ~ ptautumn + (1|Year), data = DendroClimAICSc, REML = TRUE)
ptautumnautocor <- acf(residuals(ptautumnAIC_modelREML), plot = FALSE)
ptautumnacfdf <- with(ptautumnautocor, data.frame(lag, acf))
ptautumncorplot <- ggplot(data = ptautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("Previous T"[autumn])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
twinterAIC_modelREML <- lmer(drw ~ twinter + (1|Year), data = DendroClimAICSc, REML = TRUE)
twinterautocor <- acf(residuals(twinterAIC_modelREML), plot = FALSE)
twinteracfdf <- with(twinterautocor, data.frame(lag, acf))
twintercorplot <- ggplot(data = twinteracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("T"[winter])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
tspringAIC_modelREML <- lmer(drw ~ tspring + (1|Year), data = DendroClimAICSc, REML = TRUE)
tspringautocor <- acf(residuals(tspringAIC_modelREML), plot = FALSE)
tspringacfdf <- with(tspringautocor, data.frame(lag, acf))
tspringcorplot <- ggplot(data = tspringacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("T"[spring])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
tsummerAIC_modelREML <- lmer(drw ~ tsummer + (1|Year), data = DendroClimAICSc, REML = TRUE)
tsummerautocor <- acf(residuals(tsummerAIC_modelREML), plot = FALSE)
tsummeracfdf <- with(tsummerautocor, data.frame(lag, acf))
tsummercorplot <- ggplot(data = tsummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("T"[summer])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
tautumnAIC_modelREML <- lmer(drw ~ tautumn + (1|Year), data = DendroClimAICSc, REML = TRUE)
tautumnautocor <- acf(residuals(tautumnAIC_modelREML), plot = FALSE)
tautumnacfdf <- with(tautumnautocor, data.frame(lag, acf))
tautumncorplot <- ggplot(data = tautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("T"[autumn])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ppsummerAIC_modelREML <- lmer(drw ~ ppsummer + (1|Year), data = DendroClimAICSc, REML = TRUE)
ppsummerautocor <- acf(residuals(ppsummerAIC_modelREML), plot = FALSE)
ppsummeracfdf <- with(ppsummerautocor, data.frame(lag, acf))
ppsummercorplot <- ggplot(data = ppsummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("Previous P"[summer])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ppautumnAIC_modelREML <- lmer(drw ~ ppautumn + (1|Year), data = DendroClimAICSc, REML = TRUE)
ppautumnautocor <- acf(residuals(ppautumnAIC_modelREML), plot = FALSE)
ppautumnacfdf <- with(ppautumnautocor, data.frame(lag, acf))
ppautumncorplot <- ggplot(data = ppautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("Previous P"[autumn])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pwinterAIC_modelREML <- lmer(drw ~ pwinter + (1|Year), data = DendroClimAICSc, REML = TRUE)
pwinterautocor <- acf(residuals(pwinterAIC_modelREML), plot = FALSE)
pwinteracfdf <- with(pwinterautocor, data.frame(lag, acf))
pwintercorplot <- ggplot(data = pwinteracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("P"[winter])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pspringAIC_modelREML <- lmer(drw ~ pspring + (1|Year), data = DendroClimAICSc, REML = TRUE)
pspringautocor <- acf(residuals(pspringAIC_modelREML), plot = FALSE)
pspringacfdf <- with(pspringautocor, data.frame(lag, acf))
pspringcorplot <- ggplot(data = pspringacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("P"[spring])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
psummerAIC_modelREML <- lmer(drw ~ psummer + (1|Year), data = DendroClimAICSc, REML = TRUE)
psummerautocor <- acf(residuals(psummerAIC_modelREML), plot = FALSE)
psummeracfdf <- with(psummerautocor, data.frame(lag, acf))
psummercorplot <- ggplot(data = psummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("P"[summer])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pautumnAIC_modelREML <- lmer(drw ~ pautumn + (1|Year), data = DendroClimAICSc, REML = TRUE)
pautumnautocor <- acf(residuals(pautumnAIC_modelREML), plot = FALSE)
pautumnacfdf <- with(pautumnautocor, data.frame(lag, acf))
pautumncorplot <- ggplot(data = pautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title=expression(paste("P"[autumn])), x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
NDVImodisAIC_modelREML <- lmer(drw ~ NDVImodis + (1|Year), data = DendroClimAICSc, REML = TRUE)
NDVImodisautocor <- acf(residuals(NDVImodisAIC_modelREML), plot = FALSE)
NDVImodisacfdf <- with(NDVImodisautocor, data.frame(lag, acf))
NDVImodiscorplot <- ggplot(data = NDVImodisacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="MODIS NDVI", x ="lag", y = "acf") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
min.extentAIC_modelREML <- lmer(drw ~ min.extent + (1|Year), data = DendroClimAICSc, REML = TRUE)
min.extentautocor <- acf(residuals(min.extentAIC_modelREML), plot = FALSE)
min.extentacfdf <- with(min.extentautocor, data.frame(lag, acf))
min.extentcorplot <- ggplot(data = min.extentacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="Minimum sea ice extent", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
onset.meltAIC_modelREML <- lmer(drw ~ onset.melt + (1|Year), data = DendroClimAICSc, REML = TRUE)
onset.meltautocor <- acf(residuals(onset.meltAIC_modelREML), plot = FALSE)
onset.meltacfdf <- with(onset.meltautocor, data.frame(lag, acf))
onset.meltcorplot <- ggplot(data = onset.meltacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145), colour = "#000099", linetype = 2) +
  labs(title="Sea ice melt onset date", x ="", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()

overallautocor <- grid.arrange(P2corplot, P5corplot, Pxcorplot, pPxcorplot, ptsummercorplot, ptautumncorplot, twintercorplot, tspringcorplot, tsummercorplot, tautumncorplot, ppsummercorplot, ppautumncorplot, pwintercorplot, pspringcorplot, psummercorplot, pautumncorplot, NDVImodiscorplot, min.extentcorplot, onset.meltcorplot, P1corplot)
ggsave("Autocorrelation.pdf", plot = overallautocor, width = 40, height = 40, units = "cm")

### Panel figure ----
as.data.frame(DendroClimAIC)
P2model <- lmer(drw ~ P2 + (1|Year), data = DendroClimAIC, REML = FALSE)
P5model <- lmer(drw ~ P5 + (1|Year), data = DendroClimAIC, REML = FALSE)
Pxmodel <- lmer(drw ~ Px + (1|Year), data = DendroClimAIC, REML = FALSE)
pPxmodel <- lmer(drw ~ pPx + (1|Year), data = DendroClimAIC, REML = FALSE)

P2predict <- ggpredict(P2model, terms = "P2", type = "re")
P5predict <- ggpredict(P5model, terms = "P5", type = "re")
Pxpredict <- ggpredict(Pxmodel, terms = "Px", type = "re")
pPxpredict <- ggpredict(pPxmodel, terms = "pPx", type = "re")
                               

P2plot <- ggplot(P2model, aes(P2, drw)) + geom_point(colour="#7200a3", alpha=0.5) +  
  geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
  labs(x = "Emergence (doy)", y = "Relative growth")

P5plot <- ggplot(P5model, aes(P5, drw)) + geom_point(colour="#7200a3", alpha=0.5) +  
geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
              labs(x = "Senescence (doy)", y = "Relative growth")
            
Pxplot <- ggplot(Pxmodel, aes(Px, drw)) + geom_point(colour="#7200a3", alpha=0.5) +  
geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
              labs(x = "GSL (days)", y = "Relative growth")
           
pPxplot <- ggplot(pPxmodel, aes(pPx, drw)) + geom_point(colour="#7200a3", alpha=0.5) +  
geom_smooth(colour="#7200a3", fill="#7200a3", method=lm) + theme_JB() + 
              labs(x = "Previous GSL (days)", y = "Relative growth")
            
GridPlot <- grid.arrange(P2plot, P5plot, Pxplot, pPxplot)
ggsave("PhenologyGrowthModels.pdf", width = 20, height = 20, units = "cm", GridPlot)
            
# Plot overall figure ----
groupings_order <- c("ice", "ndvi", "ice", "ice", rep("pheno", 2), rep("precip", 3), "pheno", 
                                 rep("precip", 2), rep("temp", 2), "precip", "pheno",
                                 rep("temp", 4))
            
models1Sc$alpha <- c(rep(0.3, 5), 0.6, rep(0.3, 6), 0.6, rep(0.3, 5), 0.6, 0.3)
models1Sc$type <- as.factor(groupings_order)
models1Sc$Variable <- as.factor(models1Sc$Variable)
models1Sc$estimate <- as.numeric(models1Sc$estimate)
models1Sc$std.error <- as.numeric(models1Sc$std.error)
ggplot(models1Sc, aes(colour = type, fill = type, x = Variable, y = estimate, alpha = alpha)) + 
              geom_hline(yintercept=0, alpha=0.5, linetype = 3) +
              geom_crossbar(aes(ymin = estimate - std.error, ymax = estimate + std.error,)) +
              theme_JBangled() + ylab("Standardised Effect Size") +
              theme(legend.position = "none", axis.title.x = element_blank()) +
              scale_x_discrete(limits=c("P2","P5","Px","pPx","ptsummer","ptautumn","twinter","tspring","tsummer",
                                        "tautumn","ppsummer","ppautumn","pwinter","pspring","psummer","pautumn",
                                        "NDVImodis","min.extent","onset.melt","P1"), 
                               labels =  c("Emergence", "Senescence", "GSL", "Previous GSL", expression(paste("Previous T"[summer])), 
                                           expression(paste("Previous T"[autumn])), expression(paste("T"[winter])), expression(paste("T"[spring])), 
                                           expression(paste("T"[summer])), expression(paste("T"[autumn])), expression(paste("Previous P"[summer])), 
                                           expression(paste("Previous P"[autumn])), expression(paste("P"[winter])), expression(paste("P"[spring])), 
                                           expression(paste("P"[summer])), expression(paste("P"[autumn])),
                                           "MODIS NDVI", "Minimum sea ice extent", "Sea ice melt onset date", "Date snow free")) +
              scale_fill_manual(values = c("#65c0ed","#F2AD00","#7200a3","#00A08A","#ce0000")) +
              scale_colour_manual(values = c("#65c0ed","#F2AD00","#7200a3","#00A08A","#ce0000")) +
              scale_alpha(models1Sc$alpha, range = c(0.4, 0.7))
            
ggsave("AllVariableModels.pdf", width = 20, height = 20, units = "cm")
            
#senescence tsummer correlation ----
ggplot(DendroClimAICSc, aes(P5, tsummer)) + geom_point()
cor.test(DendroClimAICSc$P5, DendroClimAICSc$tsummer)
cor.test(DendroClimAICSc$ptautumn, DendroClimAICSc$tsummer)
cor.test(DendroClimAICSc$ptautumn, DendroClimAICSc$P5)

