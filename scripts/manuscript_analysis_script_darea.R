########################################
# Phenology/growth Manuscript code   ###
# Joe Boyle                          ###
# 04/11/2020                         ###
# Adapted from Sandra Angers-Blondin ###
########################################

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(grid)
library(gridExtra)
library(gridGraphics)
library(lme4)
library(nlme)
library(effects)
library(ggeffects)
library(MuMIn)
library(brms)
library(tidybayes)
library(dplR)
library(clusterSim)
library(corrplot)
library(broom)
library(wesanderson)

# Functions ----
# Create function which ignores NA
cumsum2 <- function(x) {
  x[is.na(x)] <- 0
  cumsum(x)
}

# Create a theme for plotting
theme_JB <- function() {
  theme_classic() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
      plot.title = element_text(
        size = 16,
        vjust = 1,
        hjust = 0.5
      ),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = c(0.9, 0.9)
    )
}

theme_JBangled <- function() {
  theme_classic() +
    theme(
      axis.text.x = element_text(
        size = 12,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
      plot.title = element_text(
        size = 16,
        vjust = 1,
        hjust = 0.5
      ),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = c(0.9, 0.9)
    )
}

# Import data -------------------------------------------------------------

# Ring width data
dendro <- read.csv("data/RingWidthData.csv")

# Pith data
pithdata <- read.csv("data/PithData.csv")

# Microscope calibration data
cal_rad <-
  read.csv("data/CalData.csv") # this is the magnification of each radius
cal_table <-
  read.csv("data/CalTable.csv") # this is the conversion factor, in pixels/mm

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
qiki_modis <-
  MODIS %>% filter(site_name == "QHI") %>% na.omit() %>% group_by(year) %>% summarise(NDVI = max(NDVI))
colnames(qiki_modis)[1] <- "Year"

# Sea ice data
sea_ice_data <- read.csv("data/SeaIce.csv")
sea_ice_data <-
  sea_ice_data %>% na.omit() %>% filter (year < 2016 & year > 1990)

# determine max and 85% of sea ice extent
max_extent <- max(sea_ice_data$sea_ice_extent, na.rm = T)
min.extent <- min(sea_ice_data$sea_ice_extent, na.rm = T)
extent_85 <- max_extent * 0.85

# determine earliest date when annual minimum is reached
seaice <-
  sea_ice_data %>% group_by(year) %>% summarise(min.extent = min(sea_ice_extent, na.rm = T))
seaice$min.doy <- sapply(seaice$year, function(year) {
  min.doy <- min(sea_ice_data[sea_ice_data$year == year &
                                sea_ice_data$sea_ice_extent == seaice[seaice$year == year, ]$min.extent
                              , ]$doy, na.rm = T)
})

# find last day at which the 85 extent is reached and add 1
seaice$onset.melt <- sapply(seaice$year, function(year) {
  onset.melt <- max(sea_ice_data[sea_ice_data$year == year &
                                   sea_ice_data$doy < seaice[seaice$year == year, ]$min.doy &
                                   sea_ice_data$sea_ice_extent > extent_85, ]$doy, na.rm = T) + 1
})

colnames(seaice)[1] <- "Year"
seaice <- seaice %>% dplyr::select(-min.doy)

cal_table <- cal_table %>% mutate(Conversion = Conversion * 2)

# Combine the RW and calibration data to get actual measurements
dendro <- merge(dendro, cal_rad)
dendro <- merge(dendro, cal_table)

# Calculate the ring width in um
dendro <- mutate(dendro, rw = rwpix / Conversion * 1000) %>%
  dplyr::select(Plot, Individual, Radius, Year, rw) # ditch the columns we don't need

# Average the ring width at the individual level (mean of all radii)
dendro_av <-
  dendro %>% group_by(Plot, Individual, Year) %>% summarise(rw = mean(rw))

# Ring Count Data
dendro_av <- dendro_av %>%
  group_by(Plot, Individual) %>%
  mutate(count = seq(by = 1, length.out = n()))

# Sample Depth Plot
dendro_av$Plot <- as.factor(dendro_av$Plot)
SDall <- ggplot(AgeData, aes(RingCount, fill = Plot)) +
  geom_histogram(binwidth = 5, center = 2) +  
  theme_classic() +
  scale_fill_manual(values = wes_palette("Moonrise3")) +
  labs(x = "Years of data", y = "Count") +
  guides(fill=guide_legend(title="Transect"))
SDall
ggsave(
  "figures/SampleDepthAll.pdf",
  width = 20,
  height = 20,
  units = "cm"
)

# Create a wide rw df
dendro_wide <- dendro_av %>% ungroup() %>%
  mutate(Individual = paste(Plot, "I", Individual, sep = "")) %>%
  dplyr::select(Individual, Year, rw) %>%
  spread(Individual, rw)
dendro_wide <- as.data.frame(dendro_wide)
row.names(dendro_wide) <- dendro_wide$Year
dendro_wide <- dplyr::select(dendro_wide,-Year)
dendro_rwl <- dendro_wide %>% as.data.frame()

# Create ring area data
pithdata <-  mutate(pithdata, pw = Pithpix / Conversion * 1000) %>%
  mutate(pr = pw / 2) %>%
  dplyr::select(Plot, Individual, pr)

area_av <- merge(pithdata, dendro_av)

area_wide <- area_av %>%
  ungroup() %>%
  mutate(Individual = paste(Plot, "I", Individual)) %>%
  dplyr::select(Individual, Year, rw) %>%
  spread(Individual, rw)
row.names(area_wide) <- area_wide$Year
area_wide <- dplyr::select(area_wide,-Year)
cs_wide <- cumsum2(area_wide)
cs_wide[cs_wide == 0] <- NA
cs_long <- cs_wide %>%
  gather(na.rm = TRUE) %>%
  separate(
    key,
    c("Plot", "Individual"),
    sep = "I",
    remove = TRUE,
    convert = TRUE
  )

area_cs <- cbind(area_av, cs_long[3]) %>%
  mutate(area = (pi * (pr + value) ^ 2) - (pi * (pr + value - rw) ^ 2)) %>%
  dplyr::select(-pr,-rw,-value)

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
dendro_av <-
  filter(
    dendro_av,
    !(
      Year == 2016 & count < 7 | Year == 2015 & count < 6 |
        Year == 2014 &
        count < 5 | Year == 2013 & count < 4
    )
  )

# Filter out first two years' data for each individual
dendro_av <- filter(dendro_av, count > 2)

# Sample Depth Plot
sampledepth <- data.frame(dendro_av$Year, dendro_av$count, dendro_av$Plot)
sampledepth$dendro_av.Year <- as.numeric(as.character(sampledepth$dendro_av.Year))
sampledepth <- filter(sampledepth, dendro_av.Year > 2001, dendro_av.Year < 2016)
SDused <- ggplot(sampledepth, aes(dendro_av.Year, fill = dendro_av.Plot)) +
  geom_histogram(binwidth = 1) +
  scale_x_reverse() +
  theme_classic() +
  geom_vline(xintercept = 2001.5,
             alpha = 0.5,
             linetype = 3) +
  scale_fill_manual(values = wes_palette("Moonrise3")) +
  labs(x = "Year", y = "Number of individuals") +
  guides(fill = guide_legend(title = "Transect"))
SDused
ggsave(
  "figures/SampleDepthUsed.pdf",
  width = 20,
  height = 20,
  units = "cm"
)

legend <- function(my_ggp) {
  part1 <- ggplot_gtable(ggplot_build(my_ggp))
  part2 <- which(sapply(part1$grobs, function(x) x$name) == "guide-box")
  part3 <- part1$grobs[[part2]]
  return(part3)
}
SDall_nolegend <- SDall + theme(legend.position = "none")
SDused_nolegend <- SDused + theme(legend.position = "none")
SDlegend <- legend(SDall)
SDmulti <- grid.arrange(arrangeGrob(SDall_nolegend, SDused_nolegend, ncol = 1),
             SDlegend, nrow = 1)
ggsave(plot = SDmulti,
  "figures/FigS2_SampleDepth.pdf",
  width = 20,
  height = 20,
  units = "cm"
)

# Detrending rw
# recreate wide df
dendro_wide <- dendro_av %>% ungroup() %>%
  dplyr::select(Individual, Year, rw) %>%
  spread(Individual, rw)
dendro_wide <- as.data.frame(dendro_wide)
row.names(dendro_wide) <- dendro_wide$Year
dendro_wide <- dplyr::select(dendro_wide,-Year)
dendro_wide <- dendro_wide %>% as.data.frame()

# Test relationship shape
rwcountplot <-
  ggplot(dendro_av, aes(count, rw)) + geom_smooth(method = loess) + theme_classic()
print(rwcountplot)

detrendedrw <-
  detrend(dendro_rwl, make.plot = TRUE, method = "ModNegExp")
detrendedrwlong <- detrendedrw %>% tibble::rownames_to_column() %>%
  gather(Individual, drw, -rowname) %>% filter(drw < 4)
names(detrendedrwlong)[1] <- "Year"
dendro_av <- merge(detrendedrwlong, dendro_av)

# Detrending area
area_wide <- as.data.frame(area_wide)
row.names(area_wide) <- area_wide$Year
area_wide <- dplyr::select(area_wide,-Year)
area_rwl <- area_wide %>% as.data.frame()

# Test relationship shape
areacountplot <-
  ggplot(dendro_av, aes(count, area)) + geom_smooth(method = loess) + theme_classic()
print(areacountplot)

detrendedarea <-
  detrend(area_rwl, make.plot = TRUE, method = "Spline")
detrendedarealong <-
  detrendedarea %>% tibble::rownames_to_column() %>%
  gather(Individual, darea,-rowname) %>%
  na.omit()
names(detrendedarealong)[1] <- "Year"
dendro_av <- merge(detrendedarealong, dendro_av)

# Prepare the pheno data
pheno <- qiki_phen %>% filter(SPP == "SALARC") %>%
  group_by(Year) %>% summarise(
    P1 = mean(P1, na.rm = TRUE),
    P2 = mean(P2, na.rm = TRUE),
    P3 = mean(P3, na.rm = TRUE),
    P4 = mean(P4, na.rm = TRUE),
    P5 = mean(P5, na.rm = TRUE),
    P6 = mean(P6, na.rm = TRUE),
    P7 = mean(P7, na.rm = TRUE)
  ) %>%
  mutate(Px = P5 - P2, pPx = lag(Px))
pheno[is.na(pheno)] <- NA

#Sorting the daily data into monthly
QikTemp <- dplyr::select(QikTempStart,-X,-doy)
QikTempMonthly <- aggregate(temp ~ Month + Year , QikTemp , mean)
climQHI <- spread(QikTempMonthly, Month, temp)
colnames(climQHI) <-
  c(
    "Year",
    "tjan",
    "tfeb",
    "tmar",
    "tapr",
    "tmay",
    "tjun",
    "tjul",
    "taug",
    "tsep",
    "toct",
    "tnov",
    "tdec"
  )
climQHI <-
  climQHI %>% mutate(
    ptjun = lag(tjun, 1),
    ptjul = lag(tjul, 1),
    ptaug = lag(taug, 1),
    ptsep = lag(tsep, 1),
    ptoct = lag(toct, 1),
    ptnov = lag(tnov, 1),
    ptdec = lag(tdec, 1)
  )

#Climate data collated into seasons
climQHI$tsummer <- (climQHI$tjun + climQHI$tjul) / 2
climQHI$tautumn <- (climQHI$taug + climQHI$tsep) / 2
climQHI$ptsummer <- (climQHI$ptjun + climQHI$ptjul) / 2
climQHI$ptautumn <- (climQHI$ptaug + climQHI$ptsep) / 2
climQHI$tspring <- (climQHI$tapr + climQHI$tmay) / 2
climQHI$twinter <-
  (
    climQHI$ptoct + climQHI$ptnov + climQHI$ptdec + climQHI$tjan +
      climQHI$tfeb + climQHI$tmar
  ) / 6
climyears <- filter(climQHI, Year > 1990)

# Merge the clim data and the ring width data in one dataframe
DendroClim <-
  merge(dendro_av, climyears, by = "Year", all = TRUE) # this adds the corresponding climate variables to each ring of each year; the by.x and by.y arguments are the column name by which you wish to merge in both datasets. all = TRUE preserves the rows even when there are NAs in one of the datasets. ** You need to make sure that the joining variable (Year in this case) is the same class in both dataframes: here it is a factor.

# Add in the phenology data
DendroClim <- merge(DendroClim, pheno, by = "Year", all = TRUE)

# Add precipitation data
ERAdata <- spread(ERA, month, value)
names(ERAdata) <-
  c(
    "Year",
    "pjan",
    "pfeb",
    "pmar",
    "papr",
    "pmay",
    "pjun",
    "pjul",
    "paug",
    "psep",
    "poct",
    "pnov",
    "pdec"
  )
precip <-
  ERAdata %>% mutate(
    ppjun = lag(pjun, 1),
    ppjul = lag(pjul, 1),
    ppaug = lag(paug, 1),
    ppsep = lag(psep, 1),
    ppoct = lag(poct, 1),
    ppnov = lag(pnov, 1),
    ppdec = lag(pdec, 1)
  )
precip$psummer <- precip$pjun + precip$pjul
precip$pautumn <- precip$paug + precip$psep
precip$ppsummer <- precip$ppjun + precip$ppjul
precip$ppautumn <- precip$ppaug + precip$ppsep
precip$pspring <- precip$papr + precip$pmay
precip$pwinter <-
  precip$ppoct + precip$ppnov + precip$ppdec + precip$pjan +
  precip$pfeb + precip$pmar
precipyears <- filter(precip, Year > 1990)
DendroClim2 <-
  merge(precipyears, DendroClim, by = "Year", all = TRUE)

# Add modis data
DendroClim2 <-
  merge(qiki_modis, DendroClim2, by = "Year", all = TRUE)

# Add sea ice data
DendroClim2 <- merge(seaice, DendroClim2, by = "Year", all = TRUE) %>%
  filter(Year > 1990)

# Select key variables
DendroClim2 <-
  dplyr::select(
    DendroClim2,
    Year,
    Plot,
    Individual,
    drw,
    rw,
    area,
    darea,
    count,
    P1,
    P2,
    P5,
    pPx,
    Px,
    ptsummer,
    ptautumn,
    twinter,
    tspring,
    tsummer,
    tautumn,
    ppsummer,
    ppautumn,
    pwinter,
    pspring,
    psummer,
    pautumn,
    NDVI,
    min.extent,
    onset.melt
  )
names(DendroClim2)[26] <- "NDVImodis"

## Remove 2016 data (growth not complete) and prior to 2001 (insufficient rw)
DendroClim2 <- filter(DendroClim2, Year > 2000 & Year < 2016)

## Filter out years without full data for final AICs
DendroClimAIC <- filter(DendroClim2, Year > 2001 & Year < 2016)

# Filter out individuals with fewer than 4 years' data
DendroClimAIC <-
  filter(DendroClimAIC,
         !(Year == 2012 & count < 4 | Year == 2011 & count < 3 |
             Year == 2010 &
             count < 2))

# Scaled DF
DendroClimAICvariables1 <-
  DendroClimAIC %>% dplyr::select(-Year,-Plot,-Individual,-drw,-rw,-count,-area,-darea)
DendroClimAICvariables1Sc <-
  data.Normalization(DendroClimAICvariables1, type = "n5")
DendroClimAICvariables2 <-
  DendroClimAIC %>% dplyr::select(Year, Plot, Individual, count, drw,
                                  rw, area, darea)
DendroClimAICSc <-
  cbind(DendroClimAICvariables1Sc, DendroClimAICvariables2)
DendroClimAICScLong <-
  gather(DendroClimAICSc, Variable, Value, c(1:20))

# Data exploration ----
# Plot all individual curves
ggplot(DendroClimAIC, aes(x = count, y = area, group = Individual)) +
  geom_smooth(stat = "identity") +
  labs(x = "Year", y = "Basal area increment") +
  theme_JB()

ggplot(DendroClimAIC, aes(x = count, y = darea, group = Individual)) +
  geom_smooth(stat = "identity") +
  labs(x = "Year", y = "Detrended basal area increment") +
  theme_JB()

DendroClimAIC <- dplyr::select(DendroClimAIC,-count)

# Age Distribution across sites
ggplot(AgeData, aes(RingCount, fill = Plot)) +
  geom_histogram(binwidth = 5, center = 2) +
  theme_classic() +
  scale_fill_manual(values = wes_palette("Moonrise3")) +
  labs(x = "Years of data", y = "Count") +
  guides(fill = guide_legend(title = "Transect"))

# Ring width distribution (raw, log transformed)
DendroClimSclongnorm <- DendroClimAICSc %>%
  dplyr::select(-Plot , -Individual,-rw,-area,-count) %>%
  gather(Variable, Value, c(1:20, 22:23)) %>%
  filter(Value != "NA") %>% unique(.)

data.normality <- data.frame("Shapiro_p.value" = NA)
var_norm <- unique(DendroClimSclongnorm$Variable)
for (i in 1:length(var_norm)) {
  normdata <-
    DendroClimSclongnorm[DendroClimSclongnorm$Variable == var_norm[i], ]
  r <- diff(range(normdata$Value))
  normplot <- ggplot(normdata, aes(Value)) +
    geom_histogram(binwidth = r / 8) + theme_classic() + labs(x = var_norm[i])
  print(normplot)
  qqnorm <-
    qplot(sample = Value, data = normdata) + geom_abline(intercept = 0, slope = 1)
  print(qqnorm)
  data.normality[i, 1] <- shapiro.test(normdata$Value)[2]
}
data.normality$Shapiro_p.value <-
  as.numeric(data.normality$Shapiro_p.value)
rownames(data.normality) <- var_norm

qplot(sample = Value, data = normdata) + geom_abline(intercept = 0, slope = 1)

# Crossdating
plot_list <- unique(dendro_av$Plot)
for (i in 1:length(plot_list)) {
  plotdata <- dendro_av[dendro_av$Plot == plot_list[i], ]
  plotplot <-
    ggplot(plotdata, aes(Year, area, colour = factor(Individual))) +
    geom_line(aes(group = Individual)) + theme_classic()
  print(plotplot)
  plotplot <-
    ggplot(plotdata, aes(count, area, colour = factor(Individual))) +
    geom_line() + theme_classic()
  print(plotplot)
  
  plotplot <-
    ggplot(plotdata, aes(Year, area, colour = factor(Individual))) +
    geom_line(aes(group = Individual)) + theme_classic()
  print(plotplot)
  plotplot <-
    ggplot(plotdata, aes(count, area, colour = factor(Individual))) +
    geom_line() + theme_classic()
  print(plotplot)
}

# Mixed-model analyses ----
anova(aov(darea ~ Plot, DendroClimAIC))
plots <- c(1:5)
randoms <-
  data.frame("p-value" = NA,
             "dof" = NA,
             "F-value" = NA)[numeric(0),]
for (i in 1:length(plots)) {
  individualdata <- DendroClimAIC[DendroClimAIC$Plot == plots[i], ]
  individualsplot <-
    ggplot(individualdata, aes(Individual, area, group = Individual)) + geom_boxplot() + theme_classic()
  individualsplot <-
    boxplot(darea ~ Individual, individualdata, range = 5)
  print(individualsplot)
  randoms[i, 1] <- anova(aov(darea ~ Individual, individualdata))[1, 5]
  randoms[i, 2] <- anova(aov(darea ~ Individual, individualdata))[2, 1]
  randoms[i, 3] <- anova(aov(darea ~ Individual, individualdata))[1, 4]
}
ggplot(DendroClimAIC, aes(Year, area, group = Year)) + geom_boxplot() + theme_classic()
boxplot(darea ~ Year, DendroClimAIC, range = 5)
anova(aov(darea ~ Year, DendroClimAIC))

# Overall linear models
models1Sc <- DendroClimAICScLong %>%
  group_by(Variable) %>%
  do(estimate = unlist(coef(summary(
    lmer(darea ~ (Value) + (1 |
                            Year), data = ., REML = TRUE)
  ))[2][[1]]), std.error = unlist(coef(summary(
    lmer(darea ~ (Value) + (1 | Year), data = ., REML = TRUE)
  ))[4][[1]]))

# AIC and Distributions of resids ----
results <-
  data.frame(
    "AIC" = NA,
    "pseudo.R2m" = NA,
    "pseudo.R2c" = NA,
    "LR p" = NA,
    "LR X" = NA
  )[numeric(0),]
AIC_var <- unique(DendroClimAICScLong$Variable)
var_names <-
  c(
    "Date snow free",
    "Leaf emergence",
    "Leaf senescence",
    "Growing season",
    "Previous growing season",
    "Previous summer temperature",
    "Previous autumn temperature",
    "Winter temperature",
    "Spring temperature",
    "Summer temperature",
    "Autumn temperature",
    "Previous summer precipitation",
    "Previous autumn precipitation",
    "Winter precipitation",
    "Spring precipitation",
    "Summer precipitation",
    "Autumn precipitation",
    "MODIS max. NDVI",
    "Sea ice min. extent",
    "Sea ice melt onset"
  )
for (i in 1:length(AIC_var)) {
  modeldata <-
    DendroClimAICScLong[DendroClimAICScLong$Variable == AIC_var[i], ]
  AIC_nullREML <-
    lmer(darea ~ 1 + (1 | Year), data = modeldata, REML = TRUE)
  AIC_null <-
    lmer(darea ~ 1 + (1 | Year), data = modeldata, REML = FALSE)
  AIC_modelREML <-
    lmer(darea ~ Value + (1 | Year), data = modeldata, REML = TRUE)
  AIC_model <-
    lmer(darea ~ Value + (1 | Year), data = modeldata, REML = FALSE)
  results[i + 1, 1] <- summary(AIC_model)$AIC[1]
  results[i + 1, 2] <- r.squaredGLMM(AIC_modelREML)[1]
  results[i + 1, 3] <- r.squaredGLMM(AIC_modelREML)[2]
  results[i + 1, 4] <- anova(AIC_model, AIC_null) [2, 8]
  results[i + 1, 5] <- anova(AIC_model, AIC_null) [2, 6]
  modelplot <-
    ggplot(AIC_model, aes(Value, area)) + geom_point() +  geom_smooth(method =
                                                                       lm) +
    theme_JB() + labs(x = var_names[i])
  diagplot <-
    ggplot(AIC_model, aes(.fitted, .resid))  + geom_point() +
    geom_hline(yintercept = 0,
               alpha = 0.5,
               linetype = 3) +
    theme_JB()
  autocor <- acf(residuals(AIC_modelREML), plot = FALSE)
  acfdf <- with(autocor, data.frame(lag, acf))
  corplot <- ggplot(data = acfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_hline(aes(yintercept = 0.145)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
  grid.arrange(main = textGrob(var_names[i], gp = gpar(fontsize = 12)), modelplot, diagplot, corplot)
}
results[1, 1] <- summary(AIC_null)$AIC[1]
results[1, 2] <- r.squaredGLMM(AIC_nullREML)[1]
results[1, 3] <- r.squaredGLMM(AIC_nullREML)[2]
results[1, 4] <- "-"
results[1, 5] <- "-"
rownames(results) <- c("null", AIC_var)

write.csv(results, file = "outputs/AIC_results_table_darea.csv")

# Figure S3 Autocorrelation ----
P1AIC_modelREML <-
  lmer(darea ~ P1 + (1 | Year), data = DendroClimAICSc, REML = TRUE)
P1autocor <- acf(residuals(P1AIC_modelREML), plot = FALSE)
P1acfdf <- with(P1autocor, data.frame(lag, acf))
P1corplot <-
  ggplot(data = P1acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Snow melt date", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
P2AIC_modelREML <-
  lmer(darea ~ P2 + (1 | Year), data = DendroClimAICSc, REML = TRUE)
P2autocor <- acf(residuals(P2AIC_modelREML), plot = FALSE)
P2acfdf <- with(P2autocor, data.frame(lag, acf))
P2corplot <-
  ggplot(data = P2acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Leaf emergence", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
P5AIC_modelREML <-
  lmer(darea ~ P5 + (1 | Year), data = DendroClimAICSc, REML = TRUE)
P5autocor <- acf(residuals(P5AIC_modelREML), plot = FALSE)
P5acfdf <- with(P5autocor, data.frame(lag, acf))
P5corplot <-
  ggplot(data = P5acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Leaf senescence", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
PxAIC_modelREML <-
  lmer(darea ~ Px + (1 | Year), data = DendroClimAICSc, REML = TRUE)
Pxautocor <- acf(residuals(PxAIC_modelREML), plot = FALSE)
Pxacfdf <- with(Pxautocor, data.frame(lag, acf))
Pxcorplot <-
  ggplot(data = Pxacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "GSL", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pPxAIC_modelREML <-
  lmer(darea ~ pPx + (1 | Year), data = DendroClimAICSc, REML = TRUE)
pPxautocor <- acf(residuals(pPxAIC_modelREML), plot = FALSE)
pPxacfdf <- with(pPxautocor, data.frame(lag, acf))
pPxcorplot <-
  ggplot(data = pPxacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Prev. GSL", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ptsummerAIC_modelREML <-
  lmer(darea ~ ptsummer + (1 |
                           Year), data = DendroClimAICSc, REML = TRUE)
ptsummerautocor <-
  acf(residuals(ptsummerAIC_modelREML), plot = FALSE)
ptsummeracfdf <- with(ptsummerautocor, data.frame(lag, acf))
ptsummercorplot <-
  ggplot(data = ptsummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Prev. summer temp.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ptautumnAIC_modelREML <-
  lmer(darea ~ ptautumn + (1 |
                           Year), data = DendroClimAICSc, REML = TRUE)
ptautumnautocor <-
  acf(residuals(ptautumnAIC_modelREML), plot = FALSE)
ptautumnacfdf <- with(ptautumnautocor, data.frame(lag, acf))
ptautumncorplot <-
  ggplot(data = ptautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Prev. autumn temp.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
twinterAIC_modelREML <-
  lmer(darea ~ twinter + (1 | Year), data = DendroClimAICSc, REML = TRUE)
twinterautocor <- acf(residuals(twinterAIC_modelREML), plot = FALSE)
twinteracfdf <- with(twinterautocor, data.frame(lag, acf))
twintercorplot <-
  ggplot(data = twinteracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Winter temp.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
tspringAIC_modelREML <-
  lmer(darea ~ tspring + (1 | Year), data = DendroClimAICSc, REML = TRUE)
tspringautocor <- acf(residuals(tspringAIC_modelREML), plot = FALSE)
tspringacfdf <- with(tspringautocor, data.frame(lag, acf))
tspringcorplot <-
  ggplot(data = tspringacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Spring temp.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
tsummerAIC_modelREML <-
  lmer(darea ~ tsummer + (1 | Year), data = DendroClimAICSc, REML = TRUE)
tsummerautocor <- acf(residuals(tsummerAIC_modelREML), plot = FALSE)
tsummeracfdf <- with(tsummerautocor, data.frame(lag, acf))
tsummercorplot <-
  ggplot(data = tsummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Summer temp.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
tautumnAIC_modelREML <-
  lmer(darea ~ tautumn + (1 | Year), data = DendroClimAICSc, REML = TRUE)
tautumnautocor <- acf(residuals(tautumnAIC_modelREML), plot = FALSE)
tautumnacfdf <- with(tautumnautocor, data.frame(lag, acf))
tautumncorplot <-
  ggplot(data = tautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Autumn temp.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ppsummerAIC_modelREML <-
  lmer(darea ~ ppsummer + (1 |
                           Year), data = DendroClimAICSc, REML = TRUE)
ppsummerautocor <-
  acf(residuals(ppsummerAIC_modelREML), plot = FALSE)
ppsummeracfdf <- with(ppsummerautocor, data.frame(lag, acf))
ppsummercorplot <-
  ggplot(data = ppsummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Prev. summer precip.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
ppautumnAIC_modelREML <-
  lmer(darea ~ ppautumn + (1 |
                           Year), data = DendroClimAICSc, REML = TRUE)
ppautumnautocor <-
  acf(residuals(ppautumnAIC_modelREML), plot = FALSE)
ppautumnacfdf <- with(ppautumnautocor, data.frame(lag, acf))
ppautumncorplot <-
  ggplot(data = ppautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Prev. autumn precip.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pwinterAIC_modelREML <-
  lmer(darea ~ pwinter + (1 | Year), data = DendroClimAICSc, REML = TRUE)
pwinterautocor <- acf(residuals(pwinterAIC_modelREML), plot = FALSE)
pwinteracfdf <- with(pwinterautocor, data.frame(lag, acf))
pwintercorplot <-
  ggplot(data = pwinteracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Winter precip.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pspringAIC_modelREML <-
  lmer(darea ~ pspring + (1 | Year), data = DendroClimAICSc, REML = TRUE)
pspringautocor <- acf(residuals(pspringAIC_modelREML), plot = FALSE)
pspringacfdf <- with(pspringautocor, data.frame(lag, acf))
pspringcorplot <-
  ggplot(data = pspringacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Spring precip.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
psummerAIC_modelREML <-
  lmer(darea ~ psummer + (1 | Year), data = DendroClimAICSc, REML = TRUE)
psummerautocor <- acf(residuals(psummerAIC_modelREML), plot = FALSE)
psummeracfdf <- with(psummerautocor, data.frame(lag, acf))
psummercorplot <-
  ggplot(data = psummeracfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Summer precip.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
pautumnAIC_modelREML <-
  lmer(darea ~ pautumn + (1 | Year), data = DendroClimAICSc, REML = TRUE)
pautumnautocor <- acf(residuals(pautumnAIC_modelREML), plot = FALSE)
pautumnacfdf <- with(pautumnautocor, data.frame(lag, acf))
pautumncorplot <-
  ggplot(data = pautumnacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Autumn precip.", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
NDVImodisAIC_modelREML <-
  lmer(darea ~ NDVImodis + (1 |
                            Year), data = DendroClimAICSc, REML = TRUE)
NDVImodisautocor <-
  acf(residuals(NDVImodisAIC_modelREML), plot = FALSE)
NDVImodisacfdf <- with(NDVImodisautocor, data.frame(lag, acf))
NDVImodiscorplot <-
  ggplot(data = NDVImodisacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "MODIS max. NDVI", x = "lag", y = "acf") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
min.extentAIC_modelREML <-
  lmer(darea ~ min.extent + (1 |
                             Year), data = DendroClimAICSc, REML = TRUE)
min.extentautocor <-
  acf(residuals(min.extentAIC_modelREML), plot = FALSE)
min.extentacfdf <- with(min.extentautocor, data.frame(lag, acf))
min.extentcorplot <-
  ggplot(data = min.extentacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Min. sea ice extent", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()
onset.meltAIC_modelREML <-
  lmer(darea ~ onset.melt + (1 |
                             Year), data = DendroClimAICSc, REML = TRUE)
onset.meltautocor <-
  acf(residuals(onset.meltAIC_modelREML), plot = FALSE)
onset.meltacfdf <- with(onset.meltautocor, data.frame(lag, acf))
onset.meltcorplot <-
  ggplot(data = onset.meltacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = 0.145),
                                               colour = "#000099",
                                               linetype = 2) +
  labs(title = "Sea ice melt onset date", x = "", y = "") + geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_JB()

overallautocor <-
  grid.arrange(
    P2corplot,
    P5corplot,
    pPxcorplot,
    Pxcorplot,
    ptsummercorplot,
    ptautumncorplot,
    twintercorplot,
    tspringcorplot,
    tsummercorplot,
    tautumncorplot,
    ppsummercorplot,
    ppautumncorplot,
    pwintercorplot,
    pspringcorplot,
    psummercorplot,
    pautumncorplot,
    NDVImodiscorplot,
    min.extentcorplot,
    onset.meltcorplot,
    P1corplot
  )
ggsave(
  "figures/FigS3_Autocorrelation_darea.pdf",
  plot = overallautocor,
  width = 40,
  height = 40,
  units = "cm"
)

# AIC Comparison ----
as.data.frame(DendroClimAIC)
P2model <-
  lmer(drw ~ P2 + (1 | Year), data = DendroClimAIC, REML = FALSE)
P5model <-
  lmer(drw ~ P5 + (1 | Year), data = DendroClimAIC, REML = FALSE)
pPxmodel <-
  lmer(drw ~ pPx + (1 | Year), data = DendroClimAIC, REML = FALSE)
Pxmodel <-
  lmer(drw ~ Px + (1 | Year), data = DendroClimAIC, REML = FALSE)
tsummermodel <-
  lmer(drw ~ tsummer + (1 | Year), data = DendroClimAIC, REML = FALSE)
tautumnmodel <-
  lmer(drw ~ tautumn + (1 | Year), data = DendroClimAIC, REML = FALSE)

# Bayesian analysis ----

# Plot histograms
DendroClimAICScLong <- DendroClimAICScLong %>% mutate(Variable_full = recode(Variable, "P1" = "Snow melt", "P2" = "Leaf emergence", "P5" = "Leaf senescence", "pPx" = "Prev. GSL", "Px" = "GSL", "ptsummer" = "Prev. temp. summer", "ptautumn" = "Prev. temp. autumn", "twinter" = "Winter temp.", "tspring" = "Spring temp.", "tsummer" = "Summer temp.", "tautumn" = "Autumn temp.", "ppsummer" = "Prev. summer precip.", "ppautumn" = "Prev. autumn precip.", "pwinter" = "Winter precip.", "pspring" = "Spring precip.", "psummer" = "Summer precip.", "pautumn" = "Autumn precip.", "NDVImodis" = "MODIS max. NDVI", "min.extent" = "Sea ice min. extent", "onset.melt" = "Sea ice melt onset"), Variable_type = recode(Variable, "P1" = "ice", "P2" = "pheno", "P5" = "pheno", "pPx" = "pheno", "Px" = "pheno", "ptsummer" = "temp", "ptautumn" = "temp", "twinter" = "temp", "tspring" = "temp", "tsummer" = "temp", "tautumn" = "temp", "ppsummer" = "precip", "ppautumn" = "precip", "pwinter" = "precip", "pspring" = "precip", "psummer" = "precip", "pautumn" = "precip", "NDVImodis" = "ndvi", "min.extent" = "ice", "onset.melt" = "ice")) 

models1Sc_Bayes <- DendroClimAICScLong %>% 
  group_by(Variable) %>%
  do(ggsave(ggplot(., aes(Value)) +
              geom_density(fill = "lightgrey") +
              labs(.$Variable_full) +
              theme_JB(), 
            filename = gsub("", "", paste0("figures/histograms/", 
                                           unique(as.character(.$Variable_full)),
                                           ".pdf")), device = "pdf"))

# Overall linear models
models1Sc_Bayes <- DendroClimAICScLong %>%
  group_by(Variable) %>%
  do(estimate = unlist(summary(
    brm(darea ~ (Value) + (1 | Year), data = ., warmup = 1000, iter = 3000, 
        cores = 2, chains = 2)
    #control = list(adapt_delta = 0.95)
  )$fixed[2,])) %>% group_by(Variable) %>%
  unnest_wider(col = c(estimate)) %>% rename("Est_error" = 'Est.Error', "CI_low" = `l-95% CI`, "CI_high" = `u-95% CI`) 
  
models1Sc_Bayes <- models1Sc_Bayes %>% mutate(Variable_full = recode(Variable, "P1" = "Snow melt", "P2" = "Leaf emergence", "P5" = "Leaf senescence", "pPx" = "Prev. GSL", "Px" = "GSL", "ptsummer" = "Prev. temp. summer", "ptautumn" = "Prev. temp. autumn", "twinter" = "Winter temp.", "tspring" = "Spring temp.", "tsummer" = "Summer temp.", "tautumn" = "Autumn temp.", "ppsummer" = "Prev. summer precip.", "ppautumn" = "Prev. autumn precip.", "pwinter" = "Winter precip.", "pspring" = "Spring precip.", "psummer" = "Summer precip.", "pautumn" = "Autumn precip.", "NDVImodis" = "MODIS max. NDVI", "min.extent" = "Sea ice min. extent", "onset.melt" = "Sea ice melt onset"), Variable_type = recode(Variable, "P1" = "ice", "P2" = "pheno", "P5" = "pheno", "pPx" = "pheno", "Px" = "pheno", "ptsummer" = "temp", "ptautumn" = "temp", "twinter" = "temp", "tspring" = "temp", "tsummer" = "temp", "tautumn" = "temp", "ppsummer" = "precip", "ppautumn" = "precip", "pwinter" = "precip", "pspring" = "precip", "psummer" = "precip", "pautumn" = "precip", "NDVImodis" = "ndvi", "min.extent" = "ice", "onset.melt" = "ice")) 

write.csv(models1Sc_Bayes, file = "outputs/Bayesian_results_table_darea.csv")

# Figure 2 Growth Models ----
P2model <- brm(darea ~ P2 + (1 | Year),  
              data = DendroClimAIC, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2)

P5model <- brm(darea ~ P5 + (1 | Year),  
               data = DendroClimAIC, 
               warmup = 1000, iter = 3000, 
               cores = 2, chains = 2)

pPxmodel <- brm(darea ~ pPx + (1 | Year),  
                data = DendroClimAIC, 
                warmup = 1000, iter = 3000, 
                cores = 2, chains = 2)

Pxmodel <- brm(darea ~ Px + (1 | Year),  
               data = DendroClimAIC, 
               warmup = 1000, iter = 3000, 
               cores = 2, chains = 2)

tsummermodel <- brm(darea ~ tsummer + (1 | Year),  
                    data = DendroClimAIC, 
                    warmup = 1000, iter = 3000, 
                    cores = 2, chains = 2)

tautumnmodel <- brm(darea ~ tautumn + (1 | Year),  
                    data = DendroClimAIC, 
                    warmup = 1000, iter = 3000, 
                    cores = 2, chains = 2)

P2predict <- ggpredict(P2model, terms = "P2", type = "re")
P5predict <- ggpredict(P5model, terms = "P5", type = "re")
pPxpredict <- ggpredict(pPxmodel, terms = "pPx", type = "re")
Pxpredict <- ggpredict(Pxmodel, terms = "Px", type = "re")
tsummerpredict <- ggpredict(tsummermodel, terms = "tsummer", type = "re")
tautumnpredict <- ggpredict(tautumnmodel, terms = "tautumn", type = "re")

# Extract the prediction data frame
P2plot_data <- DendroClimAIC %>%
  add_predicted_draws(P2model)

P2plot <- ggplot() +
  stat_lineribbon(data = P2plot_data, aes(x = P2, y = .prediction), .width = c(.95, .80, .50),
                  alpha = 0.2, linetype = 0, fill = "#7200a3") +
  geom_line(data = P2predict, aes(x = x, y = predicted), colour = "#7200a3", size = 1) +
  geom_point(data = DendroClimAIC, aes(x = P2, y = darea), colour = "#7200a3", alpha = 0.5, size = 2) +
  labs(x = "\nLeaf emergence (DOY)", y = "Relative growth\n") +
  theme_JB()

P5plot_data <- DendroClimAIC %>%
  add_predicted_draws(P5model)

P5plot <- ggplot() +
  stat_lineribbon(data = P5plot_data, aes(x = P5, y = .prediction), .width = c(.95, .80, .50),
                  alpha = 0.2, linetype = 0, fill = "#7200a3") +
  geom_line(data = P5predict, aes(x = x, y = predicted), colour = "#7200a3", size = 1) +
  geom_point(data = DendroClimAIC, aes(x = P5, y = darea), colour = "#7200a3", alpha = 0.5, size = 2) +
  labs(x = "\nLeaf senescence (DOY)", y = "Relative growth\n") +
  theme_JB()

pPxplot_data <- DendroClimAIC %>%
  add_predicted_draws(pPxmodel)

pPxplot <- ggplot() +
  stat_lineribbon(data = pPxplot_data, aes(x = pPx, y = .prediction), .width = c(.95, .80, .50),
                  alpha = 0.2, linetype = 0, fill = "#7200a3") +
  geom_line(data = pPxpredict, aes(x = x, y = predicted), colour = "#7200a3", size = 1) +
  geom_point(data = DendroClimAIC, aes(x = pPx, y = darea), colour = "#7200a3", alpha = 0.5, size = 2) +
  labs(x = "\nPrev. Growing Season Length (days)", y = "Relative growth\n") +
  theme_JB()

Pxplot_data <- DendroClimAIC %>%
  add_predicted_draws(Pxmodel)

Pxplot <- ggplot() +
  stat_lineribbon(data = Pxplot_data, aes(x = Px, y = .prediction), .width = c(.95, .80, .50),
                  alpha = 0.2, linetype = 0, fill = "#7200a3") +
  geom_line(data = Pxpredict, aes(x = x, y = predicted), colour = "#7200a3", size = 1) +
  geom_point(data = DendroClimAIC, aes(x = Px, y = darea), colour = "#7200a3", alpha = 0.5, size = 2) +
  labs(x = "\nGrowing Season Length (days)", y = "Relative growth\n") +
  theme_JB()

tsummerplot_data <- DendroClimAIC %>%
  add_predicted_draws(tsummermodel)

tsummerplot <- ggplot() +
  stat_lineribbon(data = tsummerplot_data, aes(x = tsummer, y = .prediction), .width = c(.95, .80, .50),
                  alpha = 0.2, linetype = 0, fill = "#ce0000") +
  geom_line(data = tsummerpredict, aes(x = x, y = predicted), colour = "#ce0000", size = 1) +
  geom_point(data = DendroClimAIC, aes(x = tsummer, y = darea), colour = "#ce0000", alpha = 0.5, size = 2) +
  labs(x = "\nSummer Temperature (˚C)", y = "Relative growth\n") +
  theme_JB()

tautumnplot_data <- DendroClimAIC %>%
  add_predicted_draws(tautumnmodel)

tautumnplot <- ggplot() +
  stat_lineribbon(data = tautumnplot_data, aes(x = tautumn, y = .prediction), .width = c(.95, .80, .50),
                  alpha = 0.2, linetype = 0, fill = "#ce0000") +
  geom_line(data = tautumnpredict, aes(x = x, y = predicted), colour = "#ce0000", size = 1) +
  geom_point(data = DendroClimAIC, aes(x = tautumn, y = darea), colour = "#ce0000", alpha = 0.5, size = 2) +
  labs(x = "\nAutumn Temperature (˚C)", y = "Relative growth\n") +
  theme_JB()

Fig2_GrowthModels <- grid.arrange(P2plot, P5plot, pPxplot, Pxplot, tsummerplot, tautumnplot)
ggsave(plot = Fig2_GrowthModels,
  filename = "figures/Fig2_GrowthModels_darea.pdf",
  width = 25,
  height = 35,
  units = "cm"
)

# Figure 3 All Variables ----
models1Sc_Bayes$Variable_type <- as.factor(models1Sc_Bayes$Variable_type)
models1Sc_Bayes$Variable_full <- as.factor(models1Sc_Bayes$Variable_full)
models1Sc_Bayes$Estimate <- as.numeric(models1Sc_Bayes$Estimate)
models1Sc_Bayes$CI_low <- as.numeric(models1Sc_Bayes$CI_low)
models1Sc_Bayes$CI_high <- as.numeric(models1Sc_Bayes$CI_high)
Fig3_AllVariables <- ggplot(models1Sc_Bayes,
       aes(
         colour = Variable_type,
         fill = Variable_type,
         x = Variable_full,
         y = Estimate,
         alpha = 0.5
       )) +
  geom_hline(yintercept = 0,
             alpha = 0.5,
             linetype = 3) +
  geom_crossbar(aes(ymin = CI_low, ymax = CI_high, )) +
  theme_JBangled() + ylab("Standardised Effect Size\n") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_x_discrete(
    limits = c(
      "Leaf emergence", 
      "Leaf senescence", 
      "Prev. GSL", 
      "GSL", 
      "Prev. temp. summer", 
      "Prev. temp. autumn", 
      "Winter temp.", 
      "Spring temp.", 
      "Summer temp.", 
      "Autumn temp.", 
      "Prev. summer precip.", 
      "Prev. autumn precip.", 
      "Winter precip.", 
      "Spring precip.", 
      "Summer precip.", 
      "Autumn precip.", 
      "MODIS max. NDVI", 
      "Sea ice min. extent", 
      "Sea ice melt onset", 
      "Snow melt"
    ),
    labels =  c(
      "Leaf emergence",
      "Leaf senescence",
      "Prev. GSL",
      "GSL",
      "Prev. summer temp.",
      "Prev. autumn temp.",
      "Winter temp.",
      "Spring temp.",
      "Summer temp.",
      "Autumn temp.",
      "Prev. summer precip.",
      "Prev. autumn precip.",
      "Winter precip.",
      "Spring precip.",
      "Summer precip.",
      "Autumn precip.",
      "MODIS max. NDVI",
      "Min. sea ice extent",
      "Sea ice melt onset",
      "Date snow free"
    )
  ) +
  scale_fill_manual(values = c("#65c0ed", "#F2AD00", "#7200a3", "#00A08A", "#ce0000")) +
  scale_colour_manual(values = c("#65c0ed", "#F2AD00", "#7200a3", "#00A08A", "#ce0000"))

ggsave(plot = Fig3_AllVariables, 
  file = "figures/Fig3_AllVariables_darea.pdf",
  width = 20,
  height = 20,
  units = "cm"
)

# Correlation ----
ggplot(DendroClimAICSc, aes(P5, tsummer)) + geom_point()
cor.test(DendroClimAICSc$P5, DendroClimAICSc$tsummer)
cor.test(DendroClimAICSc$ptautumn, DendroClimAICSc$tsummer)
cor.test(DendroClimAICSc$ptautumn, DendroClimAICSc$P5)

corr_data <- DendroClimAICSc %>% dplyr::select(-Year, -Plot, -Individual, -count, -drw, -rw, -area, -darea) %>% distinct() %>% rename("Snow melt" = P1, "Leaf emergence" = P2, "Leaf senescence" = P5, "Prev. GSL" = pPx, GSL = Px, "Prev. temp. summer" = ptsummer, "Prev. temp. autumn" = ptautumn, "Winter temp." = twinter, "Spring temp." = tspring, "Summer temp." = tsummer, "Autumn temp." = tautumn, "Prev. summer precip." = ppsummer, "Prev. autumn precip." = ppautumn, "Winter precip." = pwinter, "Spring precip." = pspring, "Winter precip." = pwinter, "Summer precip." = psummer, "Autumn precip." = pautumn, "MODIS max. NDVI" = NDVImodis, "Sea ice min. extent" = min.extent, "Sea ice melt onset" = onset.melt)

corr_mat <- cor(corr_data, use = "all.obs")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf(file = "figures/FigS4_Correlation.pdf")

corrplot(corr_mat, method=c("color"), col=col(200),
         type="upper",
         order = 'original',
         addCoef.col = "black",
         tl.cex = 0.8,
         cl.cex = 0.8,
         number.cex = 0.5,
         number.digits = 2,
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         addgrid.col = NA)

dev.off()
