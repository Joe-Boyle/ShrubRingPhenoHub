# Qikiqtaruk Ecological Monitoring

### Description
_This repository contains code and data necessary to replicate data analysis, figures, and tables in __"Summer temperature – but not growing season length – influences growth of a dwarf willow in coastal Arctic tundra"__._

### Paper reference
Preprint pending

### Dataset DOI
https://doi.org/10.5281/zenodo.2397996

### Authors
Joseph S. Boyle, Sandra Angers-Blondin, Jakob J. Assmann, Isla H. Myers-Smith

### Contact: Joseph S. Boyle joe.boyle123@gmail.com

### Data use guidelines
Data are publicly available using a Creative Commons Attribution 4.0 International copyright (CC BY 4.0). Data are fully public but should be appropriately referenced by citing the paper. Although not mandatory, we additionally suggest that data users contact and collaborate with data contributors if this dataset will contribute a substantial proportion of observations used in a particular paper or analysis.

### Data availability & access
This dataset will be maintained at this GitHub repository (https://github.com/ShrubHub/QikiqtarukHub).

### Citation
Preprint

### Acknowledgements from the manuscript
We would like to extend our sincere gratitude to the Inuvialuit people for the opportunity to visit and conduct research on their land. We thank John Godlee and Eleanor Walker for helping with sample collection, Heather Goodare for proof-reading versions of this manuscript. We thank the Herschel Island-Qikiqtaruk Territorial Park rangers for collecting the phenology measurements and the Aurora Research Institute for logistical support in the field with particular thanks to Richard Gordon, Cameron Eckert and in particular the park rangers Edward McLeod, Sam McLeod and Ricky Joe. We thank the research group of Hugues Lantuit at the Alfred Wegener Institute and the Aurora Research Institute for logistical support. Research permits include Yukon Researcher and Explorer permits (16-48S&E) and Yukon Parks Research permits (RE-Inu-02-16).

Funding for this research was provided by NERC through the ShrubTundra standard grant (NE/M016323/1) and an equipment loan from the NERC Geophysical Equipment Facility (GEF 1063).


### Please use the following acknowledgement statement if you use these data
We thank Joseph S. Boyle, Sandra Angers-Blondin, Jakob J. Assmann, Isla H. Myers-Smith for collecting and compiling these data. We thank the Herschel Island-Qikiqtaruk Territorial Park management, Catherine Kennedy, Dorothy Cooley, and Dr. Jill F. Johnstone for establishing and maintaining the phenology and composition data from Qikiqtaruk. We thank previous rangers including LeeJohn Meyook, Jordan McLeod, Pierre Foisy, Colin Gordon, Jeremy Hansen, Albert Rufus and field assistants including Santeri Lehtonen, William Palmer, Louise Beveridge, Clara Flintrop, John Godlee, Eleanor Walker, Catherine Henry and Anika Trimble. We thank the Inuvialuit People for the opportunity to conduct research on their traditional lands.

### Other links
The phenology data are also included as a part of the 'Long-term phenology data for 47 tundra plant species at 19 high-latitude sites, 1992-2017' dataset:
https://www.polardata.ca/pdcsearch/?doi_id=12722

For more information on the Team Shrub research group see:
https://teamshrub.com/

For more information on this dataset and the perspectives of Qikiqtaruk Park Rangers please see:
https://teamshrub.com/2017/09/28/qikiqtaruk-perspectives-by-ranger-edward-mcleod/

# Data

All data for our analyses can be found in the `data` folder. The folder contains the following files:

```
+ qiki_phen.Rda
# Phenology data

+ RingWidthData.csv
# Ring width data

+ PithData.csv
# Pith data

+ CalData.csv
# Microscope calibration data

+ CalTable.csv
# Conversion factor

+ climQHI.Rdata
# Climate data

+ QikTemp.csv
# Temperature data

+ qiki_phen.Rda
# Pheno data

+ AgeData.csv
# Age Data

+ MODIS6_ShrubHub_ITEX.RData
# NDVI data

+ SeaIce.csv
# Sea ice data

+ ERA_Qik.csv
# Precipitation data

+ map_figure_data/
# Additional data used to build the map figure (Fig. 1)

+ map_figure_data/JB_GPS_salix.csv
# GPS points for transects
```

# Scripts

```
+ manuscript_analysis_script.R
# Performs statistical analyses and plots figures

+ map_figure_script.R
# Generates part of the map figure

+ tidy_GPS_points.R
# Simplifies GPS points and generates bound polygon
```

# Figures

The figures generated in `R` are stored in the `figures` folder.

```
+ map_figure/Fig1_Map
# Figure 1

+ Fig2_GrowthModels_darea.pdf
# Figure 2

+ Fig3_AllVariables_darea.pdf
# Figure 3

+ FigS1_ThinSection.jpg
# Figure S1

+ FigS2_SampleDepth.pdf 
# Figure S2

+ FigS3_Autocorrelation_darea.pdf
# Figure S3

+ FigS4_Correlation.pdf
# Figure S4

+ FigS5_Autocorrelation_drw.pdf
# Figure S5

+ FigS6_GrowthModels_drw.pdf
# Figure S6

+ FigS7_AllVariables_drw.pdf
# Figure S7

+ Histograms/
# Bayesian histograms for all modelled variables, named for each variable
```

# Model outputs

Full model outputs for all statistical analyses are stored in the `outputs` folder:

```
+ darea_AIC_results_table.csv
# Results from frequentist analysis of basal area increments

+ darea_Bayesian_results_table.csv
# Results from Bayesian analysis of basal area increments

+ drw_AIC_Results_table.csv
# Results from frequentist analysis of ring widths

+ drw_Bayesian_Results_table.csv
# Results from Bayesian analysis of ring widths
```

# Requirements

### Software
R version 3.3.3 or greater

### Packages
`brms, broom, clusterSim, corrplot, cowplot, dplR, dplyr, effects, ggeffects, ggplot2, grid, gridExtra, gridGraphics, lme4, MuMIn, nlme, rnaturalearth, rnaturalearthdata, sf, tidybayes, tidyverse, wesanderson`

_Note that for full functionality of the code, the `plyr` package has to be loaded before the `dplyr` package._
