# This script compiles FIA data from Pennsylvania: 1999 - most recent plots
# 1999 plots established under FHM program and follow a different format from FIA program plots
# This script combines these two data tables. Coordinates come from a non-public source used in previous project.
# It also fixes told taxonomic nomeclature based on _______________


options(stringsAsFactors=F)
library(dplyr)


## Define directories
data_dir = 'Data/raw'
derived_dir = 'Data/derived'
code_dir = 'Code'
working_dir = './'

## Load raw data tables
fia_lichen <- read.csv(file.path(data_dir, "PA_LICHEN_LAB.CSV"))
fia_plots <- read.csv(file.path(data_dir, "PA_PLOT.CSV"))
fhm_lichen <- read.csv(file.path(data_dir, "PA_LICHEN_ABUNDANCE_1999.CSV"))
fhm_plots <- read.csv(file.path(data_dir, "PA_PLOT_1999.CSV"))

## Load prior data with coordinates
all_fia_plots <- read.csv(file.path(data_dir, "AllLichenPlots_forJES.csv"))


######### Combine FHM and FIA Data tables #############

# Rename FHM columns to conform to FIA columns
fhm_plots <- rename(fhm_plots, 
                    STATECD = STATE,
                    COUNTYCD = COUNTY,
                    MEASYEAR = YEAR)
fhm_lichen <- rename(fhm_lichen,
                     STATECD = STATE,
                     COUNTYCD = COUNTY,
                     PLOT = P3ID,
                     LICH_SPPCD = LICHEN_SPECIES_CODE)
fhm_lichen$MEASYEAR = 1999 # because we only are using data from this year from FHM                     

# Add coordinates to FHM data
fhm_plots <- left_join(fhm_plots, all_fia_plots[, c("STATECD", "COUNTYCD", "P3ID", "FZ_LAT", "FZ_LON", "FZ_ELEV")])
fhm_plots <- rename(fhm_plots,
       LAT = FZ_LAT,
       LON = FZ_LON,
       ELEV = FZ_ELEV, 
       PLOT = P3ID)

# Add plot id to each table
fhm_plots$yrplotid <- make.yrplotid(fhm_plots)
fhm_lichen$yrplotid <- make.yrplotid(fhm_lichen)
fia_plots$yrplotid <- make.yrplotid(fia_plots)
fia_lichen$yrplotid <- make.yrplotid(fia_lichen)

# Combine data from FHM and FIA
keep_lichen_cols <- c("yrplotid", "MEASYEAR", "STATECD", "COUNTYCD", "PLOT", "LICH_SPPCD", "ABUNDANCE_CLASS")
lichen <- bind_rows(fhm_lichen[, keep_lichen_cols], fia_lichen[, keep_lichen_cols])

keep_plot_cols <- c("yrplotid", "MEASYEAR", "STATECD", "COUNTYCD", "PLOT", "LAT", "LON", "ELEV", "QA_STATUS")
plots <- bind_rows(fhm_plots[, keep_plot_cols], fia_plots[, keep_plot_cols])


# Match lichen data to all fia plot data to determine whether any missing
lichen_plots <- unique(lichen[, c("yrplotid", "STATECD")])
lichen_plots <- left_join(lichen_plots, unique(plots))

# Check whether plot coordinates missing
missing_plots <- subset(lichen_plots, is.na(LAT)) # none

# Save FIA plot data
write.csv(lichen_plots, file.path(data_dir, "fia_plots.csv"), row.names = FALSE)


######### Correct taxonomic names changes in FIA #############







