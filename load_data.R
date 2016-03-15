## This script loads the data for the FIA vs Inventory (PA-Lichens) project.
options(stringsAsFactors=F)

## Define directories
data_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Raw_data/'
derived_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Derived_data/'
git_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/'
working_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/'
fig_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/Figures/'

setwd(working_dir)

## Read in data
# Inventory data come from James Lendemer
inv_lichen_raw = read.csv(paste0(data_dir, 'ALL_LENDEMER_LICHENS_ONLY.csv'))
inv_plots = read.csv(paste0(data_dir, 'ALL_LENDEMER_LICHENS_ONLY_PLOTS.csv'))

# FIA Lichen data were generated using the following script:
# https://github.com/jescoyle/FIA-Lichens/blob/master/download_data/parse_lichen_data.txt
# after downloading from the FIA & FHM Program websites in 2011. (As of 3/14/2016, newer versions were not available online.)
fia_lichen = read.csv(paste0(data_dir, 'FIA_lichens_parsed_2012-05-13.csv'))
fia_lichen = subset(fia_lichen, STATECD==42) # Subset to PA

# Read in FIA plot data from Coyle & Hurlbert 2016
fia_plots = read.csv(paste0(data_dir, 'fia_lichen_master_data_2015-09-19.csv'))

# Add coordinates to fia lichens
fia_lichen = merge(fia_lichen, fia_plots[,c('yrplot.id','LON','LAT')], all.x=T, all.y=F)

# Add taxon info to fia lichens
fia_taxa = read.csv(paste0(data_dir, 'REF_LICHEN_SPECIES.CSV'))
fia_taxa = subset(fia_taxa, is.na(YEAREND)) # Do not use old taxa names
fia_lichen = merge(fia_lichen, fia_taxa[,c('LICH_SPPCD','GENUS','SPECIES')], all.x=T, all.y=F)







