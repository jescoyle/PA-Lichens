## This script cleans up records for PA-Lichens project so that they can be analyzed

options(stringsAsFactors=F)
library(stringr)

## Define directories
data_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Raw_data/'
derived_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Derived_data/'
git_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/'
working_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/'

setwd(working_dir)

## Load functions
source(paste0(git_dir, 'project_functions.R'))

## Read in data
# Inventory data come from James Lendemer
inv_lichen = read.csv(paste0(data_dir, 'ALL_LENDEMER_LICHENS_ONLY_PLOTS.csv'))

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

# Create species binomial that is consistant across data sets
inv_lichen$Species = sapply(strsplit(inv_lichen$DarScientificName, ' '), function(x) x[2])
inv_lichen$Binomial = str_trim(paste(inv_lichen$DarGenus, inv_lichen$Species))
inv_lichen[is.na(inv_lichen$Species),'Binomial'] = inv_lichen[is.na(inv_lichen$Species),'DarGenus'] 
fia_lichen$Binomial = str_trim(paste(fia_lichen$GENUS, fia_lichen$SPECIES))


# Add unique plot identfier based on lat/lon to inventory data
unique_locs = unique(inv_lichen[c('DarLongitude','DarLatitude')])
unique_locs$PlotID = 1:nrow(unique_locs)

inv_lichen2 = merge(inv_lichen, unique_locs)



#plot_dists = spDists(as.matrix(unique_locs), longlat=T)
#diag(plot_dists) = NA
#plot_dists = as.dist(plot_dists)
#plot_dists[order(plot_dists)][1:100]
#min(plot_dists, na.rm=T)

# Add habitat category to inventory data

# Make dictionary of words used to describe substrate
words = unlist(strsplit(inv_lichen$HabSubstrate, ' '))
words = gsub('[(),.]', '', words)
words = tolower(words)
write.csv(unique(words), 'habitat_words.csv')

# Read back in manually modified dictionary that codes substrate into: bark, wood, moss, soil, fungi, lichen
words = read.csv('habitat_words.csv')
words = subset(words, substrate!='')

# Find keywords in substrate descriptions
sub_list = sapply(inv_lichen$HabSubstrate, function(x) categorize_substrate(x, words))
num_words = sapply(sub_list, length)
max(num_words)
inv_lichen$substrate1 = sapply(sub_list, function(x) x[1])
inv_lichen$substrate2 = sapply(sub_list, function(x) x[2])
inv_lichen$substrate3 = sapply(sub_list, function(x) x[3])
write.csv(inv_lichen[,c('ecatalogue_key','Binomial','HabSubstrate','HabHabitat','substrate1','substrate2','substrate3')], 'inv_substrates.csv', row.names=F)

# Read back in manually modified substrate descritions (for those matching more than on key word)
substrate = read.csv('inv_substrates_manual.csv')
inv_lichen = merge(inv_lichen, substrate[,c('ecatalogue_key','substrate')])
inv_lichen = inv_lichen[,-which(colnames(inv_lichen) %in% c('substrate1','substrate2','substrate3'))] # Drop unnecessary columns

# Indicate whether inventory lichens were sampled from the bark of a standing tree (as in FIA)
inv_lichen$standing = inv_lichen$substrate=='bark'
inv_lichen[grep('fallen', inv_lichen$HabSubstrate, ignore.case=T),'standing'] = F

# Save derived data sets
write.csv(inv_lichen, paste0(derived_dir, 'inv_lichens.csv'), row.names=F)
write.csv(fia_lichen, paste0(derived_dir, 'fia_lichens.csv'), row.names=F)






