## This script cleans up records for PA-Lichens project so that they can be analyzed

options(stringsAsFactors=F)
library(stringr)
library(dplyr)

## Define directories
data_dir <- 'Data/raw'
derived_dir <- 'Data/derived'
code_dir <- 'Code'
working_dir <- './'

## Load functions
source(file.path(code_dir, 'project_functions.R'))

## Read in data
# Inventory data come from James Lendemer
inv_lichen <- read.csv(file.path(data_dir, 'ALL_LENDEMER_LICHENS_ONLY_PLOTS.csv'))

# FIA Lichen data were generated using the script create_fia_data.R
fia_lichen <- read.csv(file.path(derived_dir, "fia_lichens.csv"))

# Read in FIA plot data generated using the script create_fia_data.R
fia_plots <- read.csv(file.path(data_dir, "PA_lichen_plots.csv"))

# Add coordinates to fia lichens
fia_lichen <- left_join(fia_lichen, fia_plots[,c('yrplotid','LON','LAT')])

# Make table that maps taxon code to names
# Needs to use most recent entries from COMMENTS file and older entries from SPECIES file
fia_taxa <- data.frame(LICH_SPPCD = unique(fia_lichen$LICH_SPPCD))

ref_comments <- read.csv(file.path(data_dir, "REF_LICHEN_SPP_COMMENTS_PlainTextFINALUPDATES.csv"))
ref_species <- read.csv(file.path(data_dir, "REF_LICHEN_SPECIES.csv"))
ref_comments <- subset(ref_comments, is.na(YEAREND)) # Do not use old taxa names
ref_species <- subset(ref_species, is.na(YEAREND))

fia_taxa$Binomial <- sapply(fia_taxa$LICH_SPPCD, function(x) {
  if(x %in% ref_comments$LICH_SPPCD){
    this_taxon <- ref_comments[ref_comments$LICH_SPPCD == x, c("GENUS", "SPECIES")]
  } else {
    this_taxon <- ref_species[ref_species$LICH_SPPCD == x, c("GENUS", "SPECIES")]
  }
  
  str_trim(paste(this_taxon$GENUS, this_taxon$SPECIES))
})

# Add taxon names to fia lichen records
fia_lichen <- left_join(fia_lichen, fia_taxa)

# Create species binomial that is consistant across data sets
inv_lichen$Species <- sapply(strsplit(inv_lichen$DarScientificName, ' '), function(x) x[2])
inv_lichen$Binomial <- str_trim(paste(inv_lichen$DarGenus, inv_lichen$Species))
inv_lichen[is.na(inv_lichen$Species),'Binomial'] <- inv_lichen[is.na(inv_lichen$Species),'DarGenus'] 

# Write out table of names to compare to INV
fia_count <- fia_lichen %>%
  group_by(Binomial) %>%
  summarise(FIA = n())
inv_count <- inv_lichen %>%
  filter(Binomial %in% fia_taxa$Binomial) %>%
  group_by(Binomial) %>%
  summarise(INV = n())

compare_taxa <- left_join(fia_count, inv_count) %>%
  mutate(INV = ifelse(is.na(INV), 0, INV)) %>%
  rename(FIA_NAME = Binomial)
  
#write.csv(compare_taxa, file.path(derived_dir, "compare_taxa.csv"), row.names = FALSE)

## Examined taxa that were missing from INV and renamed manually in Code/FIA_TAXONOMY_conversion.csv


# Convert FIA taxa to those used in INV
taxa_mapper <- read.table(file.path(code_dir, "FIA_taxonomy_conversion.txt"),
                          header = TRUE, sep = "\t")
change_taxa = subset(taxa_mapper, FIA_NAME!=INVENTORY_NAME)
for(i in 1:nrow(change_taxa)){
	fia_lichen[fia_lichen$Binomial==change_taxa[i,'FIA_NAME'],'Binomial'] <- change_taxa[i,'INVENTORY_NAME']
}

# Add habitat category to inventory data

# Make dictionary of words used to describe substrate
# Only run this once
#words = unlist(strsplit(inv_lichen$HabSubstrate, ' '))
#words = gsub('[(),.]', '', words)
#words = tolower(words)
#write.csv(unique(words), 'habitat_words.csv')

# Read back in manually modified dictionary that codes substrate into: bark, wood, moss, soil, fungi, lichen
#words = read.csv('habitat_words.csv')
#words = subset(words, substrate!='')

# Find keywords in substrate descriptions
# Only run this once
#sub_list = sapply(inv_lichen$HabSubstrate, function(x) categorize_substrate(x, words))
#num_words = sapply(sub_list, length)
#max(num_words)
#inv_lichen$substrate1 = sapply(sub_list, function(x) x[1])
#inv_lichen$substrate2 = sapply(sub_list, function(x) x[2])
#inv_lichen$substrate3 = sapply(sub_list, function(x) x[3])
#write.csv(inv_lichen[,c('ecatalogue_key','Binomial','HabSubstrate','HabHabitat','substrate1','substrate2','substrate3')], 'inv_substrates.csv', row.names=F)

# Read back in manually modified substrate descriptions (for those matching more than one key word)
substrate = read.csv('inv_substrates_manual.csv')
inv_lichen = merge(inv_lichen, substrate[,c('ecatalogue_key','substrate')])

# Indicate whether inventory lichens were sampled from the bark of a standing tree (as in FIA)
inv_lichen$standing = inv_lichen$substrate=='bark'
inv_lichen[grep('fallen', inv_lichen$HabSubstrate, ignore.case=T),'standing'] = F

# Save derived data sets
#write.csv(inv_lichen, file.path(derived_dir, 'inv_lichens.csv'), row.names=F)
#write.csv(fia_lichen, file.path(derived_dir, 'fia_lichens_matchINV.csv'), row.names=F)

# Read in data
inv_lichen <- read.csv(file.path(derived_dir, 'inv_lichens.csv'))
fia_lichen <- read.csv(file.path(derived_dir, 'fia_lichens_matchINV.csv'))



### Make dataframes with plot-level information and environmental data
library(raster)
library(sp)
library(rgdal)

## INV
inv_plots <- unique(inv_lichen[,c('Site_Number','DarLongitude','DarLatitude','HabHabitat','DarLocality')])
#write.csv(inv_plots, file.path(derived_dir, 'inv_plots.csv'), row.names=F)

# Read back in plots after manually removing duplicates
inv_plots <- read.csv(file.path(derived_dir, 'inv_plots.csv'))
coordinates(inv_plots) <- c('DarLongitude','DarLatitude')
projection(inv_plots) <- c('+proj=longlat')

# FIA
fia_plots <- fia_plots[,c('yrplotid','LON','LAT')]
fia_plots = subset(fia_plots, yrplotid %in% unique(fia_lichen$yrplotid)) #146 plots
coordinates(fia_plots) = c('LON','LAT')
projection(fia_plots) = c('+proj=longlat')

## Environmental Data
# This data is not included in Rproject
env_dir = 'C:/Users/User/Google Drive/Research/GIS DATA/'

## PRISM 30-yr normals, 800m

# precipitation
ppt = raster(file.path(env_dir, 'PRISM/Normals81', 
	'PRISM_ppt_30yr_normal_800mM2_annual_bil', 'PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'))

# max vpd
vpdmax = raster(file.path(env_dir, 'PRISM/Normals81', 
	'PRISM_vpdmax_30yr_normal_800mM2_annual_bil', 'PRISM_vpdmax_30yr_normal_800mM2_annual_bil.bil'))

# max temp
tmax = raster(file.path(env_dir, 'PRISM/Normals81', 
	'PRISM_tmax_30yr_normal_800mM2_annual_bil', 'PRISM_tmax_30yr_normal_800mM2_annual_bil.bil'))

# mean temp
tmean = raster(file.path(env_dir, 'PRISM/Normals81', 
	'PRISM_tmean_30yr_normal_800mM2_annual_bil', 'PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'))

# total (wet+dry estimated) N and S deposition (from NTN)
tot_n = raster(file.path(env_dir, 'NTN', 'total_N_dep_2001-2010.tif'))
tot_s = raster(file.path(env_dir, 'NTN', 'total_S_dep_2001-2010.tif'))

# stack rasters 
prism = stack(ppt, vpdmax, tmax, tmean); names(prism) = c('ppt','vpdmax','tmax','tmean')
projection(prism) = c('+proj=longlat')

# Extract data and append to plots tables
fia_plots = cbind(fia_plots, extract(prism, fia_plots))
fia_plots$tot_n = extract(tot_n, fia_plots)
fia_plots$tot_s = extract(tot_s, fia_plots)

inv_plots = cbind(inv_plots, extract(prism, inv_plots))
inv_plots$tot_n = extract(tot_n, inv_plots)
inv_plots$tot_s = extract(tot_s, inv_plots)

## Save data
#write.csv(fia_plots, file.path(derived_dir, 'fia_plots.csv'), row.names=F)
#write.csv(inv_plots, file.path(derived_dir, 'inv_plots.csv'), row.names=F)

fia_plots <- read.csv(file.path(derived_dir, 'fia_plots.csv'))
inv_plots <- read.csv(file.path(derived_dir, 'inv_plots.csv'))

# Visualize

pdf(file.path(working_dir, 'Figures', 'FIA_vs_INV_environment.pdf'), height=4, width=10)
par(mfrow=c(1,3))
par(mar=c(4,5,1,1))

make_plot(xlim=c(6.5, 12.5), ylim=c(90, 140),
	xlab=expression('Mean annual temp.'~(degree~C)), ylab='Annual precip. (cm)')
points(I(ppt/10)~tmean, data=fia_plots, pch=3, col='red')
points(I(ppt/10)~tmean, data=inv_plots, pch=1, col='blue')

make_plot(xlim=c(12, 18), ylim=c(700, 1300),
	xlab=expression('Max. annual temp.'~(degree~C)), ylab='Max. annual VPD (Pa)',	ylab_loc=3.5)
points(I(100*vpdmax)~tmax, data=fia_plots, pch=3, col='red')
points(I(100*vpdmax)~tmax, data=inv_plots, pch=1, col='blue')

make_plot(xlim=c(95, 180), ylim=c(100, 400),
	xlab=expression('Tot. N deposition '~(kg~ha^-1)), ylab=expression('Tot. S deposition '~(kg~ha^-1)))
points(tot_s~tot_n, data=fia_plots, pch=3, col='red')
points(tot_s~tot_n, data=inv_plots, pch=1, col='blue')

legend('topright', c('INV','FIA'), pch=c(1,3), col=c('blue','red'))
dev.off()


# Distribution of plots across ppt gradient

fia_prec = cut(fia_plots$ppt, seq(900, 1400, 50))
inv_prec = cut(inv_plots$ppt, seq(900, 1400, 50))

sum((table(fia_prec)/length(fia_prec))[8:10])
sum((table(inv_prec)/length(inv_prec))[8:10])




