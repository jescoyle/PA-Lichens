## This script summarizes FIA and Inventory data for Pennsylvania so that species composistion and diversity can be compared.

# Load data and settings
source('C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/load_data.R')

# Load packages
library(sp)
library(rgdal)
library(rgeos)


## Convert lichen records to spatial data
coordinates(fia_lichen) = c('LON','LAT')
proj4string(fia_lichen) = CRS('+proj=longlat +ellps=WGS84')
coordinates(inv_plots) = c('DarLongitude','DarLatitude')
proj4string(inv_plots) = CRS('+proj=longlat +ellps=WGS84')

### Create grid across Pennsylvania

# Load ecoregions
eco = readOGR('pa_eco', 'pa')
eco_ll = spTransform(eco, CRS('+proj=longlat +ellps=WGS84')) # Convert to lat/lon

# Make outline of PA
pa_outline = gUnionCascaded(eco)
pa_outline_ll = gUnionCascaded(eco_ll)


# Count number of records in each ecoregion
fia_sptab = over(fia_lichen, eco_ll)
eco_ll$fiarecs = sapply(eco_ll$OBJECTID, function(x) nrow(subset(fia_sptab, OBJECTID==x)))

inv_sptab = over(inv_plots, eco_ll)
eco_ll$invrecs = sapply(eco_ll$OBJECTID, function(x) nrow(subset(inv_sptab, OBJECTID==x)))

# Add points from FIA and INV plots to maps
use_layout = list(list('sp.points', fia_lichen, which=1, col='red', pch=3, cex=.5),
	list('sp.points',inv_plots, which=2, col='red', pch=3, cex=.5))

pdf(paste0(fig_dir,'Num records inv vs fia plots.pdf'), height=6, width=4)
par(lend=1)
spplot(eco_ll, c('fiarecs','invrecs'), col='transparent', sp.layout=use_layout)
dev.off()


