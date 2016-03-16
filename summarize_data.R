## This script summarizes FIA and Inventory data for Pennsylvania so that species composistion and diversity can be compared.

# Load data and settings
source('C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/load_data.R')

# Load packages
library(sp)
library(rgdal)
library(rgeos)

# Read in color ramp
blue2red = read.csv('../blue2red_10colramp.txt')
blue2red = rgb(blue2red, maxColorValue=255)[10:1]

## Convert lichen records to spatial data
coordinates(fia_lichen) = c('LON','LAT')
proj4string(fia_lichen) = CRS('+proj=longlat +ellps=WGS84')
coordinates(inv_lichen) = c('DarLongitude','DarLatitude')
proj4string(inv_lichen) = CRS('+proj=longlat +ellps=WGS84')

### Create grid across Pennsylvania

# Load ecoregions
eco = readOGR('pa_eco', 'pa')
eco_ll = spTransform(eco, CRS('+proj=longlat +datum=NAD83')) # Convert to lat/lon but keep original datum

# Make outline of PA
pa_outline = gUnionCascaded(eco)
pa_outline_ll = gUnionCascaded(eco_ll)


## Plot number of records in each ecoregion
# Count number of records in each ecoregion
fia_sptab = over(fia_lichen, eco_ll)
eco_ll$fiarecs = sapply(eco_ll$OBJECTID, function(x) nrow(subset(fia_sptab, OBJECTID==x)))

inv_sptab = over(inv_lichen, eco_ll)
eco_ll$invrecs = sapply(eco_ll$OBJECTID, function(x) nrow(subset(inv_sptab, OBJECTID==x)))

# Add points from FIA and INV plots to maps
use_layout = list(list('sp.points', fia_lichen, which=1, col='red', pch=3, cex=.5),
	list('sp.points',inv_lichen, which=2, col='red', pch=3, cex=.5))

pdf(paste0(fig_dir,'Num records inv vs fia plots.pdf'), height=6, width=4)
par(lend=1)
spplot(eco_ll, c('fiarecs','invrecs'), col='transparent', sp.layout=use_layout)
dev.off()

## Plot number of records in grid

# Transform to Albers equal area projection
myproj = '+proj=aea +lat_1=35 +lat_2=45 +lat_0=40.8 +lon_0=-77.66 +datum=NAD83 +units=m'
eco_aea = spTransform(eco, CRS(myproj))
pa_outline_aea = gUnionCascaded(eco_aea)

fia_aea = spTransform(fia_lichen, CRS(myproj))
inv_aea = spTransform(inv_lichen,  CRS(myproj))

# Make grid
bounds = bbox(pa_outline_aea)
widths = c(25000,50000,100000)

w = 50000 # in meters

pa_grids = sapply(widths, function(w){
	ncells = ceiling(apply(bounds, 1, diff)/width)
	pa_grid_aea = SpatialGrid(GridTopology(bounds[,'min']+w/2, rep(w,2), ncells), CRS(myproj))
	pa_grid_aea = SpatialGridDataFrame(pa_grid_aea, data.frame(ID=1:prod(ncells)))

	# Count number of records in each grid cell
	fia_sptab = over(fia_aea, pa_grid_aea)$ID
	pa_grid_aea$fiarecs = tabulate(fia_sptab, prod(ncells))
	inv_sptab = over(inv_aea, pa_grid_aea)$ID
	pa_grid_aea$invrecs = tabulate(inv_sptab, prod(ncells))

	# Get species list for each grid cell
	fia_splist = sapply(pa_grid_aea$ID, function(x){
		this_data = subset(fia_aea@data, fia_sptab==x)
		get_species(this_data)
	})
	inv_splist = sapply(pa_grid_aea$ID, function(x){
		this_data = subset(inv_aea@data, inv_sptab==x)
		get_species(this_data)
	})


	# Count number of species in each grid cell
	


})





use_layout = list(list('sp.polygons', pa_outline_aea, col='grey70', first=F),
	list('sp.points', fia_aea, which=1, col='red', pch=3, cex=.5),
	list('sp.points', inv_aea, which=2, col='red', pch=3, cex=.5)
)

use_cuts = c(-1,0,10,50,100,200,500,1000,max(c(pa_grid_aea$fiarecs, pa_grid_aea$invrecs)))
use_col = colorRampPalette(blue2red[c(1:4,6:7)])(length(use_cuts))

pdf(paste0(fig_dir,'Num records inv vs fia plots grid.pdf'), height=6, width=4)
par(lend=1)
spplot(pa_grid_aea, c('fiarecs','invrecs'), col='transparent', sp.layout=use_layout, 
	at=use_cuts, col.regions=use_col)
dev.off()







