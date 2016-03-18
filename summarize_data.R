## This script summarizes FIA and Inventory data for Pennsylvania so that species composistion and diversity can be compared.

# Load data and settings
source('C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/load_data.R')

# Load packages
library(sp) # spatial objects
library(rgdal) # spatial projection
library(rgeos) # spatial union
library(lattice) # plotting
library(reshape) # converting lists, arrays, dataframes

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

# Define subsets of INV data to analyize
inv_subsets = data.frame(all=rep(T, nrow(inv_lichen)), epi=inv_lichen$standing, 
	macro=inv_lichen$Macro_vs_Micro=='Macrolichen')
inv_subsets$epi_macro = inv_subsets$epi&inv_subsets$macro

pa_grids = sapply(widths, function(w){
	ncells = ceiling(apply(bounds, 1, diff)/w)
	pa_grid = SpatialGrid(GridTopology(bounds[,'min']+w/2, rep(w,2), ncells), CRS(myproj))
	pa_grid = SpatialGridDataFrame(pa_grid, data.frame(ID=1:prod(ncells)))

	# Assign each record to a grid cell
	fia_sptab = over(fia_aea, pa_grid)$ID
	inv_sptab = over(inv_aea, pa_grid)$ID

	# Count number of records in each grid cell
	fiarecs = tabulate(fia_sptab, prod(ncells))
	invrecs = sapply(names(inv_subsets), function(s){
		tabulate(inv_sptab[inv_subsets[,s]], prod(ncells))
	})

	df = data.frame(pa_grid@data, fia=fiarecs, invrecs)
	names(df)[(ncol(df)-3):ncol(df)] = paste('inv', names(df)[(ncol(df)-3):ncol(df)] , sep='_')
	pa_grid = SpatialGridDataFrame(pa_grid, df)
	
	# Return grid
	pa_grid	
})
names(pa_grids) = widths/1000

# Tabulate lists of species names for each grid cell
pa_splists = sapply(pa_grids, function(pa_grid){

	# Assign each record to a grid cell
	fia_sptab = over(fia_aea, pa_grid)$ID
	inv_sptab = over(inv_aea, pa_grid)$ID
	
	# Get FIA species list for each grid cell
	fia_splist = sapply(pa_grid$ID, function(x){
		this_data = subset(fia_aea@data, fia_sptab==x)
		get_species(this_data$Binomial)
	})

	# Get INV species list for each grid cell: all, epiphytes, macrolichens, macrolichen epiphytes
	inv_splist = sapply(names(inv_subsets), function(s){
		sub_data = subset(inv_aea@data, inv_subsets[,s])
		sapply(pa_grid$ID, function(x){
			this_data = subset(sub_data, inv_sptab[inv_subsets[,s]]==x)
			get_species(this_data$Binomial)
		})
	})

	# Return species lists
	list(FIA=fia_splist, INV=inv_splist)
})
dimnames(pa_splists)[[2]] = widths/1000


## Calculate species richness in each grid cell
fia_rich = sapply(pa_splists['FIA',], function(splist){
	sapply(splist, length)
})

inv_rich = sapply(pa_splists['INV',], function(splist){
	apply(splist, c(1,2), function(x) length(unlist(x)))
})

for(w in names(pa_grids)){
	rich_data = data.frame(S_fia = fia_rich[[w]], inv_rich[[w]])	
	names(rich_data)[2:5] = paste('S_inv', names(rich_data)[2:5], sep='_')
	
	df = data.frame(pa_grids[[w]]@data, rich_data)
	pa_grids[[w]] = SpatialGridDataFrame(pa_grids[[w]], df)
}


## Plot number of records

# Define columns to plot
plot_cols = rev(c('fia','inv_all','inv_epi','inv_macro','inv_epi_macro'))
plot_cols_names = rev(c('FIA','INV','INV (epiphytes)','INV (macrolichens)','INV (epiphytic macrolichens)'))

# Define other elements of plot (PA outline, site locations)
use_layout = list(list('sp.polygons', pa_outline_aea, col='grey70', first=F),
	list('sp.points', fia_aea, which=5, col='red', pch=3, cex=.5),
	list('sp.points', inv_aea, which=1:4, col='red', pch=3, cex=.5)
)

pdf(paste0(fig_dir, 'Num records inv vs fia plots grid.pdf'), height=11, width=8.5)
for(w in as.character(widths/1000)){
	# Select grid
	this_grid = pa_grids[[w]]
	
	# Define color ramp for this grid
	maxrecs = max(this_grid@data[,plot_cols])
	use_cuts = 10
	while(use_cuts[length(use_cuts)] < maxrecs) use_cuts = c(use_cuts, 2*use_cuts[length(use_cuts)])
	use_cuts = c(-1,0,use_cuts)
	use_col = colorRampPalette(blue2red[c(1:4,6:7)])(length(use_cuts))

	print(
	spplot(this_grid, plot_cols, col='transparent', sp.layout=use_layout, layout=c(2,3), skip=c(F,T,F,F,F,F),
		at=use_cuts, col.regions=use_col, main=paste(w, 'km grid'),
		strip=strip.custom(factor.levels=plot_cols_names, bg='transparent')
	)
	)
}
dev.off()

## Plot number of species
plot_cols = paste('S', plot_cols, sep='_')

pdf(paste0(fig_dir, 'Richness inv vs fia plots grid.pdf'), height=11, width=8.5)
for(w in as.character(widths/1000)){
	# Select grid
	this_grid = pa_grids[[w]]
	
	# Define color ramp for this grid
	maxrecs = max(this_grid@data[,plot_cols])
	use_cuts = 5
	while(use_cuts[length(use_cuts)] < maxrecs) use_cuts = c(use_cuts, 2*use_cuts[length(use_cuts)])
	use_cuts = c(-1,0,use_cuts)
	use_col = colorRampPalette(blue2red[c(1:4,6:7)])(length(use_cuts))

	print(
	spplot(this_grid, plot_cols, col='transparent', sp.layout=use_layout, layout=c(2,3), skip=c(F,T,F,F,F,F),
		at=use_cuts, col.regions=use_col, main=paste(w, 'km grid'),
		strip=strip.custom(factor.levels=plot_cols_names, bg='transparent')
	)
	)
}
dev.off()


## Compare species lists of FIA vs epiphytic macrolichens in inventory








