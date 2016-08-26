## This script summarizes FIA and Inventory data for Pennsylvania so that species composistion and diversity can be compared.

# Load data and settings
source('C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/load_data.R')

# Load packages
library(sp) # spatial objects
library(rgdal) # spatial projection
library(rgeos) # spatial union
library(lattice) # plotting
library(reshape2) # converting lists, arrays, dataframes
library(vegan) # rarefaction
library(maptools)

# Read in color ramp
blue2red = read.csv('../blue2red_10colramp.txt')
blue2red = rgb(blue2red, maxColorValue=255)[10:1]

# Define subsets of INV data to analyize based on macrolichen and epiphytes
inv_subsets = data.frame(all=rep(T, nrow(inv_lichen)), epi=inv_lichen$standing, 
	macro=inv_lichen$Macro_vs_Micro=='Macrolichen')
inv_subsets$epi_macro = inv_subsets$epi&inv_subsets$macro
inv_lichen = cbind(inv_lichen, inv_subsets)



## Convert lichen records and plots to spatial data
coordinates(fia_lichen) = c('LON','LAT')
proj4string(fia_lichen) = CRS('+proj=longlat +ellps=WGS84')
coordinates(fia_plots) = c('LON','LAT')
proj4string(fia_plots) = CRS('+proj=longlat +ellps=WGS84')
coordinates(inv_lichen) = c('DarLongitude','DarLatitude')
proj4string(inv_lichen) = CRS('+proj=longlat +ellps=WGS84')
coordinates(inv_plots) = c('DarLongitude','DarLatitude')
proj4string(inv_plots) = CRS('+proj=longlat +ellps=WGS84')


###############################################################
### Plot based analysis to compare environmental models


## CHOOSING NEAREST PLOT TO COMPARE WON'T WORK B/C PLOTS ARE TOO FAR AND FIA COORDS INEXACT
# Set search radius distance
D = 5 

# For each FIA plot, find all INV plots within D km
dmat = spDists(fia_plots, inv_plots, longlat=T)

# Find dist between FIA plot and nearest INV plot
minDs = apply(dmat, 1, min)
minDs[order(minDs)]

# Find dist between INV plot and nearest FIA plot
minDs = apply(dmat, 2, min)
minDs[order(minDs)]

plots_within_D = sapply(1:20, function(D){
	apply(dmat, 1, function(x) sum(x<=D) )
})

## START HERE ##

## Calculate species richness in each plot
fia_splist = sapply(fia_plots$yrplot.id, function(y){
	these_sp = get_species(subset(fia_lichen, yrplot.id==y)$Binomial)
	these_sp[order(these_sp)]
})

fia_plots$S = sapply(fia_splist, length)


## Calculate species richness in each plot
inv_splist = sapply(inv_plots$Site_Number, function(y){
	these_sp = subset(inv_lichen, Site_Number==y)
	
	all = get_species(these_sp[these_sp$all,]$Binomial)
	all = all[order(all)]

	epi = get_species(these_sp[these_sp$epi,]$Binomial)
	epi = epi[order(epi)]

	macro = get_species(these_sp[these_sp$macro,]$Binomial)
	macro = macro[order(macro)]

	epi_macro = get_species(these_sp[these_sp$epi_macro,]$Binomial)
	epi_macro = epi_macro[order(epi_macro)]
	
	list(all = all, epi = epi, macro = macro, epi_macro = epi_macro)
})

inv_S = t(apply(inv_splist, 1:2, function(x) length(x[[1]])))
colnames(inv_S) = paste('S', colnames(inv_S), sep='_')

inv_plots = cbind(inv_plots, inv_S)

## Convert species lists to site X species matrices
species = unique(c(fia_lichen$Binomial,inv_lichen$Binomial))
species = species[order(species)]

fia_siteXsp = t(sapply(fia_splist, function(x) species %in% x))
colnames(fia_siteXsp) = species
inv_siteXsp = apply(inv_splist, 1:2, function(x) species %in% x[[1]])
inv_siteXsp = aperm(inv_siteXsp, c(2,3,1))
dimnames(inv_siteXsp)[[2]] = inv_plots$Site_Number
dimnames(inv_siteXsp)[[3]] = species


## Plot richness in FIA vs INV

## Richness estimators from vegan

fia_acc = specaccum(fia_siteXsp, method='exact')
inv_all_acc = specaccum(inv_siteXsp['all',,], method='exact')
inv_epi_acc = specaccum(inv_siteXsp['epi',,], method='exact')
inv_macro_acc = specaccum(inv_siteXsp['macro',,], method='exact')
inv_em_acc = specaccum(inv_siteXsp['epi_macro',,], method='exact')

pdf(file.path(fig_dir, 'compare_richness_plot_samplerarefaction.pdf'), height=6, width=9)

par(mfrow=c(1,2))
par(mar=c(4,4,1,1))

# All INV subsets
make_plot(c(0,210), c(0,525), xlab='Num. Plots', ylab='Num. Species')
plot(fia_acc, col='red', lwd=2, ci=2, ci.type='polygon', ci.col='#50000050', ci.lty=0, add=T)

text(124, 75, labels='FIA', adj=c(0, 1.1))

plot(inv_all_acc, col='blue', lwd=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
plot(inv_epi_acc, col='blue', lwd=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
plot(inv_macro_acc, col='blue', lwd=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
plot(inv_em_acc, col='blue', lwd=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)

text(204, c(513, 263, 188, 98), labels=paste('INV',names(invsets)), adj=c(1,-.2))

# Just compare epiphytc macrolichens
make_plot(c(0,210), c(0,100), xlab='Num. Plots', ylab='Num. Species')

plot(fia_acc, col='red', lwd=2, ci=2, ci.type='polygon', ci.col='#50000050', ci.lty=0, add=T)
text(124, 75, labels='FIA', pos=4)
plot(inv_em_acc, col='blue', lwd=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
text(204, 98, labels='INV Epiphytic Macrolichens', adj=c(1,-.5))

dev.off()

## Species richness estimators
comm_mats = list(fia=fia_siteXsp, inv_all=inv_siteXsp['all',,], inv_epi=inv_siteXsp['epi',,], 
	inv_macro=inv_siteXsp['macro',,], inv_epi_macro=inv_siteXsp['epi_macro',,])

rich_ests = sapply(comm_mats, specpool)
write.csv(rich_ests, 'richness_estimators_compared.csv')






################################################################
### Grid based analysis to compare species richness and composition across space

### Create grid across Pennsylvania

# Load ecoregions
eco = readOGR(file.path(working_dir,'pa_eco'), 'pa')
eco_ll = spTransform(eco, CRS('+proj=longlat')) # Convert to lat/lon 

# Aggregate by level 2 North American ecoregions
eco_L2 = unionSpatialPolygons(eco_ll, eco_ll$NA_L2NAME)

eco_L2_df = data.frame(L2_NAME = sapply(slot(eco_L2, "polygons"), function(x) slot(x, "ID")))
rownames(eco_L2_df) = eco_L2_df$L2_NAME
eco_L2_df$L2CODE = sapply(eco_L2_df$L2_NAME, function(x) unique(eco_ll$NA_L2CODE[eco_ll$NA_L2NAME==x]))
eco_L2_df$L1CODE = sapply(eco_L2_df$L2_NAME, function(x) unique(eco_ll$NA_L1CODE[eco_ll$NA_L2NAME==x]))
eco_L2_df$L1NAME = sapply(eco_L2_df$L2_NAME, function(x) unique(eco_ll$NA_L1NAME[eco_ll$NA_L2NAME==x]))
eco_L2 = SpatialPolygonsDataFrame(eco_L2, eco_L2_df)

# Make outline of PA
pa_outline = gUnionCascaded(eco)
pa_outline_ll = gUnionCascaded(eco_ll)

## Plot number of records and plots in each ecoregion
# Count number of records and plots in each ecoregion
fia_sptab = over(fia_lichen, eco_L2)
fia_sptab$yrplot.id = fia_lichen$yrplot.id
eco_L2$fiarecs = sapply(eco_L2$L2CODE, function(x) nrow(subset(fia_sptab, L2CODE==x)))
eco_L2$fiaplots = sapply(eco_L2$L2CODE, function(x) length(unique(subset(fia_sptab, L2CODE==x)$yrplot.id)))

inv_sptab = over(inv_lichen, eco_L2)
inv_sptab$Site_Number = inv_lichen$Site_Number
eco_L2$invrecs = sapply(eco_L2$L2CODE, function(x) nrow(subset(inv_sptab, L2CODE==x)))
eco_L2$invplots = sapply(eco_L2$L2CODE, function(x) length(unique(subset(inv_sptab, L2CODE==x)$Site_Number)))


# Add points from FIA and INV plots to maps
use_layout = list(list('sp.points', fia_plots, which=c(1,3), col='deeppink', pch=3, cex=.5),
	list('sp.points',inv_plots, which=c(2,4), col='deeppink', pch=3, cex=.5))


pdf(file.path(fig_dir,'Num records inv vs fia plots.pdf'), height=6, width=4)
par(lend=1)
spplot(eco_L2, c('fiarecs','invrecs'), col='black', sp.layout=use_layout,
	col.regions=blue2red, cuts=length(blue2red)-1, colorkey=list(at=seq(0,5000, 500)),
	strip=strip.custom(factor.levels=c('INV','FIA'), bg='white'))
dev.off()

pdf(file.path(fig_dir,'Num plots inv vs fia plots.pdf'), height=6, width=4)
par(lend=1)
spplot(eco_L2, c('fiaplots','invplots'), col='black', sp.layout=use_layout,
	col.regions=blue2red, cuts=length(blue2red)-1, colorkey=list(at=seq(0,200,20), tick.number=10),
	strip=strip.custom(factor.levels=c('INV','FIA'), bg='white'))
dev.off()


### Plot-based rarefaction comparing richness in ecoregions

fia_acc = by(fia_siteXsp, fia_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_em_acc = by(inv_siteXsp['epi_macro',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))

inv_all_acc = by(inv_siteXsp['all',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_epi_acc = by(inv_siteXsp['epi',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_macro_acc = by(inv_siteXsp['macro',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))


## Plot FIA and INV epiphytic macrolichens together

regions =  names(fia_acc)

pdf(file.path(fig_dir, 'compare_richness_ecoregion_samplerarefaction.pdf'), height=7, width=14)

par(lend=1)
par(mfrow=c(2,4))
par(mar=c(4,2,1,1))
par(oma=c(0,4,2,0))

# Just compare epiphytc macrolichens
for(reg in regions){
	nplots =  max(eco_L2@data[reg,c('invplots','fiaplots')])
	make_plot(c(0,nplots), c(0,100), xlab='', ylab=ifelse(reg==regions[1], 'Num. Species',''))
	plot(fia_acc[[reg]], col='red', lwd=2, lty=1, ci=2, ci.type='polygon', ci.col='#50000050', ci.lty=0, add=T)
	plot(inv_em_acc[[reg]], col='blue', lwd=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
	mtext(reg, 3, 0.5)
}

legend('topright', c('FIA','INV Epiphytic\nMacrolichens'), col=c('red','blue'), lty=1, bty='n', lwd=2)

# Compare INV lichens
for(reg in regions){
	nplots =  max(eco_L2@data[reg,c('invplots','fiaplots')])
	make_plot(c(0,nplots), c(0,420), xlab='Num. Plots', ylab=ifelse(reg==regions[1], 'Num. Species',''))

	plot(inv_all_acc[[reg]], col='blue', lwd=2, lty=1, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
	plot(inv_macro_acc[[reg]], col='blue', lwd=2, lty=2, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
	plot(inv_epi_acc[[reg]], col='blue', lwd=2, lty=3, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
	plot(inv_em_acc[[reg]], col='blue', lwd=2, lty=4, ci=2, ci.type='polygon', ci.col='#00005050', ci.lty=0, add=T)
}

legend('topright', c('All', 'Macrolichens', 'Epiphytes', 'Epiphytic\nMacrolichens'), col='blue', 
	lty=1:4, bty='n', lwd=2, seg.len=3)

dev.off()


## Species richness estimators table
comm_mats = list(fia=fia_siteXsp, inv_all=inv_siteXsp['all',,], inv_epi=inv_siteXsp['epi',,], 
	inv_macro=inv_siteXsp['macro',,], inv_epi_macro=inv_siteXsp['epi_macro',,])

fia_rich = by(fia_siteXsp, fia_plots$L2NAME, specpool)
inv_em_rich = by(inv_siteXsp['epi_macro',,], inv_plots$L2NAME, specpool)
inv_all_rich = by(inv_siteXsp['all',,], inv_plots$L2NAME, specpool)
inv_epi_rich = by(inv_siteXsp['epi',,], inv_plots$L2NAME, specpool)
inv_macro_rich = by(inv_siteXsp['macro',,], inv_plots$L2NAME, specpool)

rich_arr = sapply(regions, function(reg){
	as.matrix(rbind(inv_all = inv_all_rich[[reg]], inv_epi=inv_epi_rich[[reg]], 
		inv_macro=inv_macro_rich[[reg]], inv_em = inv_em_rich[[reg]], fia = fia_rich[[reg]]))
}, simplify='array')
names(dimnames(rich_arr)) = c('set','est','region')

rich_ests_reg = dcast(melt(rich_arr), region+set~est)

rich_ests_all = data.frame(region='ALL', set = colnames(rich_ests), t(rich_ests))

rich_ests_all = rbind(rich_ests_reg, rich_ests_all)
write.csv(as.matrix(rich_ests_all), 'richness_estimators_compared.csv', row.names=F)


## Compare species lists across ecoregions and data sets

fia_occ = by(fia_siteXsp, fia_plots$L2NAME, function(x) colSums(x)[colSums(x)>0])
inv_em_occ =  by(inv_siteXsp['epi_macro',,], inv_plots$L2NAME, function(x) colSums(x)[colSums(x)>0])
inv_macro_occ =  by(inv_siteXsp['macro',,], inv_plots$L2NAME, function(x) colSums(x)[colSums(x)>0])

comp_splists_reg = sapply(regions, function(reg){
	fia_sp = fia_occ[[reg]]
	em_sp = inv_em_occ[[reg]]
	macro_sp = inv_macro_occ[[reg]]

	# Species in both datasets
	both = names(fia_sp)[names(fia_sp) %in% names(macro_sp)]

	# Species only in FIA dataset
	fia = names(fia_sp)[!(names(fia_sp) %in% names(macro_sp))]
	
	# Macrolichen epiphytes only found in INV data
	inv = names(em_sp)[!names(em_sp) %in% names(fia_sp)]

	# Remove genus-only names that occur as species in the other dataset
	use_sp = get_species(c(fia, inv, both))
	fia = fia[fia %in% use_sp]
	inv = inv[inv %in% use_sp]

	# FIA species found in macrolichens that may not have been collected as epiphytes
	both_not_epi = both[!both %in% names(em_sp)]

	list(fia = fia_sp[fia], inv = macro_sp[inv], both = data.frame(fia=fia_sp[both], inv=macro_sp[both]), 
		both_not_epi = both_not_epi)
})


### WORKING HERE


# Which ones are found in 2 or fewer plots?

# Which ones are targeted by FIA (but not found)?





###################################
### OLD GRID-BASED ANALYSES

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
use_layout = list(list('sp.polygons', pa_outline_aea, col='black', first=F),
	list('sp.points', fia_aea, which=5, col='red', pch=3, cex=.5),
	list('sp.points', inv_aea, which=1:4, col='red', pch=3, cex=.5)
)

pdf(file.path(fig_dir, 'Num records inv vs fia plots grid.pdf'), height=11, width=8.5)
for(w in as.character(widths/1000)){
	# Select grid
	this_grid = pa_grids[[w]]
	
	# Define color ramp for this grid
	maxrecs = max(this_grid@data[,plot_cols])
	use_cuts = 10
	while(use_cuts[length(use_cuts)] < maxrecs) use_cuts = c(use_cuts, 2*use_cuts[length(use_cuts)])
	use_cuts = c(-1,0,use_cuts)
	use_col = c('#FFFFFF', colorRampPalette(blue2red[c(1:4,6:7)])(length(use_cuts)-1))

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

pdf(file.path(fig_dir, 'Richness inv vs fia plots grid.pdf'), height=11, width=8.5)
for(w in as.character(widths/1000)){
	# Select grid
	this_grid = pa_grids[[w]]
	
	# Define color ramp for this grid
	maxrecs = max(this_grid@data[,plot_cols])
	use_cuts = 5
	while(use_cuts[length(use_cuts)] < maxrecs) use_cuts = c(use_cuts, 2*use_cuts[length(use_cuts)])
	use_cuts = c(-1,0,use_cuts)
	use_col = c('#FFFFFF', colorRampPalette(blue2red[c(1:4,6:7)])(length(use_cuts)-1))

	print(
	spplot(this_grid, plot_cols, col='transparent', sp.layout=use_layout, layout=c(2,3), skip=c(F,T,F,F,F,F),
		at=use_cuts, col.regions=use_col, main=paste(w, 'km grid'),
		strip=strip.custom(factor.levels=plot_cols_names, bg='transparent')
	)
	)
}
dev.off()


## Plot richness in FIA vs INV



## Compare species lists of FIA vs epiphytic macrolichens in inventory

splists_partition = sapply(c('macro','epi_macro'), function(set){

	sapply(as.character(widths/1000), function(w){
		inv_splist = pa_splists['INV',w][[1]][,set]
		fia_splist = pa_splists['FIA',w][[1]]

		both_splist = unique(unlist(c(inv_splist, fia_splist)))
		both_splist = both_splist[order(both_splist)]

		sp_arr = array(FALSE, dim=c(length(both_splist), length(inv_splist), 2), 
			dimnames = list(species = both_splist, cell = 1:length(inv_splist), data = c('FIA','INV'))
		)

		for(i in 1:length(inv_splist)){
			fia = fia_splist[[i]]
			inv = inv_splist[[i]]
			both = unique(c(fia, inv))
			if(length(both) > 1) sp_arr[both,i,] = t(sapply(both, function(x) c(find_species(x, fia),find_species(x, inv))))
		}
		sp_arr
	})

})


# Compare number of cells each species is found in by data set
compare_ncells = apply(splists_partition['100','macro'], c(1,3), sum)
write.csv(compare_ncells, paste0(fig_dir, 'macro_species_freq_across_100km_cells.csv'), row.names=T)

sp_overlap = sapply(splists_partition['100',], function(x){
	# Tally whether each species is present in one or both data sets
	sp_cat = colSums(t(apply(x, c(1,3), sum) > 0) * c(10,1))
	tally = table(sp_cat)
	names(tally) = c('INV','FIA','BOTH')
	tally[c('INV','BOTH','FIA')]
})

pdf(paste0(fig_dir,'species_overlap_100km_cells.pdf'), height=4, width=4)
par(mar=c(3,4,1,1))
bp = barplot(sp_overlap, las=1, legend=T, args.legend = list(bty='n'), ylab='Number of Species',
	col=c('grey60','grey80','white'), names.arg = c('Macrolichens','Epiphytic Macrolichens'))

midpoints = apply(sp_overlap, 2, function(x){
	x = c(0,x)
	sapply(2:length(x), function(i) x[i]/2 + sum(x[1:(i-1)]))
})
text(bp, t(midpoints), labels=sp_overlap)
dev.off()





###########################################################################################3
#### OLD CODE TO DELETE LATER IF NOT NEEDED





## OLD RAREFACTION USING RE-SAMPLING
reps = 1000

plots_n = seq(5, nrow(fia_plots), 5)

fia_bootS = sapply(plots_n, function(N){
	replicate(reps, sum(colSums(fia_siteXsp[sample(nrow(fia_plots), N),])>0))
})

plots_n = seq(5, nrow(inv_plots), 5)

inv_bootS = sapply(plots_n, function(N){
	replicate(reps, rowSums(apply(inv_siteXsp[,sample(nrow(inv_plots), N),], c(1,3), sum)>0))
}, simplify='array')

save(fia_bootS, inv_bootS, file='bootstrap_richness.RData')

# Calculate summary stats
fia_S = apply(fia_bootS, 2, function(x) quantile(x, c(.025, .5, .975)))
inv_S = apply(inv_bootS, c(1,3), function(x) quantile(x, c(.025, .5, .975)))

# Plot rarefaction
invsets = c('all','epi','macro','epi_macro')
names(invsets) = c('All','Epiphytes','Macrolichens','Epiphytic Macrolichens')

pdf(file.path(fig_dir, 'compare_richness_plot_subsampling.pdf'), height=6, width=9)

par(mfrow=c(1,2))
par(mar=c(4,4,1,1))

# All INV subsets
make_plot(c(0,210), c(0,525), xlab='Num. Plots', ylab='Num. Species')

xvals = seq(5, nrow(fia_plots), 5)
polygon(c(xvals, rev(xvals)), c(fia_S['2.5%',],rev(fia_S['97.5%',])),col='#50000050', border=NA)
lines(xvals, fia_S['50%',], lwd=2, col='red')
text(xvals[length(xvals)], fia_S['2.5%',length(xvals)], labels='FIA', adj=c(0, 1.1))

xvals = seq(5, nrow(inv_plots), 5)

for(i in invsets){
	polygon(c(xvals, rev(xvals)), c(inv_S['2.5%',i,], rev(inv_S['97.5%',i,])),col='#00005050', border=NA)
	lines(xvals, inv_S['50%',i,], lwd=2, col='blue')
}

text(xvals[length(xvals)], inv_S['97.5%',,length(xvals)], labels=paste('INV',names(invsets)), adj=c(1,-.2))

# Just compare epiphytc macrolichens
make_plot(c(0,210), c(0,100), xlab='Num. Plots', ylab='Num. Species')

xvals = seq(5, nrow(fia_plots), 5)
polygon(c(xvals, rev(xvals)), c(fia_S['2.5%',],rev(fia_S['97.5%',])),col='#50000050', border=NA)
lines(xvals, fia_S['50%',], lwd=2, col='red')
text(xvals[length(xvals)], fia_S['50%',length(xvals)], labels='FIA', pos=4)

xvals = seq(5, nrow(inv_plots), 5)

polygon(c(xvals, rev(xvals)), c(inv_S['2.5%','epi_macro',], rev(inv_S['97.5%','epi_macro',])),col='#00005050', border=NA)
lines(xvals, inv_S['50%','epi_macro',], lwd=2, col='blue')

text(xvals[length(xvals)], inv_S['97.5%','epi_macro',length(xvals)], labels='INV Epiphytic Macrolichens', adj=c(1,-.2))

dev.off()












