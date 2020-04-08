## This script conducts analyses and makes figures and tables for the FIA and Inventory data project.

# Load data and settings
source('C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/load_data.R')

# Load functions
source(file.path(git_dir, 'project_functions.R'))

# Load packages
library(sp) # spatial objects
library(rgdal) # spatial projection
library(rgeos) # spatial union
library(raster) # raster data
library(reshape2) # converting lists, arrays, dataframes
library(vegan) # rarefaction
library(maptools) # unionSpatialPolygons
library(stringr) # str_trim
library(mgcv) # gam

# Define directory for saving final figures
paper_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Manuscript'
results_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Results'

# Define subsets of INV data to analyize based on macrolichen and epiphytes
inv_subsets = data.frame(all=rep(T, nrow(inv_lichen)), epi=inv_lichen$standing, 
	macro=inv_lichen$Macro_vs_Micro=='Macrolichen')
inv_subsets$epi_macro = inv_subsets$epi&inv_subsets$macro
inv_lichen = cbind(inv_lichen, inv_subsets)

invsets = c('all','epi','macro','epi_macro')
names(invsets) = c('All','Epiphytes','Macrolichens','Epiphytic Macrolichens')


## Convert lichen records and plots to spatial data 
coordinates(fia_lichen) = c('LON','LAT')
proj4string(fia_lichen) = CRS('+proj=longlat +ellps=WGS84')
coordinates(fia_plots) = c('LON','LAT')
proj4string(fia_plots) = CRS('+proj=longlat +ellps=WGS84')
coordinates(inv_lichen) = c('DarLongitude','DarLatitude')
proj4string(inv_lichen) = CRS('+proj=longlat +ellps=WGS84')
coordinates(inv_plots) = c('DarLongitude','DarLatitude')
proj4string(inv_plots) = CRS('+proj=longlat +ellps=WGS84')

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

# Count number of records and plots in each ecoregion
fia_sptab = over(fia_lichen, eco_L2)
fia_sptab$yrplot.id = fia_lichen$yrplot.id
eco_L2$fiarecs = sapply(eco_L2$L2CODE, function(x) nrow(subset(fia_sptab, L2CODE==x)))
eco_L2$fiaplots = sapply(eco_L2$L2CODE, function(x) length(unique(subset(fia_sptab, L2CODE==x)$yrplot.id)))
inv_sptab = over(inv_lichen, eco_L2)
inv_sptab$Site_Number = inv_lichen$Site_Number
eco_L2$invrecs = sapply(eco_L2$L2CODE, function(x) nrow(subset(inv_sptab, L2CODE==x)))
eco_L2$invplots = sapply(eco_L2$L2CODE, function(x) length(unique(subset(inv_sptab, L2CODE==x)$Site_Number)))

# Add ecoregion to plot data
fia_plots$L2NAME = over(fia_plots, eco_L2)$L2_NAME
inv_plots$L2NAME = over(inv_plots, eco_L2)$L2_NAME

# Define regions
regions =  unique(inv_plots$L2NAME); names(regions) = c('Atlantic Highlands','Ozark/Ouachita-Appalachian Forests','Mixed Wood Plains','Southeastern USA Plains')
regions = regions[c('Atlantic Highlands','Mixed Wood Plains','Ozark/Ouachita-Appalachian Forests','Southeastern USA Plains')]

# Define colors for data sets
data_cols =  c('#FF0000','#0000FF'); names(data_cols) = c('FIA','INV')


###############################################################
### Map plot level species richness

# Make outline of PA
pa_outline = gUnionCascaded(eco)
pa_outline_ll = gUnionCascaded(eco_ll)

## Load environmental data layers
env_dir = 'C:/Users/jrcoyle/Documents/Research/GIS DATA/'

# max vpd (from PRISM 30-year normals)
vpdmax = raster(file.path(env_dir, 'PRISM/Normals81', 
	'PRISM_vpdmax_30yr_normal_800mM2_annual_bil', 'PRISM_vpdmax_30yr_normal_800mM2_annual_bil.bil'))

# annual precip (from PRISM 30-year normals)
ap = raster(file.path(env_dir, 'PRISM/Normals81', 
	'PRISM_ppt_30yr_normal_800mM2_annual_bil', 'PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'))

# total (wet+dry estimated) N and S deposition (from NTN)
tot_n = raster(file.path(env_dir, 'NTN', 'total_N_dep_2001-2010.tif'))
tot_s = raster(file.path(env_dir, 'NTN', 'total_S_dep_2001-2010.tif'))

# Crop to PA
pa_outline_ll = spTransform(pa_outline_ll, proj4string(vpdmax))
totN_PA = projectRaster(tot_n, crs=CRS(proj4string(pa_outline_ll)))
totS_PA = projectRaster(tot_s, crs=CRS(proj4string(pa_outline_ll)))

use_ext = extent(pa_outline_ll)+rep(c(-0.25, .25), 2)
vpd_PA = crop(vpdmax, use_ext)
vpd_PA = mask(vpd_PA, pa_outline_ll)
ap_PA = crop(ap, use_ext)
ap_PA = mask(ap_PA, pa_outline_ll)
totN_PA = crop(totN_PA, use_ext)
totN_PA = mask(totN_PA, pa_outline_ll)
totS_PA = crop(totS_PA, use_ext)
totS_PA = mask(totS_PA, pa_outline_ll)

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


## Write out derived data tables for publication?



## Table 1: Description of Datasets

table1 = data.frame(Dataset=c('FIA','INV'), Plots=c(nrow(fia_plots), nrow(inv_plots)))
table1$Records = c(nrow(fia_lichen), nrow(inv_lichen))
table1$Total_Epiphytic_Macrolichens = c(length(get_species(fia_lichen$Binomial)), length(get_species(inv_lichen[inv_lichen$epi_macro,]$Binomial)))
table1$Total_Epiphytic_Lichens = c(length(get_species(fia_lichen$Binomial)), length(get_species(inv_lichen[inv_lichen$epi,]$Binomial)))
table1$Total_Macrolichens = c(length(get_species(fia_lichen$Binomial)), length(get_species(inv_lichen[inv_lichen$macro,]$Binomial)))
table1$Total_Lichens = c(length(get_species(fia_lichen$Binomial)), length(get_species(inv_lichen[inv_lichen$all,]$Binomial)))

# Save
write.csv(table1, file.path(paper_dir, 'Table_1_datasets.csv'), row.names=F)



## Figure 1: Map of Species Richness with Environmental Layers ##

rich_cols = colorRampPalette(c('#ffffb2','#f03b20'))(11)
env_cols = colorRampPalette(c('grey95','grey10'))(50)
reg_cols = c('white','grey80', 'grey60','grey40','grey20')
reg_accs = c('AH','MWP','OAF','SP'); names(reg_accs)= regions

# Add background environment and ecoregions from FIA and INV plots to maps
reg_lines = list('sp.polygons',eco_L2, col='black', lwd=2, first=F)
reg_labels = list('sp.text', cbind(c(-78.5, -75.5, -80.2, -76, -78.5, -76.5), c(41.7, 41.1, 41.7, 41.7, 40.5, 40)), reg_accs[c(1,1,2,2,3,4)], col='black', font=2)	
fia_S_points = list('sp.points', fia_plots, col=rich_cols[cut(fia_plots$S, seq(0,33,3), include.lowest=T)], pch=16, cex=1.1)
inv_S_em_points = list('sp.points', inv_plots, col=rich_cols[cut(inv_plots$S_epi_macro, seq(0,33,3), include.lowest=T)], pch=16, cex=1.1)	

# Layout settings
#lattice_pars = list(layout.widths=list(key.right=10))

# Order points by richness so that highest richness occurs on top
fia_plots = fia_plots[order(fia_plots$S),]
inv_plots = inv_plots[order(inv_plots$S_epi_macro),]

# Panel A: FIA + VPD
svg(file.path(paper_dir, 'Fig_1A-FIA+VPD.svg'), height=4, width=6)
use_layout = list(reg_lines, fia_S_points, reg_labels)
spplot(vpd_PA, col.regions=env_cols, cuts=49, sp.layout=use_layout,
	#par.settings=list(axis.line = list(col = 'transparent')),
	colorkey=list(height=0.8), scales=list(draw=T),
	main=list(expression(Max.~VPD~(kPa)), font=1, cex=1)
)
dev.off()

# Panel B: INV + N Deposition
svg(file.path(paper_dir, 'Fig_1B-INV+Ndep.svg'), height=4, width=6)
use_layout = list(reg_lines, inv_S_em_points, reg_labels)
spplot(totN_PA, col.regions=env_cols, cuts=49, sp.layout=use_layout,
	#par.settings=list(axis.line = list(col = 'transparent')),
	colorkey=list(height=0.8), scales=list(draw=T),
	main=list(expression(Total~N~Deposition~(kg~ha^-1)), font=1, cex=1) #, useRaster=T
)
dev.off()

# Species Richnesss key for A and B
svg(file.path(paper_dir, 'Fig_1AB-colorbar.svg'), height=6, width=6)
plot.new()
plotColorRamp(rich_cols, 11, c(.5, .1, .6, .9), labels=seq(0,33,3), title='Number of Species', horiz=F, ndig=0, side=2)
dev.off()


# Order points by richness so that highest richness occurs on top
inv_plots = inv_plots[order(inv_plots$S_macro),]
inv_S_m_points = list('sp.points', inv_plots, col=rich_cols[cut(inv_plots$S_macro, seq(0,110,10), include.lowest=T)], pch=16, cex=1.1)	

# Panel C: INV Macrolichens + Annual precip
svg(file.path(paper_dir, 'Fig_1C-INVmacro+AP.svg'), height=4, width=6)
par(mar=c(3,4,5,0))
use_layout = list(reg_lines, inv_S_m_points, reg_labels)
spplot(ap_PA, col.regions=env_cols, cuts=49, sp.layout=use_layout,
	#par.settings=list(axis.line = list(col = 'transparent')),
	colorkey=list(height=0.8), scales=list(draw=T),
	main=list(expression(Annual~Precip.~(mm)), font=1, cex=1) #, useRaster=T
)
dev.off()

# Order points by richness so that highest richness occurs on top
inv_plots = inv_plots[order(inv_plots$S_all),]
inv_S_points = list('sp.points', inv_plots, col=rich_cols[cut(inv_plots$S_all, seq(0,110,10), include.lowest=T)], pch=16, cex=1.1)	

# Panel D: INV All lichens +  S deposition
svg(file.path(paper_dir, 'Fig_1D-INVall+Sdep.svg'), height=4, width=6)
par(mar=c(3,4,5,0))
use_layout = list(reg_lines, inv_S_points, reg_labels)
spplot(totS_PA, col.regions=env_cols, cuts=49, sp.layout=use_layout,
	#par.settings=list(axis.line = list(col = 'transparent')),
	colorkey=list(height=0.8), scales=list(draw=T),
	main=list(expression(Total~S~Deposition~(kg~ha^-1)), font=1, cex=1) #, useRaster=T
)
dev.off()

# Species Richnesss key for C and D
svg(file.path(paper_dir, 'Fig_1CD-colorbar.svg'), height=6, width=6)
plot.new()
plotColorRamp(rich_cols, 11, c(.5, .1, .6, .9), labels=seq(0,110,10), title='Number of Species', horiz=F, ndig=0, side=2)
dev.off()


## Figure 2: Compare richness distributions ##

# Tabulate number of species per plot
fia_tab = table(factor(fia_plots$S, 0:33))
inv_tab = table(factor(inv_plots$S_epi_macro, 0:33))
fia_tab[fia_tab==0] = NA
inv_tab[inv_tab==0] = NA


# Distance to shift bars
jit = 0.15

# Make figure
pdf(file.path(paper_dir, 'Fig_2-richness_distributions.pdf'), height=3, width=4)
par(mar=c(4,4,1,1))
par(lend=1)
make_plot(c(0,33), c(0, 25), xlab='Number of Species', ylab='Number of Plots')
points((0:33)-jit, inv_tab, type='h', lwd=2, col=data_cols['INV'])
points((0:33)+jit, fia_tab, type='h', lwd=2, col=data_cols['FIA'])
legend('topright', names(data_cols), fill=data_cols, bty='n', border=NA)
dev.off()

svg(file.path(paper_dir, 'Fig_2-richness_distributions.svg'), height=3, width=4)
par(mar=c(4,4,1,1))
par(lend=1)
make_plot(c(0,33), c(0, 25), xlab='Number of Species', ylab='Number of Plots')
points((0:33)-jit, inv_tab, type='h', lwd=2, col=data_cols['INV'])
points((0:33)+jit, fia_tab, type='h', lwd=2, col=data_cols['FIA'])
legend('topright', names(data_cols), fill=data_cols, bty='n', border=NA)
dev.off()



###############################################################
### Compare species richness across data sets

## Convert species lists to site X species matrices
species = unique(c(fia_lichen$Binomial,inv_lichen$Binomial))
species = species[order(species)]

fia_siteXsp = t(sapply(fia_splist, function(x) species %in% x))
colnames(fia_siteXsp) = species
inv_siteXsp = apply(inv_splist, 1:2, function(x) species %in% x[[1]])
inv_siteXsp = aperm(inv_siteXsp, c(2,3,1))
dimnames(inv_siteXsp)[[2]] = inv_plots$Site_Number
dimnames(inv_siteXsp)[[3]] = species

## Rarefaction of species richness

fia_acc_tot = specaccum(fia_siteXsp, method='exact')
inv_all_acc_tot = specaccum(inv_siteXsp['all',,], method='exact')
inv_epi_acc_tot = specaccum(inv_siteXsp['epi',,], method='exact')
inv_macro_acc_tot = specaccum(inv_siteXsp['macro',,], method='exact')
inv_em_acc_tot = specaccum(inv_siteXsp['epi_macro',,], method='exact')

## Rarefaction of species richness within ecoregions

fia_acc = by(fia_siteXsp, fia_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_em_acc = by(inv_siteXsp['epi_macro',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_all_acc = by(inv_siteXsp['all',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_epi_acc = by(inv_siteXsp['epi',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))
inv_macro_acc = by(inv_siteXsp['macro',,], inv_plots$L2NAME, function(x) specaccum(x, method='exact'))

# Define line types
use_lty = 1:4; names(use_lty) = c('all','epi','macro','em')

# Define colors
use_col = data_cols

# Define CI transparency
ci_trans = '40'

## Figure 3: Rarefaction plot ##
#pdf(file.path(paper_dir, 'Fig_3.pdf'), height=9.5, width=8)
svg(file.path(paper_dir, 'Fig_3.svg'), height=9.5, width=8)

# Set up multipanel plot with two columns and 4 rows
layout(matrix(c(9,10,1,2,3,4,5,6,7,8), ncol=2, byrow=T), heights=c(.5,1,1,1,1))

# Set plot margins
par(mar=c(3,3,3,2))
par(oma=c(3,4,1,1))

# Set up graphics options
par(lend=1)

# Plot totals across all ecoregions

# Epiphytic macrolichens (INV vs FIA)
make_plot(c(0,210), c(0,105), ylab='Number of Species', xlab='Number of Plots', cex=.8)
plot(fia_acc_tot, col=use_col['FIA'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', 
	ci.col=paste0(use_col['FIA'],ci_trans), ci.lty=0, add=T)
#text(124, 75, labels='FIA', pos=4)
plot(inv_em_acc_tot, col=use_col['INV'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', 
	ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
#text(204, 98, labels='INV', pos=3)
add_panel_label(1, x=-.05, y=0.95, add_sym=') Pennsylvania')
#mtext('All Regions', 3, 2, adj=0)

# Species groups within INV
make_plot(c(0,210), c(0,525), xlab='Number of Plots', cex=.8)
plot(fia_acc_tot, col=use_col['FIA'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', 
	ci.col=paste0(use_col['FIA'],ci_trans), ci.lty=0, add=T)
#text(124, 75, labels='FIA', adj=c(0, 1.1))

plot(inv_all_acc_tot, col=use_col['INV'], lwd=2, lty=use_lty['all'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
plot(inv_epi_acc_tot, col=use_col['INV'], lwd=2, lty=use_lty['epi'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
plot(inv_macro_acc_tot, col=use_col['INV'], lwd=2, lty=use_lty['macro'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
plot(inv_em_acc_tot, col=use_col['INV'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)

#text(204, c(513, 263, 188, 98), labels=paste('INV',names(invsets)), adj=c(1,-.2))
add_panel_label(2, x=-.05, y=0.95, add_sym=') Pennsylvania')

# Plot richness within ecoregions (except southeastern plains)
which_panel=3
for(reg in regions[1:3]){
	nplots =  max(eco_L2@data[reg,c('invplots','fiaplots')])

	# Compare INV and FIA
	make_plot(c(0,nplots), c(0,100), xlab='Number of Plots', ylab='Number of Species', cex=0.8)
	plot(fia_acc[[reg]], col=use_col['FIA'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', ci.col=paste0(use_col['FIA'],ci_trans), ci.lty=0, add=T)
	plot(inv_em_acc[[reg]], col=use_col['INV'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
	add_panel_label(which_panel, x=-.05, y=0.95, add_sym=paste(')', names(regions[regions==reg])))
	which_panel = which_panel+1
	#mtext(names(regions[regions==reg]), 3, 2, adj=0)
	

	# Compare species groups within INV
	make_plot(c(0,nplots), c(0,420), 'Number of Plots', cex=0.8)
	plot(fia_acc[[reg]], col=use_col['FIA'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', ci.col=paste0(use_col['FIA'],ci_trans), ci.lty=0, add=T)
	plot(inv_all_acc[[reg]], col=use_col['INV'], lwd=2, lty=use_lty['all'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
	plot(inv_macro_acc[[reg]], col=use_col['INV'], lwd=2, lty=use_lty['macro'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
	plot(inv_epi_acc[[reg]], col=use_col['INV'], lwd=2, lty=use_lty['epi'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
	plot(inv_em_acc[[reg]], col=use_col['INV'], lwd=2, lty=use_lty['em'], ci=2, ci.type='polygon', ci.col=paste0(use_col['INV'],ci_trans), ci.lty=0, add=T)
	add_panel_label(which_panel, x=-.05, y=0.95, add_sym=paste(')', names(regions[regions==reg])))
	which_panel = which_panel+1
}

# Add legend to side
par(mar=c(0,0,0,0))
plot.new()
legend('bottom', c(paste('INV',names(invsets)[4]), 'FIA (Epiphytic Macrolichens)'), col=use_col[c(2,1)],
	lty=use_lty[4], lwd=2, seg.len=3, bty='n')

plot.new()
legend('bottom', paste('INV',names(invsets)[1:3]), col=use_col[2],
	lty=use_lty[1:3], lwd=2, seg.len=3, bty='n')

dev.off()


## Species richness estimators table

## Species richness estimators table
comm_mats = list(fia=fia_siteXsp, inv_all=inv_siteXsp['all',,], inv_epi=inv_siteXsp['epi',,], 
	inv_macro=inv_siteXsp['macro',,], inv_epi_macro=inv_siteXsp['epi_macro',,])

rich_ests = sapply(comm_mats, specpool)
colnames(rich_ests)[5] = 'inv_em'

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

# Table with all data
rich_ests_reg = dcast(melt(rich_arr), region+set~est)
rich_ests_all = data.frame(region='ALL', set = colnames(rich_ests), t(rich_ests))
rich_ests_all = rbind(rich_ests_reg, rich_ests_all)
write.csv(as.matrix(rich_ests_all), file.path(fig_dir, 'richness_estimators_compared.csv'), row.names=F)

## Table 2: Chao estimator of species richness
chao = with(rich_ests_all, paste(round(unlist(chao)), '+-', round(unlist(chao.se))))

chao_df = data.frame(rich_ests_all[,c('region', 'set')], chao)
chao_tab = dcast(melt(chao_df, id.vars=c('region','set')), region~set) 
write.csv(chao_tab, file.path(paper_dir, 'Table_2-chao_estimator.csv'), row.names=F)



#######################################################################
### Compare species lists across data sets

# Species occurence within regions
fia_occ = by(fia_siteXsp, fia_plots$L2NAME, function(x) colSums(x)[colSums(x)>0])
inv_em_occ =  by(inv_siteXsp['epi_macro',,], inv_plots$L2NAME, function(x) colSums(x)[colSums(x)>0])
inv_macro_occ =  by(inv_siteXsp['macro',,], inv_plots$L2NAME, function(x) colSums(x)[colSums(x)>0])

# Compare species lists from INV vs FIA within regions
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

# Convert to data tables and save
df_splists_reg = sapply(regions, function(reg){
	fia_sp = fia_occ[[reg]]
	em_sp = inv_em_occ[[reg]]
	macro_sp = inv_macro_occ[[reg]]

	sp = unique(c(names(fia_sp), names(em_sp)))
	
	nplots_df = data.frame(t(sapply(sp, function(x){
		fia = ifelse(x %in% names(fia_sp), fia_sp[x], 0)
		inv = ifelse(x %in% names(em_sp), em_sp[x], 0)
		c(fia, inv)
	})))
	names(nplots_df) = c('FIA','INV')

	# If the species in in INV but not counted as an epiphyte change this value
	macro_only = rownames(nplots_df)[which((rownames(nplots_df) %in% names(macro_sp))& nplots_df$INV==0)]
	nplots_df[macro_only,'INV'] = macro_sp[macro_only]

	nplots_df$Binomial = rownames(nplots_df)
	
	nplots_df[,c('Binomial','FIA','INV')]
	
})
for(reg in names(regions)){
	df = data.frame(df_splists_reg[,reg])
	regname = gsub('/', '-', reg)
	write.csv(df, file.path(results_dir, paste0('species_freqs_', regname, '.csv')), row.names=F)
}

## Figure S4: Number of species in INV and FIA across regions
pdf(file.path(paper_dir, 'SI', 'Fig_S4-compare_species_cross_region.pdf'), height=3, width=8.5)
par(mar=c(4.5,19,1,1.5))
barplot( comp_counts_reg[c('inv','both','fia'),], horiz=T, las=1, xlim=c(0,100),
	legend.text=c('INV','BOTH','FIA'), xlab='Num. Species', 
	args.legend=list(x='topright', horiz=T, bty='n'))
dev.off()


## Compare species lists from INV vs FIA across all regions
fia_sp = colSums(fia_siteXsp)[colSums(fia_siteXsp)>0]
em_sp = colSums(inv_siteXsp['epi_macro',,])[colSums(inv_siteXsp['epi_macro',,])>0]
macro_sp = colSums(inv_siteXsp['macro',,])[colSums(inv_siteXsp['macro',,])>0]

both = names(fia_sp)[names(fia_sp) %in% names(macro_sp)]
fia = names(fia_sp)[!(names(fia_sp) %in% names(macro_sp))]
inv = names(em_sp)[!names(em_sp) %in% names(fia_sp)]
use_sp = get_species(c(fia, inv, both))
fia = fia[fia %in% use_sp]
inv = inv[inv %in% use_sp]
both_not_epi = both[!both %in% names(em_sp)]

comp_splists = list(fia = fia_sp[fia], inv = macro_sp[inv], both = data.frame(fia=fia_sp[both], inv=macro_sp[both]), 
		both_not_epi = both_not_epi)

# Number of species shared across regions
comp_counts_reg = apply(comp_splists_reg, 1:2, function(x) length(unlist(x)))
comp_counts_reg['both',] = comp_counts_reg['both',]/2

# Proportion of total species pool in each data set
plot_data = comp_counts_reg[c('inv','both','fia'),]
prop_sp = t(plot_data)/colSums(plot_data)
comp_counts = sapply(comp_splists, function(x) length(unlist(x)))
comp_counts['both'] = comp_counts['both']/2
prop_sp = rbind(prop_sp, ALL=comp_counts[c('inv','both','fia')] / sum(comp_counts[c('inv','both','fia')]))

write.csv(prop_sp, file.path(working_dir, 'proportion_species_across_datasets.csv'))


## Write out species list tables with indication of which regions found in

# INV-only species
inv = data.frame(All=comp_splists$inv)
inv_reg = by(inv_siteXsp['macro',,rownames(inv)], inv_plots$L2NAME, colSums)
inv = cbind(inv, sapply(regions, function(reg) inv_reg[[reg]]))		

# Which ones are targeted by FIA (but not found)?
fia_ref = read.csv(file.path(working_dir, 'REF_LICHEN_SPECIES.csv'))
fia_ref$Binomial = str_trim(paste(fia_ref$GENUS, fia_ref$SPECIES))

# Convert species names that differ between FIA and INV
taxa_mapper = read.csv('C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/FIA_TAXONOMY_conversion.csv')
change_taxa = subset(taxa_mapper, FIA_NAME!=INVENTORY_NAME)
for(i in 1:nrow(change_taxa)){
	fia_ref[fia_ref$Binomial==change_taxa[i,'FIA_NAME'],'Binomial'] = change_taxa[i,'INVENTORY_NAME']
}
inv$FIA_targeted = rownames(inv) %in% fia_ref$Binomial

# FIA-only species
fia = data.frame(All=comp_splists$fia)
fia_reg = by(fia_siteXsp[,rownames(fia)], fia_plots$L2NAME, colSums)
fia = cbind(fia, sapply(regions, function(reg) fia_reg[[reg]]))

# Calculate Chi-sq statistic for species found in both FIA and INV
comp_chisq = function(Nfia, Ninv, Tfia, Tinv){
	tab = matrix(c(Tfia-Nfia, Nfia, Tinv-Ninv, Ninv), nrow=2) 

	sr <- rowSums(tab)
      sc <- colSums(tab)
      E <- outer(sr, sc, "*")/sum(tab)
	any(E < 5)

	chitest = chisq.test(tab)
	c(Chi_sq=as.numeric(chitest$statistic), P=as.numeric(chitest$p.value), use_P=all(E >= 5))
}

both = comp_splists$both
both = cbind(both, t(sapply(1:nrow(both), function(i) comp_chisq(both[i,1], both[i,2], nrow(fia_siteXsp), dim(inv_siteXsp)[2]))))

# Write out tables used in Table S1
write.csv(inv, file.path(results_dir, 'INV_only_species.csv'))
write.csv(fia, file.path(results_dir, 'FIA_only_species.csv'))
write.csv(both, file.path(results_dir, 'INV_and_FIA_species.csv'))




#######################################################################
### Environmental Models

xvars = c('ppt','vpdmax','tmax','tmean','tot_n','tot_s')
xvarnames = expression('Annual precipitation (mm)', 'Max. VPD (kPa)', 'Max. temperature'~(degree~C),
	'Mean temperature'~(degree~C), 'Total N deposition'~(kg~ha^-1), 'Total S deposition'~(kg~ha^-1))
names(xvarnames) = xvars

fia_data = fia_plots@data
inv_data = inv_plots@data
rownames(inv_data) = inv_data$Site_Number


## Figure S1: Location of plots in climate space

pdf(file.path(paper_dir, 'Fig_S1-FIA_vs_INV_environment.pdf'), height=3.5, width=9)

par(mar=c(4.5, 4.5, 1, 1))
par(mfrow=c(1, 3))

# Panel A: annual precip vs. mean annual temp
make_plot(c(900, 1400), c(6, 12.5), xlab=xvarnames['ppt'], ylab=xvarnames['tmean'], xlab_loc=3, cex=0.8)
points(tmean~ppt, data=fia_plots, col=data_cols['FIA'], pch=3)
points(tmean~ppt, data=inv_plots, col=data_cols['INV'], pch=1)

# Panel B: max temp vs. vpd
make_plot(c(12,18), c(7, 13), xlab=xvarnames['tmax'], ylab=xvarnames['vpdmax'], xlab_loc=3, cex=0.8)
points(vpdmax~tmax, data=fia_plots, col=data_cols['FIA'], pch=3)
points(vpdmax~tmax, data=inv_plots, col=data_cols['INV'], pch=1)

# Panel C: N vs S deposition
make_plot(c(90, 180), c(100, 400), xlab=xvarnames['tot_n'], ylab=xvarnames['tot_s'], xlab_loc=3, cex=0.8)
points(tot_s~tot_n, data=fia_plots, col=data_cols['FIA'], pch=3)
points(tot_s~tot_n, data=inv_plots, col=data_cols['INV'], pch=1)

legend('topright', c('FIA','INV'), col=data_cols[c('FIA','INV')], pch=c(3, 1), bty='n')


dev.off()


## Figure S2: Species richness vs environment ##
pdf(file.path(paper_dir, 'Fig_S2-richness_env.pdf'), height=8, width=16)
layout(matrix(1:(3*length(xvars)), nrow=3, byrow=F))
par(mar=c(1.5,1.5,1,1))
par(oma=c(3,4.5,0,0))
for(x in xvars){
	use_xlim = range(c(fia_data[,x], inv_data[,x]))
	
	# Plot FIA richnes vs env
	plot(fia_data[,x], fia_data$S, las=1, xlim=use_xlim, xlab='', ylab='')
	if(x==xvars[1]) mtext('FIA Num. Species', 2, 3)

	# Add GAM lines and 2*SE lines
	xvals = data.frame(seq(min(fia_data[,x]), max(fia_data[,x]), length.out=100))
	names(xvals)=x
	mod = gam(bquote(S~s(.(as.name(x)))), data=fia_data, family=nb(link='log'))
	pred = predict(mod, newdata=xvals, type='response', se.fit=T)
	lines(xvals[,1], pred$fit, col=2)
	lines(xvals[,1], pred$fit + pred$se.fit*2, col=2, lty=2)
	lines(xvals[,1], pred$fit - pred$se.fit*2, col=2, lty=2)

	# Plot INV epiphytic macrolichen richness vs env
	plot(inv_data[,x], inv_data$S_epi_macro, las=1, xlim=use_xlim, xlab='', ylab='')
	if(x==xvars[1]) mtext('INV Epiphytic Macroclichen\nNum. Species', 2, 3)

	# Add GAM lines and 2*SE lines
	xvals = data.frame(seq(min(inv_data[,x]), max(inv_data[,x]), length.out=100))
	names(xvals)=x
	mod = gam(bquote(S_epi_macro~s(.(as.name(x)))), data=inv_data, family=nb(link='log'))
	pred = predict(mod, newdata=xvals, type='response', se.fit=T)
	lines(xvals[,1], pred$fit, col=2)
	lines(xvals[,1], pred$fit + pred$se.fit*2, col=2, lty=2)
	lines(xvals[,1], pred$fit - pred$se.fit*2, col=2, lty=2)	

	# Plot INV all species richness vs env
	plot(inv_data[,x], inv_data$S_epi, las=1, xlim=use_xlim)
	mtext(xvarnames[x], 1, 3)
	if(x==xvars[1]) mtext('INV Epiphyte\nNum. Species', 2, 3)

	# Add GAM lines and 2*SE lines
	mod = gam(bquote(S_epi~s(.(as.name(x)))), data=inv_data, family=nb(link='log'))
	pred = predict(mod, newdata=xvals, type='response', se.fit=T)
	lines(xvals[,1], pred$fit, col=2)
	lines(xvals[,1], pred$fit + pred$se.fit*2, col=2, lty=2)
	lines(xvals[,1], pred$fit - pred$se.fit*2, col=2, lty=2)

}
dev.off()

## Figure 4: Individual Species Distributions ##

# Focus on macrolichen species occuring on at least 6 sites in both datasets
inv_m_comm = inv_siteXsp['macro',,]
fia_comm = fia_siteXsp
focal_sp = c('Hypogymnia physodes','Physcia millegrana','Usnocetraria oakesiana','Flavoparmelia caperata')
focal_names = sapply(focal_sp, function(x) paste(unlist(strsplit(x, ' ')), collapse='_'))

use_xvars = c('vpdmax','tmean','tot_n','tot_s')

# Add species occurences to dataframes for modeling
fia_sp = fia_comm[,focal_sp]; colnames(fia_sp) = focal_names
fia_data = cbind(fia_data, fia_sp)
inv_sp = inv_m_comm[,focal_sp]; colnames(inv_sp) = focal_names
inv_data = cbind(inv_data, inv_sp)

# Define tansparency for CI and data points
ci_trans = '50'
rug_trans = '90'

#pdf(file.path(paper_dir, 'Fig_4-species_env_distributions.pdf'), height=8.5, width=10)
svg(file.path(paper_dir, 'Fig_4-species_env_distributions.svg'), height=8.5, width=10)

par(lend=1)
par(mfrow=c(length(focal_names), length(use_xvars)))
par(mar=c(3, 3, 2.5, 1))
par(oma=c(2, 8, 0, 0))

for(sp in focal_names){
for(x in use_xvars){
	counter = which(focal_names==sp)

	fia_mod = gam(bquote(.(as.name(sp)) ~ s(.(as.name(x)), k=4)), data=fia_data, family=binomial(link='logit'))
	inv_mod = gam(bquote(.(as.name(sp)) ~ s(.(as.name(x)), k=4)), data=inv_data, family=binomial(link='logit'))

	xrange = range(c(fia_data[,x], inv_data[,x]))
	
	# FIA
	make_plot(xrange, c(0,1), xlab=ifelse(counter%%length(use_xvars)==0, xvarnames[x], ''), cex=.8)
	points(fia_data[,x], fia_data[,sp], pch='|', col=paste0(data_cols['FIA'], rug_trans))
	xvals = data.frame(seq(min(fia_data[,x]), max(fia_data[,x]), length.out=100))
	names(xvals)=x
	fia_pred = predict(fia_mod, newdata=xvals, type='response', se.fit=T)
	polygon(c(xvals[,1], rev(xvals[,1])), c(fia_pred$fit-fia_pred$se.fit, rev(fia_pred$fit+fia_pred$se.fit)),
		col=paste0(data_cols['FIA'], ci_trans), border=NA)
	use_lty = ifelse(summary(fia_mod)$s.table[1,4] < 0.05, 1, 3)
	lines(xvals[,1], fia_pred$fit, col=data_cols['FIA'], lwd=2, lty=use_lty)

	# INV
	points(inv_data[,x], inv_data[,sp], pch='|', col=paste0(data_cols['INV'], rug_trans))
	xvals = data.frame(seq(min(inv_data[,x]), max(inv_data[,x]), length.out=100))
	names(xvals)=x
	inv_pred = predict(inv_mod, newdata=xvals, type='response', se.fit=T)
	polygon(c(xvals[,1], rev(xvals[,1])), c(inv_pred$fit-inv_pred$se.fit, rev(inv_pred$fit+inv_pred$se.fit)),
		col=paste0(data_cols['INV'], ci_trans), border=NA)
	use_lty = ifelse(summary(inv_mod)$s.table[1,4] < 0.05, 1, 3)
	lines(xvals[,1],inv_pred$fit, col=data_cols['INV'], lwd=2, lty=use_lty)

	mtext(format(summary(fia_mod)$dev.expl, digits=2), 3, 0, col=data_cols['FIA'], adj=0, cex=0.8)
	mtext(format(summary(inv_mod)$dev.expl, digits=2), 3, 0, col=data_cols['INV'], adj=1, cex=0.8)

	if(x==use_xvars[1]){
		mtext('Prob. occurence', 2, 2.5, cex=0.8)	
		par(xpd=NA)
		add_panel_label(counter, add_sym=paste(')', focal_sp[counter]), x=-.1, y=.99)
		par(xpd=F)
	}
}
}
dev.off()


# Code for Figure S3 containing all species with > 6 occurances can be found in summarize_data.R script.


### NMDS
### From Will-wolf et al 2006:
### Unconstrained NMS with Sorenson dissimilarity 
### Exclude species found at < 3 plots
### Exclude outlier plots with >2.5 std dev away from avg pairwise Sorensen distance
### Standardized by total abundance
### 3 axes
### Pearson and kendall correlations with ordination acxes

# Make community data matrices
fia_comm = fia_siteXsp[,colSums(fia_siteXsp)>0]
inv_em_comm = inv_siteXsp['epi_macro',,]
inv_em_comm = inv_em_comm[, colSums(inv_em_comm)>0]
inv_e_comm = inv_siteXsp['epi',,]
inv_e_comm = inv_e_comm[,colSums(inv_e_comm)>0]

# Remove species occuring fewer than 3 times
dim(fia_comm) #75
fia_comm = fia_comm[,colSums(fia_comm)>2]; dim(fia_comm) # 37 species, lost 38

dim(inv_e_comm) # 263
inv_e_comm = inv_e_comm[,colSums(inv_e_comm)>2]; dim(inv_e_comm) # 146 species, lost 137
inv_e_comm = inv_e_comm[rowSums(inv_e_comm)>0,]; dim(inv_e_comm) # 193 plots, lost 12

dim(inv_em_comm) # 98
inv_em_comm = inv_em_comm[,colSums(inv_em_comm)>2]; dim(inv_em_comm) # 55 species, lost 43
inv_em_comm = inv_em_comm[rowSums(inv_em_comm)>0,]; dim(inv_em_comm) # 185 plots, lost 20

# Run analyses with data subset to the same size as FIA
# These analyses are reported in the SI
# DONT RUN UNLESS INTERESTED IN SUBSETTING INV DATA` 
# NOTE: when saving files, names have _subsample added to them
nplots = nrow(fia_comm)
#inv_em_comm = inv_em_comm[sample(nrow(inv_em_comm), nplots),]
#inv_e_comm = inv_e_comm[sample(nrow(inv_e_comm), nplots),]

# Calculate Sorenson dissimilarity to view outliers
# None need to be excluded
fact = 2.5

fia_dist = 1-vegdist(fia_dat, method='bray')
which(rowMeans(as.matrix(fia_dist)) > mean(fia_dist) + fact*sqrt(var(fia_dist)))

inv_e_dist = 1-vegdist(inv_e_dat, method='bray')
which(rowMeans(as.matrix(inv_e_dist)) > mean(inv_e_dist) + fact*sqrt(var(inv_e_dist)))

inv_em_dist = 1-vegdist(inv_em_dat, method='bray')
which(rowMeans(as.matrix(inv_em_dist)) > mean(inv_em_dist) + fact*sqrt(var(inv_em_dist)))

# Perform NMDS and fit environmental variables to FIA data
fia_mds = metaMDS(fia_comm, distance='bray', k=4, autotransform=F, trymax=300, maxit=500)
fia_ef = envfit(fia_mds, fia_data[,c(xvars,'L2NAME')], choices=1:4)
fia_surf = sapply(xvars, function(x){
	mod12 = eval(bquote(with(fia_data, ordisurf(fia_mds, .(as.name(x)), choices=1:2))))
	mod34 = eval(bquote(with(fia_data, ordisurf(fia_mds, .(as.name(x)), choices=3:4))))
	list(mod12=mod12, mod34=mod34)
})
fia_ss = scores(fia_mds, 'sites', choices=1:4) # Oksanen suggests not doing this.
fia_cor = cor(fia_comm, fia_ss, method='kendall')

# Perform NMDS and fit environmental variables to INV epiphytes
inv_e_mds = metaMDS(inv_e_comm, distance='bray', k=4, autotransform=F, trymax=300, maxit=500)
inv_e_ef = envfit(inv_e_mds, inv_data[rownames(inv_e_comm),c(xvars,'L2NAME')], choices=1:4)
inv_e_surf = sapply(xvars, function(x){
	mod12 = eval(bquote(with(inv_data[rownames(inv_e_comm),], ordisurf(inv_e_mds, .(as.name(x)), choices=1:2))))
	mod34 = eval(bquote(with(inv_data[rownames(inv_e_comm),], ordisurf(inv_e_mds, .(as.name(x)), choices=3:4))))
	list(mod12=mod12, mod34=mod34)
})
inv_e_ss = scores(inv_e_mds, 'sites', choices=1:4) # Oksanen suggests not doing this.
inv_e_cor = cor(inv_e_comm, inv_e_ss, method='kendall')

# Perform NMDS and fit environmental variables to INV epiphytic macrolichens
inv_em_mds = metaMDS(inv_em_comm, distance='bray', k=4, autotransform=F, trymax=300, maxit=500)
inv_em_ef = envfit(inv_em_mds, inv_data[rownames(inv_em_comm),c(xvars,'L2NAME')], choices=1:4)
inv_em_surf = sapply(xvars, function(x){
	mod12 = eval(bquote(with(inv_data[rownames(inv_em_comm),], ordisurf(inv_em_mds, .(as.name(x)), choices=1:2))))
	mod34 = eval(bquote(with(inv_data[rownames(inv_em_comm),], ordisurf(inv_em_mds, .(as.name(x)), choices=3:4))))
	list(mod12=mod12, mod34=mod34)
})
inv_em_ss = scores(inv_em_mds, 'sites', choices=1:4) # Oksanen suggests not doing this.
inv_em_cor = cor(inv_em_comm, inv_em_ss, method='kendall')

# Save NMDS objects for future analyses
#save(fia_mds, inv_e_mds, inv_em_mds, file=file.path(working_dir, 'NMDS_results_subsample.RData'))
load(file.path(working_dir, 'NMDS_results.RData'))

## Compare environmental fits to NMDS across data sets

## Table 2 and S2: Environmental correlations with NMS axes
compare_ef = data.frame(FIA_r2 = c(fia_ef$vectors$r, fia_ef$factors$r), FIA_P = c(fia_ef$vectors$pvals, fia_ef$factors$pvals),
	INV_epi_r2 = c(inv_e_ef$vectors$r, inv_e_ef$factors$r), INF_epi_P = c(inv_e_ef$vectors$pvals, inv_e_ef$factors$pvals),
	INV_em_r2 = c(inv_em_ef$vectors$r, inv_em_ef$factors$r), INF_em_P = c(inv_em_ef$vectors$pvals, inv_em_ef$factors$pvals)
)

# Save the table depending on whether NMDS was performed on a reduced INV data set (subsampled) or not
#write.csv(compare_ef, file.path(working_dir, 'Table_3-compare_envfit_nmds.csv'))
#write.csv(compare_ef, file.path(working_dir, 'Table_S2-compare_envfit_nmds_subsample.csv'))



###########################################
# Not included in manuscript:
# Explore non-linearity through plots
pdf(file.path(fig_dir, 'env_surfaces_nmds_subsample.pdf'), height=9, width=6)
par(mfrow=c(3, 2))
par(mar=c(4,4,1,1))
par(oma=c(0,3,3,0))
for(x in xvars){
	
	# FIA
	plot(fia_mds, display='sites',choices=1:2, las=1)
	plot(fia_surf['mod12',x][[1]], add=T)
	mtext('FIA', 2, 5)
	plot(fia_mds, display='sites',choices=3:4, las=1)
	plot(fia_surf['mod34',x][[1]], add=T)

	# INV EM
	plot(inv_em_mds, display='sites',choices=1:2, las=1)
	plot(inv_em_surf['mod12',x][[1]], add=T)
	mtext('INV Epiphytic Macrolichens', 2, 5)
	plot(inv_em_mds, display='sites',choices=3:4, las=1)
	plot(inv_em_surf['mod34',x][[1]], add=T)

	# INV E
	plot(inv_e_mds, display='sites',choices=1:2, las=1)
	plot(inv_e_surf['mod12',x][[1]], add=T)
	mtext('INV Epiphytic Lichens', 2, 5)
	plot(inv_e_mds, display='sites',choices=3:4, las=1)
	plot(inv_e_surf['mod34',x][[1]], add=T)

	mtext(xvarnames[x], 3, 0, outer=T)
}

dev.off()

# Ordination plots to see how species map on axes
plot(fia_mds, choices=1:2, las=1)

plot(inv_em_mds, display='sites', choices=1:2, las=1)
plot(inv_em_surf['mod12','ppt'][[1]], add=T)
text(inv_em_mds, display='species', col='blue', cex=0.7)

# Anaolsis of variance on distance matrices
fia_ad = sapply(c(xvars, 'L2NAME'), function(x){
	eval(bquote(adonis(fia_comm ~ .(as.name(x)), fia_data, method='bray')))
})
inv_e_ad = sapply(c(xvars, 'L2NAME'), function(x){
	eval(bquote(adonis(inv_e_comm ~ .(as.name(x)), inv_data[rownames(inv_e_comm),], method='bray')))
})
inv_em_ad = sapply(c(xvars, 'L2NAME'), function(x){
	eval(bquote(adonis(inv_em_comm ~ .(as.name(x)), inv_data[rownames(inv_em_comm),], method='bray')))
})



# Multivariate ANOVA
fia_aov = t(sapply(c(xvars, 'L2NAME'), function(x) fia_ad['aov.tab',x][[1]][1,]))
inv_em_aov = t(sapply(c(xvars, 'L2NAME'), function(x) inv_em_ad['aov.tab',x][[1]][1,]))
inv_e_aov = t(sapply(c(xvars, 'L2NAME'), function(x) inv_e_ad['aov.tab',x][[1]][1,]))

cbind(fia_aov[,'R2'], inv_em_aov[,'R2'], inv_e_aov[,'R2'])

#write.csv(fia_aov, file.path(working_dir, 'fia_adonis.csv'))
#write.csv(inv_em_aov, file.path(working_dir, 'inv_em_adonis_subsample.csv'))
#write.csv(inv_e_aov, file.path(working_dir, 'inv_e_adonis_subsample.csv'))


