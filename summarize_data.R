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
library(maptools) # unionSpatialPolygons
library(stringr) # str_trim
library(mgcv) # gam

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


##################################
### Environmental Models


xvars = c('ppt','vpdmax','tmax','tmean','tot_n','tot_s')
xvarnames = expression('Annual precipitation (mm)', 'Max. VPD (kPa)', 'Max. temperature'~(degree~C),
	'Mean temperature'~(degree~C), 'Total N deposition'~(kg~ha^-1), 'Total S deposition'~(kg~ha^-1))
names(xvarnames) = xvars

fia_data = fia_plots@data
inv_data = inv_plots@data
rownames(inv_data) = inv_data$Site_Number

pdf(file.path(fig_dir, 'richness_env.pdf'), height=9, width=18)
layout(matrix(1:(3*length(xvars)), nrow=3, byrow=F))
par(mar=c(1.5,1.5,1,1))
par(oma=c(3,4.5,0,0))
for(x in xvars){
	use_xlim = range(c(fia_data[,x], inv_data[,x]))
	
	plot(fia_data[,x], fia_data$S, las=1, xlim=use_xlim, xlab='', ylab='')
	if(x==xvars[1]) mtext('FIA Num. Species', 2, 3)

	xvals = data.frame(seq(min(fia_data[,x]), max(fia_data[,x]), length.out=100))
	names(xvals)=x
	mod = gam(bquote(S~s(.(as.name(x)))), data=fia_data, family=nb(link='log'))
	lines(xvals[,1], predict(mod, newdata=xvals, type='response'), col=2)

	plot(inv_data[,x], inv_data$S_epi_macro, las=1, xlim=use_xlim, xlab='', ylab='')
	if(x==xvars[1]) mtext('INV Epiphytic Macroclichen\nNum. Species', 2, 3)

	xvals = data.frame(seq(min(inv_data[,x]), max(inv_data[,x]), length.out=100))
	names(xvals)=x
	mod = gam(bquote(S_epi_macro~s(.(as.name(x)))), data=inv_data, family=nb(link='log'))
	lines(xvals[,1], predict(mod, newdata=xvals, type='response'), col=2)	

	plot(inv_data[,x], inv_data$S_epi, las=1, xlim=use_xlim)
	mtext(xvarnames[x], 1, 3)
	if(x==xvars[1]) mtext('INV Epiphyte\nNum. Species', 2, 3)


	mod = gam(bquote(S_epi~s(.(as.name(x)))), data=inv_data, family=nb(link='log'))
	lines(xvals[,1], predict(mod, newdata=xvals, type='response'), col=2)

}
dev.off()


### RDA

## Things to evaluate in RDA
# 1. Explanatory power of each env variable and ecoregion
# 2. Variance partitioning


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

# Remove species occuring fewer than 3 time
dim(fia_comm) #75
fia_comm = fia_comm[,colSums(fia_comm)>2]; dim(fia_comm) # 37 species, lost 38

dim(inv_e_comm) # 263
inv_e_comm = inv_e_comm[,colSums(inv_e_comm)>2]; dim(inv_e_comm) # 146 species, lost 137
inv_e_comm = inv_e_comm[rowSums(inv_e_comm)>0,]; dim(inv_e_comm) # 193 plots, lost 12

dim(inv_em_comm) # 98
inv_em_comm = inv_em_comm[,colSums(inv_em_comm)>2]; dim(inv_em_comm) # 55 species, lost 43
inv_em_comm = inv_em_comm[rowSums(inv_em_comm)>0,]; dim(inv_em_comm) # 185 plots, lost 20

# Run analyses with data subset to the same size as FIA
# NOTE: when saving files, names have _subsample added to them
nplots = nrow(fia_comm)
inv_em_comm = inv_em_comm[sample(nrow(inv_em_comm), nplots),]
inv_e_comm = inv_e_comm[sample(nrow(inv_e_comm), nplots),]

# Calculate Sorenson dissimilarity to view outliers
# None need to be excluded
fact = 2.5

fia_dist = 1-vegdist(fia_dat, method='bray')
which(rowMeans(as.matrix(fia_dist)) > mean(fia_dist) + fact*sqrt(var(fia_dist)))

inv_e_dist = 1-vegdist(inv_e_dat, method='bray')
which(rowMeans(as.matrix(inv_e_dist)) > mean(inv_e_dist) + fact*sqrt(var(inv_e_dist)))

inv_em_dist = 1-vegdist(inv_em_dat, method='bray')
which(rowMeans(as.matrix(inv_em_dist)) > mean(inv_em_dist) + fact*sqrt(var(inv_em_dist)))


fia_mds = metaMDS(fia_comm, distance='bray', k=4, autotransform=F, trymax=300, maxit=500)
fia_ef = envfit(fia_mds, fia_data[,c(xvars,'L2NAME')], choices=1:4)
fia_surf = sapply(xvars, function(x){
	mod12 = eval(bquote(with(fia_data, ordisurf(fia_mds, .(as.name(x)), choices=1:2))))
	mod34 = eval(bquote(with(fia_data, ordisurf(fia_mds, .(as.name(x)), choices=3:4))))
	list(mod12=mod12, mod34=mod34)
})
fia_ad = sapply(c(xvars, 'L2NAME'), function(x){
	eval(bquote(adonis(fia_comm ~ .(as.name(x)), fia_data, method='bray')))
})

fia_ss = scores(fia_mds, 'sites', choices=1:4) # Oksanen suggests not doing this.
fia_cor = cor(fia_comm, fia_ss, method='kendall')


inv_e_mds = metaMDS(inv_e_comm, distance='bray', k=4, autotransform=F, trymax=300, maxit=500)
inv_e_ef = envfit(inv_e_mds, inv_data[rownames(inv_e_comm),c(xvars,'L2NAME')], choices=1:4)
inv_e_surf = sapply(xvars, function(x){
	mod12 = eval(bquote(with(inv_data[rownames(inv_e_comm),], ordisurf(inv_e_mds, .(as.name(x)), choices=1:2))))
	mod34 = eval(bquote(with(inv_data[rownames(inv_e_comm),], ordisurf(inv_e_mds, .(as.name(x)), choices=3:4))))
	list(mod12=mod12, mod34=mod34)
})
inv_e_ad = sapply(c(xvars, 'L2NAME'), function(x){
	eval(bquote(adonis(inv_e_comm ~ .(as.name(x)), inv_data[rownames(inv_e_comm),], method='bray')))
})
inv_e_ss = scores(inv_e_mds, 'sites', choices=1:4) # Oksanen suggests not doing this.
inv_e_cor = cor(inv_e_comm, inv_e_ss, method='kendall')


inv_em_mds = metaMDS(inv_em_comm, distance='bray', k=4, autotransform=F, trymax=300, maxit=500)
inv_em_ef = envfit(inv_em_mds, inv_data[rownames(inv_em_comm),c(xvars,'L2NAME')], choices=1:4)
inv_em_surf = sapply(xvars, function(x){
	mod12 = eval(bquote(with(inv_data[rownames(inv_em_comm),], ordisurf(inv_em_mds, .(as.name(x)), choices=1:2))))
	mod34 = eval(bquote(with(inv_data[rownames(inv_em_comm),], ordisurf(inv_em_mds, .(as.name(x)), choices=3:4))))
	list(mod12=mod12, mod34=mod34)
})
inv_em_ad = sapply(c(xvars, 'L2NAME'), function(x){
	eval(bquote(adonis(inv_em_comm ~ .(as.name(x)), inv_data[rownames(inv_em_comm),], method='bray')))
})
inv_em_ss = scores(inv_em_mds, 'sites', choices=1:4) # Oksanen suggests not doing this.
inv_em_cor = cor(inv_em_comm, inv_em_ss, method='kendall')

# Save NMDS objects for future analyses
#save(fia_mds, inv_e_mds, inv_em_mds, file=file.path(working_dir, 'NMDS_results_subsample.RData'))

## Compare data sets

# Environmental correlations
compare_ef = data.frame(FIA_r2 = c(fia_ef$vectors$r, fia_ef$factors$r), FIA_P = c(fia_ef$vectors$pvals, fia_ef$factors$pvals),
	INV_epi_r2 = c(inv_e_ef$vectors$r, inv_e_ef$factors$r), INF_epi_P = c(inv_e_ef$vectors$pvals, inv_e_ef$factors$pvals),
	INV_em_r2 = c(inv_em_ef$vectors$r, inv_em_ef$factors$r), INF_em_P = c(inv_em_ef$vectors$pvals, inv_em_ef$factors$pvals)
)
#write.csv(compare_ef, file.path(working_dir, 'compare_envfit_nmds_subsample.csv'))

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

# Multivariate ANOVA
fia_aov = t(sapply(c(xvars, 'L2NAME'), function(x) fia_ad['aov.tab',x][[1]][1,]))
inv_em_aov = t(sapply(c(xvars, 'L2NAME'), function(x) inv_em_ad['aov.tab',x][[1]][1,]))
inv_e_aov = t(sapply(c(xvars, 'L2NAME'), function(x) inv_e_ad['aov.tab',x][[1]][1,]))

cbind(fia_aov[,'R2'], inv_em_aov[,'R2'], inv_e_aov[,'R2'])

#write.csv(fia_aov, file.path(working_dir, 'fia_adonis.csv'))
#write.csv(inv_em_aov, file.path(working_dir, 'inv_em_adonis_subsample.csv'))
#write.csv(inv_e_aov, file.path(working_dir, 'inv_e_adonis_subsample.csv'))


#### Species Distributions

# Focus on macrolichen species occuring on at least 6 sites in both datasets
inv_m_comm = inv_siteXsp['macro',,]
inv_m_comm = inv_m_comm[,colSums(inv_m_comm)>=6]
fia_comm = fia_siteXsp[,colSums(fia_siteXsp)>=6]
focal_sp = colnames(fia_comm)[colnames(fia_comm) %in% colnames(inv_m_comm)] # 15 species
focal_names = sapply(focal_sp, function(x) paste(unlist(strsplit(x, ' ')), collapse='_'))

# Add species occurences to dataframes for modeling
fia_sp = fia_comm[,focal_sp]; colnames(fia_sp) = focal_names
fia_data = cbind(fia_data, fia_sp)
inv_sp = inv_m_comm[,focal_sp]; colnames(inv_sp) = focal_names
inv_data = cbind(inv_data, inv_sp)

# Change species names to better variable names

pdf(file.path(fig_dir, 'species_distributions_singlevars.pdf'), height=18, width=18)
par(lend=1)
par(mfrow=c(5, length(xvars)))
par(mar=c(3, 3, 3, 1))
par(oma=c(2, 15, 0, 0))

for(sp in focal_names){
for(x in xvars){
	counter = which(focal_names==sp)

	fia_mod = gam(bquote(.(as.name(sp)) ~ s(.(as.name(x)), k=4)), data=fia_data, family=binomial(link='logit'))
	inv_mod = gam(bquote(.(as.name(sp)) ~ s(.(as.name(x)), k=4)), data=inv_data, family=binomial(link='logit'))

	xrange = range(c(fia_data[,x], inv_data[,x]))
	
	# FIA
	make_plot(xrange, c(0,1), xlab=ifelse(counter%%5==0, xvarnames[x], ''))
	points(inv_data[,x], inv_data[,sp], pch='|', col='#0000ff90')
	xvals = data.frame(seq(min(inv_data[,x]), max(inv_data[,x]), length.out=100))
	names(xvals)=x
	fia_pred = predict(fia_mod, newdata=xvals, type='response', se.fit=T)
	polygon(c(xvals[,1], rev(xvals[,1])), c(fia_pred$fit-fia_pred$se.fit, rev(fia_pred$fit+fia_pred$se.fit)),
		col='#0000ff50', border=NA)
	use_lty = ifelse(summary(fia_mod)$s.table[1,4] < 0.05, 1, 3)
	lines(xvals[,1], fia_pred$fit, col='blue', lwd=2, lty=use_lty)

	# INV
	points(fia_data[,x], fia_data[,sp], pch='|', col='#ff000090')
	xvals = data.frame(seq(min(fia_data[,x]), max(fia_data[,x]), length.out=100))
	names(xvals)=x
	inv_pred = predict(inv_mod, newdata=xvals, type='response', se.fit=T)
	polygon(c(xvals[,1], rev(xvals[,1])), c(inv_pred$fit-inv_pred$se.fit, rev(inv_pred$fit+inv_pred$se.fit)),
		col='#ff000050', border=NA)
	use_lty = ifelse(summary(inv_mod)$s.table[1,4] < 0.05, 1, 3)
	lines(xvals[,1],inv_pred$fit, col='red', lwd=2, lty=use_lty)

	mtext(paste('FIA prop. explained dev:', format(summary(fia_mod)$dev.expl, digits=2)), 3, 1, col='blue', adj=0, cex=0.8)
	mtext(paste('INV prop. explained dev:', format(summary(inv_mod)$dev.expl, digits=2)), 3, 0, col='red', adj=0, cex=0.8)

	if(x==xvars[1]){
		mtext('Prob. occurence', 2, 2.5, cex=0.8)	
		par(xpd=NA)
		text(750, 0.5, labels=focal_sp[counter], adj=1, cex=1.2)
		par(xpd=F)
	}
}
}
dev.off()


################################################################
### Comparison of species richness and composition across space

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

# Assign plots to ecoregions
fia_plots$L2NAME = over(fia_plots, eco_L2)$L2_NAME
inv_plots$L2NAME = over(inv_plots, eco_L2)$L2_NAME

# Add points from FIA and INV plots to maps
use_layout = list(list('sp.points', fia_plots, which=c(1,3), col='deeppink', pch=3, cex=.5),
	list('sp.points',inv_plots, which=c(2,4), col='deeppink', pch=3, cex=.5))


pdf(file.path(fig_dir,'Num records inv vs fia plots.pdf'), height=6, width=4)
par(lend=1)
spplot(eco_L2, c('fiarecs','invrecs'), col='black', sp.layout=use_layout,
	col.regions=blue2red, cuts=length(blue2red)-1, colorkey=list(at=seq(0,5000, 500)),
	strip=strip.custom(factor.levels=c('FIA','INV'), bg='white'))
dev.off()

pdf(file.path(fig_dir,'Num plots inv vs fia plots.pdf'), height=6, width=4)
par(lend=1)
spplot(eco_L2, c('fiaplots','invplots'), col='black', sp.layout=use_layout,
	col.regions=blue2red, cuts=length(blue2red)-1, colorkey=list(at=seq(0,200,20), tick.number=10),
	strip=strip.custom(factor.levels=c('FIA','INV'), bg='white'))
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

## Compare species lists for the entire data set
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


## Barplot across regions
comp_counts_reg = apply(comp_splists_reg, 1:2, function(x) length(unlist(x)))
comp_counts_reg['both',] = comp_counts_reg['both',]/2

# Proportion of species pool in each data set
plot_data = comp_counts_reg[c('inv','both','fia'),]
prop_sp = t(plot_data)/colSums(plot_data)
comp_counts = sapply(comp_splists, function(x) length(unlist(x)))
comp_counts['both'] = comp_counts['both']/2
prop_sp = rbind(prop_sp, ALL=comp_counts[c('inv','both','fia')] / sum(comp_counts[c('inv','both','fia')]))

write.csv(prop_sp, file.path(working_dir, 'proportion_species_across_datasets.csv'))

pdf(file.path(fig_dir, 'compare_species_cross_region.pdf'), height=3, width=8.5)
par(mar=c(4.5,19,1,1.5))
barplot( comp_counts_reg[c('inv','both','fia'),], horiz=T, las=1, xlim=c(0,100),
	legend.text=c('INV','BOTH','FIA'), xlab='Num. Species', 
	args.legend=list(x='topright', horiz=T, bty='n'))
dev.off()

## Write out species list tables with indication of which regions found in

# INV-only species
inv = data.frame(All=comp_splists$inv)
inv_reg = by(inv_siteXsp['macro',,rownames(inv)], inv_plots$L2NAME, colSums)
inv = cbind(inv, sapply(regions, function(reg) inv_reg[[reg]]))		

# Which ones are targeted by FIA (but not found)?
setwd(working_dir)
fia_ref = read.csv('REF_LICHEN_SPECIES.csv')
fia_ref$Binomial = str_trim(paste(fia_ref$GENUS, fia_ref$SPECIES))

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

# Found in both lists, include Chi-sq

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

# Write out tables
write.csv(inv, file.path(working_dir, 'INV_only_species.csv'))
write.csv(fia, file.path(working_dir, 'FIA_only_species.csv'))
write.csv(both, file.path(working_dir, 'INV_and_FIA_species.csv'))


# Which ones are found in 2 or fewer plots?
# annotated in saved tables



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


## Plotting species lists across regions NOT DISPLAYING WELL. 
par(mfrow=c(2,2))
par(mar=c(2,4,1,1))

for(reg in regions){
	plot_data = comp_counts_reg[c('inv','both','fia'),reg, drop=F]
	bp = barplot(plot_data, las=1, xlim=c(0,6), ylab='Num. Species',
		ylim=c(0,sum(plot_data)+6), names.arg='')
	abline(v=par('usr')[1])
	abline(h=par('usr')[3])

	inv = unlist(comp_splists_reg['inv',reg])
	print_inv = names(inv)[names(inv) %in% names(comp_splists[['inv']])]
	print_inv = paste(print_inv, ' (', inv[print_inv], ')', sep='')
	
	fia = unlist(comp_splists_reg['fia',reg])
	print_fia = names(fia)[names(fia) %in% names(comp_splists[['fia']])]
	print_fia = paste(print_fia, ' (', fia[print_fia], ')', sep='')

	print_fia = paste(print_fia, collapse='\n')
	text(1.5, sum(plot_data[c('inv','both'),]) + 0.5*plot_data['fia',], print_fia, pos=4) 

	if(length(print_inv) > plot_data['inv',]/2){
		print_inv1 = print_inv[1:((length(print_inv)+1)/2)]
		print_inv2 = print_inv[(length(print_inv1)+1):length(print_inv)]
		print_inv1 = paste(print_inv1, collapse='\n')
		print_inv2 = paste(print_inv2, collapse='\n')

		text(1.5, 0.5*plot_data['inv',], print_inv1, pos=4)
		text(4, 0.5*plot_data['inv',], print_inv2, pos=4)	

	} else {
		print_inv = paste(print_inv, collapse='\n')
		text(1.5, 0.5*plot_data['inv',], print_inv, pos=4)
	}

	mtext(reg, 3, 1)
}



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












