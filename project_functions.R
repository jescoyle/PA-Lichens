### This script contains functions used by the PA-Lichens project


## A function that uses a dictionary to assign a substrate to a record
categorize_substrate = function(x, words){
	found_words = words[sapply(1:nrow(words), function(i) length(grep(words[i,1], x, ignore.case=T))>0),]
	unique(found_words$substrate)
} 

## A function that returns a list of unique species and drops genera already present in the species list
# recs	= a vector species binomials: Genus species
get_species = function(recs){
	recs_split = strsplit(recs, ' ')
	
	# Find records with only a genus records
	just_genera = unique(recs[sapply(recs_split, length)==1])
	
	# Find records with species epithet
	species = unique(recs[sapply(recs_split, length)==2])

	# Keep genera in records with genus only if not already represented in species list
	sp_genera = unique(sapply(strsplit(species, ' '), function(sp) sp[1]))
	keep_genera = just_genera[!(just_genera %in% sp_genera)]
	
	# Return unique species list
	c(species, keep_genera)
}


## A function that determines whether a taxon is present in a vector of species names.
# If the taxon is just a genus, then the function returns FALSE only if there are no representative species in that genus
find_species = function(x, splist){

	if(length(splist) > 0){
		x_split = strsplit(x, ' ')[[1]]
	
		if(length(x_split)==1){
			splist_gen = unique(sapply(strsplit(splist, ' '), function(sp) sp[1]))
			found = x_split %in% splist_gen
		} else {
			found = x %in% splist
		}
	} else {
		found = F
	}
	found
}



## A function for making a pretty plot
make_plot = function(xlim, ylim, xlab=NULL, ylab=NULL, cex=1, xlab_loc=2.5, ylab_loc=3){	
	plot.new()
	plot.window(xlim=xlim, ylim=ylim)
	axis(1)
	abline(h=par('usr')[3], lwd=3)
	axis(2, las=1)
	abline(v=par('usr')[1], lwd=3)
	if(!is.null(xlab)) mtext(xlab, 1, xlab_loc, cex=cex)
	if(!is.null(ylab)) mtext(ylab, 2, ylab_loc, cex=cex)

}

# Capitalize firt letter of each word
# Function from R help
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# Adds a letter labeling a panel to the current figure in the upper left corner
add_panel_label = function(which_panel, x=0, y=1, cap=TRUE, add_sym='', cex=1.5){
	use_letter = letters[which_panel]
	if(cap) use_letter = toupper(use_letter)
	
	use_letter = paste0(use_letter, add_sym)
	
	# Get figure dimensions
	usr = par('usr')
	deltaX = usr[2]-usr[1]
	deltaY = usr[4]-usr[3]
	plt = par('plt')
	fracX = plt[2]-plt[1]
	fracY = plt[4]-plt[3]
	omd = par('omd')
	
	xpos = (deltaX/fracX)*(x - plt[1]) + usr[1]
	ypos = (deltaY/fracY)*(y - plt[3]) + usr[3]
	
	# Add text	
	xpd_old = par('xpd')
	par(xpd=NA)
	text(xpos, ypos, use_letter, adj=c(0,1), cex=cex)
	par(xpd=xpd_old)
}

# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
# mycex   : size of label and title text
# uneven.lab : TRUE when numeric labels are not evenly spaced along length of color bar
# labrange : use when uneven.lab=T to specify minimum and maximum values of color range c(min, max)
# ndig    : number of digits beyond decimal place to print for numeric labels 
plotColorRamp = function(cols, n, barends, labels=NULL, title=NA, mycex=1.5, uneven.lab = F, labrange = NA, ndig=1, horiz=FALSE, side=0){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	dx = dX/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	lend.old = par('lend')
	par(lend=1)

	usecols = colorRampPalette(cols)(n)

	if(horiz){
		for(i in 1:n){
			rect(barends[1]+dx*(i-1), barends[2], barends[1]+dx*i, barends[4], col=usecols[i], border=usecols[i])
		}	
	} else {
		for(i in 1:n){
			rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=usecols[i])
		}
	}

	if(exists('labels')){
		# Determine side to label
		if(side==0 & horiz) side=1
		if(side==0 & !horiz) side=4
	
		if(is.numeric(labels)){
			labels.round = format(round(labels, ndig), nsmall=ndig, trim=F)

			if(uneven.lab){
				dz = ifelse(horiz, dX/diff(labrange), dY/diff(labrange))
				Yposition = barends[2] + dz*(labels-labrange[1])
				Xposition = barends[1] + dz*(labels-labrange[1])
			} else {
				dZ = labels[length(labels)]-labels[1]
				dz = ifelse(horiz, dX/dZ, dY/dZ)
				Yposition = barends[2] + dz*(labels-labels[1])
				Xposition = barends[1] + dz*(labels-labels[1])
			}

		} else {
			labels.round = labels
			dz = dY/length(labels)
			Yposition = barends[2] + dz*(0:(length(labels)-1))
			Xposition = barends[1] + dz*(0:(length(labels)-1))
		}
		
		if(horiz){
			text(Xposition, ifelse(side==1, barends[2]-dY*0.5, barends[4]+dY*0.5), labels.round, pos=side, cex=mycex)
			segments(Xposition, ifelse(side==1, barends[2], barends[4]), Xposition, ifelse(side==1, barends[2]-dY*0.5, barends[4]+dY*0.5))
		} else {
			text(ifelse(side==4, barends[3]+dX*0.5, barends[1]-dX*0.5), Yposition, labels.round, pos=side, cex=mycex)
			segments(ifelse(side==4, barends[3], barends[1]), Yposition, ifelse(side==4, barends[3]+dX*0.5, barends[1]-dX*0.5), Yposition)
		}
			
	}
	if(!is.na(title)){
			
		## Determine how many characters away to place title
		digits = max(nchar(labels.round)) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		
		small.chars = grep('[-.]', largest) # Does the largest label have a small character?
		
		if(length(small.chars)==length(largest)) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		if(horiz){
			text(barends[1]+dX*0.5, ifelse(side==1, barends[2]-dY*0.5-par('cxy')[2]*mycex*2, barends[4]+dY*0.5+par('cxy')[2]*mycex*2), labels=title, cex=mycex)
		} else {
			text(ifelse(side==4, barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[1]-dX*0.5-par('cxy')[1]*mycex*(digits+.5)), barends[2]+0.5*dY, labels=title, srt=ifelse(side==4, -90, 90), cex=mycex)
		}
	}
	par(xpd=xpd.old)
	par(lend=lend.old)
}
	