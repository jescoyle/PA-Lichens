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
