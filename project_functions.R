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


