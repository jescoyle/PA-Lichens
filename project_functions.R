### This script contains functions used by the PA-Lichens project


## A function that uses a dictionary to assign a substrate to a record
categorize_substrate = function(x, words){
	found_words = words[sapply(1:nrow(words), function(i) length(grep(words[i,1], x, ignore.case=T))>0),]
	unique(found_words$substrate)
} 
