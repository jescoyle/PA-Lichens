# This script compiles FIA data from Pennsylvania: 1999 - most recent plots
# 1999 plots established under FHM program and follow a different format from FIA program plots
# This script combines these two data tables. Coordinates come from a non-public source used in previous project.
# It also fixes told taxonomic nomeclature based on _______________


options(stringsAsFactors=F)
library(dplyr)


## Define directories
data_dir = 'Data/raw'
derived_dir = 'Data/derived'
code_dir = 'Code'
working_dir = './'

## Load functions
source(file.path(code_dir, "project_functions.R"))

## Load raw data tables
fia_lichen <- read.csv(file.path(data_dir, "PA_LICHEN_LAB.CSV"))
fia_plots <- read.csv(file.path(data_dir, "PA_PLOT.CSV"))
fhm_lichen <- read.csv(file.path(data_dir, "PA_LICHEN_ABUNDANCE_1999.CSV"))
fhm_plots <- read.csv(file.path(data_dir, "PA_PLOT_1999.CSV"))

## Load prior data with coordinates
all_fia_plots <- read.csv(file.path(data_dir, "AllLichenPlots_forJES.csv"))


######### Combine FHM and FIA Data tables #############

# Rename FHM columns to conform to FIA columns
fhm_plots <- rename(fhm_plots, 
                    STATECD = STATE,
                    COUNTYCD = COUNTY,
                    MEASYEAR = YEAR)
fhm_lichen <- rename(fhm_lichen,
                     STATECD = STATE,
                     COUNTYCD = COUNTY,
                     PLOT = P3ID,
                     LICH_SPPCD = LICHEN_SPECIES_CODE)
fhm_lichen$MEASYEAR = 1999 # because we only are using data from this year from FHM                     

# Add coordinates to FHM data
fhm_plots <- left_join(fhm_plots, all_fia_plots[, c("STATECD", "COUNTYCD", "P3ID", "FZ_LAT", "FZ_LON", "FZ_ELEV")])
fhm_plots <- rename(fhm_plots,
       LAT = FZ_LAT,
       LON = FZ_LON,
       ELEV = FZ_ELEV, 
       PLOT = P3ID)

# Add plot id to each table
fhm_plots$yrplotid <- make.yrplotid(fhm_plots)
fhm_lichen$yrplotid <- make.yrplotid(fhm_lichen)
fia_plots$yrplotid <- make.yrplotid(fia_plots)
fia_lichen$yrplotid <- make.yrplotid(fia_lichen)

# Combine data from FHM and FIA
keep_lichen_cols <- c("yrplotid", "MEASYEAR", "STATECD", "COUNTYCD", "PLOT", "LICH_SPPCD", "ABUNDANCE_CLASS")
lichen <- bind_rows(fhm_lichen[, keep_lichen_cols], fia_lichen[, keep_lichen_cols])

keep_plot_cols <- c("yrplotid", "MEASYEAR", "STATECD", "COUNTYCD", "PLOT", "LAT", "LON", "ELEV", "QA_STATUS")
plots <- bind_rows(fhm_plots[, keep_plot_cols], fia_plots[, keep_plot_cols])


# Match lichen data to all fia plot data to determine whether any missing
lichen_plots <- unique(lichen[, c("yrplotid", "STATECD")])
lichen_plots <- left_join(lichen_plots, unique(plots))

# Check whether plot coordinates missing
missing_plots <- subset(lichen_plots, is.na(LAT)) # none

# Save FIA plot data
#write.csv(lichen_plots, file.path(data_dir, "PA_lichen_plots.csv"), row.names = FALSE)


######### Correct taxonomic names changes in FIA #############

# Add column to indicate taxonomy change
lichen$taxon_change <- NA

# Read in REF_LICHEN_SPP_COMMENTS file explaining FIA name changes
ref_comments <- read.csv(file.path(data_dir, "REF_LICHEN_SPP_COMMENTS_PlainTextFINALUPDATES.csv"))

# Which species in this data set need to be fixed?
these_sp <- unique(lichen$LICH_SPPCD)
fix_sp <- unique(ref_comments$LICH_SPPCD)

these_sp[these_sp %in% fix_sp]

## 601: ACTION 2: 601 Bryoria abbreviata is a synonym and should be combined with 4551 Nodobryoria abbreviata for all analyses.
lichen[lichen$LICH_SPPCD == 601, 'taxon_change'] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 601, 'LICH_SPPCD'] <- 4551

## 1002: ACTION 0: 1002 Cetraria aurescens is a synonym of Tuckermannopsis aurescens. 
#        FIA never adopted the latter name due to lack of scientific consensus. 
#        This species is now named 1002 Nephromopsis aurescens (Divakar et al. 2017).
# do nothing

## 1012: ACTION 0: 1012 Cetraria oakesiana was re-named to 1012 Usnocetraria oakesiana (Thell et al. 2009), 
#        then re-named back to Cetraria oakesiana (Divakar et al. 2017).
# do nothing

## 1200: ACTION 0: If a specialist finds a small Cladina-like specimen unidentifiable to species, it is placed in the 1180 Cladina-form category. 
#        Ahti and DePriest (2001) shifted all Cladina species to Cladonia.
# do nothing

## 1203: ACTION 0: Cladonia bacillaris = Cladonia macilenta var. bacillaris. The former name is retained for brevity.
# lump to match INV data
lichen[lichen$LICH_SPPCD == 1203, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 1203, "LICH_SPPCD"] <- 1225

## 1210: ACTION 0: Species in the C. chlorophaea morphological group that can't be distinguished without TLC are recorded as 1210 C. chlorophaea.
# do nothing

## 1211 & 1228: ACTION 5/ACTION 2: WEST - 1211 Cladonia coniocraea should be combined into 1228 C. ochrochlora. 
#        EAST - 1228 C. ochrochlora should be combined into 1211 C. coniocraea, which is the more common taxon in the east. 
#       See (Pino-Bodas et al. 2011) - C. ochrochlora & C. coniocraea may be conspecific.
lichen[lichen$LICH_SPPCD == 1228, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 1228, "LICH_SPPCD"] <- 1211

## 2702 & 2704: ACTION 0: Some material is intermediate between Flavopunctelia flaventior and F. soredica 
#       (i.e. with only marginal soralia but also numerous pseudocyphellae).  
#       Starting with 1998 data, we follow Hale's heavier weighting of pseudocyphellae and class these intermediates with F. flaventior.  
#       For pre-1998 data, species names were assigned with heavier weighting on soralia characteristics.
# do nothing

## 2901: ACTION 0: For analyzing data from multiple years crossing 2012, 2904 H. confusa should be combined into 2901 H. adglutinata. 
#        2901 H. adglutinata s. lat. has been shown to include two distinct species with broad overlap - H. adglutinata s. str. and H. confusa (Esslinger et al. 2012). 
#        Records for H. adglutinata from 2012 and prior may include some H. confusa.
# do nothing: already matches INV data

## 3000: ACTION 1: Exclude for most analyses. Squamulose/crustose growth forms were not consistently collected.
lichen <- lichen[lichen$LICH_SPPCD != 3000, ]

## 3218: ACTION 0: The species concept for Hypotrachyna showmanii changed in 2006 (Lendemer and Harris, 2006). 
#        Older specimens identified as 5103 Parmelinopsis spumosa were reexamined in 2006-07 and re-assigned to H. showmanii as appropriate.
# do nothing

## 4000: ACTION 0: Starting in 2005, specimens identified to species were assigned to either Melanoelixia or Melanohalea (Blanco et al. 2004). 
#        All specimens too small to identify to species will remain 4000 Melanelia sp.
# do nothing

## 4015: ACTION 0: 4015 Melanelia subaurifera was re-named to 4015 Melanelixia subaurifera (Blanco et al. 2004).
# do nothing    

## 4101: ACTION 3: For analyzing data for multiple years crossing 2009, 4101 Menegazzia terebrata should be combined into 4102 M. subsimilis. 
#        Specimens identified as M. terebrata before 2009 were most likely M. subsimilis, which is the more common species (Bjerke 2003).
lichen[lichen$LICH_SPPCD == 4101, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 4101, "LICH_SPPCD"] <- 4102

## 4806: ACTION 3: For analyzing data from multiple years crossing 2010, 4808 Parmelia barrenoae should be combined into 4806 P. sulcata. 
#        In 2010, P. barrenoae was split from P. sulcata  (Hodkinson et al. 2010).
# do nothing

## 5101: ACTION 0: 5101 Parmelinopsis horrescens was re-named to 5101 Hypotrachyna horrescens (Divakar et al. 2013).
# do nothing

## 5300: ACTION 0: 770 Canomaculina should be combined into 5300 Parmotrema. Canomaculina was renamed to Parmotrema (Blanco et al. 2005).
# do nothing


## 5303: ACTION 0: 5303 Parmotrema chinense was re-named to 5303 Parmotrema perlatum (Hawksworth 2004).
# do nothing

## 5702: ACTION 0: The name Physcia aipolia is applied sensu lato. 
#       ACTION 0: For analyzing data for multiple years crossing 2014, 
#       5728 P. alnophila should be combined into 5702 P. aipolia. 
#       Molecular data justify recognition of P. aipolia var. alnophila as a 
#       distinct species, P. alnophila (Lohtander et al. 2009).  
#       ACTION 5/ACTION 0: ALASKA -  the name P. alnophila is applied for 
#       FIA data from Alaska without requiring TLC. 
#       ACTION 5/ACTION 0: OUTSIDE ALASKA - P. alnophila is applied only if 
#       identified using TLC, not normally done for FIA specimens. 
#       P. alnophila has a more northerly distribution than P. aipolia although
#       there is consider geographic overlap in the lower 48 states and intermediate 
#       forms of these species can not be reliably separated without TLC (Brodo et al. 2013). 
#       ACTION 0: P. aipolia and P. stellaris are distinguished solely on the K reaction of the medulla; 
#       other characters are not reliably correlated with the K reaction. 
#       ACTION 0: SOUTHEAST U.S.- Based on New York Botanical Garden herbarium records (http://sweetgum.nybg.org/science/vh/), 
#       Physcia aipolia may not occur in the most southeastern USA states. Re-evaluation of
#       specimens in other herbaria will be required to confirm this. 
#       See comment for Physcia pumilior.
# do nothing

## 5723: ACTION 0: Physcia aipolia and P. stellaris are distinguished solely on the K reaction of the medulla; 
#        other characters are not reliably correlated with the K reaction. 
#        P. biziana intergrades with P. stellaris, both having K- medulla. 
#        Moderately pruinose specimens with short, rounded, scalloped lobes were named P. biziana; 
#        moderately pruinose specimens with narrower lobes of P. aipolia type were named P. stellaris.
#        Specimens with little or no pruinosity were named P. stellaris regardless of lobe size.
# do nothing

## 5801: ACTION 0: In arid W habitats, small specimens of 5605 Phaeophyscia hirsuta, 5611 P. nigricans, 5711 Physcia dubia, 
#        and 5801 Physciella chloantha may be morphologically indistinguishable. 
#        FIA examines cells of the lower cortex to a limited extent in order to assign accurate 
#        abundance codes.
# do nothing

## 5901: ACTION 5/ACTION 3: WEST- 5907 Physconia isidiigera, 5906 P. perisidiosa, and 5911 P. leucoleiptes 
#        were recognized for the west starting in 1998. Prior to the split these species were identified as 
#       5901 P. detersa, which is rare in the west (Brodo et al. 2016). For analysis of data crossing 1998, 
#       5901 P. detersa should be combined into 5906 P. perisidiosa, the most common species in the group. 
#       ACTION 5/ACTION 0: EAST - P. detersa, P. isidiigera, P. leucoleiptes, and P. perisidiosa have been 
#       distinguished in all inventory years.
# do nothing

## 6705: ACTION 0: 6712 Punctelia punctilla should be combined into 6705 P. missouriensis for all analyses. 
#        P. missouriensis is the correct name for this taxon (Aptroot 2003, Lendemer & Hodkinson 2010).
lichen[lichen$LICH_SPPCD == 6712,] # none
# do nothing

## 6706: ACTION 0: Lendemer & Hodkinson (2010) narrowed the species concept for P. perreticulata. 
#       ACTION 5/ACTION 2: WEST- 6706 P. perreticulata should be combined into 6714 P. jeckeri 
#       for all analyses. ACTION 5/ACTION 4: EAST -for data collected before 2010, combine 
#       6706 P. perreticulata into 6713 P. caseana for all analyses. It is likely that eastern 
#       specimens collected outside the Ozarks are P. caseana.
lichen[lichen$LICH_SPPCD == 6706, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 6706, "LICH_SPPCD"] <- 6713

## 6711: ACTION 0: Punctelia subrudecta is no longer a valid name for USA specimens (Lendemer & Hodkinson 2010). 
#        ACTION 5/ACTION 2: EAST- 6711 Punctelia subrudecta should be combined into 6713 P. caseana for all analyses.  
#        ACTION 5/ACTION 2: WEST - 6711 P. subrudecta should be combined into 6714 P. jeckeri for all analyses.
lichen[lichen$LICH_SPPCD == 6711, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 6711, "LICH_SPPCD"] <- 6713

## 6803: ACTION 5/ACTION 0: Southern (S) FIA Region  -  The names 6801 Pyxine albovirens, 6803 P. caesiopruinosa,
#        and 6809 P. subcinerea may have been misapplied in early years in the SE. ACTION 5/ACTION 4:  For any analysis 
#       including Southern data collected in 1999 or prior, P. albovirens and P. caesiopruinosa should be combined into 
#       P. subcinerea.
lichen[lichen$LICH_SPPCD == 6803, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 6803, "LICH_SPPCD"] <- 6809

## 7100: ACTION 0: 7100 Rimelia was re-named to 5300 Parmotrema (Blanco et al. 2005).
lichen[lichen$LICH_SPPCD == 7100, "taxon_change"] <- "REF_FIA"
lichen[lichen$LICH_SPPCD == 7100, "LICH_SPPCD"] <- 5300

## 8000: do nothing

## 8301: ACTION 0: When analyzing data from multiple years crossing 2002, 8303 Candelaria pacifica should be combined into 8301 C. concolor.  
#        This is a distinct species segregated from 8301 C. concolor in 2002 (Westberg & Nash 2002)
#        but formally described by Westberg & Arup (2011). ACTION 0: C. concolor data 
#        collected before 2002 likely includes some C. pacifica.
# do nothing


## Save
write.csv(lichen, file.path(derived_dir, "fia_lichens.csv"), row.names = FALSE)

