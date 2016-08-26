## This script loads the data for the FIA vs Inventory (PA-Lichens) project.

options(stringsAsFactors=F)

## Define directories
data_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Raw_data/'
derived_dir = 'C:/Users/jrcoyle/Dropbox/Pennsylvania_FIA_vs_inventory/Derived_data/'
git_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/GitHub/PA-Lichens/'
working_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/'
fig_dir = 'C:/Users/jrcoyle/Documents/Research/PA-Lichens/Figures/'

setwd(working_dir)

## Load functions
source(file.path(git_dir, 'project_functions.R'))

## Read in data

inv_lichen = read.csv(file.path(derived_dir, 'inv_lichens.csv'))
fia_lichen = read.csv(file.path(derived_dir, 'fia_lichens.csv'))

inv_plots = read.csv(file.path(derived_dir, 'inv_plots.csv'))
fia_plots = read.csv(file.path(derived_dir, 'fia_plots.csv'))
