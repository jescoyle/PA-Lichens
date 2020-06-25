## This script loads the data for the FIA vs Inventory (PA-Lichens) project.

options(stringsAsFactors=F)

## Define directories
data_dir <- 'Data/raw'
derived_dir <- 'Data/derived'
code_dir <- 'Code'
working_dir <- './'
fig_dir <- 'Figures'

## Load functions
source(file.path(code_dir, 'project_functions.R'))

## Read in data

inv_lichen = read.csv(file.path(derived_dir, 'inv_lichens.csv'))
fia_lichen = read.csv(file.path(derived_dir, 'fia_lichens_matchINV.csv'))

inv_plots = read.csv(file.path(derived_dir, 'inv_plots.csv'))
fia_plots = read.csv(file.path(derived_dir, 'fia_plots.csv'))
