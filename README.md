# fishsuit
Framework to model climate change impacts on the geographical ranges of freshwater fish species

*Please, refer to the publication for additional details:* 
*Barbarossa, V., Bosmans, J., Wanders, N. et al. Threats of global warming to the world’s freshwater fishes. Nat Commun 12, 1701 (2021). https://doi.org/10.1038/s41467-021-21655-w*

# R Scripts

* config.R/config_local.R: global parameters used in the R scripts

* exec.R: generates batch file for streamlined execution to be run on a linux slurm cluster

## Preprocessing

* scripts/R/pre/pcrglobwb_hist.R: preprocess pcrglobwb historical runs to obtain one long-term layer

* scripts/R/pre/pcrglobwb_fut.R: preprocess pcrglobwb future scenario runs to obtain one long-term layer at each warming target-RCP-GCM combination

* scripts/R/preprocess/species2points.R: this script transforms the polygon ranges of the species to spatial points based on the resolution of the environmental variables, to use a matrix-like framework to model range changes

* scripts/R/preprocess/species2points_dispersal.R: this script transforms the polygon ranges of the species to spatial points based on the resolution of the environmental variables, to use a matrix-like framework to model range changes

## Modelling

* scripts/R/map_thresholds.R: mapping of species-specific thresholds

* scripts/R/model_occurrence.R: modelling of species-specific percentage of threatened range under no dispersal assumption

* scripts/R/model_occurrence_dispersal.R: DEPRECATED 

* scripts/R/model_occurrence_dispersal2.R: modelling of species-specific percentage of threatened range under maximal dispersal assumption

## Phyloenetic regression

* scripts/R/res/phyloreg_stochastic_array.R: runs phylogenetic regression based on stochastically generated trees

* scripts/R/res/phyloreg_stochastic.R: tabulates results of coefficients and variable importance for the phylogenetic regression as described in

* scripts/R/res/phyloreg.R: DEPRECATED

* scripts/R/res/phyloreg_array.R: DEPRECATED

* scripts/R/res/phyloreg_explore.R: DEPRECATED

## Extras

* scripts/R/xtr/retrieve_habitat_iucn.R: uses the API from IUCN to retrieve info on species habitat

* scripts/R/xtr/analyze_habitat_iucn.R: compiles info from retrieve_habitat_iucn into a table

* scripts/R/xtr/calc_global_mean_air_temp.R: NEED TO CHECK IF THIS IS THE ONE USED IN THE FINAL MODEL

* scripts/R/xtr/compile_species_ranges_dataset: merge datasets of fish ranges from IUCN, Barbarossa et al. 2020 (PNAS), AmazonFISH into one lobal dataset with polygon vertices based on a 5 arcmin grid

* scripts/R/xtr/compile_traits.R: retrieve info on species' traits, range area, climate zone and compile them in a table

* scripts/R/xtr/thresholds_selection.R: analyze thresholds mapped for the species and check multicollinearity

* scripts/R/xtr/extract_climate_zones.R: assign a climate zone to each species

* scripts/R/xtr/reference_species_ranges2hybas12.R: DEPRECATED

## Functions

* scripts/R/fun/generic.R: generic functions

* scripts/R/fun/HighstatLibV10.R: function to calculate VIFs from Highland Statistics LTD

* scripts/R/fun/map_variable2species.R: function to retrieve quantiles of an environmental variable for a species range

* scripts/R/fun/range2table.R: function that convert polygon ranges to single spatial points

* scripts/R/fun/variable_importance.R: function to calculate variable importance

## Figures

*Figure numbers refer to Barbarossa et al., 2021. Threats of global warming to the world's freshwater fishes. Nature Communications.*

* scripts/R/fig/compare_CTmax.R: Figure S7

* scripts/R/fig/compare_SR_basin_w_tedesco.R: Figure S10

* scripts/R/fig/pcrglobwb_anomaly.R: Figure S5

* scripts/R/fig/radarplot_order.R: Figure S6

* scripts/R/fig/RC_basins_barplots.R: Figures 3 and S4

* scripts/R/fig/RC_boxviolins.R: Figures 1 and S1

* scripts/R/fig/RC_overall_maps: Figures 2, 4, S2, S3, S9

* scripts/R/fig/thresholds_violin_plots.R: Figure S8

# Additional files

* thresholYears_4targets.csv: year which the specified warming target is reached at for the different GCM-RCP combinations

* several batch files used for slurm execturion of the scripts
