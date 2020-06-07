# fishsuit
Framework to model climate change and dam-driven fragmentation impacts on the geographical ranges of freshwater fish species

This is a beta release. More details on input data and streamline of processes will be made available upon acceptance of the related manuscript. For further information, please contact Valerio Barbarossa (v.barbarossa@fnwi.ru.nl)

# R Scripts

## Preprocessing

* scripts/R/pre/pcrglobwb_hist.R: preprocess pcrglobwb historical runs to obtain one long-term layer

* scripts/R/pre/pcrglobwb_fut.R: preprocess pcrglobwb future scenario runs to obtain one long-term layer at each warming target-RCP-GCM combination

* scripts/R/preprocess/species2points.R: this script transforms the polygon ranges of the species to spatial points based on the resolution of the environmental variables, to use a matrix-like framework to model range changes

## Modelling

* scripts/R/mod/map_thresholds.R

* cripts/R/mod/model_occurrence.R

## Extras

* scripts/R/xtr/retrieve_habitat_iucn.R: uses the API from IUCN to retrieve info on species habitat

* scripts/R/xtr/analyze_habitat_iucn.R: compiles info from retrieve_habitat_iucn into a table

* scripts/R/xtr/calc_global_mean_air_temp.R: NEED TO CHECK IF THIS IS THE ONE USED IN THE FINAL MODEL

* scripts/R/xtr/compile_species_ranges_dataset: merge datasets of fish ranges from IUCN, Barbarossa et al. 2020 (PNAS), AmazonFISH into one lobal dataset with polygon vertices based on a 5 arcmin grid

* scripts/R/xtr/compile_traits.R: retrieve info on species' traits, range area, climate zone and compile them in a table

* scripts/R/xtr/thresholds_selection.R: analyze thresholds mapped for the species and check multicollinearity

## Functions

* scripts/R/fun/generic.R: generic functions

* scripts/R/fun/map_variable2species.R: function to retrieve quantiles of an environmental variable for a species range

* scripts/R/fun/range2table.R: function that convert polygon ranges to single spatial points

