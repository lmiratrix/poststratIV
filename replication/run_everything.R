
# Run all the scripts to generate all results in paper and paper appendix


# Load libraries

library( tidyverse )
library( poststratIV )
library( cli )

# Set to true to run the simulations
RUN_SIMULATIONS = FALSE

# Number of iterations per scenario
R = 1000
M_CHUNK = 10

# Run all the simulations
if ( RUN_SIMULATIONS ) {
    ONE_SIDED_SIMULATION = FALSE
    source("simulation_full_twosided.R")
    source("simulation_pred_C_twosided.R")

    ONE_SIDED_SIMULATION = TRUE
    source("simulation_full_twosided.R")
    source("simulation_pred_C_twosided.R")

    cat( "Finished with the double simulation scripts\n" )

    source("simulation_exclusion_restriction.R")
    source("simulation_random_stratification.R")
}

# Now make the figures and results
cat( "\n\nMoving to analysis files\n\n" )

source("analysis_main_simulation.R")
source("analysis_predictC_simulation.R")
source("analysis_exclusion.R")
source("analysis_appendix_se_instability_plot.R")
source("analysis_random_stratification.R")
