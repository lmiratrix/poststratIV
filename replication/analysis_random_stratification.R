


##### Analyze the random stratification simulation results #####

cat( "\n\nAnalyzing the 'random stratification' simulation\n" )

library( tidyverse )


# Where to load results from
results_file = here::here( "results/random_strat_sim.rds" )

dtT = readRDS( file = results_file )

THRESHOLD <- 100

dtT

# Calculate how many simulation trials there were
R = dtT %>%
    dplyr::filter( num_strata==1, nt_shift==-1, one_sided==TRUE, method=="UNSTRAT" ) %>%
    nrow()
R

figures_dir = here::here( "figures/" )
if ( !dir.exists( figures_dir ) ) {
    dir.create( figures_dir )
}


mean( is.na( dtT$LATE.hat ) )

sum <- dtT %>%
    filter( !is.na( LATE.hat ) ) %>%
    mutate( LATE.hat = pmin( pmax( LATE.hat, -THRESHOLD ), THRESHOLD ) ) %>%
    group_by(  num_strata, nt_shift, one_sided, method ) %>%
    summarise( pextreme = mean( abs( LATE.hat ) == THRESHOLD ),
               Eest = mean( LATE.hat, na.rm = TRUE ),
               medEst = median( LATE.hat, na.rm = TRUE ),
               SE = sd( LATE.hat, na.rm=TRUE ),
               RMSE = sqrt( mean( (LATE.hat - LATE)^2, na.rm = TRUE ) ),
               bias = mean( LATE.hat - LATE ),
               LATE = mean( LATE ),
               IQR = quantile( LATE.hat, 0.75, na.rm = TRUE ) - quantile( LATE.hat, 0.25, na.rm = TRUE ),
               SE_ITT = sd( ITT.hat ),
               SE_f = sd( pi.hat ),
               pdrop = sum( empty > 0 ) / R,
               pNA = 1 - n() / R, #mean( is.na( LATE.hat ) ),
               n = n(), .groups = "drop" )
sum

summary( sum$pdrop )

# Density of LATE estimates
dtT %>%
    filter( !is.na( LATE.hat ) ) %>%
    mutate( LATE.hat = pmin( pmax( LATE.hat, -THRESHOLD ), THRESHOLD ) ) %>%
    ggplot( aes( LATE.hat, col = method ) ) +
    facet_grid( nt_shift ~ num_strata ) +
    geom_density()


# Dropped strata plot
summary( sum$pdrop )
ggplot( sum, aes( num_strata, pdrop ) ) +
    facet_grid( one_sided ~ nt_shift, labeller = label_both ) +
    geom_line() + geom_point() +
    geom_hline( yintercept = 0 ) +
    scale_x_continuous( breaks = unique( sum$num_strata ) )



# Bias plot
ggplot( sum, aes( num_strata, bias, col=method ) ) +
    facet_grid( one_sided ~ nt_shift, labeller = label_both ) +
    geom_line() + geom_point() +
    geom_hline( yintercept = 0 ) +
    scale_x_continuous( breaks = unique( sum$num_strata ) )



# SE plot
ggplot( sum, aes( num_strata, SE, pch = method, col=method ) ) +
    facet_grid( one_sided ~ nt_shift, labeller = label_both ) +
    geom_line() + geom_point( size = 4, alpha=0.5 ) +
    geom_hline( yintercept = 0 ) +
    scale_x_continuous( breaks = unique( sum$num_strata ) ) +
    theme_minimal()


# RMSE plot (one-sided only)
sum %>% filter( one_sided == TRUE ) %>%
    ggplot( aes( num_strata, RMSE, pch = method, col=method ) ) +
    facet_grid( one_sided ~ nt_shift ) +
    geom_line() + geom_point( size = 4, alpha=0.5 ) +
    geom_hline( yintercept = 0 ) +
    scale_x_continuous( breaks = unique( sum$num_strata ) )



# One-sided plot, bias, SE, and RMSE all together
sumOS <- sum %>% filter( one_sided == TRUE ) %>%
    pivot_longer( cols=c( "bias", "SE", "RMSE" ),
                  names_to = "metric",
                  values_to = "value" ) %>%
    mutate( metric = factor( metric, levels=c("bias","SE", "RMSE" )))

ggplot( sumOS, aes( num_strata, value, col=method ) ) +
    facet_grid(  nt_shift ~ metric ) +
    geom_line() + geom_point( size = 2, alpha=0.5 ) +
    geom_hline( yintercept = 0 ) +
    scale_x_continuous( breaks = unique( sum$num_strata ) ) +
    scale_y_continuous( limits = c(-2,4), breaks = c( -2, 0, 2, 4 ) ) +
    theme_minimal() +
    labs( x = "number of random strata", y = "standard deviations" ) +
    theme(plot.margin=margin(0,10,0,5),
          panel.spacing = unit(1, "lines") ) # Adjust space between facets



ggsave( file=paste0( figures_dir, "random_strat_nohet.pdf" ),
        width = 6, height = 4 )



# Two-sided plot, bias, SE, and RMSE all together
sumTS <- sum %>% filter( one_sided == FALSE ) %>%
    pivot_longer( cols=c( "bias", "SE", "RMSE" ),
                  names_to = "metric",
                  values_to = "value" ) %>%
    mutate( metric = factor( metric, levels=c("bias","SE", "RMSE" )))

ggplot( sumTS, aes( num_strata, value, col=method ) ) +
    facet_grid(  nt_shift ~ metric ) +
    geom_line() + geom_point( size = 4, alpha=0.5 ) +
    geom_hline( yintercept = 0 ) +
    scale_x_continuous( breaks = unique( sum$num_strata ) ) +
    theme_minimal()



#### Integrity check: What is going on with single strata differences? #####


# This verifies that in a single strata we do not have disagreement between
# the three estimators.  They can all be undefined together in some cases.

sss = filter( dtT, num_strata == 1 )

head( sss )

sw <- sss %>% dplyr::select(method, nt_shift, LATE.hat ) %>%
    mutate( #LATE.hat = ifelse( is.na(LATE.hat), 0, LATE.hat ),
        id = rep( 1:(n()/3), each=3 ) ) %>%
    pivot_wider( names_from="method", values_from="LATE.hat" )
sw

sw_bad <- filter( sw, is.na( UNSTRAT ) | UNSTRAT != IV_a )
sw_bad
all( is.na( sw_bad$IV_w ) )
all( is.na( sw_bad$IV_a ) )
all( is.na( sw_bad$UNSTRAT ) )


