

# These functions allow analysis using factorial MCAFE estimates from
# Blackwell & Pashley (2023)
# See R package factiv for standard implementation for MCAFE estimation
#
# 2024, by Pashley & Miratrix




#' Estimate MCAFE for given dataset.
#'
#' Same as poststratIV's IV.est() method, but updated to calculate MCAFE
mcafe_est <- function( data,  Yobs="Yobs", S="S", Z="Z", W = "W", W2=NULL){
  data$Yobs = data[[Yobs]]
  data$S = data[[S]]
  data$Z = data[[Z]]
  data$W = data[[W]]
  if(!is.null(W2)){
    K = 3
    dataG = data %>%
      group_by( Z, W, W2 ) %>%
      summarise(Y.bar = mean( Yobs ),
                p.S = mean(S),
                var.Y = var( Yobs ),
                var.S = var( S ),
                cov.Y.s = cov( Yobs, S ),
                n = n(),
                .groups = "drop" )
  }else{
  K = 2
  dataG = data %>%
    group_by( Z, W ) %>%
    summarise(Y.bar = mean( Yobs ),
            p.S = mean(S),
            var.Y = var( Yobs ),
            var.S = var( S ),
            cov.Y.s = cov( Yobs, S ),
            n = n(),
            .groups = "drop" )
  }
  pi.hat <- mean(dataG$p.S[dataG$Z == 1]) - mean(dataG$p.S[dataG$Z == 0])
  if ( pi.hat < 0 ) {
    browser()
  }

  ITT.hat = mean(dataG$Y.bar[dataG$Z == 1]) - mean(dataG$Y.bar[dataG$Z == 0])

  # Store estimates
  res = data.frame( ITT.hat=ITT.hat,
                    pi.hat = pi.hat,
                    LATE.hat = ITT.hat / pi.hat )

  # Calc primary SEs of the direct estimands
  res = mutate( res,
                SE.ITT = with( dataG, sqrt( sum( var.Y / n )/2^(2*(K-1)) ) ),
                SE.pi = with( dataG, sqrt( sum( var.S / n ) /2^(2*(K-1))) ),
                cov.ITT.pi = with( dataG,  sum( cov.Y.s / n )/2^(2*(K-1))  ) )

  # Calculate SEs, including delta method SE (dropping correlation term)
  res = mutate( res,
                SE.wald = SE.ITT / pi.hat,
                SE.delta = sqrt( SE.ITT^2 / pi.hat^2 + ITT.hat^2 * SE.pi^2 / pi.hat^4 - 2*ITT.hat * cov.ITT.pi / pi.hat^3 )
  )

  res$n = nrow( data )

  # clean up if pi.hat = 0
  if ( res$pi.hat == 0 ) {
    res$LATE.hat = NA
    res$SE.wald = NA
    res$SE.delta = NA
  }

  # and done!
  res



}

#' Estimate post-stratified MCAFE
#'
#' Same as IV.est.strat() but updated to use MCAFE estimate instead of
#' pure marginal IV estimates
IV.est.strat.fac = function( data,  Yobs="Yobs", S="S", Z="Z", W = "W", W2 = NULL, strat_var ="X",
                         DSS_cutoff = 0.02,
                         drop_without_warning = FALSE,
                         include_blocks = TRUE,
                         include_FSS = TRUE ) {

  stopifnot( "S" != Yobs && "S" != Z )
  stopifnot( "Yobs" != S && "Yobs" != Z )
  stopifnot( "Z" != Yobs && "Z" != S )

  data$Yobs = data[[Yobs]]
  data$S = data[[S]]
  data$Z = data[[Z]]
  data$W = data[[W]]
  if(!is.null(W2)){
    data$W2 = data[[W2]]
  }

  # give canonical name to our passed variable
  data$Xblk = as.character( data[[strat_var]] )
  stopifnot( !is.null(data$Xblk) )

  dataG = data %>%
    group_by( Xblk ) %>%
    nest()

  gsum = set_names( dataG$data, nm=dataG$Xblk ) %>%
    map_df( mcafe_est, .id="Xblk" )

  # Zero out empty estimated strata
  empty = gsum$pi.hat <= 0 | !is.finite( gsum$LATE.hat )
  gsum$LATE.hat[ empty  ] = 0
  gsum$SE.wald[ empty ] = 0
  gsum$SE.delta[ empty ] = 0

  # Estimate the number of compliers in each strata
  gsum = mutate( gsum,
                 n_comp = pi.hat*n ) %>%
    arrange( Xblk )
  gsum

  # Drop all strata without at least 2 Tx and Co units
  nr = nrow(gsum)
  gsum = filter( gsum, !is.na( SE.ITT ) )
  nr_post = nrow(gsum)
  if ( !drop_without_warning && nr > nr_post ) {
    warning( glue::glue( "{nr - nr_post} blocks dropped (out of {nr}) due to missing Tx or Co units" ) )
  }


  # IV_w: The weighted average of the strata-level estimates
  IV_w <- gsum %>%
    filter( pi.hat > 0 ) %>%
    aggregate_strata()

  # IV_a: The other version of stratification: First calculate
  # post-stratified ITT.hat, SE.ITT, pi.hat, and SE.pi, and then
  # calculate LATE with the post-stratified estimates of ITT and pi.

  W = sum(gsum$n)^2
  gsum2 <- gsum %>% summarise( ITT.hat.ps = weighted.mean( ITT.hat, w=.$n ),
                               pi.hat.ps = weighted.mean( pi.hat, w=.$n ),
                               SE.ITT.ps = sqrt(sum(SE.ITT^2*n^2)/W),
                               SE.pi.ps = sqrt(sum(SE.pi^2*n^2)/W),
                               cov.ITT.pi.ps = sum(cov.ITT.pi*n^2)/W,
                               n = sum(n),
                               n_comp = sum(n_comp),
                               .groups = "drop")
  stopifnot( nrow(gsum2)==1 )
  if ( is.na( gsum2$ITT.hat.ps ) ) {
    browser()
  }
  IV_a = gsum2 %>% mutate( Xblk = "IV_a",
                           LATE.hat = ITT.hat.ps/pi.hat.ps,
                           SE.wald = SE.ITT.ps / abs(pi.hat.ps),
                           SE.delta = sqrt( SE.ITT.ps^2 / pi.hat.ps^2 +
                                              (ITT.hat.ps/pi.hat.ps)^2 * SE.pi.ps^2 / pi.hat.ps^2 -
                                              2 * (ITT.hat.ps/pi.hat.ps) * cov.ITT.pi.ps / pi.hat.ps^2 ),
                           n = n,
                           n_comp = n_comp,
                           pi.hat = pi.hat.ps,
                           ITT.hat = ITT.hat.ps,
                           SE.ITT = SE.ITT.ps,
                           SE.pi = SE.pi.ps,
                           cov.ITT.pi=cov.ITT.pi.ps)

  if ( is.infinite( IV_a$LATE.hat ) ) {
    IV_a$LATE.hat = NA
    IV_a$SE.wald = NA
    IV_a$SE.delta = NA
  }


  # The "drop small strata" estimator
  res_DSS = gsum %>%
    filter( pi.hat > DSS_cutoff ) %>%
    aggregate_strata()
  res_DSS$Xblk = "DSS"


  # "Drop Small F", or "Testing-F" estimator
  res_DSF = NULL
  if ( include_FSS ) {
    Fvals <- data %>% nest_by( Xblk ) %>%
      summarise( Fval = anova( lm( S ~ Z, data=data ) )$`F value`[[1]],
                 .groups = "drop" )
    res_DSF <- gsum %>%
      left_join( Fvals, by="Xblk" ) %>%
      filter( Fval > 10 ) %>%
      aggregate_strata()
    res_DSF$Xblk = "DSF"

  }

  # The "double weight" estimator
  # res_DDW = gsum %>%
  #     filter( pi.hat > 0 ) %>%
  #     aggregate_strata( weighting = "double" )
  # res_DDW$Xblk = "DDW"

  # The "double weight" (precision-weighted, actually) estimator
  res_DDW = gsum %>%
    filter( pi.hat > 0 ) %>%
    aggregate_strata( weighting = "precision" )
  res_DDW$Xblk = "PWIV"

  # No stratification as baseline
  res_unstrat = mcafe_est( data )
  res_unstrat = mutate( res_unstrat,
                        Xblk = "UNSTRAT",
                        n_comp = n * pi.hat )

  gsum$k = 1
  res_unstrat$k = nrow(gsum)
  IV_a$k = nrow(gsum)

  if ( !include_blocks ) {
    gsum = NULL
  }

  res = bind_rows( gsum, res_unstrat, IV_w, IV_a, res_DSS, res_DDW, res_DSF ) %>%
    relocate( Xblk )

  res$empty = sum( empty )

  if ( any( is.nan( res$LATE.hat ) ) || any( is.infinite( res$LATE.hat ) ) ) {
    stop( "Uneexpected NaN or Infinite value" )
  }

  res
}






