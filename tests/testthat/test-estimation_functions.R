
# Testing the IV estimation code
library( poststratIV )
library( testthat )


test_that("do the core estimation functions run?", {

    dat = make.dat.tuned( N = 1000,
                          pi_c = c( 0.01, 0.1, 0.35, 0.80 ) )
    dat = rand.exp(dat)

    expect_true( is.data.frame(dat) )
    expect_true( nrow(dat) == 1000 )
    expect_true( all( c( "Y0", "Y1", "X", "complier", "S0", "S1", "Y1", "S", "Yobs" ) %in% names(dat) ) )


    # Just vanilla IV
    res <- IV.est( dat )
    expect_true( is.data.frame(res) )

    head( dat )
    res <- IV.est.strat( dat )
    res
    expect_true( is.data.frame(res) )


    res = filter( res, Xblk %in% c( "IV_w", "IV_a" ) )
    res
    round( diff( res$LATE.hat ), digits=3 )
    round( diff( res$SE.wald ), digits= 3 )
    round( diff( res$SE.delta ), digits= 3 )


})








test_that("basic check of estimation functions", {


    set.seed( 40404 )
    dat = make.dat.tuned( N = 1000,
                          pi_c = c( 0.01, 0.1, 0.35, 0.80 ) )
    dat = rand.exp(dat)
    expect_true( is.data.frame(dat) )

    head( dat )
    table( dat$X, dat$complier )

    # Just vanilla IV
    res <- IV.est( dat )
    expect_true( nrow(res) == 1 )

    head( dat )
    res <- IV.est.strat( dat )
    res
    expect_true( all( !is.na( res$ITT.hat ) ) )


    res = filter( res, Xblk %in% c( "IV_w", "IV_a" ) )
    res

    a = round( diff( res$LATE.hat ), digits=3 )
    b = round( diff( res$SE.wald ), digits= 3 )
    expect_equal( a, b )

    round( diff( res$SE.delta ), digits= 3 )

    # Checking that custom names of variables works
    dat2 = rename( dat,
                   XX = X,
                   YY = Yobs,
                   SS = S,
                   ZZ = Z )
    r2 <- IV.est.strat( dat2, Yobs = "YY", S = "SS", Z = "ZZ", strat_var = "XX" ) %>%
        tidy_IV_table()
    r1 = IV.est.strat(dat) %>%
        tidy_IV_table()
    expect_equal( r2$ITT.hat, r1$ITT.hat )
} )





test_that( "missing data handled", {

    # Looking at many many strata: drop automatically?
    set.seed( 33303 )
    dat = make.dat.tuned( N = 200,
                          pi_c = c( 0.01, 0.1, 0.35, 0.80 ) )
    dat = rand.exp(dat)

    dat
    dat$X = paste0( dat$X, "-", sample(1:7, size=nrow(dat), replace=TRUE ) )
    table( dat$X, dat$Z )

    dd = filter( dat, X == "X4-1" ) #X3-6" ) #X4-2" )
    dd
    IV.est( dd )

    # Drop automatically, give warning
    expect_warning( res0 <- IV.est.strat( dat, include_blocks = FALSE, include_FSS = FALSE ) )

    res <- IV.est.strat( dat, include_blocks = FALSE,
                         drop_without_warning = TRUE, include_FSS = FALSE )
    res



    # Looking at missing data on stratification covariate
    dat = make.dat.tuned( N = 100,
                          pi_c = c( 0.01, 0.1, 0.35, 0.80 ) )
    dat = rand.exp(dat)

    dat$X[1:10] = NA
    # Just vanilla IV
    IV.est( dat )

    expect_warning( res <- IV.est.strat( dat ) )
    res
    expect_true( res$Xblk[[1]] == "NA" )

    # should throw an error.
    dat$X[11:13] = "NA"
    expect_error(
        res <- IV.est.strat( dat )
    )

})





test_that( "aggegate_strata() works right", {


    set.seed( 333034 )
    dat = make.dat.tuned( N = 1000,
                          pi_c = c( 0.01, 0.1, 0.35, 0.80 ) )
    dat = rand.exp(dat)

    gg = IV.est.strat( dat, return_strata_only = TRUE )
    gg
    gg$LATE.hat = 1:4
    gg$n = c( 100, 200, 300, 400 )
    gg$pi.hat = c( 0, 0.2, 0.3, 0.9 )
    gg$SE.wald = c( 1, 5, 10, 100 )

    emp <- poststratIV:::aggregate_strata( gg %>% filter( Xblk == "none" ), weighting = "normal" )
    expect_true( nrow( emp ) == 1 )
    expect_true( is.na( emp$LATE.hat ) )

    nm <- poststratIV:::aggregate_strata( gg, weighting = "normal" )
    expect_equal( nm$LATE.hat, weighted.mean( gg$LATE.hat, gg$n_comp ) )


    nm2 <- poststratIV:::aggregate_strata( gg, weighting = "double" )
    nm2
    expect_equal( nm2$LATE.hat, weighted.mean( gg$LATE.hat, gg$n_comp^2 ) )

    nm3 <- poststratIV:::aggregate_strata( gg, weighting = "precision" )
    nm3
    expect_equal( nm3$LATE.hat, weighted.mean( gg$LATE.hat, 1 / gg$SE.wald^2 ) )

    expect_equal( nm3$ITT.hat, weighted.mean( gg$ITT.hat, gg$n ) )
    expect_equal( nm3$pi.hat, weighted.mean( gg$pi.hat, gg$n ) )


})



test_that( "IVa works right", {

    set.seed( 333034 )
    pps <- c( 0.2, 0.7, 0.05, 0.05 )
    dat = make.dat.tuned( N = 10000,
                          pstrat = pps,
                          pi_c = c( 0.01, 0.1, 0.35, 0.80 ) )
    dat = rand.exp(dat)

    t = table( dat$X )
    expect_true( all( abs( as.numeric(t) - 10000 * pps ) < 100 ) )

    gg = IV.est.strat( dat, return_strata_only = TRUE )

    raw <- poststratIV:::calc_IVa( gg )
    expect_true( raw$SE.delta > 0 )

    gg
    gg$n = c( 100, 200, 300, 400 )
    gg$ITT.hat = 4:1
    gg$pi.hat = c(0.1,0.2,0.9,0)
    p2 <- poststratIV:::calc_IVa( gg )
    expect_equal( p2$ITT.hat, weighted.mean( gg$ITT.hat, gg$n ) )
    expect_equal( p2$pi.hat, weighted.mean( gg$pi.hat, gg$n ) )
})



