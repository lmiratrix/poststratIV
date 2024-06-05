

test_that("multiplication works", {

    set.seed( 40404 )


    #### Checking make.dat.tuned() ####

    pc = c( 0.05, 0.15, 0.35, 0.80 )
    pn = c( 0.1, 0.1, 0.2, 0.2 )
    dat = make.dat.tuned( N = 1000,
                          pi_c = pc,
                          pi_n = pn )
    dat = rand.exp(dat)
    head( dat )

    expect_true( is.data.frame(dat) )
    expect_true( nrow(dat) == 1000 )
    expect_true( all( c( "Y0", "Y1", "X", "complier", "S0", "S1", "Y1", "S", "Yobs" ) %in% names(dat) ) )

    tb <- table( S0=dat$S0, S1=dat$S1 )
    tb
    expect_true( sum(tb==0) == 1 )

    #debug( describe_sim_data )
    dd <- describe_sim_data( dat )
    dd
    dd$nt[1:4] / pn
    dd$at[1:4] / (1-pc-pn)
    expect_true( is.data.frame(dd) )


    # Check residual variation: sd0 = 1 unless there are shifts for
    # the never-takers or always-takers.
    pc = c( 0.05, 0.15, 0.35, 0.80 )
    pn = 1 - pc #                           c( 0.1, 0.1, 0.2, 0.2 ),

    dat = make.dat.tuned( N = 100000, nt_shift = 0,
                          pi_c = pc,
                          pi_n = pn,
                          sd0 = c( 1, 1, 1, 100 ) )
    dat = rand.exp(dat)
    head( dat )

    table( S0=dat$S0, S1=dat$S1 )

    #debug( describe_sim_data )
    dd <- describe_sim_data( dat )
    expect_equal( dd$sd0[1:4] / c( 1,1,1,100), c(1,1,1,1), tolerance=0.1 )



    #### make_data() and two sided DGP ####

    dat = make_dat( N = 3000, pi_c = 0.15, het_tx = "yes",
                    pred_comp = "yes", pred_Y = "yes",
                    Ybar0 = c(0.5, 0.8, 1.2, 1.6),
                    pstrat = c(0.4, 0.3, 0.1, 0.2),
                    scaled_C = 0.5, one_sided = FALSE )
    dat
    table( Comp=dat$complier, X=dat$X )
    mean( dat$complier )

    # Three types of folks?
    tb <- table( dat$S0, dat$S1 )
    expect_true( length(tb) == 4 )
    expect_true( sum( tb == 0 ) == 1 )

    # Check the observed lines up
    dat %>% group_by( X, S0, S1 ) %>%
        summarise( tau = mean( Y1 - Y0 ),
                   tau_hat = mean(Yobs[Z==1] - mean(Yobs[Z==0] ) ),
                   Z_hat = mean( S[Z==1] ) - mean( S[Z==0] ),
                   n = n(), .groups = "drop" ) %>%
        group_by( X ) %>%
        mutate( prop = n / sum(n) )


    dat = make_dat( N = 3000, pi_c = 0.15, het_tx = "no",
                    pred_comp = "yes", pred_Y = "no",
                    pstrat = c(0.4,0.3,0.1,0.2),
                    scaled_C = 1, one_sided = TRUE )
    dat2 = make_dat( N = 3000, pi_c = 0.15, het_tx = "no",
                     pred_comp = "yes", pred_Y = "no",
                     pstrat = c(0.4,0.3,0.1,0.2),
                     scaled_C = 1, one_sided = FALSE )
    prop.table( table( dat$complier, dat$X ), margin = 2 )
    prop.table( table( dat2$complier, dat2$X ), margin = 2 )

    expect_equal( mean( dat$complier ), 0.15, tolerance = 0.1 )
    expect_equal( mean( dat2$complier ), 0.15, tolerance= 0.1 )

    head( dat2 )


})
