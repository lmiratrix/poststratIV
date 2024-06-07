

# Plot setup utils to get consistant colors, etc

library( scales )

# Make factor right order
fix_method_fct <- function( method ) {
    method = recode_factor(method,
                           `UNSTRAT` = "Unstrat",
                           `IV_a` = "IV[a]",
                           `IV_w` = "IV[w]",
                           `2SLS` = "`2SLS`" )
    fct_relevel( method,
                 rev( c( "Oracle",
                         "DSF",
                         "PWIV",
                         "DSS0",
                         "IV[w]",
                         "`2SLS`",
                         "IV[a]",
                         "Unstrat" )  ) )
}



# Add color for method
add_color_line_thing <- function( plt ) {
    plt +
        scale_color_manual( values = c( "UNSTRAT" = "darkgrey",
                                        "Unstrat" = "darkgrey",
                                        "2SLS" = "orange",
                                        "`2SLS`" = "orange",
                                        "IV_a" = "brown",
                                        "IV_w" = "blue",
                                        "IV[a]" = "brown",
                                        "IV[w]" = "blue",
                                        "DSS0" = "blue",
                                        "DSF" = "darkgreen",
                                        "PWIV" = "darkgreen",
                                        "Oracle" = "darkgrey" ), 
                            labels = parse_format() ) +
        scale_linetype_manual( values = c( "UNSTRAT" = "solid",
                                           "Unstrat" = "solid",
                                           "2SLS" = "dotted",
                                           "`2SLS`" = "dotted",
                                           "IV_a" = "solid",
                                           "IV_w" = "solid",
                                           "IV[a]" = "solid",
                                           "IV[w]" = "solid",
                                           "DSS0" = "dotted",
                                           "DSF" = "dotted",
                                           "PWIV" = "solid",
                                           "Oracle" = "twodash" ),
                               labels = parse_format() )
}