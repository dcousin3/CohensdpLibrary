##############################################################################
# DEFINITIONS of the three display methods: 
#   print: The usual, just show a list of three numbers, (lowerCI, dp, upperCI)
#   summarize (or summary): provides a human-readable output
#   explain : provides all sorts of information regarding the computations
##############################################################################

# AS a reminder, the CohensdpOjbect has all these keys
#   type       = "dp"
#   estimate   = res[1],
#   interval   = c(res[2], res[3])
#   statistics = statistics, 
#   design     = design, 
#   gamma      = gamma



 
#' @export 
print.CohensdpObject <- function(x, ...) {
    switch (x$type, 
        "dp" = { 
            print( c(x$interval[1], x$estimate, x$interval[2] ) )
        },
        "gp" = {
            print( x$estimate )
        },
        "J"  = {
            print ( x$estimate )
        },
        print ("Unknown type of object...")
    )
    return(invisible(x))
}


#' @method summarize CohensdpObject 
#' @export 
summarize.CohensdpObject <-  function(x, ...) {
    format = getOption("CohensdpLibrary.FORMAT")
    sts <- x$statistics
    nu = switch(x$design,
            "single"  = sts$n-1,
            "within"  = 2*(sts$n-1),
            "between" = sts$n1+sts$n2-2
        )
    switch(x$type,
        "dp" = { cat(
            sprintf(paste("Cohen's dp         = ",format,"\n", sep=""), x$estimate),
            sprintf(paste("%5.1f%% Confidence interval = [",format,", ",format,"]\n",sep=""), 
                          x$gamma*100, x$interval[1], x$interval[2])             
            ) },
        "gp" = { cat(
            sprintf(paste("Hedges' gp         = ",format,"\n", sep=""), x$estimate),
            sprintf("No confidence interval known for Hedges's gp\n") 
            ) },
        "J"  = { cat(
            switch(x$design,
                "single"  = sprintf(paste("Correction factor J(%d) = ",format,"\n", sep=""), nu, x$estimate ),
                "within"  = sprintf(paste("Correction factor J(%d, ",format,") = ",format,"\n", sep=""), nu, sts$rho, x$estimate ),
                "between" = sprintf(paste("Correction factor J(%d) = ",format,"\n", sep=""), nu, x$estimate  )
            )
            ) },
        print("Unknown type of object to summarize...")
    )
    return(invisible(x))
}


#' @method explain CohensdpObject 
#' @export
explain.CohensdpObject <- function(x, ...) {
    format = getOption("CohensdpLibrary.FORMAT")
    sts <- x$statistics
    nu = switch(x$design,
            "single"  = sts$n-1,
            "within"  = 2*(sts$n-1),
            "between" = sts$n1+sts$n2-2
        )

    switch(x$type,
        "dp" = { cat(
            sprintf(paste("Cohen's dp         = ",format,"\n",sep=""), x$estimate),
            switch(x$design, 
                "single"  = sprintf(paste("\t sample mean ",format," is compared to assumed mean ",format,"\n", sep=""), sts$m, sts$m0 ),
                "within"  = sprintf(paste("\t first measure's sample mean ",format," is compared to second measure's sample mean ",format,"\n", sep=""), sts$m1, sts$m2 ),
                "between" = sprintf(paste("\t first group's sample mean ",format," is compared to second group's sample mean ",format,"\n", sep=""), sts$m1, sts$m2 )
            ),
            switch(x$design, 
                "single"  = sprintf(paste("\t sample standard deviation ",format," is the denominator\n", sep=""), sts$s ),
                "within"  = sprintf(paste("\t pooled standard deviation ",format," is the denominator\n", sep=""), sqrt( (sts$s1^2 + sts$s2^2) / 2) ),
                "between" = sprintf(paste("\t pooled standard deviation ",format," is the denominator\n", sep=""), sqrt( ( (sts$n1-1)*sts$s1^2 + (sts$n2-1)*sts$s2^2) / (sts$n1+sts$n2-2) ) )
            ),
            sprintf(paste("%5.1f%% Confidence interval = [",format,", ",format,"]\n",sep=""), x$gamma*100, x$interval[1], x$interval[2]),
            switch(x$design, 
                "single"  = sprintf("\t*: confidence interval obtained from the method lambda-prime with %d degrees of freedom (Lecoutre, 2007, Journal of Modern Applied Statistical Methods)\n", sts$n-1),
                "within"  = if ("rho" %in% names(sts)) 
                                sprintf("\t*: confidence interval obtained from the method lambda-second with %d degrees of freedom (Cousineau, 2022, The Quantitative Methods for Psychology)\n", 2*(sts$n-1))
                            else
                                sprintf("\t*: confidence interval obtained from the Adjusted lambda' method with %d degrees of freedom (Cousineau & Goulet-Pelletier, 2021, The Quantitative Methods for Psychology)\n", 2*(sts$n-1))
                ,
                "between" = sprintf("\t*: confidence interval obtained from the lambda-prime method with %d degrees of freedom (Lecoutre, 2007, Journal of Modern Applied Statistical Methods)\n", sts$n1+sts$n2-2) 
            )
            ) },
        "gp" = { cat(
            sprintf(paste("Hedges' gp         = ",format,"\n",sep=""), x$estimate),
            switch(x$design, 
                "single"  = sprintf(paste("\t sample mean ",format," is compared to assumed mean ",format,"\n", sep=""), sts$m, sts$m0 ),
                "within"  = sprintf(paste("\t first measure's sample mean ",format," is compared to second measure's sample mean ",format,"\n", sep=""), sts$m1, sts$m2 ),
                "between" = sprintf(paste("\t first group's sample mean ",format," is compared to second group's sample mean ",format,"\n", sep=""), sts$m1, sts$m2 )
            ),
            switch(x$design, 
                "single"  = sprintf(paste("\t sample standard deviation ",format," is the denominator\n", sep=""), sts$s ),
                "within"  = sprintf(paste("\t pooled standard deviation ",format," is the denominator\n", sep=""), sqrt( (sts$s1^2 + sts$s2^2) / 2) ),
                "between" = sprintf(paste("\t pooled standard deviation ",format," is the denominator\n", sep=""), sqrt( ( (sts$n1-1)*sts$s1^2 + (sts$n2-1)*sts$s2^2) / (sts$n1+sts$n2-2) ) )
            ),
            sprintf("\t*: there is no known confidence interval for an (unbiased) Hedges' gp.\n")
            ) },
        "J" = { cat(
            switch(x$design,
                "single"  = sprintf(paste("Correction factor J(%d) = ",format,"\n", sep=""), nu, x$estimate  ),
                "within"  = sprintf(paste("Correction factor J(%d, ",format,") = ",format,"\n", sep=""), nu, sts$rho, x$estimate ),
                "between" = sprintf(paste("Correction factor J(%d) = ",format,"\n", sep=""), nu, x$estimate   )
            ),
            switch(x$design, 
                "single"  = sprintf(paste("\t*: degrees of freedom nu = n-1 = %d-1; Hedges, 1981, Journal of Educational Statistics.\n", sep=""), sts$n),
                "within"  = sprintf(paste("\t*: degrees of freedom nu = 2(n-1) = 2(%d-1); Cousineau, 2022, The Quantitative Methods for Psychology.\n", sep=""), sts$n),
                "between" = sprintf(paste("\t*: degrees of freedom nu = n1 + n2 -2 = %d+%d-2; Hedges, 1981, Jjournal of Educational Statistics.\n", sep=""), sts$n1,sts$n2) 
            )
            ) },
        print("Unknown type of object to explain...")
    )
    return(invisible(x))
}



#' @export
summary.CohensdpObject <- function(object, ...) {
    summarize(object)
}

