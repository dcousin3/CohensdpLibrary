# Cohen's dp library: Getting the Cohen's $d_p$ and its confidence interval in any design

<!-- badges: start -->
[![CRAN Status](https://www.r-pkg.org/badges/version/CohensdpLibrary)](https://cran.r-project.org/package=CohensdpLibrary)
<!-- badges: end -->

This repository contains source code and a package suitable for R 3.0 and above. It computes 
the Cohen's $d_p$ and its confidence intervals in within-subject and between-subject designs.

See :
Cousineau, D. (under preparation) The Exact confidence interval of the Cohen's dp in repeated-measure designs.
for the formal mathematical derivations of the results.

Note that for between-subject designs, the library MBESS already has a potent function to get confidence intervals (Kelley, 2019).

You can install this library on you computer from CRAN (note the uppercase C and uppercase L)
```{r}
install.packages("CohensdpLibrary")
```

or if the library devtools is installed with:
```{r}
devtools::install_github("dcousin3/CohensdpLibrary")<br>
```

and before using it:
```{r}
library(CohendspLibrary)
```

The main function is Cohensdp, which returns the Cohen's $d_p$ and its 
confidence intervals under various designs.

```{r}
# this returns the lower 95% confidence interval bound, d_p then the upper bound
Cohensdp( statistics = list(m1=76, m2=72, n=20, s1=14.8, s2=18.8, rho=0.2),
          design = "within"
)
#  -0.3221494  0.2364258  0.7918157
```

You get a more readable output with ``summarize``, e.g.,

```{r}
summarize(Cohensdp( statistics = list(m1=76, m2=72, n=20, s1=14.8, s2=18.8, rho=0.2),
                    design = "within")
)
# Cohen's dp         = 0.236
#   95.0% Confidence interval = [-0.322, 0.792]
```

Check the web site for more on the functions.
 ``help(CohensdpLibrary)`` will get you started.


