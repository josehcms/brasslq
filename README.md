# brasslq: Brass-Log-Quadratic Integrated Method for Under-5 Mortality Estimation

<!-- badges: start -->
<!-- badges: end -->

## Overview

**brasslq** integrates the Brass indirect method for estimating child mortality with the
log-quadratic model of under-5 mortality (Guillot et al. 2022) to reconstruct detailed age
patterns of under-5 mortality from **summary birth histories** (children ever born and
children dead by age group of the mother), as collected in censuses and surveys.

Compared with the classic Brass method (Trussell variant), the package:

- replaces fixed model life tables with an **extended version of the log-quadratic model**
  (up to age 50), indexed by the level of under-5 mortality, q(5), and the shape parameter *k*;
- computes the Brass multipliers and the Coale-Trussell time-location coefficients
  **from user-chosen fertility schedules** (Coale-Trussell model schedules or UN WPP 2024
  schedules), instead of fixed tables;
- **iterates** between the Brass method and the log-quadratic model until the estimated
  q(5) converges;
- returns full **under-5 mortality schedules** (22 cumulative probabilities of dying, from
  the first week of life to age 5) using the [logquad5q0](https://github.com/josehcms/logquad5q0) package.

To report issues or suggest improvements, please open an
[issue](https://github.com/josehcms/brasslq/issues) or contact <josehcms@sas.upenn.edu>.

#### Author

- **José H C Monteiro da Silva** - University of Pennsylvania

## Installation

``` r
install.packages("remotes")

remotes::install_github("josehcms/brasslq")

library(brasslq)
```

## Main function: `brass_logquad`

| Argument | Description |
|---|---|
| `PDx5` | Proportion of children dead by 5-year age group of mother (15-19, ..., 45-49) |
| `P1`, `P2`, `P3` | Mean number of children ever born to women aged 15-19, 20-24 and 25-29 |
| `tsurv` | Reference date of the survey/census (decimal year) |
| `fertility_schedules` | Fertility schedules used to derive the coefficients (`fert_scheds_ct`, `fert_scheds_un`, or a user-supplied data frame with columns `id`, `x`, `fx`, `tx_subset`) |
| `k` | Shape parameter of the log-quadratic model (default = 0) |
| `r` | Population growth rate (default = 0.02) |
| `Rd` | Annual rate of mortality decline used for time location (default = 0.03) |
| `q5_init` | Initial value of q(5) for the iteration (default = 0.150) |
| `tol`, `max_iter` | Convergence tolerance (default = 1e-4) and maximum number of iterations (default = 20) |
| `n_cores` | Number of cores for looping over fertility schedules (default = 1; parallel computing is not available on Windows) |

The function returns a list with:

- `brasslq_iterative`: Brass estimates by age group of mother: q(a), reference dates and coefficients;
- `u5mr_schedule`: detailed under-5 mortality schedules from the log-quadratic model, fitted to
  q(2), q(3) and q(5), each one with its reference date.

## Example

Data on children ever born and children surviving by age of mother from the **2008 Census
of Malawi** (reference date: 18 June 2008, i.e., 2008.46), as used in the worked example of the
IUSSP *Tools for Demographic Estimation* chapter on the indirect estimation of child mortality
(<https://demographicestimation.iussp.org/content/indirect-estimation-child-mortality>):

``` r
library(brasslq)

# Malawi, 2008 Census
P1    = 0.2833; P2 = 1.5316; P3 = 2.8487
PDx5  = c( 0.1034, 0.1146, 0.1330, 0.1597, 0.1848, 0.2186, 0.2458 )
tsurv = 2008.46

# Coale-Trussell fertility schedules (fast)
res_ct = brass_logquad( PDx5 = PDx5,
                        P1 = P1, P2 = P2, P3 = P3,
                        tsurv = tsurv,
                        fertility_schedules = fert_scheds_ct,
                        k = 0, Rd = 0.02, r = 0.02,
                        q5_init = PDx5[4],
                        n_cores = 1 )

# Brass estimates by age group of mother
res_ct$brasslq_iterative

# Full under-5 mortality schedules
res_ct$u5mr_schedule
```

The UN WPP 2024 schedules (`fert_scheds_un`) contain one schedule per country-year (17,464
schedules). Using all of them is slower, so a subset (for example, selected years, a region or
a set of countries selected with `country_code`) and/or `n_cores > 1` is recommended. Here, only
schedules for years after 1990 are used:

``` r
# UN WPP 2024 fertility schedules, years after 1990 only
fert_un_1990 = fert_scheds_un[ fert_scheds_un$year > 1990, ]

res_un = brass_logquad( PDx5 = PDx5,
                        P1 = P1, P2 = P2, P3 = P3,
                        tsurv = tsurv,
                        fertility_schedules = fert_un_1990,
                        k = 0, Rd = 0.02, r = 0.02,
                        q5_init = PDx5[4],
                        n_cores = 1 )

res_un$brasslq_iterative
res_un$u5mr_schedule
```

## Data included

- `fert_scheds_ct`: Coale-Trussell model fertility schedules;
- `fert_scheds_un`: UN World Population Prospects 2024 fertility schedules, by country and year (1950-2023);
- `extlogquad_coeffs`: coefficients of the extended log-quadratic model (up to age 50).

## References

Brass W. 1968. Methods of analysis and estimation. In Brass et al., *The Demography of Tropical Africa*. Princeton University Press.

Coale A. J., and J. Trussell. 1977. Annex I: Estimating the time to which Brass estimates apply. *Population Bulletin of the United Nations* 10: 87-89.

Guillot M., J. Romero Prieto, A. Verhulst, and P. Gerland. 2022. Modeling age patterns of under-5 mortality: Results from a log-quadratic model applied to high-quality vital registration data. *Demography* 59(1): 321-347. <https://doi.org/10.1215/00703370-9709538>

Trussell J. 1975. A re-estimation of the multiplying factors for the Brass technique for determining childhood survivorship rates. *Population Studies* 29(1): 97-107.

United Nations. 1983. *Manual X: Indirect Techniques for Demographic Estimation*. New York: United Nations.
