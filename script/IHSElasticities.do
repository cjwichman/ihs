
*****************************************************************
* Replication file for:
* "Elasticities and the inverse hyperbolic sine transformation"
* by Marc F. Bellemare & Casey J. Wichman (2019)
*****************************************************************


clear
drop _all
set obs 1000
set seed 12345

gen x = rnormal(10,2)
gen y = 5 + 0.5*x + rnormal(0,2)


egen xbar = mean(x)
egen ybar = mean(y)

gen ihs_y = asinh(y)
gen ihs_x = asinh(x)

* Note: nlcom relies on the delta method by default for standard errors


* Standard Elasticity: Regression of y on x *

reg y x
nlcom _b[x]*xbar/ybar


* Linear-arcinsh Case: Regression of y on arcsinh(x) *

reg y ihs_x
nlcom (_b[ihs_x]*xbar)/(ybar*sqrt(xbar^2 + 1))


* arcsinh-Linear Case: Regression of arcsinh(y) on x *


reg ihs_y x
nlcom _b[x]*xbar*((sqrt(ybar^2 + 1))/ybar)


* arcsinh-arcsinh Case: Regression of arcsinh(y) on arcsinh(x) *

reg ihs_y ihs_x
nlcom (_b[ihs_x]*xbar*(sqrt(ybar^2 + 1)))/(ybar*sqrt(xbar^2 + 1))

