
*****************************************************************
* Replication file for:
* "Elasticities and the inverse hyperbolic sine transformation"
* by Marc F. Bellemare & Casey J. Wichman (2019)
*****************************************************************


clear
drop _all
scalar drop _all
set obs 1000
set seed 12345

* create sample data
gen x = rnormal(10,2)
gen y = 5 + 0.5*x + rnormal(0,2)

* generate transformed variables
egen xbar = mean(x)
egen ybar = mean(y)
gen ihs_y = asinh(y)
gen ihs_x = asinh(x)
qui foreach i in 1 5 10 50 100 1000 10000 {
	gen ihs_y_`i' = asinh(y*`i')
	gen ihs_x_`i' = asinh(x*`i')
	gen y_`i' = y*`i'
	gen x_`i' = x*`i'
	gen ybar_`i' = ybar*`i'
	gen xbar_`i' = xbar*`i'
}





** Model 1 **
* Standard Elasticities: y on x *
qui foreach i in 1 5 10 50 100 1000 10000 {
	reg y_`i' x_`i'
	nlcom _b[x_`i']*xbar_`i'/ybar_`i', post
	scalar m1_`i' = _b[_nl_1]
}

** Model 2 **
* y on IHS(x) *
qui foreach i in 1 5 10 50 100 1000 10000 {
	reg y_`i' ihs_x_`i'
	nlcom (_b[ihs_x_`i']*xbar_`i')/(ybar_`i'*sqrt(xbar_`i'^2 + 1)), post
	scalar m2_`i' = _b[_nl_1]
}

** Model 3 **
* IHS(y) on x *
qui foreach i in 1 5 10 50 100 1000 10000 {
	reg ihs_y_`i' x_`i'
	nlcom _b[x_`i']*xbar_`i'*((sqrt(ybar_`i'^2 + 1))/ybar_`i'), post
	scalar m3_`i' = _b[_nl_1]
}

** Model 4 **
* IHS(y) on IHS(x) *
qui foreach i in 1 5 10 50 100 1000 10000 {
	reg ihs_y_`i' ihs_x_`i'
	nlcom (_b[ihs_x_`i']*xbar_`i'*(sqrt(ybar_`i'^2 + 1)))/(ybar_`i'*sqrt(xbar_`i'^2 + 1)), post
	scalar m4_`i' = _b[_nl_1]
}






* post statitics to scalars for inclusing in estout table
# delimit;
eststo clear;
eststo;
qui estadd scalar ss1 = m1_1, replace;
qui estadd scalar ss5 = m1_5, replace;
qui estadd scalar ss10 = m1_10, replace;
qui estadd scalar ss50 = m1_50, replace;
qui estadd scalar ss100 = m1_100, replace;
qui estadd scalar ss1000 = m1_1000, replace;
qui estadd scalar ss10000 = m1_10000, replace;

eststo;
qui estadd scalar ss1 = m2_1, replace;
qui estadd scalar ss5 = m2_5, replace;
qui estadd scalar ss10 = m2_10, replace;
qui estadd scalar ss50 = m2_50, replace;
qui estadd scalar ss100 = m2_100, replace;
qui estadd scalar ss1000 = m2_1000, replace;
qui estadd scalar ss10000 = m2_10000, replace;

eststo;
qui estadd scalar ss1 = m3_1, replace;
qui estadd scalar ss5 = m3_5, replace;
qui estadd scalar ss10 = m3_10, replace;
qui estadd scalar ss50 = m3_50, replace;
qui estadd scalar ss100 = m3_100, replace;
qui estadd scalar ss1000 = m3_1000, replace;
qui estadd scalar ss10000 = m3_10000, replace;

eststo;
qui estadd scalar ss1 = m4_1, replace;
qui estadd scalar ss5 = m4_5, replace;
qui estadd scalar ss10 = m4_10, replace;
qui estadd scalar ss50 = m4_50, replace;
qui estadd scalar ss100 = m4_100, replace;
qui estadd scalar ss1000 = m4_1000, replace;
qui estadd scalar ss10000 = m4_10000, replace;
#delimit cr






* write table to tex file
#delimit;
esttab using "results/mcmctable.tex", replace
		nonotes noobs nonumber label substitute(\_ _)
		drop(_nl_1) nogaps compress
		mgroups("Empirical specification", pattern(1 0 0 0)
		prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
		mtitle("Linear-linear" "Linear-arcsinh" 
			   "arcsinh-linear" "arcsinh-arcsinh")
		scalars("x Values of \textit{k}:"
		"ss1 1"
		"ss5 5"
		"ss10 10"
		"ss50 50"
		"ss100 100"
		"ss1000 1000"
		"ss10000 10,000")
		sfmt(a6)
		;
#delimit cr



