
*****************************************************************
* Replication file for:
* "Elasticities and the inverse hyperbolic sine transformation"
* by Marc F. Bellemare & Casey J. Wichman (2020)
*****************************************************************

clear all


* Data downloaded from: https://economics.mit.edu/faculty/angrist/data1/mhe/dehejia
* file = https://economics.mit.edu/files/3828
use "data/nswre74.dta", clear




* generate IHS transformations transformations [and log(x+1) for comparison]
** check -help(trig_functions)- for trigonometric functions in Stata
g ihsre75 = asinh(re75)
g ihsre78 = asinh(re78)

g lnre75 = ln(re75+1)
g lnre78 = ln(re78+1)

label var re78 "$\text{Earnings}_{78}$"
label var re75 "$\text{Earnings}_{75}$"
label var treat "NSW"
label var ihsre78 "$\text{arcsinh(Earnings}_{78}\text{)}$"
label var ihsre75 "$\text{arcsinh(Earnings}_{75}\text{)}$"
label var lnre75 "$\ln(\text{Earnings}_{75}+1)$"
label var lnre78 "$\ln(\text{Earnings}_{78}+1)$"


***************************
*** generate results table
***************************

* Run primary regressions + save estimates
#delimit ;
eststo clear;
eststo: reg re78 treat re75; qui estimates save "results/lin_lin", replace;
qui: eststo: reg re78 treat ihsre75; qui estimates save "results/lin_ihs", replace;
qui: eststo: reg ihsre78 treat re75; qui estimates save "results/ihs_lin", replace;
qui: eststo: reg ihsre78 treat ihsre75; qui estimates save "results/ihs_ihs", replace;
qui: eststo: reg lnre78 treat re75; qui estimates save "results/log_lin", replace;
qui: eststo: reg lnre78 treat lnre75; qui estimates save "results/log_log", replace;
#delimit cr




eststo clear
scalar drop _all





* calculate sample means, store as locals for calculating elasticities.
* calculate sample means for re75 and re78
	qui sum re75, meanonly
	scalar mean_re75 = r(mean)
	qui sum re78, meanonly
	scalar mean_re78 = r(mean)
	
* calculate sample means for re78 & re75 by treatment status
	qui sum re75 if treat==1, meanonly
	scalar mean_re75t = r(mean)
	qui sum re78 if treat==1, meanonly
	scalar mean_re78t = r(mean)
	qui sum re75 if treat==0, meanonly
	scalar mean_re75c = r(mean)
	qui sum re78 if treat==0, meanonly
	scalar mean_re78c = r(mean)
	
* calculate sample mean treatment status
	qui sum treat, meanonly
	scalar mean_treat = r(mean)


*****************************************************
** Column 1 -- linear re78, dummy treat, linear re75
*****************************************************

	* elasticity b/w re78 and re75
	est use "results/lin_lin"
	nlcom (_b[re75] * (mean_re75 / mean_re78)), post
	scalar e_re78_re75 = _b[_nl_1]
	scalar e_re78_re75_SE = _se[_nl_1]
	
	* semi-elasticity for treat
	est use "results/lin_lin"		
	nlcom (_b[treat] / mean_re78), post		
	scalar e_re78_treat = _b[_nl_1]
	scalar e_re78_treat_SE = _se[_nl_1]
	
	* store variables
	est use "results/lin_lin"		
	qui estadd scalar elas = e_re78_re75
	qui estadd scalar elas_SE = e_re78_re75_SE
	qui estadd scalar semielas2 = e_re78_treat
	qui estadd scalar semielas2_SE = e_re78_treat_SE
	qui eststo

*****************************************************
** Column 2 -- linear re78, dummy treat, IHS(re75)
*****************************************************

	* elasticity b/w re78 and re75
	est use "results/lin_ihs"
	nlcom (_b[ihsre75] / sqrt(mean_re75^2 +1)) * (mean_re75 / mean_re78), post
	scalar e_re78_ihsre75 = _b[_nl_1]
	scalar e_re78_ihsre75_SE = _se[_nl_1]
	
	* semi-elasticity for treat
	est use "results/lin_ihs"
	nlcom (_b[treat] / mean_re78), post
	scalar e_re78_treat = _b[_nl_1]
	scalar e_re78_treat_SE = _se[_nl_1]

	* store variables
	est use "results/lin_ihs"
	qui estadd scalar elas = e_re78_ihsre75
	qui estadd scalar elas_SE = e_re78_ihsre75_SE
	qui estadd scalar semielas2 = e_re78_treat
	qui estadd scalar semielas2_SE = e_re78_treat_SE
	qui eststo


*****************************************************
** Column 3 -- asinh(re78), dummy treat, linear re75
*****************************************************
	
	* save predicted asinh(y)
	est use "results/ihs_lin"
	capture drop yhat
	qui predict yhat, xb
	sum yhat, meanonly
	scalar mean_yhat = r(mean)

	* calculate elasticity b/w re78 and re75
	est use "results/ihs_lin"
	nlcom (_b[re75]*mean_re75*(sqrt(mean_re78^2+1)/mean_re78)), post
	scalar e_ihsre78_re75 = _b[_nl_1]
	scalar e_ihsre78_re75_SE = _se[_nl_1]
	
	* semi-elasticity for treat using exponential approximation
	est use "results/ihs_lin"
	nlcom (exp(_b[treat] - 0.5*_se[treat]^2) - 1), post
	scalar e_ihsre78_treat_exp = _b[_nl_1]
	scalar e_ihsre78_treat_exp_SE = _se[_nl_1]

	* semi-elasticity using exact calculations
	est use "results/ihs_lin"
	nlcom ([sinh(mean_yhat) / sinh(mean_yhat - _b[treat])] - 1), post
	scalar e_ihsre78_treat_exact = _b[_nl_1]
	scalar e_ihsre78_treat_exact_SE = _se[_nl_1]
	
	* store variables
	est use "results/ihs_lin"
	qui estadd scalar elas = e_ihsre78_re75
	qui estadd scalar elas_SE = e_ihsre78_re75_SE
	qui estadd scalar semielas =  e_ihsre78_treat_exp
	qui estadd scalar semielas_SE =  e_ihsre78_treat_exp_SE
	qui estadd scalar semielas_ihs = e_ihsre78_treat_exact
	qui estadd scalar semielas_ihs_SE = e_ihsre78_treat_exact_SE
	qui eststo

	
*****************************************************
** Column 4 -- asinh(re78), dummy treat, asinh(re75)
*****************************************************

	* save predicted asinh(y)*
	est use "results/ihs_ihs"
	capture drop yhat
	qui predict yhat, xb
	sum yhat, meanonly
	scalar mean_yhat = r(mean)

	* calculate elasticity b/w re78 and re75
	est use "results/ihs_ihs"
	nlcom (_b[ihsre75] * (sqrt(mean_re78^2+1)/mean_re78) * (mean_re75 / sqrt(mean_re75^2+1))), post
	scalar e_ihsre78_ihsre75 = _b[_nl_1]
	scalar e_ihsre78_ihsre75_SE = _se[_nl_1]
	
	* semi-elasticity for treat using exponential approximation
	est use "results/ihs_ihs"
	nlcom (exp(_b[treat] - 0.5*_se[treat]^2) - 1), post
	scalar e_ihsre78_treat_exp = _b[_nl_1]
	scalar e_ihsre78_treat_exp_SE = _se[_nl_1]

	* semi-elasticity using exact calculations
	est use "results/ihs_ihs"
	nlcom ([sinh(mean_yhat) / sinh(mean_yhat - _b[treat])]-1), post
	scalar e_ihsre78_treat_exact = _b[_nl_1]
	scalar e_ihsre78_treat_exact_SE = _se[_nl_1]
	
	* store variables
	est use "results/ihs_ihs"	
	qui estadd scalar elas = e_ihsre78_ihsre75
	qui estadd scalar elas_SE = e_ihsre78_ihsre75_SE
	qui estadd scalar semielas = e_ihsre78_treat_exp
	qui estadd scalar semielas_SE = e_ihsre78_treat_exp_SE
	qui estadd scalar semielas_ihs = e_ihsre78_treat_exact
	qui estadd scalar semielas_ihs_SE = e_ihsre78_treat_exact_SE
	qui eststo

	
	
*****************************************************
** Column 5 -- log(re78+1), dummy treat, linear re75
*****************************************************
		
	* calculate elasticity b/w re78 and re75
	est use "results/log_lin"
	nlcom _b[re75] * mean_re75, post
	scalar e_lnre78_re75 = _b[_nl_1]
	scalar e_lnre78_re75_SE = _se[_nl_1]
	
	* semi-elasticity for treat using exponential approximation
	est use "results/log_lin"
	nlcom (exp(_b[treat] - 0.5*_se[treat]^2) - 1), post
	scalar e_lnre78_treat_exp = _b[_nl_1]
	scalar e_lnre78_treat_exp_SE = _se[_nl_1]
	
	* store variables
	est use "results/log_lin"	
	qui estadd scalar elas = e_lnre78_re75
	qui estadd scalar elas_SE = e_lnre78_re75_SE
	qui estadd scalar semielas =  e_lnre78_treat_exp
	qui estadd scalar semielas_SE =  e_lnre78_treat_exp_SE
	qui eststo
	
	
	
*****************************************************
** Column 6 -- log(re78+1), dummy treat, log(re75)
*****************************************************
	
	* calculate elasticity b/w re78 and re75
	est use "results/log_log"
	nlcom _b[lnre75], post
	scalar e_lnre78_re75 = _b[_nl_1] 
	scalar e_lnre78_re75_SE = _se[_nl_1] 
	
	* semi-elasticity for treat using exponential approximation
	est use "results/log_log"	
	nlcom (exp(_b[treat] - 0.5*_se[treat]^2) - 1), post
	scalar e_lnre78_treat_exp = _b[_nl_1]
	scalar e_lnre78_treat_exp_SE = _se[_nl_1]
	
	* store variables
	est use "results/log_log"	
	qui estadd scalar elas = e_lnre78_re75
	qui estadd scalar elas_SE = e_lnre78_re75_SE
	qui estadd scalar semielas =  e_lnre78_treat_exp
	qui estadd scalar semielas_SE =  e_lnre78_treat_exp_SE
	qui eststo

	
	
******************************************	
* write table to file	
******************************************
#delimit ;
esttab using "results/reg1.tex", replace nonotes substitute(\_ _)
nogaps compress b(a3) se(a3) ar2(a3) star(* 0.10 ** 0.05 *** 0.01) 
eqlabel(none) nodep label nomtitles
mgroups("$\text{Earnings}_{78}$" "$\text{arcsinh(Earnings}_{78}\text{)}$" 
"$\ln\text{(Earnings}_{78}+1\text{)}$", pattern(1 0 1 0 1 0) 
prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
scalars("aic AIC"
		"bic BIC"
		"ll Log-likelihood"
		"a \hline"
		"b $\text{\textbf{Calculated (semi-)elasticities:}}$"
		"semielas2 $\xi(\text{Earnings}_{78},\text{ NSW})$"
		"semielas2_SE \hspace{5mm} Std. Err."
		"elas $\xi(\text{Earnings}_{78},\text{ Earnings}_{75})$"
		"elas_SE \hspace{5mm} Std. Err."
		"semielas $\tilde{P}(\text{Earnings}_{78},\text{ NSW})/100$"
		"semielas_SE \hspace{5mm} Std. Err."
		"semielas_ihs $\bar{P}(\text{Earnings}_{78},\text{ NSW})/100$"
		"semielas_ihs_SE \hspace{5mm} Std. Err.")
		sfmt(a3);
#delimit cr


#delimit ;
esttab using "results/reg1.tex", replace nonotes substitute(\_ _)
	nogaps compress b(a3) se(a3) ar2(a3) star(* 0.10 ** 0.05 *** 0.01) 
	eqlabel(none) nodep label nomtitles
	mgroups("$\text{Earnings}_{78}$" "$\text{arcsinh(Earnings}_{78}\text{)}$" 
	"$\ln\text{(Earnings}_{78}+1\text{)}$", pattern(1 0 1 0 1 0) 
	prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
	stats(N r2_a aic bic ll a b semielas2 semielas2_SE elas elas_SE 
		  semielas semielas_SE semielas_ihs semielas_ihs_SE, 
	labels(	"Obs"
			"Adj. R$^2$"
			"AIC" 
			"BIC" 
			"Log-likelihood"
			"\hline"
			"$\text{\textbf{Calculated (semi-)elasticities:}}$"
			"$\xi(\text{Earnings}_{78},\text{ NSW})$"
			" "
			"$\xi(\text{Earnings}_{78},\text{ Earnings}_{75})$"
			" "
			"$\tilde{P}(\text{Earnings}_{78},\text{ NSW})/100$"
			" "
			"$\bar{P}(\text{Earnings}_{78},\text{ NSW})/100$"
			" "
			)
	layout("@" "@" "@" "@" "@" "@" "@" "@" "(@)" "@" "(@)" "@" "(@)" "@" "(@)"))
	;

