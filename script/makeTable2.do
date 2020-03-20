
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


* save number of zero obs per var/treatment group
qui count if re78==0 & treat==1
local re78_tz = r(N)
qui count if re75==0 & treat==1
local re75_tz = r(N)
qui count if re78==0 & treat==0
local re78_cz = r(N)
qui count if re75==0 & treat==0
local re75_cz = r(N)


** generate summary statistics table
#delimit ;
eststo clear ;
eststo treat: quietly estpost sum re78 re75 if treat==1;
	qui estadd local re78z = `re78_tz';
	qui estadd local re75z = `re75_tz';
	
eststo ctrl: quietly estpost sum re78 re75 if treat==0;
	qui estadd local re78z = `re78_cz';
	qui estadd local re75z = `re75_cz';

esttab treat ctrl using "results/sumstats.tex", replace nogaps label substitute(\_ _)
	cells("mean(pattern(1 1) fmt(1)) sd(pattern(1 1) fmt(1))") 
	mgroups("Treatment Group" "Control Group", pattern(0 1)               
	prefix(\multicolumn{@span}{c}{) suffix(})   
	span erepeat(\cmidrule(lr){@span})) 
	nonumbers collabels("Mean" "SD")
	scalars("re78z $\text{Earnings}_{78}=0$" "re75z $\text{Earnings}_{75}=0$")
	;
#delimit cr

