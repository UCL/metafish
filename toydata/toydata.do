/*
toydata.do
IW 5/9/2025
for github 30/1/2026
minor revision 9apr2026
*/

// User-specific settings
cd C:\ian\git\metafish\toydata
adopath ++ C:\ian\git\metafish\ado
set scheme mrc

prog drop _all
version 18
cap log close
set linesize 100
log using toydata, text replace

version
which metafish

// Two-stage analyses
use toydata, clear
reshape wide d py, i(study) j(z)
gen est = log(d1/py1) - log(d0/py0)
gen se = sqrt(1/d1+1/d0)

* perpare to store results in other frame
frame change default
cap frame drop examples
frame create examples str8 method b se

* Poisson-approximation method
metafish est se, d(d1 d0) 
frame post examples ("Poisson") (r(eff)) (r(se_eff))

* Normal-approximation method
metan est se, nograph model(fixed)
frame post examples ("Normal") (r(eff)) (r(se_eff))

* Normal-approximation method with bias correction
gen estbc = log((d1+0.5)/py1) - log((d0+0.5)/py0)
gen sebc1 = sqrt(1/(d1+0.5)+1/(d0+0.5))
gen sebc2 = sqrt(d1/(d1+0.5)^2+d0/(d0+0.5)^2)
metan estbc sebc1, nograph model(fixed)
frame post examples ("N-BC-1") (r(eff)) (r(se_eff))
* both SEs decrease, but that of study 1 decreases more (sparser)

metan estbc sebc2, nograph model(fixed)
* both SEs decrease, but that of study 1 decreases more (sparser)
frame post examples ("N-BC-2") (r(eff)) (r(se_eff))

frame examples {
	gen lci = b - invnorm(.975)*se
	gen uci = b + invnorm(.975)*se
	gen rr = exp(b)
	gen rrlci = exp(lci)
	gen rruci = exp(uci)
	l method b lci uci rr*
	format b se lci uci %6.3f
	format rr* %6.2f
	l method b lci uci rr*
}

// Compute PLLFs: a bit slow
use toydata, clear
forvalues i=1/2 {
	poisson d z if study==`i', exposure(py) 
	local b`i'=_b[z]
	local se`i'=_se[z]
	pllf, profile(z) range(-3 2) nograph diff norm n(500) ///
		gen(beta`i' pll`i' norm`i'): ///
		poisson d z if study==`i', exposure(py) 
	label var pll`i' "Study `i', exact PLLF"
	label var norm`i' "Study `i', normal approx"
	label var beta`i' "beta: log rate ratio for treatment"
}
assert beta1==beta2
drop study z d py beta2
rename beta1 beta
save toydata_pllf, replace

// Rearrange data
use toydata_pllf, clear
gen pll12=pll1+pll2
gen norm12=norm1+norm2
summ pll12
replace pll12=pll12-r(max)
summ norm12
replace norm12=norm12-r(max)

reshape long pll norm, i(beta) j(method)
label var pll "Exact profile log-likelihood function"
label var norm "Normal approximation"
label def method 1 "Study 1" 2 "Study 2" 12 "Meta-analysis"
label val method method

// Draw graph
foreach var in pll norm {
	replace `var'=. if `var'<-15
}
keep if inrange(beta,-3,1)
line pll norm beta, lc(black =) lp(solid dash) lw(*2 =) ///
	ytitle(Relative profile log-likelihood,size(small)) ///
	legend(col(1) size(small)) ///
	xsize(6) ysize(7) xtitle(,size(small)) ytitle(,size(small)) ///
	ylab(,labsize(small)) xlab(,labsize(small)) ///
	by(method, note("") col(1) legend(pos(6))) ///
	subtitle(,size(small)) xline(-0.355, lcol(gray)) ///
	name(toydata_pllf,replace) saving(toydata_pllf,replace)

log close
