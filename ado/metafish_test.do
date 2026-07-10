/*
Test script for metafish
IW 17feb2026
Added more tests using 25 siloes data, 10jul2026
	and use mytest routine
*/

cd C:\ian\git\metafish\ado

prog drop _all
prog def mytest // prints two numbers and checks they are very close
args first second tol
di `first', `second'
if mi("`tol'") local tol 1E-7
assert abs(`first' - `second') < `tol'
end

set linesize 100
cap log close
log using metafish_test, replace text nomsg
which metafish

// Toy data: show metafish reproduces 1-stage Poisson

clear
input study z d py
1 0 40 1000 
1 1 10 1000
2 0 60 1000 
2 1 60 1000
end
poisson d z i.study, exp(py)
local true = _b[z]
local truese = _se[z]

gen lograte = log(d/py)
gen var = 1/d
reshape wide d py lograte var, i(study) j(z)
gen est = lograte1-lograte0
gen se = sqrt(var1+var0)
metan est se, nograph model(fixed)

metafish est se, d(d1 d0) 
mytest r(eff) `true'
mytest r(se_eff) `truese'

// Reversing groups gives same answer (reversed)

gen negest = -est
metafish negest se, d(d0 d1) 
mytest r(eff) -`true'
mytest r(se_eff) `truese'

// Check that options run

metafish est se, d(d1 d0) study(study) py(py1 py0) wt centre verb list irr eform re poissonopt(noheader intpoints(5)) wttol(20)

metafish est se, d(d1 d0) study(study) py(py1 py0) wt centre verb list irr eform poissonopt(robust) wttol(50)
* robust isn't allowed with mepoisson

// se not required if no wt
metafish est, d(d1 d0) 
cap noi metafish est, d(d1 d0) wt
assert _rc==198

// Now test on a larger data set
use ..\brcancer\siloes25_12_meta, clear
metafish b se, d(d1 d0) study(silo) py(p1 p0) wt centre re
local b = r(eff)
local se = r(se_eff)

* test same answers with different sort order
sort b
metafish b se, d(d1 d0) study(silo) py(p1 p0) wt centre re
mytest `b' r(eff)
mytest `se' r(se_eff)

* test same answers with different variable names
rename (*) (z*)
metafish zb zse, d(zd1 zd0) study(zsilo) py(zp1 zp0) wt centre re
mytest `b' r(eff)
mytest `se' r(se_eff)

use ..\brcancer\siloes25_12_meta, clear
metafish b se, d(d1 d0) study(silo) py(p1 p0) 
local b = r(eff)
local se = r(se_eff)

* test analogous answers with data multiplied by 4 (CE model only)
for var d? p?: replace X=X*4
replace se=se/2
metafish b se, d(d1 d0) study(silo) py(p1 p0)
mytest `b' r(eff)
mytest `se'/2 r(se_eff)

log close
