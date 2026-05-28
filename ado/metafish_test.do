/*
very basic test script for metafish
IW 17feb2026
*/

set linesize 100
cap log close
prog drop _all
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
assert abs(r(eff)-`true')<1E-10
assert abs(r(se_eff)-`truese')<1E-10

// Reversing groups gives same answer (reversed)

gen negest = -est
metafish negest se, d(d0 d1) 
assert abs(-r(eff)-`true')<1E-10
assert abs(r(se_eff)-`truese')<1E-10

// Check that options run

metafish est se, d(d1 d0) study(study) py(py1 py0) wt centre verb list irr eform re poissonopt(noheader intpoints(5)) wttol(20)
metafish est se, d(d1 d0) study(study) py(py1 py0) wt centre verb list irr eform poissonopt(robust) wttol(50)
* robust isn't allowed with mepoisson

// se not required if no wt
metafish est, d(d1 d0) 
cap noi metafish est, d(d1 d0) wt
assert _rc==198

log close
