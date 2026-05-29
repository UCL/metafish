/*
siloes.do
Multiply split the German breast cancer data set into siloes, 
	analyse by IPD MA,
	and compare metafish with metan for one of the splits
IW 23jan2026
minor revision 9apr2026, 24apr2026
bug fix to make the data splits reproducible, 28may2026
*/

// User-specific settings
cd C:\ian\git\metafish\brcancer
adopath ++ C:\ian\git\metafish\ado
set scheme mrc

version 18
cap log close
set linesize 100
log using siloes, text replace

version
which metafish

local reps 50
set seed 461860

webuse brcancer, clear
stset rectime, fail(censrec)
stcox hormon, nohr
global b = _b[hormon]
global se = _se[hormon]

if c(frame)=="results" {
	di as error "Can't start in frame results"
	exit 499
}
cap frame drop results
frame create results siloes rep b se
foreach siloes in 2 5 10 25 {
	frame post results (`siloes') (0) ($b) ($se)
	forvalues rep=1/`reps' {
		if `rep'==1 _dots 0, title("`siloes' siloes: simulation running (`reps' repetitions)")
		_dots `rep' 0
		cap drop random silo
		sort id // restore original sort order
		gen double random = runiform() // double is included to make ties much less likely
		sort random id // id is included in order to sort consistently not arbitrarily, in the unlikely event of a tie on random
		gen silo = int(`siloes'*(_n-1)/_N)+1
		qui ipdmetan, study(silo) nograph: stcox hormon
		frame post results (`siloes') (`rep') (r(eff)) (r(se_eff))
		if `siloes'==25 & `rep'==12 qui save siloes`siloes'_`rep', replace
		* a fairly typical data set
	}
	frame results: label def siloes `siloes' "`siloes' siloes", modify	
}

// Table 4 for paper: average logHR and SE across splits
frame results {
	sort siloes rep
	label val siloes siloes
	save siloes, replace
	gen siloes2 = cond(rep,siloes,1)
	table (siloes2), stat(mean b) stat(semean b) stat(mean se) stat(semean se) nototal nformat(%6.3f)
}

// Figure 2 for paper: SE vs logHR across splits
frame results : ///
	scatter se b if rep==0, ms(X) msize(*2) ///
	|| scatter se b if rep>0, /*mlab(rep)*/ ms(oh) ///
	|| scatter se b if rep==12 & siloes==25, /*mlab(rep)*/ ms(O) msize(*1.5) ///
	by(siloes, note("")) legend(order(1 "Without splitting" 2 "2-stage Normal approx")) ///
	xline($b) yline($se) ///
	xtitle(Log hazard ratio) ytitle(Standard error) ///
	name(siloes, replace) $PPT


// explore a typical data set

use siloes25_12, clear

cap frame drop res25_12
frame create res25_12 str4 method b se
// one-stage
stcox hormon, nohr
frame post res25_12 ("1Sun") (_b[hormon]) (_se[hormon])

// one-stage
stcox hormon i.silo, nohr
frame post res25_12 ("1Sadj") (_b[hormon]) (_se[hormon])

tab silo hormon if _d
byvar silo, b(hormon) se(hormon) unique gen: stcox hormon
egen d1 = sum(_d * hormon), by(silo)
egen d0 = sum(_d * (1-hormon)), by(silo)
egen p1 = sum(hormon), by(silo)
egen p0 = sum((1-hormon)), by(silo)
keep if !mi(B)
keep silo B S d? p?
rename B b
rename S se
format b se %6.3f

// Table 5 for paper: the selected data set
l, noo clean

// two-stage Normal
metan b se
frame post res25_12 ("2SN") (r(eff)) (r(se_eff))
// two-stage Poisson
metafish b se, d(d1 d0) py(p1 p0)
frame post res25_12 ("2SP") (r(eff)) (r(se_eff))
metafish b se, d(d1 d0) py(p1 p0) wt
frame post res25_12 ("2SPW") (r(eff)) (r(se_eff))

frame res25_12: l, noo clean

log close
