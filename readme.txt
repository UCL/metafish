readme.txt for metafish
=======================

Folder ado
Description
	adofiles for other programs
Contents (packages)
	pllf - compute profile log likelihood
	metafish (unpublished) - two-stage Poisson meta-analysis
	simrun (unpublished) - utility to run simulation studies
	simsum - compute performance measures for simulation studies
	siman - tabulate and graph performance measures for simulation studies

Folder toydata
Description
	data and analyses for the hypothetical data of two studies with grouped Poisson outcome
Inputs
	toydata.dta - the data in Stata format
	toydata.do - analyse by 1-stage and 2-stage Normal. Draw PLLFs with their Normal approximations.
Outputs	
	toydata.log
	toydata_pllf.dta
	toydata_pllf.gph

Folder brcancer
Description
	Analyses of German breast cancer data to illustrate bias due to 2-stage Normal approximation and how 2-stage Poisson analysis can fix it
Inputs
	siloes.do - repeatedly split the data into siloes and analyse by 1-stage and 2-stage Normal. Extract one typical data set and also analyse it by 2-stage Poisson.
Outputs
	siloes.log
	siloes25_12.dta - the one typical data set

Folder simulation
Description
	Simulation program and results
Inputs
	twoSPsimprog.do - defines programs twoSPsimgen to generate data and twoSPsimana to analyse them
	twoSPsimrun.do - uses simrun package to run simulation and compute graphs & tables of performance measures
	twoSPsimrun2.do - as twoSPsimrun.do, but different methods
	twoSPcombine.do - takes both sets of results and combines into twoSPcombine.dta
	twoSPsimresults.do - produce all graphs and tables from twoSPcombine.dta
Outputs
	simrun_results (folder) - estimates datasets and random number generator states
	twoSPsimrun.log - main results
	twoSPcombine.dta - combined data set
	twoSPsimresults.log
Notes
	The command metafish had the earlier name of meta2p when the simulation was run. Similarly, the repository metafish was previously called TwoStagePoisson. Both names have been changed by editing the do files and the log files.
	"twoSP" in the file names reflects an earlier concept of the Poisson approximation method as a two-stage meta-analysis method, while we now present it as an aggregate data meta-analysis method.

Graph names from twoSPsimresults.do
	Fig 3 - discrep_ref
	Fig 4 - discrep_sd
	Fig 5 - homCE_1_pctbias
	Fig 6 - homCE_1_empse
	Fig 7 - homCE_1_cover
	Fig 8 - homCE_1_power
	Fig 9 - het_1_pctbias
	Fig 10 - het_1_empse
	Fig 11 - het_1_cover
	Fig 12 - het_1_power
	Fig 13 - [not from simulation]
	Fig 14 - weightingCE_1_pctbias
	Fig 15 - weightingCE_1_modelse
	Fig 16 - weightingRE_1_pctbias
	Fig 17 - weightingRE_1_modelse
	Fig 18 - centringPU_1_pctbias
	Fig 19 - centringPU_1_cover
	Fig 20 - centringWei_1_pctbias
	Fig 21 - centringWei_1_cover
	Fig 22 - compareC_1_pctbias
	Fig 23 - compareC_1_modelse
	Fig 24 - tauhet_1_mean