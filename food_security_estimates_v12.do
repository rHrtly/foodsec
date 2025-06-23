
clear all
macro drop _all
global version = 12

set more off, permanently
set varabbrev off, permanently

cd "G:\My Drive\Interchange\Research\FoodSec\"

local dt = `"`=subinstr(subinstr("${S_DATE}_${S_TIME}"," ","",.),":","_",.)'"'
capture log close FSEC_est
log using "LogFiles\FSEC_estimates_`dt'", replace text name(FSEC_est)

/*******************************************************************************
********************************************************************************

  "The Persistence of Food Security Status Across Generations"
  
	Robert Paul Hartley
	Jaehyun Nam
	Christopher Wimer
	
	June 2025
	
********************************************************************************
*******************************************************************************/

use "food_security_v${version}_estimation_sample.dta", clear

**# Global macro setup
********************************************************************************
{
/* Family fixed effects */
tab Siblings_Cousins, g(S_C_)
unab fam: S_C_1-S_C_`r(r)'
global FAM `fam'

/* State fixed effects */
tab STATEALPHA, g(SD)
unab sd: SD1-SD`r(r)'
global SD `sd'

/* Year fixed effects */
tab YEAR, g(YD)
unab yd: YD1-YD`r(r)'
global YD `yd'

/* Age profiles */
global age AGE AGE2 PAGE PAGE2 AGE_h AGE_h2

/* Control variables */
global controls spm_historical_under100_S unemployment_rate_S FEMALE		 ///
	BlackNonHispanic OtherNonHispanic Hispanic RaceEthnicity_imputed		 ///
	NC1 NC2 NC3 NC4

/* Mean childhood family earnings-to-poverty ratio and net wealth/equity */
global pov EarnPoverty_c lnequity_c

/* Main covariates, not including family fixed effects */
global X ${age} ${controls} ${pov} ${SD} ${YD}

/* Main 2SLS covariates, not including state and family fixed effects */
global ZX ${age} ${controls} ${pov} ${YD}

/* Aggregated covariates: Means across adult observation years */
foreach g in age controls {
	global m_`g'
	foreach x of global `g' {
		global m_`g' ${m_`g'} m_`x'
	}
}
global m_X ${m_age} ${m_controls} ${pov} ${SD} ${YD} ${FAM}

/* SNAP policy index instruments; mean childhood and contemporaneous adult */
global IV snap_eli_cit_usda_S_c snap_tra_oap_usda_S_c snap_tra_crt_usda_S_c
global iv snap_eli_cit_usda_S snap_tra_oap_usda_S snap_tra_crt_usda_S

/* Core instrument set: SNAP and TANF benefit generosity, SNAP policy index */
global IV1 real_snap_bene_c real_tanf_bene_c ${IV}
global iv1 real_snap_bene real_tanf_bene ${iv}

/* Extended 7-instrument set with additional SNAP policy index variables */
global IV2 ${IV1} snap_eli_bbc_usda_S_c snap_tra_rep_usda_S_c
global iv2 ${iv1} snap_eli_bbc_usda_S snap_tra_rep_usda_S

/* Extended 9-instrument set with additional SNAP policy index variables */
global IV3 ${IV2} snap_eli_veh_usda_S_c snap_stg_fpr_usda_S_c
global iv3 ${iv2} snap_eli_veh_usda_S snap_stg_fpr_usda_S

/* Extended 12-instrument set with additional SNAP policy index variables */
global IV4 ${IV3} outreach_SPD_c snap_eli_vex_usda_S_c snap_stg_ebt_usda_S_c
global iv4 ${iv3} outreach_SPD snap_eli_vex_usda_S snap_stg_ebt_usda_S

/* Instruments for Lee-Solon (2009) age adjustment */
global CAIV ${IV1}
global caiv ${iv1}
forvalues i = 1/4 {
	foreach x of global IV1 {
		global CAIV ${CAIV} C`i'_`x'
	}
	foreach x of global iv1 {
		global caiv ${caiv} C`i'_`x'
	}
}

/* OLS controls for Lee-Solon (2009) age adjustment */
local all_X ${X} ${FAM}
local except ${age}
local leesolon: list all_X - except
global leesolon `leesolon' CA1 CA2 CA3 CA4 MA1 MA2 MA3 MA4

/* 2SLS controls for Lee-Solon (2009) age adjustment */
local all_ZX ${ZX}
local except ${age}
local leesolonZ: list all_ZX - except
global leesolonZ `leesolonZ' CA1 CA2 CA3 CA4 MA1 MA2 MA3 MA4

/* Food module item outcomes, recoded as binary indicator variables */
global FOOD FOOD01 FOOD02 FOOD03 FOOD04 FOOD05 FOOD06 FOOD07 FOOD08 FOOD09	 ///
	FOOD10 FOOD11 FOOD12 FOOD13 FOOD14 FOOD15 FOOD16 FOOD17 FOOD18

/* Continuous food-security-related outcome variables */
global YC FSRAW latentFS latentFS_pr foodspend_pr
}

**# Programs for organizing estimation results
********************************************************************************
{
/* Correlations with Gelbach (2016) decompositions */
cap program drop fb1x2
program fb1x2, rclass
	quietly {
		tempname b0 V0 b1 V1 D E W S P R
		matrix `D' = e(Delta)
		matrix `E' = `D'[1,`=colsof(`D')']
		matrix `D' = `D'[1,1..`=colsof(`D')-1']
		matrix `W' = e(Covdelta)
		matrix `b0' = e(b1base)
		matrix `V0' = e(V1base)
		matrix `b1' = e(b1full)
		matrix `V1' = e(V1full)
		forvalues i = 1/`=colsof(`D')' {
			matrix `S' = nullmat(`S') \ `D'[1,`i'] \ `W'[`i',`i']^.5
			matrix `P' = nullmat(`P') \ `D'[1,`i'] / (`b0'[1,1] - `b1'[1,1]) ///
				\ 2 * (1 - normal(abs(`D'[1,`i'] / `W'[`i',`i']^.5)))
		}
		local Es = `W'[`=colsof(`W')',`=colsof(`W')']^.5
		if (`e(N)' == 6476) {
			local n = 1701
			local N = 6476
		}
		else {
			xtsum `e(depvar)' if e(sample)
			local n = `r(n)'
			local N = `r(N)'
		}
		matrix `R' = (`b0'[1,1], `b1'[1,1]) \  (`V0'[1,1]^.5, `V1'[1,1]^.5)	 ///
			\ (J(`=colsof(`D')*2',1,.), `S') \ (J(2,1,.),					 ///
			(`E'[1,1] \ `Es')) \ J(1,2,`n') \ J(1,2,`N')					 ///
			\ 2 * (1 - t(`=e(N_clust)-1', abs(`b0'[1,1] / `V0'[1,1]^.5)),	 ///
			1 - t(`=e(N_clust)-1', abs(`b1'[1,1] / `V1'[1,1]^.5)))			 ///
			\ (J(`=colsof(`D')*2',1,.), `P')
		matrix `R' = `R'[1..2,1] \ `R'[1...,2]
		local rnames
		forvalues i = 1/`=colsof(`D')' {
			local rnames = "`rnames' `: word `i' of `e(groupnames)'' se"
		}
		local rnames2
		forvalues i = 1/`=colsof(`D')' {
			local rnames2 = "`rnames2' pctd-`: word `i' of `e(groupnames)'' pv"
		}
		matrix rownames `R' = uncond se cond se `rnames' Total se n N pvalue ///
			`rnames2'
		matrix coleq `R' = `e(depvar)'
		matrix colnames `R' = `: word 1 of `: coleq `D'''
		noisily matrix list `R', format(%6.3f) title("Results")
		return matrix results = `R'
	}
end

/* IV estimates including 2SLS, weak IV tests, and LIML */
capture program drop iv_est
program iv_est, eclass
	noisily {
		tempname b V r F R
		matrix `R' = 1
		
		local first
		if (length(`"=subinword("`0'", "first", "", 1)"') == length(`"`0'"')) {
			local first = " first"
		}
		if (!regexm(`"`0'"', "savefirst")) {
			local first = " savefirst"
		}
		
		/* IV estimate with first stage */
		ivreg2 `0'`first'
		matrix `b' = e(b)
		matrix `V' = e(V)
		matrix `R' = _b[`e(instd)'] \ _se[`e(instd)'] \ `e(idstat)' \		 ///
			`e(idp)' \ `e(j)' \ `e(jp)' \ `e(widstat)' \ e(first)
		matrix coleq `R' = `: word 1 of `0''
		matrix colnames `R' = `e(instd)'
		matrix rownames `R' = b se id idp j jp wid rmse sheapr2 pr2 F df	 ///
			df_r pvalue SWF SWFdf1 SWFdf2 SWFp SWchi2 SWchi2p SWr2 APF		 ///
			APFdf1 APFdf2 APFp APchi2 APchi2p APr2
		matrix roweq `R' = IV

		/* First-stage results */
		estimates replay _ivreg2_`e(instd)'
		matrix `r' = r(table)
		matrix `F' = J(24, 1, .)
		forvalues i = 1/`=colsof(`r')' {
			matrix `F'[`=(`i'-1)*2+1', 1] = `r'[1..2, `i']
		}
		matrix rownames `F' = snap snapse tanf tanfse cit citse oap oapse	 ///
			crt crtse bbc bbcse rep repse veh vehse fpr fprse out outse vex	 ///
			vexse ebt ebtse
		matrix roweq `F' = FirstStage
		matrix `R' = `R' \ `F'
		
		/* Weak IV test (95% confidence) */
		ivreg2 `=regexreplace("`0'", "partial\(.*?\)", "")'
		weakivtest, level(.05)
		matrix `r' = `r(F_eff)' \ `r(c_TSLS_5)' \ `r(c_TSLS_10)' \			 ///
			`r(c_TSLS_20)' \ `r(c_TSLS_30)' \ `r(c_LIML_5)' \				 ///
			`r(c_LIML_10)' \ `r(c_LIML_20)' \ `r(c_LIML_30)'
		matrix rownames `r' = F_eff TSLS_5_95 TSLS_10_95 TSLS_20_95			 ///
			TSLS_30_95 LIML_5_95 LIML_10_95 LIML_20_95 LIML_30_95
		matrix roweq `r' = WeakIVtest
		matrix `R' = `R' \ `r'
		
		/* Weak IV test (90% confidence) */
		weakivtest, level(.1)
		matrix `r' = `r(c_TSLS_5)' \ `r(c_TSLS_10)' \ `r(c_TSLS_20)' \		 ///
			`r(c_TSLS_30)' \ `r(c_LIML_5)' \ `r(c_LIML_10)' \				 ///
			`r(c_LIML_20)' \ `r(c_LIML_30)'
		matrix rownames `r' = TSLS_5_90 TSLS_10_90 TSLS_20_90 TSLS_30_90	 ///
			LIML_5_90 LIML_10_90 LIML_20_90 LIML_30_90
		matrix roweq `r' = WeakIVtest
		matrix `R' = `R' \ `r'
		
		/* Weak-IV-robust estimates */
		weakiv ivreg2 `0'
		local s1 = strpos("`e(clr_cset)'", "[")
		local s2 = strpos("`e(clr_cset)'", "]")
		local c = strpos("`e(clr_cset)'", ",")
		local s1a = strpos("`e(ar_cset)'", "[")
		local s2a = strpos("`e(ar_cset)'", "]")
		local ca = strpos("`e(ar_cset)'", ",")
		matrix `r' = `e(clr_p)' \ `e(k_p)' \ `e(k_df)' \ `e(j_p)' \			 ///
			`e(j_df)' \ `e(kj_p)' \ `e(ar_p)' \ `e(ar_df)' \				 ///
			`=real(substr("`e(clr_cset)'", `=`s1'+1', `=`c'-`s1'-1'))' \	 ///
			`=real(substr("`e(clr_cset)'", `=`c'+1', `=`s2'-`c'-1'))' \		 ///
			`=real(substr("`e(ar_cset)'", `=`s1a'+1', `=`ca'-`s1a'-1'))' \	 ///
			`=real(substr("`e(ar_cset)'", `=`ca'+1', `=`s2a'-`ca'-1'))'
		weakiv ivreg2 `0' level(90)
		local s1 = strpos("`e(clr_cset)'", "[")
		local s2 = strpos("`e(clr_cset)'", "]")
		local c = strpos("`e(clr_cset)'", ",")
		local s1a = strpos("`e(ar_cset)'", "[")
		local s2a = strpos("`e(ar_cset)'", "]")
		local ca = strpos("`e(ar_cset)'", ",")
		matrix `r' = `r' \ 													 ///
			`=real(substr("`e(clr_cset)'", `=`s1'+1', `=`c'-`s1'-1'))' \	 ///
			`=real(substr("`e(clr_cset)'", `=`c'+1', `=`s2'-`c'-1'))' \		 ///
			`=real(substr("`e(ar_cset)'", `=`s1a'+1', `=`ca'-`s1a'-1'))' \	 ///
			`=real(substr("`e(ar_cset)'", `=`ca'+1', `=`s2a'-`ca'-1'))'
		matrix rownames `r' = clr_p k_p k_df j_p j_df kj_p ar_p ar_df		 ///
			clr_ll_95 clr_ul_95 ar_ll_95 ar_ul_95 clr_ll_90 clr_ul_90		 ///
			ar_ll_90 ar_ul_90
		matrix roweq `r' = WeakIVrobust
		matrix `R' = `R' \ `r'
		
		/* LIML estimates */
		ivreg2 `0'`first' liml
		matrix `r' = _b[`e(instd)'] \ _se[`e(instd)']
		matrix rownames `r' = b_liml se_liml
		matrix roweq `r' = LIML
		matrix `R' = `R' \ `r'
		
		/* Robust/non-clustered */
		local rreg = `"`=regexreplace("`0'", "cluster\(.*?\)", "rob")'"'
		ivreg2 `rreg'
		local wid = `e(widstat)'
		local rreg = `"`=regexreplace(`"`rreg'"', "partial\(.*?\)", "")'"'
		ivreg2 `rreg'
		weakivtest
		matrix `r' = `wid' \ `r(F_eff)'
		matrix rownames `r' = wid_rob F_eff_rob
		matrix roweq `r' = Robust_only
		matrix `R' = `R' \ `r'
		
		/* Number of individuals and observations */
		xtsum CLWT if e(sample)
		matrix `r' = `r(n)' \ `r(N)'
		matrix rownames `r' = n N
		matrix roweq `r' = Observations
		matrix `R' = `R' \ `r'
		
		/* Return results */
		ereturn post `b' `V'
		ereturn matrix results = `R'
		ereturn local cmdline `"iv_est `0'"'
	}
end
}

********************************************************************************
**# Manuscript Tables
********************************************************************************

**# Table 1
// 	Intergenerational Correlations of Food Security Status
********************************************************************************
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop Table1`i'
	foreach y in FSRAW mglsecure lowsecure vlwsecure {
		tempname R
		foreach x in mglsecure lowsecure vlwsecure childinsecure {
			tempname results
			b1x2 `y' if touse_sc == 1 & `s' == 1 [aw = CLWT],				 ///
				x1all(Any_`x'_c) x2all(${X} ${FAM})							 ///
				x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :		 ///
				earn = EarnPoverty_c : equity = lnequity_c :				 ///
				family = ${FAM}) x1only(Any_`x'_c) cluster(Siblings_Cousins)
			fb1x2
			matrix `results' = (`=r(results)[1,1]' \ `=r(results)[2,1]'),	 ///
				(`=r(results)[3,1]' \ `=r(results)[4,1]')
			matrix coleq `results' = `y'
			matrix colnames `results' = uncond cond
			matrix rownames `results' = `x' se
			matrix `R' = nullmat(`R') \ `results'
			matrix drop `results'
		}
		matrix Table1`i' = nullmat(Table1`i'), `R'
		matrix drop `R'
	}
}
matrix list Table1A, format(%6.3f) title("Full sample")
matrix list Table1B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("Table1") modify
putexcel A1 = matrix(Table1A), names
putexcel A12 = matrix(Table1B), names

**# Table 2
// 	2SLS Estimates of Childhood Marginal, Low, or Very Low Food Security
// 	Effects on Early Adulthood Food Security Status
********************************************************************************
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop Table2`i'
	foreach y in FSRAW mglsecure lowsecure vlwsecure {
		tempname results
		iv_est `y' (Any_mglsecure_c = ${IV1}) ${ZX} ${iv1}					 ///
			if touse_sc == 1 & `s' == 1 [aw = CLWT], partial(${ZX} ${iv1})	 ///
			cluster(Siblings_Cousins)
		matrix `results' = e(results)
		matrix Table2`i' = nullmat(Table2`i'), (							 ///
				`results'["IV:b".."IV:se",1] \								 ///
				`results'["Robust_only:wid_rob",1] \						 ///
				`results'["WeakIVtest:F_eff",1] \							 ///
				`results'["WeakIVrobust:k_p",1] \							 ///
				`results'["IV:id".."IV:jp",1]								 ///
			)
		matrix drop `results'
	}
}
matrix list Table2A, format(%6.3f) title("Full sample")
matrix list Table2B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("Table2") modify
putexcel A1 = matrix(Table2A), names
putexcel A13 = matrix(Table2B), names

**# Table 3.
// 	Intergenerational Conditional Correlations and 2SLS Estimates for
// 	Childhood Food Security Status and Early Adulthood Outcomes Potentially
// 	Related to Adult Food Insecurity
********************************************************************************
capture matrix drop Table3A
capture matrix drop Table3B
foreach y in EARNPOV_pr LNREALWAGE foodspend_pr GOODHEALTH ANYCOLLEGE FSSNAP {
	tempname R
	foreach x in FSRAW mglsecure lowsecure vlwsecure childinsecure {
		tempname results
		if ("`x'" == "FSRAW") local a = "Mean"
		else local a = "Any"
		b1x2 `y' if touse_sc == 1 & all == 1 [aw = CLWT],					 ///
			x1all(`a'_`x'_c) x2all(${X} ${FAM})								 ///
			x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :			 ///
			earn = EarnPoverty_c : equity = lnequity_c :					 ///
			family = ${FAM}) x1only(`a'_`x'_c) cluster(Siblings_Cousins)
		fb1x2
		matrix `results' = `=r(results)[3,1]' \ `=r(results)[4,1]'
		matrix roweq `results' = OLS
		matrix rownames `results' = `x' se
		matrix coleq `results' = `y'
		matrix colnames `results' = cond
		matrix `R' = nullmat(`R') \ `results'
		matrix drop `results'
	}
	matrix Table3A = nullmat(Table3A), `R'
	matrix drop `R'
	tempname results
	iv_est `y' (Any_mglsecure_c = ${IV1}) ${ZX} ${iv1}						 ///
		if touse_sc == 1 & all == 1 [aw = CLWT], partial(${ZX} ${iv1})		 ///
		cluster(Siblings_Cousins)
	matrix `results' = e(results)
	matrix Table3B = nullmat(Table3B), (									 ///
			`results'["IV:b".."IV:se",1] \									 ///
			`results'["Robust_only:wid_rob",1] \							 ///
			`results'["WeakIVtest:F_eff",1] \								 ///
			`results'["WeakIVrobust:k_p",1] \								 ///
			`results'["IV:id".."IV:jp",1] \									 ///
			`results'["Observations:n".."Observations:N",1]					 ///
		)
	matrix drop `results'
}
matrix list Table3A, format(%6.3f) title("Conditional correlations")
matrix list Table3B, format(%6.3f) title("2SLS: Any marginal or food insecure")
putexcel set FoodSec_Results_v${version}, sheet("Table3") modify
putexcel A1 = matrix(Table3A), names
putexcel A14 = matrix(Table3B), names

********************************************************************************
**# Manuscript Figures
********************************************************************************

**# Figure 1 (see below for estimates with extended dataset)
// 	Percent of Children by Food Security and Poverty Status
********************************************************************************

**# Figure 2 (see below for estimates with extended dataset)
// 	Within-Family Persistence in Food Security Status,
// 	Earnings, and Food Spending Relative to the 1997 Survey Year
********************************************************************************

**# Figure 3
// 	Intergenerational Correlations of Marginal, Low, or Very Low Food Security,
// 	by Samples for Varying Earnings-to-Needs Status Thresholds in Childhood
********************************************************************************
capture matrix drop Figure3A
capture matrix drop Figure3B
foreach i of numlist 1/100 {
	tempname R_A mA R_B mB
	
	matrix `R_A' = `i'
	foreach m in Any Mean {
		reg mglsecure `m'_mglsecure_c ${age} ${controls} ${pov} ${SD} ${YD}	 ///
			if touse_sc == 1 & EP_c <= `i' [aw = CLWT],						 ///
			cluster(Siblings_Cousins)
		matrix `R_A' = nullmat(`R_A'),										 ///
			`=r(table)[1,1]', `=r(table)[5,1]', `=r(table)[6,1]'
	}
	tabstat Any_mglsecure_c Mean_mglsecure_c [aw = CLWT]					 ///
		if touse_sc == 1 & EP_c <= `i' & EP_c != ., save s(mean sd)
	tabstatmat `mA'	
	count if touse_sc == 1 & EP_c <= `i' & EP_c != .
	matrix Figure3A = nullmat(Figure3A) \									 ///
		(`R_A', `mA'[1...,1]', `mA'[1...,2]', r(N))
	matrix drop `R_A' `mA'
	
	matrix `R_B' = `i'
	foreach m in Any Mean {
		reg mglsecure `m'_mglsecure_c ${age} ${controls} ${pov} ${SD} ${YD}	 ///
			if touse_sc == 1 & EP_c >= `i' [aw = CLWT],						 ///
			cluster(Siblings_Cousins)
		matrix `R_B' = nullmat(`R_B'),										 ///
			`=r(table)[1,1]', `=r(table)[5,1]', `=r(table)[6,1]'
	}
	tabstat Any_mglsecure_c Mean_mglsecure_c [aw = CLWT]					 ///
		if touse_sc == 1 & EP_c >= `i' & EP_c != ., save s(mean sd)
	tabstatmat `mB'	
	count if touse_sc == 1 & EP_c >= `i' & EP_c != .
	matrix Figure3B = nullmat(Figure3B) \									 ///
		(`R_B', `mB'[1...,1]', `mB'[1...,2]', r(N))
	matrix drop `R_B' `mB'
	
}
local c p any_b any_ll any_ul mean_b mean_ll mean_ul any any_s mean mean_s N
matrix colnames Figure3A = `c'
matrix colnames Figure3B = `c'
matrix list Figure3A, format(%6.3f)
matrix list Figure3B, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure3A") modify
putexcel A1 = matrix(Figure3A), colnames
putexcel set FoodSec_Results_v${version}, sheet("Figure3B") modify
putexcel A1 = matrix(Figure3B), colnames

/* Percentiles corresponding to earnings-to-needs ratios of 1, 2, and 4 */
capture matrix drop Figure3p
sum CLWT if touse_sc == 1 & EP_c != .
local T = `r(sum)'
foreach i of numlist 1/2 4 {
	sum CLWT if touse_sc == 1 & EARNPOV_c < `i' & EP_c != .
	local p = `r(sum)'
	matrix Figure3p = nullmat(Figure3p) \ 100 * (`i', `p' / `T')
}
matrix colnames Figure3p = "EARNPOV" "percentile"
matrix list Figure3p, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure3p") modify
putexcel A1 = matrix(Figure3p), colnames

**# Figure 4
// 	Intergenerational Correlation in Early Adult Marginal, Low, or
// 	Very Low Food Insecurity by Childhood Measurement Timing and
// 	Severity of Insecurity
********************************************************************************
capture matrix drop Figure4
foreach x in mglsecure lowsecure vlwsecure {
	foreach a in "c" "c0_5" "c6_11" "c12_17" {
		b1x2 mglsecure if touse_sc == 1 & all == 1 [aw = CLWT],				 ///
			x1all(Any_`x'_`a') x2all(${X} ${FAM})							 ///
			x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :			 ///
			earn = EarnPoverty_c : equity = lnequity_c :					 ///
			family = ${FAM}) x1only(Any_`x'_`a') cluster(Siblings_Cousins)
		fb1x2
		matrix Figure4 = nullmat(Figure4), r(results)
	}
}
matrix coleq Figure4 = ""
matrix list Figure4, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure4") modify
putexcel A1 = matrix(Figure4), names
putexcel A1 = "stat"

**# Figure 5
// 	Conditional Correlations and 2SLS Estimates of Childhood Marginal,
// 	Low, or Very Low Food Security Effects on Food Security Module
// 	Questionnaire Responses
********************************************************************************
capture matrix drop Figure5
foreach s in all children {
	foreach y of global FOOD {
		tempname results
		reg `y' Any_mglsecure_c ${X} ${FAM} if touse_sc == 1 & `s' == 1		 ///
			[aw = CLWT], cluster(Siblings_Cousins)
		matrix `results' = _b[Any_mglsecure_c] \ _se[Any_mglsecure_c]
		matrix roweq `results' = "OLS"
		matrix rownames `results' = "cond" "se"
		matrix colnames `results' = `y'`=substr("`s'", 1, 1)'
		iv_est `y' (Any_mglsecure_c = ${IV1}) ${ZX} ${iv1}					 ///
			if touse_sc == 1 & `s' == 1 [aw = CLWT], partial(${ZX} ${iv1})	 ///
			cluster(Siblings_Cousins)
		matrix Figure5 = nullmat(Figure5), (`results' \ e(results))
		matrix drop `results'
	}
}
matrix list Figure5, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure5") modify
putexcel A1 = matrix(Figure5), names
putexcel A1 = "model"
putexcel B1 = "stat"

********************************************************************************
**# Abstract Tables
********************************************************************************

**# Table A1 (see below for estimates with extended dataset)
// 	Descriptive Statistics by Generational Life Stage
********************************************************************************

**# Table A2
// 	Intergenerational Food Security Correlations by Childhood Proportion
// 	of Years with Insecurity
********************************************************************************
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop TableA2`i'
	foreach y in FSRAW mglsecure lowsecure vlwsecure {
		tempname R
		foreach x in FSRAW mglsecure lowsecure vlwsecure childinsecure {
			tempname results
			b1x2 `y' if touse_sc == 1 & `s' == 1 [aw = CLWT],				 ///
				x1all(Mean_`x'_c) x2all(${X} ${FAM})						 ///
				x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :		 ///
				earn = EarnPoverty_c : equity = lnequity_c :				 ///
				family = ${FAM}) x1only(Mean_`x'_c)							 ///
				cluster(Siblings_Cousins)
			fb1x2
			matrix `results' = (`=r(results)[1,1]' \ `=r(results)[2,1]'),	 ///
				(`=r(results)[3,1]' \ `=r(results)[4,1]')
			matrix coleq `results' = `y'
			matrix colnames `results' = uncond cond
			matrix rownames `results' = `x' se
			matrix `R' = nullmat(`R') \ `results'
			matrix drop `results'
		}
		matrix TableA2`i' = nullmat(TableA2`i'), `R'
		matrix drop `R'
	}
}
matrix list TableA2A, format(%6.3f) title("Full sample")
matrix list TableA2B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("TableA2") modify
putexcel A1 = matrix(TableA2A), names
putexcel A14 = matrix(TableA2B), names

**# Table A3
// 	Intergenerational Correlations of Aggregated Measures of
// 	Food Security Status in Adulthood and Childhood
********************************************************************************
foreach m in Any Mean {
	foreach s in all children {
		if ("`m'" == "Any" & "`s'" == "all")		local i = "A"
		if ("`m'" == "Any" & "`s'" == "children")	local i = "B"
		if ("`m'" == "Mean" & "`s'" == "all")		local i = "C"
		if ("`m'" == "Mean" & "`s'" == "children")	local i = "D"
		if ("`s'" == "children") local a = "c"
		else local a
		capture matrix drop TableA3`i'
		foreach y in FSRAW mglsecure lowsecure vlwsecure {
			tempname R
			if ("`y'" == "FSRAW") local k = "Mean"
			else local k = "`m'"
			foreach x in FSRAW mglsecure lowsecure vlwsecure childinsecure {
				if ("`x'" == "FSRAW" & "`m'" == "Any") continue
				tempname results
				b1x2 `k'_`y'_a if touse_sc_a`a' == 1 & `s' == 1 [aw = m_WT], ///
					x1all(`m'_`x'_c) x2all(${m_X}) x2delta(Age = ${m_age} :	 ///
					X = ${m_controls} ${SD} ${YD} : earn = EarnPoverty_c :	 ///
					equity = lnequity_c : family = ${FAM}) x1only(`m'_`x'_c) ///
					cluster(Siblings_Cousins)
				fb1x2
				matrix `results' =											 ///
					(`=r(results)[1,1]' \ `=r(results)[2,1]'),				 ///
					(`=r(results)[3,1]' \ `=r(results)[4,1]')
				matrix coleq `results' = `y'
				matrix colnames `results' = uncond cond
				matrix rownames `results' = `x' se
				matrix `R' = nullmat(`R') \ `results'
				matrix drop `results'
			}
			matrix TableA3`i' = nullmat(TableA3`i'), `R'
			matrix drop `R'
		}
	}
}
matrix list TableA3A, format(%6.3f) title("Full sample: Any")
matrix list TableA3B, format(%6.3f) title("Children present: Any")
matrix list TableA3C, format(%6.3f) title("Full sample: Mean")
matrix list TableA3D, format(%6.3f) title("Children present: Mean")
putexcel set FoodSec_Results_v${version}, sheet("TableA3") modify
putexcel A1 = matrix(TableA3A), names
putexcel A12 = matrix(TableA3B), names
putexcel A23 = matrix(TableA3C), names
putexcel A36 = matrix(TableA3D), names

**# Table A4
// 	Intensity of Childhood Food Insecurity Conditional on Levels of
// 	Any Childhood Exposure
********************************************************************************
capture matrix drop TableA4
local y mglsecure lowsecure vlwsecure childinsecure insufficient_screen		 ///
	insufficient_q123
foreach x of local y {
	tempname R
	sum Any_`x'_c if touse_sc == 1 [aw = CLWT]
	matrix `R' = nullmat(`R'), (`r(mean)' \ `r(sd)')
	sum Mean_`x'_c if touse_sc == 1 & Any_`x'_c == 1 [aw = CLWT]
	matrix `R' = nullmat(`R'), (`r(mean)' \ `r(sd)')
	sum Mean_FSRAW_c if touse_sc == 1 & Any_`x'_c == 1 [aw = CLWT]
	matrix `R' = nullmat(`R'), (`r(mean)' \ `r(sd)')
	matrix TableA4 = nullmat(TableA4) \ `R'
	matrix drop `R'
}
matrix coleq TableA4 = "" "Conditional" "Conditional"
matrix colnames TableA4 = "Any" "Mean" "Score"
matrix rownames TableA4 = "Marginal" "sd" "Insecure"  "sd" "Very low"  "sd"	 ///
	"Child insecure"  "sd" "Screener" "sd" "Questions1-3" "sd"
matrix list TableA4, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("TableA4") modify
putexcel A1 = matrix(TableA4), names

**# Table A5
// 	Intergenerational Correlations of Food Insufficiency in Childhood
********************************************************************************
capture matrix drop TableA5
foreach y of varlist FSRAW mglsecure lowsecure vlwsecure {
	tempname R
	foreach x of varlist Any_insufficient_screen_c Any_insufficient_q123_c {
		tempname results
		b1x2 `y' if touse_sc == 1 & all == 1 [aw = CLWT],					 ///
			x1all(`x') x2all(${X} ${FAM})									 ///
			x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :			 ///
			earn = EarnPoverty_c : equity = lnequity_c :					 ///
			family = ${FAM}) x1only(`x') cluster(Siblings_Cousins)
		fb1x2
		matrix `results' = (`=r(results)[1,1]' \ `=r(results)[2,1]'),	 ///
			(`=r(results)[3,1]' \ `=r(results)[4,1]')
		matrix coleq `results' = `y'
		matrix colnames `results' = uncond cond
		matrix rownames `results' = `x' se
		matrix `R' = nullmat(`R') \ `results'
		matrix drop `results'
	}
	matrix TableA5 = nullmat(TableA5), `R'
	matrix drop `R'
}
matrix list TableA5, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("TableA5") modify
putexcel A1 = matrix(TableA5), names

**# Table A6
// 	Intergenerational Correlations of Early Adulthood Food Security Status
// 	and Partial or Complete Childhood Exposure to Marginal, Low, or
// 	Very Low Food Security
********************************************************************************
local x0
local x1 ${X} ${FAM}
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop TableA6`i'
	foreach y of varlist FSRAW mglsecure lowsecure vlwsecure {
		forvalues x = 0/1 {
			reg `y' i.Cat_mglsecure_c `x`x'' if touse_sc == 1 & `s' == 1	 ///
				[aw = CLWT], cluster(Siblings_Cousins)
			matrix TableA6`i' = nullmat(TableA6`i'),						 ///
				(_b[1.Cat_mglsecure_c] \ _se[1.Cat_mglsecure_c] \			 ///
				_b[2.Cat_mglsecure_c] \ _se[2.Cat_mglsecure_c] \ e(N))
		}
	}
	matrix coleq TableA6`i' = FSRAW FSRAW mgl mgl low low vlw vlw
	matrix colnames TableA6`i' = uncond cond uncond cond uncond cond uncond cond
	matrix rownames TableA6`i' = Partial se Complete se N
}
matrix list TableA6A, format(%6.3f) title("Full sample")
matrix list TableA6B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("TableA6") modify
putexcel A1 = matrix(TableA6A), names
putexcel A9 = matrix(TableA6B), names

**# Table A7
// 	Intergenerational Correlations of Continuous Measures of Food Security
// 	or Food Spending Relative to Needs
********************************************************************************
capture matrix drop TableA7A
capture matrix drop TableA7B
capture matrix drop TableA7C
capture matrix drop TableA7D
local x0
foreach s in all children {
	foreach y of varlist FSRAW latentFS latentFS_pr foodspend_pr {
		if ("`s'" == "all") local i = "A"
		if ("`s'" == "children") local i = "B"
		local x1 ${X} ${FAM}
		forvalues x = 0/1 {
			reg `y' Mean_`y'_c `x`x'' if touse_sc == 1 & `s' == 1			 ///
				[aw = CLWT], cluster(Siblings_Cousins)
			local b = _b[Mean_`y'_c]
			local se = _se[Mean_`y'_c]
			xtsum CLWT if e(sample) == 1
			local n = r(n)
			local N = r(N)
			margins, eyex(Mean_`y'_c) post atmeans
			matrix TableA7`i' = nullmat(TableA7`i'),						 ///
				(`b' \ `se' \ _b[Mean_`y'_c] \ `=r(table)[4,1]' \ `n' \ `N')
		}
		if ("`s'" == "all") local i = "C"
		if ("`s'" == "children") local i = "D"
		if ("`s'" == "children") local a = "c"
		else local a
		local x1 ${m_X}
		forvalues x = 0/1 {
			reg Mean_`y'_a Mean_`y'_c `x`x'' if touse_sc_a`a' == 1 &		 ///
				`s' == 1 [aw = m_WT], cluster(Siblings_Cousins)
			local b = _b[Mean_`y'_c]
			local se = _se[Mean_`y'_c]
			xtsum m_WT if e(sample) == 1
			local n = r(n)
			local N = r(N)
			margins, eyex(Mean_`y'_c) post atmeans
			matrix TableA7`i' = nullmat(TableA7`i'),						 ///
				(`b' \ `se' \ _b[Mean_`y'_c] \ `=r(table)[4,1]' \ `n' \ `N')
		}
	}
}
foreach i in A B C D {
	matrix rownames TableA7`i' = b se elast pvalue n N
	matrix colnames TableA7`i' = uncond cond uncond cond uncond cond uncond cond
	matrix coleq TableA7`i' = FSRAW FSRAW latent latent latent_pr latent_pr	 ///
		fspend_pr fspend_pr
}
matrix list TableA7A, format(%6.3f) title("Full sample: Panel")
matrix list TableA7B, format(%6.3f) title("Children present: Panel")
matrix list TableA7C, format(%6.3f) title("Full sample: Aggregate")
matrix list TableA7D, format(%6.3f) title("Children present: Aggregate")
putexcel set FoodSec_Results_v${version}, sheet("TableA7") modify
putexcel A1 = matrix(TableA7A), names
putexcel A10 = matrix(TableA7B), names
putexcel A19 = matrix(TableA7C), names
putexcel A28 = matrix(TableA7D), names

**# Table B1
// 	Decompositions of Unconditional Intergenerational Correlations of
// 	Food Security Status
********************************************************************************
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop TableB1`i'
	tempname R
	foreach y in FSRAW mglsecure lowsecure vlwsecure {
		if ("`y'" == "FSRAW") local a = "Mean"
		else local a = "Any"
		b1x2 `y' if touse_sc == 1 & `s' == 1 [aw = CLWT],				 ///
			x1all(`a'_`y'_c) x2all(${X} ${FAM})							 ///
			x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :		 ///
			earn = EarnPoverty_c : equity = lnequity_c :				 ///
			family = ${FAM}) x1only(`a'_`y'_c) cluster(Siblings_Cousins)
		fb1x2
		matrix `R' = nullmat(`R'), r(results)
	}
	matrix TableB1`i' = `R'[1..4,1...] \ `R'[15..16,1...] \ `R'[5..14,1...]
	matrix drop `R'
}
matrix list TableB1A, format(%6.3f) title("Full sample")
matrix list TableB1B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("TableB1") modify
putexcel A1 = matrix(TableB1A), names
putexcel A20 = matrix(TableB1B), names

**# Table B2
// 	Intergenerational Correlations of Early Adulthood Marginal, Low, or
// 	Very Low Food Security and Childhood Food Security Status by
// 	Age of Childhood Exposure
********************************************************************************
foreach m in Any Mean {
	if ("`m'" == "Any") local i = "A"
	if ("`m'" == "Mean") local i = "B"
	capture matrix drop TableB2`i'
	foreach a in "c0_5" "c6_11" "c12_17" "c" {
		tempname R
		foreach x in mglsecure lowsecure vlwsecure childinsecure {
			tempname results
			b1x2 mglsecure if touse_sc == 1 & all == 1 [aw = CLWT],			 ///
				x1all(`m'_`x'_`a') x2all(${X} ${FAM})						 ///
				x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :		 ///
				earn = EarnPoverty_c : equity = lnequity_c :				 ///
				family = ${FAM}) x1only(`m'_`x'_`a') cluster(Siblings_Cousins)
			fb1x2
			matrix `results' = (`=r(results)[1,1]' \ `=r(results)[2,1]'),	 ///
				(`=r(results)[3,1]' \ `=r(results)[4,1]')
			matrix coleq `results' = y_mglsecure_`a'
			matrix colnames `results' = uncond cond
			matrix rownames `results' = `x' se
			matrix `R' = nullmat(`R') \ `results'
			matrix drop `results'
		}
		xtsum CLWT if e(sample) == 1
		local n = r(n)
		local N = r(N)
		matrix TableB2`i' = nullmat(TableB2`i'), (`R' \ J(1,2,`n') \ J(1,2,`N'))
		matrix drop `R'
	}
}
matrix list TableB2A, format(%6.3f) title("Any exposure")
matrix list TableB2B, format(%6.3f) title("Mean exposure")
putexcel set FoodSec_Results_v${version}, sheet("TableB2") modify
putexcel A1 = matrix(TableB2A), names
putexcel A14 = matrix(TableB2B), names

**# Table C1
// 	LIML Estimates of Childhood Marginal, Low, or Very Low Food Security
// 	Effects on Early Adulthood Food Security Status
********************************************************************************
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop TableC1`i'
	foreach y in FSRAW mglsecure lowsecure vlwsecure {
		tempname results
		ivreg2 `y' (Any_mglsecure_c = ${IV1}) ${ZX} ${iv1}					 ///
			if touse_sc == 1 & `s' == 1 [aw = CLWT], partial(${ZX} ${iv1})	 ///
			cluster(Siblings_Cousins) liml
		matrix `results' = _b[Any_mglsecure_c] \ _se[Any_mglsecure_c]
		weakiv ivreg2 `y' (Any_mglsecure_c = ${IV1}) ${ZX} ${iv1}			 ///
			if touse_sc == 1 & `s' == 1 [aw = CLWT], partial(${ZX} ${iv1})	 ///
			cluster(Siblings_Cousins) liml
		matrix TableC1`i' = nullmat(TableC1`i'), (`results' \ e(clr_p))
		matrix drop `results'
	}
	matrix rownames TableC1`i' = b se p
	matrix colnames TableC1`i' = FSRAW mglsecure lowsecure vlwsecure
}
matrix list TableC1A, format(%6.3f) title("Full sample")
matrix list TableC1B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("TableC1") modify
putexcel A1 = matrix(TableC1A), names
putexcel A6 = matrix(TableC1B), names

**# Table D1
// 	Life-Cycle-Adjusted Intergenerational Conditional Correlations and
// 	2SLS Estimates for Childhood Food Security Status and Early Adulthood
// 	Outcomes Potentially Related to Adult Food Insecurity
********************************************************************************
capture matrix drop TableD1A
capture matrix drop TableD1B
foreach y in FSRAW mglsecure lowsecure vlwsecure {
	tempname R
	foreach x in mglsecure lowsecure vlwsecure childinsecure {
		tempname cm
		reg `y' Any_`x'_c C1_`x' C2_`x' C3_`x' C4_`x' ${leesolon}			 ///
			if touse_sc == 1 & all == 1 [aw = CLWT], cluster(Siblings_Cousins)
		tabstat CA1 CA2 CA3 CA4 if e(sample) == 1 [aw = CLWT], save
		tabstatmat `cm'
		nlcom (F_LS: _b[Any_`x'_c] +										 ///
			_b[C1_`x'] * `cm'[1,1] +										 ///
			_b[C2_`x'] * `cm'[1,2] +										 ///
			_b[C3_`x'] * `cm'[1,3] +										 ///
			_b[C4_`x'] * `cm'[1,4])
		matrix `R' = nullmat(`R') \ `=r(b)[1,1]' \ sqrt(`=r(V)[1,1]')
		matrix drop `cm'
	}
	matrix TableD1A = nullmat(TableD1A), `R'
	matrix drop `R'
	tempname cm
	ivreg2 `y' (Any_mglsecure_c C1_mglsecure C2_mglsecure C3_mglsecure 		 ///
		C4_mglsecure = ${CAIV}) ${leesolonZ} ${caiv}						 ///
		if touse_sc == 1 & all == 1 [aw = CLWT],							 ///
		partial(${leesolonZ} ${caiv}) cluster(Siblings_Cousins)	
	tabstat CA1 CA2 CA3 CA4 if e(sample) == 1 [aw = CLWT], save
	tabstatmat `cm'
	nlcom (F_LS: _b[Any_mglsecure_c] +										 ///
		_b[C1_mglsecure] * `cm'[1,1] +										 ///
		_b[C2_mglsecure] * `cm'[1,2] +										 ///
		_b[C3_mglsecure] * `cm'[1,3] +										 ///
		_b[C4_mglsecure] * `cm'[1,4])
	matrix TableD1B = nullmat(TableD1B), (`=r(b)[1,1]' \ sqrt(`=r(V)[1,1]'))
	matrix drop `cm'
}
matrix colnames TableD1A = FSRAW mglsecure lowsecure vlwsecure
matrix rownames TableD1A = mglsecure se lowsecure se vlwsecure se			 ///
	childinsecure se
matrix colnames TableD1B = FSRAW mglsecure lowsecure vlwsecure
matrix rownames TableD1B = mglsecure se
matrix list TableD1A, format(%6.3f) title("Conditional correlations")
matrix list TableD1B, format(%6.3f) title("2SLS: Any marginal or insecure")
putexcel set FoodSec_Results_v${version}, sheet("TableD1") modify
putexcel A1 = matrix(TableD1A), names
putexcel A11 = matrix(TableD1B), names

**# Table E1
// 	Intergenerational Correlations of Childhood Food Security Module
// 	Questionnaire Responses and Early Adulthood Food Security Status
********************************************************************************
foreach s in all children {
	if ("`s'" == "all") local i = "A"
	if ("`s'" == "children") local i = "B"
	capture matrix drop TableE1`i'
	foreach y of varlist FSRAW mglsecure lowsecure vlwsecure {
		tempname R
		foreach x of global FOOD {
			tempname results
			b1x2 `y' if touse_sc == 1 & `s' == 1 [aw = CLWT],				 ///
				x1all(Any_`x'_c) x2all(${X} ${FAM})							 ///
				x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :		 ///
				earn = EarnPoverty_c : equity = lnequity_c :				 ///
				family = ${FAM}) x1only(Any_`x'_c) cluster(Siblings_Cousins)
			fb1x2
			matrix `results' = (`=r(results)[1,1]' \ `=r(results)[2,1]'),	 ///
				(`=r(results)[3,1]' \ `=r(results)[4,1]')
			matrix coleq `results' = `y'
			matrix colnames `results' = uncond cond
			matrix rownames `results' = `x' se
			matrix `R' = nullmat(`R') \ `results'
			matrix drop `results'
		}
		matrix TableE1`i' = nullmat(TableE1`i'), `R'
		matrix drop `R'
	}
}
matrix list TableE1A, format(%6.3f) title("Full sample")
matrix list TableE1B, format(%6.3f) title("Children present")
putexcel set FoodSec_Results_v${version}, sheet("TableE1") modify
putexcel A1 = matrix(TableE1A), names
putexcel A40 = matrix(TableE1B), names

********************************************************************************
**# Abstract Figures
********************************************************************************

**# Figure A1 (tabulated from CPS December Supplement, survey years 1998-2004)
//	Frequency of Food Hardship relative to the Food Security Raw Score
********************************************************************************

**# Figure C1
// 	Confidence Regions for 2SLS Intergenerational Estimates
// 	of Marginal, Low, or Very Low Food Security Effects
********************************************************************************
capture matrix drop FigureC1
foreach y of varlist FSRAW mglsecure lowsecure vlwsecure {
	tempname results
	weakiv ivreg2 `y' (Any_mglsecure_c = ${IV1}) ${ZX} ${iv1}				 ///
		if touse_sc == 1 & all == 1 [aw = CLWT], partial(${ZX} ${iv1})		 ///
		cluster(Siblings_Cousins)
	matrix `results' = e(citable)
	matrix `results' = `results'[1...,"null1"],								 ///
		`results'[1...,"wald_p".."clr_p"]
	matrix coleq `results' = "`y'"
	matrix FigureC1 = nullmat(FigureC1), `results'
	matrix drop `results'
}
matrix list FigureC1, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("FigureC1") modify
putexcel A1 = matrix(FigureC1), colnames

**# Figure C2
// 	Sensitivity Analysis for 2SLS Estimates of Childhood Marginal, Low, or
// 	Very Low Food Security Effects on Early Adulthood Food Security Status
********************************************************************************
capture matrix drop FigureC2
local SD0
local SD1 ${SD}
foreach y of varlist FSRAW mglsecure lowsecure vlwsecure {
	foreach s of numlist 0/1 {
		foreach z of numlist 1/4 {
			tempname results
			iv_est `y' (Any_mglsecure_c = ${IV`z'}) ${ZX} ${iv`z'} `SD`s''	 ///
				if touse_sc == 1 & all == 1 [aw = CLWT],					 ///
				partial(${ZX} ${iv`z'} `SD`s'') cluster(Siblings_Cousins)
			matrix `results' = e(results)
			matrix coleq `results' = ""
			matrix colnames `results' = "`y'_S`s'_Z`z'"
			matrix FigureC2 = nullmat(FigureC2), `results'
			matrix drop `results'
		}
	}
}
matrix list FigureC2, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("FigureC2") modify
putexcel A1 = matrix(FigureC2), names
putexcel A1 = "model"
putexcel B1 = "stat"

********************************************************************************
**# Supplementary Estimates
********************************************************************************

**# Gelbach (2016) decomposition replication
********************************************************************************
quietly {
	tempname R B P1 P2 P3 P4 P5 D

	/* Table B1, panel A, column (1) */
	b1x2 FSRAW if touse_sc == 1 & all == 1 [aw = CLWT],						 ///
		x1all(Mean_FSRAW_c) x2all(${X} ${FAM})								 ///
		x2delta(Age = ${age} : X = ${controls} ${SD} ${YD} :				 ///
		earn = EarnPoverty_c : equity = lnequity_c : family = ${FAM})		 ///
		x1only(Mean_FSRAW_c) cluster(Siblings_Cousins)
	matrix `R' = e(Delta)
	fb1x2

	/* For model Y = F\beta + H\rho + \epsilon (see Appendix Section B1): */

	/* For each group, collect estimates for H_k on F: (F'F)^-1 F'H_k */
	global G1 ${age}
	global G2 ${controls} ${SD} ${YD}
	global G3 EarnPoverty_c
	global G4 lnequity_c
	global G5 ${FAM}
	foreach g of numlist 1/5 {
		tempname G`g'
		foreach x of global G`g' {
			reg `x' Mean_FSRAW_c if touse_sc == 1 & all == 1 [aw = CLWT]
			matrix `G`g'' = nullmat(`G`g'') \ _b[Mean_FSRAW_c]
		}
	}

	/* Estimate base model for comparison, excluding H */
	reg FSRAW Mean_FSRAW_c if touse_sc == 1 & all == 1 [aw = CLWT]
	local b0 = _b[Mean_FSRAW_c]

	/* For the fully specified model, collect estimates for \rho\hat */
	reg FSRAW Mean_FSRAW_c ${G1} ${G2} ${G3} ${G4} ${G5}					 ///
		if touse_sc == 1 & all == 1 [aw = CLWT]
	local b = _b[Mean_FSRAW_c]
	matrix `B' = e(b)

	/* Sum product of coefficients by group */
	matrix `P1' =															 ///
		`B'[1,"`: word 1 of ${G1}'".."`: word `: word count ${G1}' of ${G1}'"]
	matrix `P2' =															 ///
		`B'[1,"`: word 1 of ${G2}'".."`: word `: word count ${G2}' of ${G2}'"]
	matrix `P3' = _b[EarnPoverty_c]
	matrix `P4' = _b[lnequity_c]
	matrix `P5' =															 ///
		`B'[1,"`: word 1 of ${G5}'".."`: word `: word count ${G5}' of ${G5}'"]
	matrix `D' = `P1' * `G1' \ `P2' * `G2' \ `P3' * `G3'  \ `P4' * `G4' \	 ///
		`P5' * `G5'
	matrix `D' = `D' \ J(1,5,1) * `D'

	/* Compare Gelbach decomposition estimates from b1x2 and regressions */
	noisily display "Unconditional:" %8.3f `b0'
	noisily display "Conditional:" %10.3f `b'
	noisily display "Difference:" %11.3f `=`b0'-`b''
	matrix Gelbach = `R'', `D', `R'' - `D'
	matrix colnames Gelbach = b1x2 manual difference
	matrix drop `R' `G1' `G2' `G3' `G4' `G5' `B' `P1' `P2' `P3' `P4' `P5' `D'
	noisily matrix list Gelbach, format(%12.10f)
	putexcel set FoodSec_Results_v${version}, sheet("Gelbach") modify
	putexcel A1 = matrix(Gelbach), names
	putexcel A1 = "Replication"
}

/* Lee-Solon (2009) adjustment sensitivity to recentering in early adulthood */
********************************************************************************
quietly {
	capture matrix drop LSxA
	capture matrix drop LSxB
	foreach age of numlist 19/35 {
		forvalues i = 1/4 {
			replace CA`i' = (AGE - `age')^`i'
			replace C`i'_mglsecure = CA`i' * Any_mglsecure_c
			foreach x of global IV1 {
				replace C`i'_`x' = CA`i' * `x'
			}
			foreach x of global iv1 {
				replace C`i'_`x' = CA`i' * `x'
			}
		}
		tempname cm
		reg FSRAW Any_mglsecure_c C1_mglsecure C2_mglsecure C3_mglsecure	 ///
			C4_mglsecure ${leesolon} if touse_sc == 1 & all == 1			 ///
			[aw = CLWT], cluster(Siblings_Cousins)
		tabstat CA1 CA2 CA3 CA4 if e(sample) == 1 [aw = CLWT], save
		tabstatmat `cm'
		nlcom (F_LS: _b[Any_mglsecure_c] +									 ///
			_b[C1_mglsecure] * `cm'[1,1] +									 ///
			_b[C2_mglsecure] * `cm'[1,2] +									 ///
			_b[C3_mglsecure] * `cm'[1,3] +									 ///
			_b[C4_mglsecure] * `cm'[1,4])
		matrix LSxA = nullmat(LSxA) \ (`=r(b)[1,1]', sqrt(`=r(V)[1,1]'))
		matrix drop `cm'
		
		tempname cm
		ivreg2 FSRAW (Any_mglsecure_c C1_mglsecure C2_mglsecure  			 ///
			C3_mglsecure C4_mglsecure = ${CAIV}) ${leesolonZ} ${caiv}		 ///
			if touse_sc == 1 & all == 1 [aw = CLWT],						 ///
			partial(${leesolonZ} ${caiv}) cluster(Siblings_Cousins)
		tabstat CA1 CA2 CA3 CA4 if e(sample) == 1 [aw = CLWT], save
		tabstatmat `cm'
		nlcom (F_LS: _b[Any_mglsecure_c] +									 ///
			_b[C1_mglsecure] * `cm'[1,1] +									 ///
			_b[C2_mglsecure] * `cm'[1,2] +									 ///
			_b[C3_mglsecure] * `cm'[1,3] +									 ///
			_b[C4_mglsecure] * `cm'[1,4])
		matrix LSxB = nullmat(LSxB) \ (`=r(b)[1,1]', sqrt(`=r(V)[1,1]'))
		matrix drop `cm'
	}
	forvalues i = 1/4 {
		replace CA`i' = (AGE - 25)^`i'
		replace C`i'_mglsecure = CA`i' * Any_mglsecure_c
		foreach x of global IV1 {
			replace C`i'_`x' = CA`i' * `x'
		}
		foreach x of global iv1 {
			replace C`i'_`x' = CA`i' * `x'
		}
	}
	numlist "19/35"
	matrix rownames LSxA = `r(numlist)'
	matrix rownames LSxB = `r(numlist)'
	matrix colnames LSxA = OLS_b OLS_se
	matrix colnames LSxB = IV_b IV_se
	noisily matrix list LSxA, title("Lee-Solon OLS sensitivity")
	noisily matrix list LSxB, title("Lee-Solon IV sensitivity")
	putexcel set FoodSec_Results_v${version}, sheet("Lee-Solon") modify
	putexcel A1 = matrix(LSxA), names
	putexcel D1 = matrix(LSxB), colnames
	putexcel A1 = "Age"
}


********************************************************************************
**# Estimates using extended dataset
********************************************************************************
use "food_security_v${version}_extended_sample.dta", clear

/* State fixed effects */
tab STATEALPHA, g(SD)
unab sd: SD1-SD`r(r)'
global SD `sd'

/* Age profiles */
global age AGE AGE2 PAGE PAGE2 AGE_h AGE_h2

/* Control variables */
global controls spm_historical_under100_S unemployment_rate_S FEMALE		 ///
	BlackNonHispanic OtherNonHispanic Hispanic RaceEthnicity_imputed		 ///
	NC1 NC2 NC3 NC4

**# Figure 1
// 	Percent of Children by Food Security and Poverty Status
********************************************************************************
tempname insecure childinsecure
tabstat lowsecure if AGE >= 0 & AGE < 18 & FamUnit > 0 &					 ///
	inlist(yearcats, 1, 3) & inlist(SURVEY, 1, 2) [aw = CLWT],				 ///
	by(SURVEYYEAR) s(mean semean) nototal save
tabstatmat `insecure'
tabstat childinsecure if AGE >= 0 & AGE < 18 & FamUnit > 0 &				 ///
	inlist(yearcats, 1, 3) & inlist(SURVEY, 1, 2) [aw = CLWT],				 ///
	by(SURVEYYEAR) s(mean semean) nototal save
tabstatmat `childinsecure'
matrix Figure1 = `insecure', `childinsecure'
matrix drop `insecure' `childinsecure'
matrix colnames Figure1 = "childHH" "childHH_se" "childinsecure"			 ///
	"childinsecure_se"
matrix rownames Figure1 = `: roweq Figure1'
matrix roweq Figure1 = ""
matrix list Figure1, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure1") modify
putexcel A1 = matrix(Figure1), names
putexcel A1 = "Year"

**# Figure 2
// 	Within-Family Persistence in Food Security Status,
// 	Earnings, and Food Spending Relative to the 1997 Survey Year
********************************************************************************
capture matrix drop Figure2A
foreach x in mglsecure lowsecure vlwsecure {
	tempname results
	reg `x' i.FS_`x'_1997##i.SURVEYYEAR ${age} ${controls} ${SD} [aw = CLWT] ///
		if touse_cohort == 1 & FamUnit > 0 & SURVEYYEAR > 1997,				 ///
		cluster(Siblings_Cousins)
	margins, dydx(FS_`x'_1997) at(SURVEYYEAR=(1999 2001 2003 2014 2015 2017	 ///
		2019 2020 2021)) post
	matrix `results' = r(table)
	
	numlist "1999 2001 2003 2014 2015 2017 2019 2020 2021"
	matrix rownames `results' = b_`x' s_`x' t_`x' p_`x' ll_`x' ul_`x' ""
	matrix colnames `results' = `r(numlist)' `r(numlist)'
	matrix coleq `results' = ""
	
	matrix Figure2A = nullmat(Figure2A), `results'[1,10..18]',				 ///
		`results'[5,10..18]', `results'[6,10..18]'
	matrix drop `results'
}
matrix list Figure2A, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure2A") modify
putexcel A1 = matrix(Figure2A), names
putexcel A1 = "SurveyYear"

capture matrix drop Figure2B
foreach x in foodspend lnEARNPOV lowFN {
	tempname results
	local a = "c."
	local j = "1..12"
	if ("`x'" == "lowFN") {
		local a = "i."
		local j = "13..24"
	}
	reg `x' `a'`x'3_1997##i.SURVEYYEAR ${age} ${controls} ${SD} [aw = CLWT]	 ///
		if touse_cohort == 1 & FamUnit > 0 & SURVEYYEAR > 1997,				 ///
		cluster(Siblings_Cousins)
	margins, dydx(`x'3_1997) at(SURVEYYEAR = (1999(2)2021)) post
	matrix `results' = r(table)
	
	numlist "1999(2)2021"
	matrix rownames `results' = b_`x' s_`x' t_`x' p_`x' ll_`x' ul_`x' ""
	matrix colnames `results' = `r(numlist)'
	matrix coleq `results' = ""
	
	matrix Figure2B = nullmat(Figure2B), `results'[1,`j']',				 ///
		`results'[5,`j']', `results'[6,`j']'
	matrix drop `results'
}
matrix list Figure2B, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("Figure2B") modify
putexcel A1 = matrix(Figure2B), names
putexcel A1 = "SurveyYear"

**# Table A1
// 	Descriptive Statistics by Generational Life Stage
********************************************************************************
tempvar C1 C2 C3 C4
mark `C1' if TOUSE_SC == 1 & AGE >= 0 & AGE < 18 & FamUnit > 0 &			 ///
	FScats != . & yearcats == 1
mark `C2' if touse_sc == 1
mark `C3' if touse_sc == 1 & children == 1
mark `C4' if nouse == 1
capture matrix drop TableA1
forvalues c = 1/4 {
	tempname x1 x2 x3 results
	tabstat FINC000 EARN000 FOOD000 [aw = CLWT] if `C`c'' == 1, s(med iqr n) ///
		save
	tabstatmat `x1'
	tabstat SNAP000 [aw = CLWT] if FSSNAP == 1 & `C`c'' == 1, s(med iqr n) save
	tabstatmat `x2'
	tabstat FSSNAP foodspend lowFN FS_status_1 FS_status_2 FS_status_3		 ///
		FS_status_4 foodinsecure_children POV LTEHS AGE_h AnyMarried		 ///
		FamUnit_NCU18 BlackNonHispanic WhiteNonHispanic OtherNonHispanic	 ///
		Hispanic [aw = CLWT] if `C`c'' == 1, s(mean sd n) save
	tabstatmat `x3'
	forvalues i = 1/3 {
		forvalues j = 1/`=colsof(`x`i'')' {
			matrix `results' = nullmat(`results') \ `x`i''[1..2,`j']
		}
	}
	xtsum CLWT if `C`c'' == 1
	matrix TableA1 = nullmat(TableA1), (`results' \ `r(n)' \ `r(N)')
	matrix drop `x1' `x2' `x3' `results'
}
drop `C1' `C2' `C3' `C4'
matrix coleq TableA1 = "Childhood, ages 0-17" "Early adulthood, ages 18-34"	 ///
	"Early adulthood with children present" "X-sample"
matrix colnames TableA1 = "1997-2003" "2014-2019" "2014-2019" "2014-2019"
matrix rownames TableA1 = "Family income" "-" "Family earnings" "-"			 ///
	"Food expenditure" "-" "Food stamps SNAP value" "-"						 ///
	"Receives food stamps SNAP?" "-" "Food spending per TFP" "-"			 ///
	"Food spending below TFP?" "-" "Food secure?" "-"						 ///
	"Marginal food secure?" "-" "Low food secure?" "-"						 ///
	"Very low food secure?" "-" "Food-insecure children?" "-"				 ///
	"Poverty status?" "-" "Education high school or less?" "-"				 ///
	"Age of head of family" "-" "Married couple in family?" "-"				 ///
	"Number of children in family" "-" "Black, non-Hispanic?" "-"			 ///
	"White, non-Hispanic?" "-" "Other, non-Hispanic?" "-" "Hispanic" "-"	 ///
	"Individuals" "Observations"
matrix list TableA1, format(%6.3f)
putexcel set FoodSec_Results_v${version}, sheet("TableA1") modify
putexcel A1 = matrix(TableA1), names

**# Figure A1 (tabulated from CPS December Supplement, survey years 1998-2004)
//	Frequency of Food Hardship relative to the Food Security Raw Score
********************************************************************************


********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
capture log close FSEC_est
exit
