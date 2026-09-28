clear all
set more off
* Move to the project folder (works whether Stata starts there or in Code/)
capture confirm file "Code/globals.do"
if _rc cd ..
do "Code/globals.do"
use "$data", clear

* 1. Panel structure and data
xtset ccode year

gen regtrans_alt = 0
replace regtrans_alt = 1  if trans_democ == 1
replace regtrans_alt = -1 if trans_autoc == 1

capture drop d_polity2 lag_polity2
gen d_polity2 = D.polity2
gen lag_polity2 = L.polity2

gen valid_obs = (missing(d_polity2) == 0 & ///
                 missing(petro_nx_index_growth) == 0 & ///
                 missing(lag_polity2) == 0)

bysort ccode: egen total_valid = sum(valid_obs)

summarize total_valid, meanonly
local max_T = r(max)

gen in_balanced_sample = (total_valid == `max_T')

local balance_vars polity2 growth petro_nx_index_growth trans_democ trans_autoc

estpost ttest `balance_vars', by(in_balanced_sample)
esttab using "$out_ext/balance_table.rtf", replace ///
    cells("mu_2(fmt(3)) mu_1(fmt(3)) b(star fmt(3)) p(fmt(3))") ///
    label collabels("Balanced" "Unbalanced" "Difference" "p-value") ///
    title("Table A1: Balance Test (Balanced vs. Full Sample)") ///
    nonumbers mtitle("Sample Comparison")

* Print summary
di "--- SUMMARY ---"
count if in_balanced_sample == 1
di "Observations in Balanced Panel: " r(N)
count if in_balanced_sample == 0
di "Dropped Observations: " r(N)

keep if in_balanced_sample == 1
drop if missing(d_polity2, petro_nx_index_growth, lag_polity2)

save "$balanced", replace
di "Strictly balanced dataset successfully saved."

* 2. Unit Root Test (stationarity of the shock is key)
xtunitroot fisher petro_nx_index_growth, dfuller lags(1)

* 3. Serial Correlation Test
xtserial polity2 petro_nx_index_growth

* 4. Check variation of the Threshold Variable
histogram lag_polity2, discrete width(1) color(navy%80) ///
    xtitle("Lagged Polity2 Score (State Variable)") ///
    title("Distribution of the Threshold Variable") ///
    note("Ensure mass is not exclusively at boundaries (-10 and 10).") ///
    addplot(pci 0 0 .15 0, lcolor(red) lpattern(dash) lwidth(medthick) legend(off))

* Hansen threshold
* 1. Generate clean Time Dummies for the specific balanced timeframe
capture drop td_*
quietly tab year, gen(td_)

* 2. Expand dummy list
unab all_dummies : td_*
tokenize "`all_dummies'"
macro shift 
local final_dummies "`*'"

* 3. Run the Endogenous Threshold Model
di "Executing xthreg. This will take a few minutes..."

xthreg d_polity2 `final_dummies', ///
    rx(petro_nx_index_growth) ///
    qx(lag_polity2) ///
    thnum(1) trim(0.15) bs(100)
	
* 4. Get a list of all variables to export for visualizaiton 
ds
describe

* Robustness checks

* 1. Define the alternative dependent variables from Table 15
local alt_vars polright fh_ind reg ev_ps1

foreach var of local alt_vars {
    preserve
    
    capture drop d_`var'
    gen d_`var' = D.`var'
    
    drop if missing(d_`var', petro_nx_index_growth, lag_polity2)
    
    bysort ccode: gen T_count = _N
    summarize T_count, meanonly
    keep if T_count == r(max)
    local current_T = r(max)
    
    capture drop td_*
    quietly tab year, gen(td_)
    unab all_dummies : td_*
    tokenize "`all_dummies'"
    macro shift 
    local final_dummies "`*'"
    
    di ""
    di "====================================================================="
    di "ROBUSTNESS CHECK: Dependent Variable = Change in `var'"
    di "Balanced Panel Time Periods (T) = `current_T'"
    di "NOTE: If 'var' is polright, fh_ind, or reg, signs should be FLIPPED!"
    di "====================================================================="
    
    xthreg d_`var' `final_dummies', ///
        rx(petro_nx_index_growth) ///
        qx(lag_polity2) ///
        thnum(1) trim(0.15) bs(100)
        
    restore
}


* 2. Reload full unbalanced dataset 
use "$data", clear
xtset ccode year

* 3. Define asymmetric shocks (Booms vs. Busts)
gen shock_boom = max(0, petro_nx_index_growth)
gen shock_bust = abs(min(0, petro_nx_index_growth))

* 4. Define the regimes (kagged state variables)
gen autoc_0 = (L.polity2 <= 0)
gen democ_0 = (L.polity2 > 0)

gen autoc_6 = (L.polity2 <= -6)
gen democ_6 = (L.polity2 > -6)

gen boom_autoc6 = shock_boom * autoc_6
gen boom_democ6 = shock_boom * democ_6

gen bust_autoc6 = shock_bust * autoc_6
gen bust_democ6 = shock_bust * democ_6

gen boom_autoc0 = shock_boom * autoc_0
gen boom_democ0 = shock_boom * democ_0
gen bust_autoc0 = shock_bust * autoc_0
gen bust_democ0 = shock_bust * democ_0

* 5. The Local Projections Loop

local max_h = 5

forvalues h = 0/`max_h' {
    
    * Generate the forward-looking dependent variable: (D_{t+h} - D_{t-1})
    capture drop fwd_change_`h'
    if `h' == 0 {
        gen fwd_change_`h' = polity2 - L.polity2
    }
    else {
        gen fwd_change_`h' = F`h'.polity2 - L.polity2
    }
    
    * Estimate using High-Dimensional FE (Country and Year), clustering at country level
    quietly reghdfe fwd_change_`h' boom_autoc6 boom_democ6 bust_autoc6 bust_democ6, ///
        absorb(ccode year) vce(cluster ccode)
        
    * Store the estimates for the IRF plots
    estimates store lp_h`h'
}
capture which coefplot
if _rc ssc install coefplot

capture drop d_polity2
gen d_polity2 = D.polity2

local max_h = 5

forvalues h = 0/`max_h' {
    capture drop fwd_change_`h'
    if `h' == 0 {
        gen fwd_change_`h' = polity2 - L.polity2
    }
    else {
        gen fwd_change_`h' = F`h'.polity2 - L.polity2
    }
    quietly reghdfe fwd_change_`h' boom_autoc6 boom_democ6 bust_autoc6 bust_democ6 ///
        autoc_6 L.polity2 L.d_polity2 L.shock_boom L.shock_bust, ///
        absorb(ccode year) vce(cluster ccode)
        
    estimates store lp_final_h`h'
}

* 6. Plot IRFs

* Plot A: The Autocratic Upgrading Hypothesis (Booms)
coefplot lp_final_h0 lp_final_h1 lp_final_h2 lp_final_h3 lp_final_h4 lp_final_h5, ///
    keep(boom_autoc6) vertical recast(connected) ciopts(recast(rline) lpattern(dash)) ///
    yline(0, lcolor(red)) ///
    title("BULLETPROOF IRF: Oil Booms in Hard Autocracies (Polity <= -6)") ///
    xtitle("Years after shock (Horizon h)") ytitle("Cumulative Change in Polity2") ///
    rename(lp_final_h0="0" lp_final_h1="1" lp_final_h2="2" lp_final_h3="3" lp_final_h4="4" lp_final_h5="5")

* Plot B: The Dictator's Crackdown (Busts)
coefplot lp_final_h0 lp_final_h1 lp_final_h2 lp_final_h3 lp_final_h4 lp_final_h5, ///
    keep(bust_autoc6) vertical recast(connected) ciopts(recast(rline) lpattern(dash)) ///
    yline(0, lcolor(red)) ///
    title("BULLETPROOF IRF: Oil Busts in Hard Autocracies (Polity <= -6)") ///
    xtitle("Years after shock (Horizon h)") ytitle("Cumulative Change in Polity2") ///
    rename(lp_final_h0="0" lp_final_h1="1" lp_final_h2="2" lp_final_h3="3" lp_final_h4="4" lp_final_h5="5")
	
capture drop boom_autoc0 boom_democ0 bust_autoc0 bust_democ0
gen boom_autoc0 = shock_boom * autoc_0
gen boom_democ0 = shock_boom * democ_0
gen bust_autoc0 = shock_bust * autoc_0
gen bust_democ0 = shock_bust * democ_0

local max_h = 5

forvalues h = 0/`max_h' {
    capture drop fwd_change_`h'
    if `h' == 0 {
        gen fwd_change_`h' = polity2 - L.polity2
    }
    else {
        gen fwd_change_`h' = F`h'.polity2 - L.polity2
    }
    quietly reghdfe fwd_change_`h' boom_autoc0 boom_democ0 bust_autoc0 bust_democ0 ///
        autoc_0 L.polity2 L.d_polity2 L.shock_boom L.shock_bust, ///
        absorb(ccode year) vce(cluster ccode)
        
    estimates store lp_bullet0_h`h'
}


* Plot A: The Immunity Hypothesis (Booms)
coefplot lp_bullet0_h0 lp_bullet0_h1 lp_bullet0_h2 lp_bullet0_h3 lp_bullet0_h4 lp_bullet0_h5, ///
    keep(boom_democ0) vertical recast(connected) ciopts(recast(rline) lpattern(dash)) ///
    yline(0, lcolor(red)) ///
    title("BULLETPROOF IRF: Oil Booms in Democracies (Polity > 0)") ///
    xtitle("Years after shock (Horizon h)") ytitle("Cumulative Change in Polity2") ///
    rename(lp_bullet0_h0="0" lp_bullet0_h1="1" lp_bullet0_h2="2" lp_bullet0_h3="3" lp_bullet0_h4="4" lp_bullet0_h5="5") ///
    note("Note: Base state, L.polity2, and lagged dynamics controlled. 95% CI with clustered SEs.")

* Plot B: The Immunity Hypothesis (Busts)
coefplot lp_bullet0_h0 lp_bullet0_h1 lp_bullet0_h2 lp_bullet0_h3 lp_bullet0_h4 lp_bullet0_h5, ///
    keep(bust_democ0) vertical recast(connected) ciopts(recast(rline) lpattern(dash)) ///
    yline(0, lcolor(red)) ///
    title("BULLETPROOF IRF: Oil Busts in Democracies (Polity > 0)") ///
    xtitle("Years after shock (Horizon h)") ytitle("Cumulative Change in Polity2") ///
    rename(lp_bullet0_h0="0" lp_bullet0_h1="1" lp_bullet0_h2="2" lp_bullet0_h3="3" lp_bullet0_h4="4" lp_bullet0_h5="5") ///
    note("Note: Base state, L.polity2, and lagged dynamics controlled. 95% CI with clustered SEs.")
