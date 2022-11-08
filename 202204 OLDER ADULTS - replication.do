********************************************************************************
* CALPROTECT - CDCR DATA
* Created by:	ADA KWAN
* Created on:	26 JUL 2021
* Modified on:	08 APR 2022

/********************************************************************************

	REPLICATION FILE FOR: 
	
	The Impact Of COVID-19 On The
	Health Of Incarcerated Older
	Adults In California State Prisons
	
	by kwan, garcia-grossman, sears, bertozzi, williams
	
	https://www.healthaffairs.org/doi/pdf/10.1377/hlthaff.2022.00132

********************************************************************************/
clear all
set more off
** DIRECTORY, GLOBAL (paths, datasets)

*	directory
	glo directory "/Users/adak/Documents/AK/"
	cd "/Users/adak/Documents/AK/"

*	global paths 

	glo ado "dofiles/adofiles/"
	glo constructed "constructed/"
	glo data "data/"
		glo inter "data/intermediate/" // intermediate data files are stored here
		glo tmp "data/intermediate/tmp/" // intermediate data files are stored here
		glo logs "data/logs/"

			
	glo dofiles "dofiles/"
		glo do_analyze "dofiles/analysis"

			
/* 	loading ado files
	local adoFiles : dir `"${ado}"' files "*.ado"
	local adoFiles = subinstr(`" `adoFiles' "', `"""' , "" , .)
	foreach adoFile in `adoFiles' {
		qui do "${ado}`adoFile'"
		}
		*/
		do "dofiles/adofiles/chartable.do"
		do "dofiles/adofiles/chartable_fw.do"
		do "dofiles/adofiles/forest.do"

	net from http://www.stata.com/users/vwiggins
	*net install grc1leg
	*ssc install missings
	
*	other
	local fmt eps	
	
	local graph_opts ///
		graphregion(color(white)) ylab(,angle(0) nogrid) yscale(noline) xscale(noline) xsize(7)
	local title_opts color(black) span pos(12)

	*colorpalette viridis
	
	
	global pct `" 0 "0.0%" .002 "0.2%" .004 "0.4%" .006 "0.6%" .008 "0.8%" .01 "1.0%" "'
********************************************************************************

**************************************************
** EXHIBIT 1 // TABLE 1 - DESCRIPTIVES
**************************************************

global covars ///
age male total_inst total_days always_cohort recent_cohort sec_level_first sec_level_last cust_score_first cust_score_last disability_mh_any disability_cog_any ever_disability_hear ever_disability_mob ever_disability_speech ever_disability_vision last_covid last_cancer_high last_copd last_copd_high last_immuno last_dialysis last_adv_liver ever_flulike ever_hosp ever_quarantined n_days_quarantine ever_isolated n_days_isolated ever_tested n_times_tested ever_infected ever_hosp_covid ever_treated_bam ever_died ever_vax ///
tot_roomtype_180c tot_roomtype_270c tot_roomtype_c tot_roomtype_270d tot_roomtype_d tot_roomtype_r tot_roomtype_closed
				
			tabstat $covars, stat(n mean sd semean) column(statistics) save 
			tabstat $covars if older == 0, stat(n mean sd semean) column(statistics) save 
			tabstat $covars if older == 1, stat(n mean sd semean) column(statistics) save 
			foreach var of global covars {
				ta `var' older,m
				ttest `var', by(older)	
			}
			*

			
			foreach var of varlist race3 security {
				ta `var',m
				ta `var' if older == 0,m
				ta `var' if older == 1,m
			}
			
			
			
**************************************************
** EXHIBIT 2 // FIGURE 1 - WEEKLY CASES BY OLDER
**************************************************
use "${inter}/analytic_house_demog_infect_20_21.dta", clear

	keep residentid date not_older older ViralTestStatus agegrp age
	
	
	epiweek date, epiw(epiweek) epiy(epiyear)
	egen wy = concat(epiyear epiweek), p(w)
	epiweek2 wy, s(from) e(to)
			format from %tdMon_dd,_CCYY
			format to %tdMon_dd,_CCYY
			
		save "${inter}/analytic_house_demog_infect_20_21_tmp.dta", replace
	use "${inter}/analytic_house_demog_infect_20_21_tmp.dta", clear

		keep if epiyear == 2020 | epiyear == 2021
			
				duplicates drop residentid from, force
				
				collapse (sum) tot_old = older tot_notold = not_older, by(from) 
				
				replace tot_old = tot_old
				replace tot_notold = tot_notold
				
				bysort from: gen tot_pop = tot_old + tot_notold
				
				save "${tmp}/tot.dta", replace
	
	use "${inter}/analytic_house_demog_infect_20_21_tmp.dta", clear

		keep if epiyear == 2020 | epiyear == 2021
		
			gen positive = 1 if ViralTestStatus == "Positive"
				replace positive = 0 if ViralTestStatus == "Negative"
			
			bysort residentid from: egen pos = max(positive)
			
				bysort residentid (from): egen res_pos_count = sum(positive)
			
				gen last_pos = 1 if !missing(positive) & positive == 1
				bysort residentid (from): replace last_pos = last_pos[_n-1] if missing(last_pos)

			duplicates drop residentid from pos, force
			drop positive
			
			collapse (sum) positive = pos cum_pos = last_pos, by(from older) 
							merge m:1 from using "${tmp}/tot.dta"
							
							drop if tot_pop == 0
		
			sort from older
			
					
				gen cum_prop_notold = cum_pos / tot_notold if older == 0
				gen cum_prop_old = cum_pos / tot_old if older == 1
				gen cum_prop_unknown = cum_pos / (tot_pop - tot_old - tot_notold) if older == .
					
			gen label = string(from, "%tdMon_dd,_CCYY")
			gen from2 = from - 2
					colorpalette viridis, n(20) nograph
						return list
						local col1 = r(p1)
						local col2 = r(p7)
			
			sort from
		
		
		tw ///
			bar positive from if older == 0, barw(2) yaxis(1) color("ebblue*0.45") ///
		|| ///
			bar positive from2 if older == 1, barw(2) yaxis(1) color("red*0.45") ///
		|| ///
			line cum_prop_notold from,  yaxis(2) lcolor("ebblue*0.75") lpattern(dash) lwidth(thin) ///
		|| ///
			line cum_prop_old from,  yaxis(2) lcolor("red") lpattern(dash) lwidth(thin) ///
			ytitle("Weekly new COVID-19 cases by age group", axis(1) size(vsmall)) ///
			ytitle("Cumulative weekly % of residents infected by age group", axis(2) size(vsmall)) ///
			xtitle(" ") ///
			ylabel(, axis(1) labsize(vsmall)) ylabel(0 "0%" 0.2 "20%" 0.4 "40%" 0.6 "60%", axis(2) labsize(vsmall)) ///
			xlabel(#9, labsize(vsmall) angle(45)) legend(col(1) pos(6) symysize(vsmall) size(vsmall)  ///
						label(1 "Younger than 55 Years") ///
						label(2 "55 Years and Older") ///
						label(3 "Cumulative % of <55y infected (right axis)") ///
						label(4 "Cumulative % of 55y and older infected (right axis)") ///
						order(4 2 3 1))
					graph save Graph "${goutputs}/fig2_epicurve_cumulative.gph", replace
					graph export "${goutputs}/fig2_epicurve_cumulative.png", replace
					
					

**************************************************
** EXHIBIT 3 // FIGURE 2 - ODDS RATIOS
**************************************************

	*** table with agegrp2
	glo outputs  "outputs/202204 OLDER ADULTS/outputs/"
	glo outcomes ///
		not_recent_cohort always_cohort ever_quarantined ever_isolated ever_tested ever_positive_new ever_c19hosp_er ever_c19hosp_hosp ever_c19hosp_icu ever_died


	** frequencies and percents with older, with agegrp2
	
		glo main_outcomes /// 
			ever_tested ever_quarantined ever_isolated ever_positive_new ///
			pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died always_cohort not_recent_cohort
			
		foreach var of glo main_outcomes {
			ta `var' if older == 0
			ta `var' if older == 1
			ta `var' if agegrp2 == 2
			ta `var' if agegrp2 == 3
			ta `var' if agegrp2 == 4
		}
		*
		
		glo main_outcomes ///
			pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died day_part_vax day_fully_vax day_fully_vax_boost not_recent_cohort
			
		foreach var of glo main_outcomes {
			ta `var' 
			ta `var' if older == 0
			ta `var' if older == 1
			ta `var' if agegrp2 == 2
			ta `var' if agegrp2 == 3
			ta `var' if agegrp2 == 4
			ttest `var', by(older)
		}

		global main_outcomes ///
			ever_tested ever_quarantined ever_isolated ever_positive_new ///
			pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died ///
			day_part_vax day_fully_vax day_fully_vax_boost not_recent_cohort
		
	
	** unadjusted and adjusted ORs models with agegrp2
	
		 foreach var of glo main_outcomes {
			logistic `var' i.agegrp2 [fw=daily_active_cases]
				estimates store `var'
				local N = e(N)	
				local N = e(N)
				summ `var' if older == 0
				local `var'_bar= r(mean)
				summ `var' if agegrp2 == 1
				local `var'_bar_age3= r(mean)	
				summ `var' if male == 0
				local `var'_bar_m= r(mean)	
				summ `var' if race3 == 1
				local `var'_bar_r= r(mean)	
				outreg2 [`var'] using "${outputs}/regressions_or_202003_202110.xlsx", ///
					eform stats(coef ci pval) addstat(Younger than 55 years (Group Mean), ``var'_bar', Female (Group Mean), ``var'_bar_m', White (Group Mean), ``var'_bar_r', Observations, `N') label dec(3) noobs nor2  noaster paren(ci) 
			logistic `var' i.agegrp2 male i.race3 perc_anycell perc_anydorm [fw=daily_active_cases]
				estimates store `var'
				local N = e(N)	
				local N = e(N)
				summ `var' if older == 0
				local `var'_bar= r(mean)
				summ `var' if agegrp2 == 1
				local `var'_bar_age3= r(mean)	
				summ `var' if male == 0
				local `var'_bar_m= r(mean)	
				summ `var' if race3 == 1
				local `var'_bar_r= r(mean)	
				outreg2 [`var'] using "${outputs}/regressions_or_202003_202110.xlsx", ///
					eform stats(coef ci pval) addstat(Younger than 55 years (Group Mean), ``var'_bar', Female (Group Mean), ``var'_bar_m', White (Group Mean), ``var'_bar_r', Observations, `N') label dec(3) noobs nor2  noaster paren(ci) 
			}
			*
			
		

**************************************************
** EXHIBIT 4 // FIGURE 3 - % RELEASED BY AGE GROUP
**************************************************


	** SETUP **
	use "data/intermediate/20211010_NightlyHousing2_Demog.dta", clear

	keep residentid date agegrp agegrp2
	
	
	ren date dates
	gen date = date(dates, "YMD")
	format date %td
	
	merge 1:1 residentid date using "data/intermediate/20211010_Deaths.dta"
	
		epiweek date, epiw(epiweek) epiy(epiyear)
		egen wy = concat(epiyear epiweek), p(w)
		epiweek2 wy, s(from) e(to)
		
	save "data/intermediate/20211010_NightlyHousing2_Demog_tmp.dta", replace

		
	use "data/intermediate/20211010_NightlyHousing2_Demog_tmp.dta", clear
		keep residentid date agegrp from c19_death
		bys residentid: egen ever_c19death = max(c19_death)
		
		gen age1 = (agegrp == 13 | agegrp == 14 | agegrp == 15) // 13 "75 to 79" 14 "80 to 84" 15 "85 and older"
		gen age2 = (agegrp == 11 | agegrp == 12) // 11 "65 to 69" 12 "70 to 74" 
		gen age3 = (agegrp == 9 | agegrp == 10) // 9 "55 to 59" 10 "60 to 64" ///
		gen age4 = (agegrp == 7 | agegrp == 8) // 7 "45 to 49" 8 "50 to 54" ///
		gen age5 = (agegrp == 5 | agegrp == 6) // 5 "35 to 39" 6 "40 to 44"
		gen age6 = (agegrp == 1 | agegrp == 2 | agegrp == 3 | agegrp == 4) // 1 "19 and younger" 2 "20 to 24" 3 "25 to 29" 4 "30 to 34"
		drop agegrp
			format from %tdMon_dd,_CCYY

		
		save "data/intermediate/20211010_NightlyHousing2_Demog_tmp2.dta", replace
		
		sort residentid from
		duplicates tag, gen(tag)
		drop if tag > 0
		drop tag
		
		gen older = (age >=55)
		codebook residentid if dates == "2020-03-01",c
		 
		save "data/intermediate/20211010_NightlyHousing2_Demog_tmp3.dta", replace
	use "data/intermediate/20211010_NightlyHousing2_Demog_tmp3.dta", clear
		sort residentid date 
		tostring date, gen(dates)
		bysort residentid (date): gen comers = _n
				replace comers = 0 if dates == "21550" // Jan 1, 2019
				replace comers = 0 if comers != 1
			gen stayers = 1 if (comers == 0)
		
		
		gen negdate = - date
		sort negdate
			bysort residentid (negdate): gen leavers = _n
				replace leavers = 0 if dates == "22562" // Oct 9, 2021
				replace leavers = 0 if leavers != 1
				
				replace leavers = 0 if ever_c19death == 1
			
		sort date
		drop negdate
		
		gen agegrp6 = 1 if age6 == 1
			replace agegrp6 = 2 if age5 == 1
			replace agegrp6 = 3 if age4 == 1
			replace agegrp6 = 4 if age3 == 1
			replace agegrp6 = 5 if age2 == 1
			replace agegrp6 = 6 if age1 == 1
			lab def agegrp6 1 "Younger than 35" 2 "35 to 44" 3 "45 to 54" 4 "55 to 64" 5 "65 to 74" 6 "75 or older" 
			la val agegrp6 agegrp6
			la var agegrp6 "Age Group"

			
			
		preserve
			collapse (sum) stayers_w = stayers comers_w = comers leavers_w = leavers (count) residents = residentid, by(agegrp6 from)
			ren from date
			
			gen perc_leavers_w = leavers_w / (stayers_w + comers_w + leavers_w)
			gen perc_comers_w = comers_w / (stayers_w + comers_w + leavers_w)
			
			sort date
				tostring date, gen(date2)
				destring date2, replace
		gen month=month(date)
					tostring month, replace
					replace month = string(real(month),"%02.0f")
				gen year=year(date)
					tostring year, replace

			gen months = month + "-" + year
			
			
			drop if date2 == 22556
			
			save "data/tmp/cdcr_residents_201901_202110_weeks.dta", replace
		restore
		

				
		collapse (sum) stayers = stayers comers = comers leavers = leavers (count) residents = residentid, by(agegrp6 date)
			save "data/tmp/cdcr_residents_201901_202110.dta", replace
				

	
	** TOTALS **
	use  "data/tmp/cdcr_residents_201901_202110_weeks.dta", clear
	

		*total leavers only
		* line
		tw line leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			line leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			line leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			line leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			line leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			line leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Total Weekly Releases by Age Group")  ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6))
		graph save "Graph" "${goutputs}weekly_comers_agegrp6_tmp.gph"
		

		
**************************************************
** S3 leavers and comers: percent 
**************************************************

			
	** PERCENT **
	
		*percent leavers only
		* line
		tw line perc_leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			line perc_leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			line perc_leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			line perc_leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			line perc_leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			line perc_leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents Released by Age Group")  ///
			ylabel(0 "0%" 0.001 "0.10%" 0.002 "0.20%" 0.003 "0.30%" 0.004 "0.40%" 0.005 "0.50%", labsize(small)) ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(i_leavers_0114)
		graph save "Graph" "${goutputs}perc_weekly_leavers_agegrp6_tmp.gph"

		
		*percent leavers only
		* lpoly
		tw lpoly perc_leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			lpoly perc_leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			lpoly perc_leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents Released")  ///
			ylabel(0 "0%" 0.0003 "0.03%" 0.0006 "0.06%" 0.0009 "0.09%" 0.0012 "0.12%" 0.0015 "0.15%") ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(lpoly_leavers_0114)
		graph save "Graph" "${goutputs}lp_perc_weekly_comers_agegrp6_tmp.gph", replace
		
			
**************************************************
** SUPPLEMENT
**************************************************

**************************************************
** S1 chartables comparing older vs. not older
**************************************************
	
		chartable ///
			 not_recent_cohort always_cohort ever_quarantined ever_isolated ever_tested ever_positive_new ever_c19hosp_er ever_c19hosp_hosp ever_c19hosp_icu ever_died ///
			[fw=daily_active_cases] ,	${graph_opts} command(logit) or ///
			rhs(older male i.race2 perc_anycell perc_anydorm perc_anyother) ///
			case0(Less than 55) case1(55 or older) ///
				xsize(8)
				graph export "${goutputs}chartable_weight_older_notolder_race2.png", as(png) name("Graph") replace
		*/

	
	* with race3, "inst FE"


		* unweighted, WITH cell / dorm percent, WO inst FE
		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died not_recent_cohort day_fully_vax ///
			,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_unweight_older_notolder_race3_raw.gph", replace
				graph export "${goutputs}chartable_unweight_older_notolder_race3_raw.png", as(png) name("Graph") replace
				
				
		* unweighted, WITH cell / dorm percent, WITH inst FE
		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died not_recent_cohort day_fully_vax ///
			,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_unweight_older_notolder_race3_raw.gph", replace
				graph export "${goutputs}chartable_unweight_older_notolder_race3_raw.png", as(png) name("Graph") replace

		* weighted, WITH cell / dorm percent, WITH inst FE

		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died not_recent_cohort day_fully_vax ///
			[fw=daily_active_cases] ,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_weight_older_notolder_race3_raw.gph", replace
				graph export "${goutputs}chartable_weight_older_notolder_race3_raw.png", as(png) name("Graph") replace

				

use "data/intermediate/2021-1013 cdcr_residents_collapsed_master_master_weights.dta", clear
				
		preserve
		* reviewer 1
		keep if always_cohort == 1 | ever_died == 1
		* same as above but always_cohort == 1 & ever_died == 1
				
		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died day_fully_vax /// not_recent_cohort
			[fw=daily_active_cases] ,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_weight_older_notolder_race3_raw_alwayscohort_ordied.gph", replace
				graph export "${goutputs}chartable_weight_older_notolder_race3_raw_alwayscohort_ordied.png", as(png) name("Graph") replace
			
			
		keep if always_cohort == 1 
		* same as above but always_cohort == 1
				
		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu  day_fully_vax /// not_recent_cohort
			[fw=daily_active_cases] ,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_weight_older_notolder_race3_raw_alwayscohort.gph", replace
				graph export "${goutputs}chartable_weight_older_notolder_race3_raw_alwayscohort.png", as(png) name("Graph") replace
				
				
		restore
use "data/intermediate/2021-1013 cdcr_residents_collapsed_master_master_weights.dta", clear
		preserve
		keep if always_cohort == 0
		* same as above but always_cohort == 0
				
		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu not_recent_cohort day_fully_vax /// not_recent_cohort pos_ever_died
			[fw=daily_active_cases] ,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_weight_older_notolder_race3_raw_notalwayscohort.gph", replace
				graph export "${goutputs}chartable_weight_older_notolder_race3_raw_notalwayscohort.png", as(png) name("Graph") replace
			
					
		chartable ///
			 ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu not_recent_cohort day_fully_vax /// not_recent_cohort pos_ever_died
			,	${graph_opts} command(logit) or ///
			rhs(older male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
			case0(Less than 55) case1(55 or older) /// 		 regopts(cl(inst_id)) ///
				xsize(8)
				graph save "Graph" "${goutputs}chartable_weightno_older_notolder_race3_raw_notalwayscohort.gph", replace
				graph export "${goutputs}chartable_weightno_older_notolder_race3_raw_notalwayscohort.png", as(png) name("Graph") replace
			
				
		restore

**************************************************
** S2 leavers and comers: totals
**************************************************
	
	** TOTALS **
	use  "data/tmp/cdcr_residents_201901_202110_weeks.dta", clear
	
		*total comers and leavers
		tw line comers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) || line leavers_w date if agegrp6 == 1, lcolor(ebblue*0.8) lp(dash) || ///
			line comers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) || line leavers_w date if agegrp6 == 2, lcolor(ebblue*0.5) lp(dash) || ///
			line comers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || line leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(dash) || ///
			line comers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) || line leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(dash) || ///
			line comers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || line leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(dash) || ///
			line comers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid) || line leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(dash) ///
				legend(col(1) pos(3) symysize(small))  xline(21975)
		 graph save "Graph" "${goutputs}weekly_comers_leavers_agegrp6_tmp.gph"
		 
		*total comers only
		tw line comers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			line comers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			line comers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			line comers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			line comers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			line comers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Total Weekly Intake to CDCR by Age Group")  ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6))
		graph save "Graph" "${goutputs}weekly_comers_agegrp6_tmp.gph"
		
		*total leavers only
		tw line leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			line leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			line leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			line leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			line leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			line leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Total Weekly Releases by Age Group")  ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6))
		graph save "Graph" "${goutputs}weekly_comers_agegrp6_tmp.gph"
		

		
**************************************************
** S2 leavers and comers: percent 
**************************************************

		*percent comers only
		tw line perc_comers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			line perc_comers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			line perc_comers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			line perc_comers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			line perc_comers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			line perc_comers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents who were New Intakes by Age Group")  ///
			ylabel(0 "0%" 0.0005 "0.05%" 0.001 "0.10%" 0.0015 "0.15%") ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(i_comers_0114)
		graph save "Graph" "${goutputs}perc_weekly_comers_agegrp6_tmp.gph"
		
		*percent leavers only
		tw line perc_leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			line perc_leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			line perc_leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			line perc_leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			line perc_leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			line perc_leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents Released by Age Group")  ///
			ylabel(0 "0%" 0.001 "0.10%" 0.002 "0.20%" 0.003 "0.30%" 0.004 "0.40%" 0.005 "0.50%", labsize(small)) ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(i_leavers_0114)
		graph save "Graph" "${goutputs}perc_weekly_comers_agegrp6_tmp.gph"

		
		gr combine i_leavers_0114 i_comers_0114 ///
	,	cols(1)
	
	graph save "Graph" "${goutputs}perc_weekly_comers_leavers_raw.gph", replace
	graph export "${goutputs}perc_weekly_comers_leavers_raw.png", as(png) name("Graph")
		
	** PERCENT **
		
		*percent comers only
		tw lpoly perc_comers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			lpoly perc_comers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			lpoly perc_comers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			lpoly perc_comers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			lpoly perc_comers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			lpoly perc_comers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents as New Intakes")  ///
			ylabel(0 "0%" 0.0003 "0.03%" 0.0006 "0.06%" 0.0009 "0.09%" 0.0012 "0.12%" 0.0015 "0.15%") ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(lpoly_comers_0114)
		graph save "Graph" "${goutputs}lp_perc_weekly_comers_agegrp6_tmp.gph", replace
		
		*percent leavers only
		tw lpoly perc_leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			lpoly perc_leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			lpoly perc_leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents Released")  ///
			ylabel(0 "0%" 0.0003 "0.03%" 0.0006 "0.06%" 0.0009 "0.09%" 0.0012 "0.12%" 0.0015 "0.15%") ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(lpoly_leavers_0114)
		graph save "Graph" "${goutputs}lp_perc_weekly_comers_agegrp6_tmp.gph", replace

		
		gr combine lpoly_leavers_0114 lpoly_comers_0114 ///
	,	cols(1)
	
	graph save "Graph" "${goutputs}lp_perc_weekly_comers_leavers_raw.gph", replace
	graph export "${goutputs}lp_perc_weekly_comers_leavers_raw.png", as(png) name("Graph")
		

		
		
**************************************************
** other exhibit 1 -- comers minus leavers: percent
**************************************************

gen neg_perc_comers_w = -perc_comers_w
gen perc_comers_minus_leavers = perc_comers_w - perc_leavers_w
		
		*percent comers (dash) + leavers (solid) only
		tw lpoly perc_comers_minus_leavers date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			lpoly perc_comers_minus_leavers date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			lpoly perc_comers_minus_leavers date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			lpoly perc_comers_minus_leavers date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			lpoly perc_comers_minus_leavers date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			lpoly perc_comers_minus_leavers date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly %Δ of All Residents (Coming - Leaving)")  ///
			ylabel(-0.0018 "-0.18%" -0.0015 "-0.15%" -0.0012 "-0.12%" -0.0009 "-0.09%" -0.0006 "-0.06%" -0.0003 "-0.03%" 0 "0%" 0.0003 "0.03%" 0.0006 "0.06%" ) ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(lpoly8)
		graph save "Graph" "${goutputs}lp_perc_weekly_leavers_comers_agegrp6.gph", replace
			* saved as lp_perc_weekly_leavers_comers_agegrp6_edited.gph
			* saved as Exhibit5 lp_perc_weekly_leavers_comers_agegrp6.gph
		
		*percent leavers only
		tw lpoly perc_leavers_w date if agegrp6 == 1, lcolor(ebblue*.8) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 2, lcolor(ebblue*.5) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 3, lcolor(ebblue*0.2) lp(solid) || ///
			lpoly perc_leavers_w date if agegrp6 == 4, lcolor(red*0.2) lp(solid) ||  ///
			lpoly perc_leavers_w date if agegrp6 == 5, lcolor(red*0.5) lp(solid) || ///
			lpoly perc_leavers_w date if agegrp6 == 6, lcolor(red*0.8) lp(solid)  ///
			xlab(, angle(45))  ytit("Weekly % of All Residents Released")  ///
			ylabel(0 "0%" 0.0003 "0.03%" 0.0006 "0.06%" 0.0009 "0.09%" 0.0012 "0.12%" 0.0015 "0.15%") ///
				xline(21975) legend(col(1) pos(3) symysize(small)  ///
						label(1 "<35 years") ///
						label(2 "35 to 44") ///
						label(3 "45 to 54") ///
						label(4 "55 to 64") ///
						label(5 "65 to 74") ///
						label(6 "75 or older") ///
						order(1 2 3 4 5 6)) name(lpoly_leavers_0114)
		graph save "Graph" "${goutputs}lp_perc_weekly_comers_agegrp6_tmp.gph", replace
		
		

	
**************************************************
** other exhibit 2 -- forest plots
**************************************************
	
	
	*forest logit (check2 foreign) , t(check) or
	
		forest logit ///
			( ever_tested ever_quarantined ever_isolated ever_positive_new pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died always_cohort not_recent_cohort) [fw=daily_active_cases]  ///
		,	treatment(older) or controls(male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) b critical(0.05) graph($graph_opts)
			graph save "Graph" "${goutputs}forest_weight_older_notolder_race3_instFE.gph", replace
				
				
**************************************************
** other exhibit 3 -- iemargins plots
**************************************************
						
	** iemargins

	/*colorpalette, vertical n(10): HCL heat
	colorpalette viridis, n(10) nograph  reverse
	grstyle set color plasma, n(4) reverse
	colorpalette viridis
	grstyle set color Greens, n(8)
	*/
	grstyle set color ebblue
	

	* without cell/dorm
		foreach var of varlist ever_positive_new {
	iemargins `var' [fw=daily_active_cases] ///
		, treatment(agegrp2) ///
		controls(male i.race3 ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
		graph( title(, justification(left) color(black) span pos(11)) ///
		 graphregion(color(white) lc(white) lw(med)) bgcolor(white) ///
		xlab(,angle(45)) ///
		ylab(0 "0%" 0.05 "5%" 0.1 "10%" 0.15 "15%" 0.2 "20%" 0.25 "25%" 0.3 "30%" 0.35 "35%" 0.4 "40%" 0.45 "45%" 0.5 "50%" 0.55 "55%",angle(0) nogrid) xtit(,placement(left)  justification(left)) ///
		yscale(noline) xscale(noline) legend(region(lc(none) fc(none))) ///
		plotopts(lc(white) lw(thin) la(center) fi(100)) ciopts(lc(black)) name(i_`var'_woh)) 
		}
		
	foreach var of varlist pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died {
	iemargins `var' [fw=daily_active_cases] ///
		, treatment(agegrp2) ///
		controls(male i.race3 ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
		graph( title(, justification(left) color(black) span pos(11)) ///
		 graphregion(color(white) lc(white) lw(med)) bgcolor(white) ///
		xlab(,angle(45)) ///
		ylab(0 "0%" 0.05 "5%" 0.1 "10%" 0.15 "15%" 0.2 "20%",angle(0) nogrid) xtit(,placement(left)  justification(left)) ///
		yscale(noline) xscale(noline) legend(region(lc(none) fc(none))) ///
		plotopts(lc(white) lw(thin) la(center) fi(100)) ciopts(lc(black)) name(i_`var'_woh)) 
		}
	
	gr combine i_ever_positive_new_woh i_pos_ever_c19hosp_er_woh i_pos_ever_c19hosp_hosp_woh i_pos_ever_c19hosp_icu_woh i_pos_ever_died_woh ///
	, 	rows(1)
	
	graph save "Graph" "${goutputs}iemargins_5outcomes_woh_raw.gph", replace
	graph export "${goutputs}iemargins_5outcomes_woh_raw.png", as(png) name("Graph")


	


	* with cell / dorm
	foreach var of varlist ever_positive_new {
	iemargins `var' [fw=daily_active_cases] ///
		, treatment(agegrp2) ///
		controls(male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
		graph( title(, justification(left) color(black) span pos(11)) ///
		 graphregion(color(white) lc(white) lw(med)) bgcolor(white) ///
		xlab(,angle(45)) ///
		ylab(0 "0%" 0.05 "5%" 0.1 "10%" 0.15 "15%" 0.2 "20%" 0.25 "25%" 0.3 "30%" 0.35 "35%" 0.4 "40%" 0.45 "45%" 0.5 "50%" 0.55 "55%",angle(0) nogrid) xtit(,placement(left)  justification(left)) ///
		yscale(noline) xscale(noline) legend(region(lc(none) fc(none))) ///
		plotopts(lc(white) lw(thin) la(center) fi(100)) ciopts(lc(black)) name(i_`var'_wh)) 
		}
		
	foreach var of varlist pos_ever_c19hosp_er pos_ever_c19hosp_hosp pos_ever_c19hosp_icu pos_ever_died {
	iemargins `var' [fw=daily_active_cases] ///
		, treatment(agegrp2) ///
		controls(male i.race3 perc_anycell perc_anydorm ///
ever_ASP ever_CAC ever_CAL ever_CCC ever_CCI ever_CCWF ever_CEN ever_CHCF ever_CIM ever_CIW ever_CMC ever_CMF ever_COR ever_CRC ever_CTF ever_CVSP ever_DVI ever_FSP ever_HDSP ever_ISP ever_KVSP ever_LAC ever_MCSP ever_NKSP ever_PBSP ever_PVSP ever_RJD ever_SAC ever_SATF ever_SCC ever_SOL ever_SQ ever_SVSP ever_VSP ever_WSP) ///
		graph( title(, justification(left) color(black) span pos(11)) ///
		 graphregion(color(white) lc(white) lw(med)) bgcolor(white) ///
		xlab(,angle(45)) ///
		ylab(0 "0%" 0.05 "5%" 0.1 "10%" 0.15 "15%" 0.2 "20%",angle(0) nogrid) xtit(,placement(left)  justification(left)) ///
		yscale(noline) xscale(noline) legend(region(lc(none) fc(none))) ///
		plotopts(lc(white) lw(thin) la(center) fi(100)) ciopts(lc(black)) name(i_`var'_wh)) 
		}
	
	gr combine i_ever_positive_new_wh i_pos_ever_c19hosp_er_wh i_pos_ever_c19hosp_hosp_wh i_pos_ever_c19hosp_icu_wh i_pos_ever_died_wh ///
	, 	rows(1)
	
	graph save "Graph" "${goutputs}iemargins_5outcomes_wh_raw.gph", replace
	graph export "${goutputs}iemargins_5outcomes_wh_raw.png", as(png) name("Graph")

		
