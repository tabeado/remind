*** |  (C) 2006-2023 Potsdam Institute for Climate Impact Research (PIK)
*** |  authors, and contributors see CITATION.cff file. This file is part
*** |  of REMIND and licensed under AGPL-3.0-or-later. Under Section 7 of
*** |  AGPL-3.0, you are granted additional permissions described in the
*** |  REMIND License Exception, version 1.0 (see LICENSE file).
*** |  Contact: remind@pik-potsdam.de
*** SOF ./modules/45_carbonprice/functionalForm/postsolve.gms

***-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*** Part 0 (Actual CO2 budget): If iterative_target_adj = 0, 7, 9 or 10, compute actual CO2 peak budget in current iteration.
***  If iterative_target_adj = 5, compute actual CO2 end-of-century budget in current iteration.
***  If iterative_target_adj = 10, compute actual CO2 end-of century budget in current iteration additionally to the peak budget 
***-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(cm_iterative_target_adj = 5,  !! End-of-century budget
  s45_actualbudgetco2 = sum(t$(t.val eq 2100),pm_actualbudgetco2(t)); 
else !! Peak budget
  s45_actualbudgetco2 = smax(t$(t.val le cm_peakBudgYr AND t.val le 2100),pm_actualbudgetco2(t));
  o45_peakBudgYr_Itr(iteration) = cm_peakBudgYr;
);

if((cm_iterative_target_adj ge 10) OR ((cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj gt 0)),  !! End-of-century budget (additionally to the peak budget)
  if(cm_CDRtypesInTarget eq 1, 
    s45_actualbudgetco2_2100 = sum(t$(t.val eq 2100), pm_actualbudgetco2(t)); 
  elseif (cm_CDRtypesInTarget eq 2),
    s45_actualbudgetco2_2100 = sum(t$(t.val eq 2100), pm_actualbudgetco2_noLULUCF(t)); 
  elseif (cm_CDRtypesInTarget eq 3),
    s45_actualbudgetco2_2100 = sum(t$(t.val eq 2100), pm_actualbudgetco2_noNetNegLU(t)); 
    );
  display s45_actualbudgetco2_2100;
);

display pm_actualbudgetco2, s45_actualbudgetco2;

*** Copied from postsolve algorithm for cm_iterative_target_adj = 5. TODO: Check where cm_emiscen eq 6 is used and if this should be kept.
if ((cm_emiscen eq 6) AND (cm_iterative_target_adj eq 5), 
	if(o_modelstat eq 2 AND ord(iteration)<cm_iteration_max ,   !!only for optimal iterations, and not after the last one
	display sm_budgetCO2eqGlob;		
		sm_budgetCO2eqGlob = sm_budgetCO2eqGlob * (cm_budgetCO2from2020/s45_actualbudgetco2);
		pm_budgetCO2eq(regi) = pm_budgetCO2eq(regi) * (cm_budgetCO2from2020/s45_actualbudgetco2);
	else
		sm_budgetCO2eqGlob = sm_budgetCO2eqGlob;
	);
	display sm_budgetCO2eqGlob;
);

*** Only run adjustment of carbon price trajectory if cm_emiscen eq 9 and if cm_iterative_target_adj is equal to 5,7, 9 or larger than 10 (i.e. separate CDR and emi prices post peak or throughout).
if((cm_emiscen eq 9) AND ((cm_iterative_target_adj eq 5) OR (cm_iterative_target_adj eq 7) OR (cm_iterative_target_adj ge 9)),

*** Save pm_taxCO2eq and p45_taxCO2eq_anchor over iterations for debugging
pm_taxCO2eq_iter(iteration,ttot,regi) = pm_taxCO2eq(ttot,regi);
p45_taxCO2eq_anchor_iter(iteration,ttot) = p45_taxCO2eq_anchor(ttot);
if (cm_iterative_target_adj ge 10,
  pm_taxCDR_iter(iteration,ttot,regi) = pm_taxCDR(ttot,regi); !! this is currently happening in two places, also postsolve of core - why? Is it necessary?
  p45_taxCDR_anchor_iter(iteration,t) = p45_taxCDR_anchor(t);
);

*** Compute absolute deviation of actual budget from target budget
sm_globalBudget_absDev = s45_actualbudgetco2 - cm_budgetCO2from2020;

!! End-of-century budget (additionally to the peak budget)
if ((cm_iterative_target_adj ge 10) OR (cm_iterative_target_adj eq 9 AND (cm_postPeakCpAdj gt 0)),  
  !! Compute absolute deviation of actual 2100 budget from target 2100 budget
  sm_globalBudget2100_absDev = s45_actualbudgetco2_2100 - cm_addbudgetCO2from2020to2100;
  p45_globalBudget2100_absDev_iter(iteration) = sm_globalBudget2100_absDev;
);

*** ---------------------------------------------------------------------------------------------------------------------
*** --------ALGORITHM for cm_iterative_target_adj eq 11: Separate CDR price throughout ----------------------------------
*** ---------------------------------------------------------------------------------------------------------------------
if(cm_iterative_target_adj eq 11,
 
 !! 0) save General Information about the change of the budget and CDR tax change (as of iter 2)
    !! (Possible in all iterations, even if price was not changed, because budget change is the denominator - but not always meaningful)
        if((iteration.val ge 2),
        s45_taxCdrChange =  (p45_taxCDR_anchor_iter(iteration,"2100") - p45_taxCDR_anchor_iter(iteration-1,"2100"))$(cm_CDRpriceShape ne 4)
                            + (sum(ttot2$(ttot2.val eq cm_startYear), p45_taxCDR_anchor_iter(iteration,ttot2))
                                - sum(ttot2$(ttot2.val eq cm_startYear), p45_taxCDR_anchor_iter(iteration-1,ttot2)))$(cm_CDRpriceShape eq 4);
        p45_taxCdrChange_iter(iteration) =  s45_taxCdrChange;

        !!s45_EOCbudgetChange = p45_globalBudget2100_absDev_iter(iteration) -  p45_globalBudget2100_absDev_iter(iteration - 1);
        s45_EOCbudgetChange = p45_globalBudget2100_absDev_iter(iteration) -  p45_globalBudget2100_absDev_iter(iteration - 1);
        p45_EOCbudgetChange_iter(iteration) = s45_EOCbudgetChange;

        s45_TaxBudget_ChangeSlope = s45_taxCdrChange / s45_EOCbudgetChange; 
        p45_TaxBudget_ChangeSlope_iter(iteration) = s45_TaxBudget_ChangeSlope; 
        ); !! iteration > 2 for saving information

 !! 1) Check if within tolerance.
  if (abs(p45_globalBudget2100_absDev_iter(iteration)) le cm_budget2100CO2_absDevTol,
    !! if deviation within tolerance -> no rescaling needed
    s45_factorRescale_CDRtax = 1; 
    p45_factorRescale_CDRtax_iter(iteration) = s45_factorRescale_CDRtax;

  else
 !! 2) Calculate CDRtax rescaling factor if not yet within tolerance 
   !!_____ 2a) if tax change leads to budget change in the right direction:_____
    if ((iteration.val ge 3) AND (s45_TaxBudget_ChangeSlope lt 0), 
      !! needed budget change * (taxChange / budgetChange)
      s45_neededCDRtaxChange2100 =  - sm_globalBudget2100_absDev * s45_TaxBudget_ChangeSlope ;
      p45_neededCDRtaxChange2100_iter(iteration) = s45_neededCDRtaxChange2100;

      p45_newtaxCDR_anchor2100(iteration) =  p45_taxCDR_anchor_iter(iteration,"2100")$(cm_CDRpriceShape ne 4) 
                                              + sum(ttot2$(ttot2.val eq cm_startYear), p45_taxCDR_anchor_iter(iteration,ttot2))$(cm_CDRpriceShape eq 4) 
                                              + s45_neededCDRtaxChange2100;
      
      s45_factorRescale_CDRtax =  (p45_newtaxCDR_anchor2100(iteration) / p45_taxCDR_anchor("2100"))$(cm_CDRpriceShape ne 4)
                                  +  (p45_newtaxCDR_anchor2100(iteration) / sum(ttot2$(ttot2.val eq cm_startYear), p45_taxCDR_anchor(ttot2)))$(cm_CDRpriceShape eq 4);
      p45_factorRescale_CDRtax_iter(iteration) = s45_factorRescale_CDRtax;
  
   !!_____ 2b) No meaningful price-budget information available => get ratio (actual / target budget) 
                !! NOTE: The current implementation only works for positive targets! (see EOC Work)
     else
      s45_factorRescale_CDRtax = s45_actualbudgetco2_2100 / cm_addbudgetCO2from2020to2100;
      p45_factorRescale_CDRtax_iter(iteration) = s45_factorRescale_CDRtax; 
    ); !! calculate CDR rescaling factor end
 
  ); !! Check if you need to rescale the CDRtax

 !! 3) Funnel the CDRtax rescaling factor 
  p45_factorRescale_CDRtax_Funneled(iteration)                                        
        = max(min( cm_funnelFactor * EXP( -cm_funnelExponent * iteration.val) + 1 + cm_funnelLower,   !! a) a maximum adjustment value which decreases with the number of iterations
                    s45_factorRescale_CDRtax),          
              1/ ( cm_funnelFactor * EXP( -cm_funnelExponent * iteration.val) + 1 + cm_funnelLower)   !! b) a minimum adjustment value which increases with the number of iterations (0.95 for iter 25)
          );
); !! (cm_iterative_target_adj eq 11),

*** ---------------------------------------------------------------------------------------------------------------------
*** --------Preparation for endogenous adjustment of the NNE tax  ----------------------------------
*** ---------------------------------------------------------------------------------------------------------------------

if((cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 1),
  s45_TaxBudget_ChangeSlope =  sum(t2$(t2.val eq 2100), p45_taxCO2eq_anchor(t2) - p45_taxCO2eq_anchor_iter(iteration - 1, t2))
                            / (p45_globalBudget2100_absDev_iter(iteration) -  p45_globalBudget2100_absDev_iter(iteration - 1));
  p45_TaxBudget_ChangeSlope_iter(iteration) = s45_TaxBudget_ChangeSlope; 

  !! if the slope is negative, then update the used slope. 
  if(s45_TaxBudget_ChangeSlope < 0,
    s45_TaxBudget_ChangeSlopeBest = s45_TaxBudget_ChangeSlope;
  );
  p45_TaxBudgetSlopeBest_iter(iteration) = s45_TaxBudget_ChangeSlopeBest;
); 


if((cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 2),
  !! 0) save information about the effective NNE tax change
  sm_effectiveNNEtax = p45_taxCO2eq_anchor("2100") * (1-cm_frac_NetNegEmi);
  pm_effectiveNNEtax_iter(iteration) = sm_effectiveNNEtax;
  
  s45_TaxBudget_ChangeSlope = (pm_effectiveNNEtax_iter(iteration) - pm_effectiveNNEtax_iter(iteration-1))
                            / (p45_globalBudget2100_absDev_iter(iteration) -  p45_globalBudget2100_absDev_iter(iteration - 1));
  p45_TaxBudget_ChangeSlope_iter(iteration) = s45_TaxBudget_ChangeSlope; 

  !! if the slope is negative, then update the used slope. 
  if(s45_TaxBudget_ChangeSlope < 0,
    s45_TaxBudget_ChangeSlopeBest = s45_TaxBudget_ChangeSlope;
  );
  p45_TaxBudgetSlopeBest_iter(iteration) = s45_TaxBudget_ChangeSlopeBest;
 
  !! 1) Check if within Tolerance
  if (abs(p45_globalBudget2100_absDev_iter(iteration)) le cm_budget2100CO2_absDevTol,
    !! if deviation within tolerance -> no rescaling needed, keep the same effective tax
    p45_neweffectiveNNEtax(iteration) = sm_effectiveNNEtax;
  
  !! 2) Calculate rescaling factor for effective tax level
  else 
  !!_____ 2a) if there is information on reaction to change: __________
    if ((iteration.val ge 3) AND (s45_TaxBudget_ChangeSlopeBest lt 0), 
      !! needed budget change * (taxChange / budgetChange)
      s45_neededCDRtaxChange2100 =  - sm_globalBudget2100_absDev * s45_TaxBudget_ChangeSlopeBest;
      p45_neededCDRtaxChange2100_iter(iteration) = s45_neededCDRtaxChange2100;

      p45_neweffectiveNNEtax(iteration) =  sm_effectiveNNEtax + s45_neededCDRtaxChange2100;
      
    !!_____ 2b) No meaningful price-budget information available => get ratio (actual / target budget) and scale the previously effective tax
                  !! NOTE: This current implementation only works for positive targets! (see EOC Work)
     else
      p45_neweffectiveNNEtax(iteration)  = (s45_actualbudgetco2_2100 / cm_addbudgetCO2from2020to2100) * sm_effectiveNNEtax;
    ); !! calculate new NNE tax 
   ); !! tolerance Check
); !! (cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 2)


***--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*** Part I and II (Global anchor trajectory and post-peak behaviour): Adjustment of global anchor trajectory to meet (peak or end-of-century) CO2 budget target prescribed via cm_budgetCO2from2020.
***    If iterative_target_adj = 7, 9 or 10, cm_peakBudgYr automatically adjusted (within the time window 2040--2100)
***--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*** --------ALGORITHM for cm_iterative_target_adj eq 10 --------------------------------------------------------------------------------------------
***  Calculate the new annual CO2 and CDR price increase after the peak year
!! GEKÜRZTE KOMMENTARVARIANTE 
if (cm_iterative_target_adj eq 10, 
  p45_ord_iteration(iteration) = ord(iteration)
  display s45_messageX, p45_ord_iteration;

******* (I): Determine CDR price slope change
  !! Do not update CDR price slope in early iterations.
  if (ord(iteration) lt 15,
    o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = 0;

  !! Case 1) Initialize the difference in the post-peak CDR slope in iteration 15,
   elseif (ord(iteration) eq 15),
    display s45_message01, p45_ord_iteration;
    !! Calculate difference of CDR price slope to the next iteration from the
    !! budget difference of the current iteration.
    o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = max(min(p45_globalBudget2100_absDev_iter(iteration) / 15, 5), -5);

  !! Case 2) iteration number is high enough + not yet in the last iteration, 
  elseif (ord(iteration) < cm_iteration_max AND ord(iteration) gt 15),
    display s45_message02, p45_ord_iteration, cm_iteration_max;
    
    !! 2.1) Iteration is feasible:
    if ((o_modelstat eq 2) OR (o_modelstat eq 7),
      display s45_message03, o_modelstat;
    
      !! 2.11) Close enough to target --> no slope update
      if (abs(p45_globalBudget2100_absDev_iter(iteration)) le cm_budget2100CO2_absDevTol,
        display s45_message04, p45_globalBudget2100_absDev_iter, cm_budget2100CO2_absDevTol;
        o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = 0;     
        !! Save the ratio between cumulative emissions and slope differences
        !! to be used in the next iteration missing the target.
        if (s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good eq 0, !! Only update if the last iteration had a slope change
          display s45_message05, s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good;
          s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good = 
            (p45_globalBudget2100_absDev_iter(iteration) - p45_globalBudget2100_absDev_iter(iteration-1))
            / (o45_taxCDR_IncAfterPeakBudgYr_iter(iteration) - o45_taxCDR_IncAfterPeakBudgYr_iter(iteration-1));
        );
      
      !! 2.12) Not close enough to the target --> derive a new slope
      else
        display s45_message06, p45_globalBudget2100_absDev_iter, cm_budget2100CO2_absDevTol;

        !! 2.121) Slope and budget change information from the last two iterations is available:
             !! Slope Change = needed budget change (current deviation) * (PreviousSlopeChange / PreviousBudgetChange)
              !! NOTE: may need to introduce a check if the change was counterintuitive, e.g. EOC decreased despite a steeper slope (i.e. PreviousSlopeChange / PreviousBudgetChange > 0)
        if (s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good eq 0, !! i.e. there was a change of slope last iteration
          display s45_message07, s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good;
          o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = p45_globalBudget2100_absDev_iter(iteration)  * o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration)  
                                                              / (p45_globalBudget2100_absDev_iter(iteration-1) - p45_globalBudget2100_absDev_iter(iteration));
            
        !! 2.122) Slope and budget change info *not* available (because close enough to target or infeasible):
        !! use the latest ratio between cumulative emission differences and CDR slope differences.
        !! Afterwards, set this helper parameter back to zero so that it can be adjusted again in case 2.11 or 2.2
        else
          display s45_message08, s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good;
          o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = - p45_globalBudget2100_absDev_iter(iteration) / s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good;
          s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good = 0;
        );
        !! applying to each 2.12 case: Limit the jumps in CDR slope differences between iterations
        o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = max(min(o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1), 2), -2);
      ); !! if (abs(p45_globalBudget2100_absDev_iter(iteration)) le cm_budget2100CO2_absDevTol,
    
    !! 2.2) Iteration is infeasible -> do not update the CO2 and CDR price slopes
    else 
      display s45_message09, o_modelstat;
      o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1) = 0;
      !! 2.21) Ratio between cumulative 2100 emission differences and CDR slope differences 
      !!      does *not* yet exist --> save information from previous runs
        if (s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good eq 0,
        display s45_message10, s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good;
        
        !! 2.211) It is iteration 16 (i.e. no o45_taxCDR_IncAfterPeakBudgYr_iter(iteration-2))
        if (ord(iteration) eq 16,    display s45_message11, p45_ord_iteration;
          s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good = -15; !! very simple approximation
        
        !! 2.212) Save information from the two previous iterations:
        else
          display s45_message12, p45_ord_iteration;
          s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good = 
            (p45_globalBudget2100_absDev_iter(iteration-1) - p45_globalBudget2100_absDev_iter(iteration-2))
            / (o45_taxCDR_IncAfterPeakBudgYr_iter(iteration-1) - o45_taxCDR_IncAfterPeakBudgYr_iter(iteration-2));
        
        );
      !! 2.22) Ratio between cumulative 2100 emission differences and CDR slope differences 
      !!     *does* exist --> no need to save additional information
      else 
        display s45_message13, s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good;
      ); !! (s45_ratio_emi_IterDiff_taxCDR_IterDiff_last_good eq 0,
    ); !! if ((o_modelstat eq 2) OR (o_modelstat eq 7),
 
  !! Case 3): Last iteration: do nothing
  else  
    display s45_message14, p45_ord_iteration, cm_iteration_max;
  ); !! if (ord(iteration) eq 15, elseif (ord(iteration)<cm_iteration_max AND ord(iteration) gt 15),
******* (II): Update the CDR and possibly CO2 tax slope based on the above derived changes.
  !! Update CDR tax slope
  s45_taxCDR_IncAfterPeakBudgYr = s45_taxCDR_IncAfterPeakBudgYr + o45_taxCDR_IncAfterPeakBudgYr_iterDiff(iteration+1);
  
  !! Update the CO2 tax slope; only has an effect if the derived slope is positive, i.e. Prices need to increase post peak to reach the EOC target
  s45_taxCO2_IncAfterPeakBudgYr_current = s45_taxCO2_IncAfterPeakBudgYr;
  s45_taxCO2_IncAfterPeakBudgYr = max(0, s45_taxCDR_IncAfterPeakBudgYr);
  o45_taxCO2_IncAfterPeakBudgYr_iterDiff(iteration+1) = s45_taxCO2_IncAfterPeakBudgYr - s45_taxCO2_IncAfterPeakBudgYr_current;

  !! Save slope of next iteration for diagnostics.
  o45_taxCDR_IncAfterPeakBudgYr_iter(iteration+1) = s45_taxCDR_IncAfterPeakBudgYr;
  o45_taxCO2_IncAfterPeakBudgYr_iter(iteration+1) = s45_taxCO2_IncAfterPeakBudgYr;

  !! Display target deviation and current levels of CDR and CO2 price slopes.
  display p45_globalBudget2100_absDev_iter, o45_taxCDR_IncAfterPeakBudgYr_iterDiff, o45_taxCDR_IncAfterPeakBudgYr_iter, o45_taxCO2_IncAfterPeakBudgYr_iter;
  p45_ord_iteration(iteration) = 0;
  display s45_messageY;
); !! if(cm_iterative_target_adj eq 10,

*** --------ALGORITHM for cm_iterative_target_adj eq 5 , 9, 10 or 11 ----------------------------------------------------------------------------------------
*** --------A: calculate the new CO2 price path, beginning with the CO2 tax rescale factor----------------------------------------------------------
*** --------   this step applies for peak budget and end-of-century budget targets -----------------------------------------------------------------
if((cm_iterative_target_adj eq 5) OR (cm_iterative_target_adj eq 9) OR (cm_iterative_target_adj eq 10) OR (cm_iterative_target_adj eq 11),

  if((cm_iterative_target_adj eq 9) OR (cm_iterative_target_adj eq 10) or (cm_iterative_target_adj eq 11), !! stronger sensitivity of CO2 price adjustment to CO2 budget deviation for peak budget targets
    s45_factorRescale_taxCO2_exponent_before10 = 3;
    s45_factorRescale_taxCO2_exponent_from10 = 2;
  else !! less sensitivity of CO2 price adjustment to CO2 budget deviation for peak budget targets
    s45_factorRescale_taxCO2_exponent_before10 = 2;
    s45_factorRescale_taxCO2_exponent_from10 = 1;
  );

  if( (o_modelstat ne 2) OR (abs(sm_globalBudget_absDev) le cm_budgetCO2_absDevTol) OR (ord(iteration) = cm_iteration_max), 
    !! keep CO2 tax constant if model was not optimal, if maximal number of iterations is reached, or if budget already reached
    p45_factorRescale_taxCO2(iteration)          = 1;
    p45_factorRescale_taxCO2_Funneled(iteration) = p45_factorRescale_taxCO2(iteration);
  else !! adjust CO2 tax 
    if (s45_actualbudgetco2 > 0, !! if budget positive

      !! if end-of-century budget is higher than budget at peak point, AND end-of-century budget is already in the range of the target budget (+/- 50 GtC), treat as end-of-century budget 
      !! for this iteration. Only do this rough approach (jump to 2100) for the first iterations - at later iterations the slower adjustment of the peaking time should work better
      if( ((cm_iterative_target_adj eq 9) OR (cm_iterative_target_adj ge 10)) AND ( pm_actualbudgetco2("2100") > 1.1 * s45_actualbudgetco2 ) AND ( abs(cm_budgetCO2from2020 - s45_actualbudgetco2) < 50 ) AND (iteration.val < 12), 
        display iteration;
        display "this is likely an end-of-century budget with no net negative emissions at all. Shift cm_peakBudgYr to 2100";
        cm_peakBudgYr = 2100;
        !! due to the potential strong jump in cm_peakBudgYr, which implies that the CO2 price will increase over a longer time horizon,
        !! take the average of the budget at the old peak time and the new peak time
        s45_actualbudgetco2 = 0.5 * (pm_actualbudgetco2("2100") + s45_actualbudgetco2); 
      );

      !! CO2 tax rescale factor
      if(iteration.val lt 10,
        p45_factorRescale_taxCO2(iteration) = max(0.1, (s45_actualbudgetco2/cm_budgetCO2from2020) ) ** s45_factorRescale_taxCO2_exponent_before10;
      else
        p45_factorRescale_taxCO2(iteration) = max(0.1, (s45_actualbudgetco2/cm_budgetCO2from2020) ) ** s45_factorRescale_taxCO2_exponent_from10;
      );
      p45_factorRescale_taxCO2_Funneled(iteration) 
        = max(min( 2 * EXP( -0.15 * iteration.val ) + 1.01 ,p45_factorRescale_taxCO2(iteration)),
              1/ ( 2 * EXP( -0.15 * iteration.val ) + 1.01)
          );
    else !! if budget has turned negative, reduce CO2 price by 20%
      !! CO2 tax rescale factor
      p45_factorRescale_taxCO2(iteration) = 0.8;
      p45_factorRescale_taxCO2_Funneled(iteration) = p45_factorRescale_taxCO2(iteration);
    );
    display p45_taxCO2eq_anchor, p45_taxCO2eq_anchor_until2150, p45_factorRescale_taxCO2, p45_factorRescale_taxCO2_Funneled;

    !! Apply CO2 tax rescale factor
    p45_taxCO2eq_anchor_until2150(ttot)$(ttot.val ge 2005) = p45_taxCO2eq_anchor_until2150(ttot) * p45_factorRescale_taxCO2_Funneled(iteration);
    display p45_taxCO2eq_anchor_until2150;

    !! If functionalForm is linear, re-adjust global anchor trajectory to go through the point (cm_taxCO2_historicalYr, cm_taxCO2_historical) 
$ifThen.taxCO2functionalForm4 "%cm_taxCO2_functionalForm%" == "linear"
    p45_taxCO2eq_anchor_until2150(ttot)$(ttot.val ge s45_taxCO2_historicalYr) = s45_taxCO2_historical 
        + (sum(t2$(t2.val eq cm_peakBudgYr), p45_taxCO2eq_anchor_until2150(t2)) - s45_taxCO2_historical) / (cm_peakBudgYr - s45_taxCO2_historicalYr) !! Yearly increase of CO2 price that interpolates between cm_taxCO2_historical in cm_taxCO2_historicalYr and p45_taxCO2eq_anchor_until2150 in peak year
                                      * (ttot.val - s45_taxCO2_historicalYr) ;
    display p45_taxCO2eq_anchor_until2150;
$endIf.taxCO2functionalForm4 

    !! Use rescaled p45_taxCO2eq_anchor_until2150 as starting point for re-defining p45_taxCO2eq_anchor
    p45_taxCO2eq_anchor(ttot)$(ttot.val ge 2005) = p45_taxCO2eq_anchor_until2150(ttot);
    
    if(cm_iterative_target_adj = 9 or cm_iterative_target_adj = 11, !! After cm_peakBudgYr, the global anchor trajectory increases linearly with fixed annual increase given by cm_taxCO2_IncAfterPeakBudgYr
      p45_taxCO2eq_anchor(t)$(t.val gt cm_peakBudgYr) = sum(t2$(t2.val eq cm_peakBudgYr), p45_taxCO2eq_anchor_until2150(t2)) !! CO2 tax in peak budget year
                                                  + (t.val - cm_peakBudgYr) * cm_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;  !! increase by cm_taxCO2inc_after_peakBudgYr per year 
    );  
    !! Always set carbon price constant after 2100 to prevent huge taxes after 2100 and the resulting convergence problems
    p45_taxCO2eq_anchor(t)$(t.val gt 2100) = p45_taxCO2eq_anchor("2100");

    !! Compute difference for debugging
    pm_taxCO2eq_anchor_iterationdiff(t) = p45_taxCO2eq_anchor(t) - p45_taxCO2eq_anchor_iter(iteration,t);
    o45_taxCO2eq_anchor_iterDiff_Itr(iteration) = pm_taxCO2eq_anchor_iterationdiff("2100");

    display p45_taxCO2eq_anchor, pm_taxCO2eq_anchor_iterationdiff, o45_taxCO2eq_anchor_iterDiff_Itr;

  ); !! if( (o_modelstat ne 2) OR (abs(sm_globalBudget_absDev) le cm_budgetCO2_absDevTol) OR (ord(iteration) = cm_iteration_max), 
); !! if((cm_iterative_target_adj eq 5) OR (cm_iterative_target_adj eq 9) OR (cm_iterative_target_adj eq 10) OR (cm_iterative_target_adj = 11),


*** -------B: checking the peak timing, if cm_peakBudgYr is still correct or needs to be shifted-----------------------
*** --------  this step only applies for peak budget targets-----------------------------------------------------------
if((cm_iterative_target_adj eq 9) OR (cm_iterative_target_adj eq 10) OR (cm_iterative_target_adj eq 11),
  o45_diff_to_Budg(iteration) = (cm_budgetCO2from2020 - s45_actualbudgetco2);
  o45_totCO2emi_peakBudgYr(iteration) = sum(t$(t.val = cm_peakBudgYr), sum(regi2, vm_emiAll.l(t,regi2,"co2")) );
  o45_totCO2emi_allYrs(t,iteration) = sum(regi2, vm_emiAll.l(t,regi2,"co2") );
	
  !! calculate how fast emissions are changing around the peaking time to get an idea how close it is possible to get to 0 due to the 5(10) year time steps 	
  o45_change_totCO2emi_peakBudgYr(iteration) = sum(ttot$(ttot.val = cm_peakBudgYr), (o45_totCO2emi_allYrs(ttot-1,iteration) - o45_totCO2emi_allYrs(ttot+1,iteration) )/4 );  !! Only gives a tolerance range, exact value not important. Division by 4 somewhat arbitrary - could be 3 or 5 as well. 

  display cm_peakBudgYr, o45_diff_to_Budg, o45_peakBudgYr_Itr, o45_totCO2emi_allYrs, o45_totCO2emi_peakBudgYr, o45_change_totCO2emi_peakBudgYr;


  !!----B1: check if cm_peakBudgYr should be shifted left or right: 
  if( abs(o45_diff_to_Budg(iteration)) < 20, !! only think about shifting peakBudgYr if the budget is close enough to target budget
    display "close enough to target budget to check timing of peak year";
	 
	  !!  check if the target year was just shifted back left after being shifted right before
	  if ( (iteration.val > 2) AND ( o45_peakBudgYr_Itr(iteration - 1) > o45_peakBudgYr_Itr(iteration) ) AND ( o45_peakBudgYr_Itr(iteration - 2) = o45_peakBudgYr_Itr(iteration) ),
	    o45_pkBudgYr_flipflop(iteration) = 1; 
        display "flipflop observed (before loop)";
	  );
	 
    loop(ttot$(ttot.val = cm_peakBudgYr), !! look at the peak timing
      if(  ( (o45_totCO2emi_peakBudgYr(iteration) < -(0.1 + o45_change_totCO2emi_peakBudgYr(iteration)) ) AND (cm_peakBudgYr > 2040) ), !! no peaking time before 2040
        display "shift peakBudgYr left";
	      o45_peakBudgYr_Itr(iteration+1) =  pm_ttot_val(ttot - 1);                
        p45_taxCO2eq_anchor(t)$(t.val gt pm_ttot_val(ttot - 1)) = p45_taxCO2eq_anchor_until2150(ttot-1) + (t.val - pm_ttot_val(ttot - 1)) * cm_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;  !! increase by cm_taxCO2_IncAfterPeakBudgYr per year after peakBudgYr
       
	    elseif (( o45_totCO2emi_peakBudgYr(iteration) > (0.1 + o45_change_totCO2emi_peakBudgYr(iteration)) ) AND (cm_peakBudgYr < 2100)), !! if peaking time would be after 2100, keep 2100 budget year
        if(  (o45_pkBudgYr_flipflop(iteration) eq 1), !! if the target year was just shifted left after being shifted right, and would now be shifted right again
          display "peakBudgYr was left, right, left and is now supposed to be shifted right again -> flipflop, thus go into separate loop";
          o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration); !! don't shift right again immediately, but go into a different loop:
          o45_delay_increase_peakBudgYear(iteration) = 1;
	      elseif ( o45_delay_increase_peakBudgYear(iteration) eq 1 ),
	        display "still in separate loop trying to resolve flip-flop behavior";
	      	o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration); !! keep current peakBudgYr,
        else
	        display "shift peakBudgYr right";
          o45_peakBudgYr_Itr(iteration+1) =  pm_ttot_val(ttot + 1);  !! ttot+1 is the new peakBudgYr
	    	  loop(t$(t.val ge pm_ttot_val(ttot + 1)),
            p45_taxCO2eq_anchor(t) = p45_taxCO2eq_anchor_until2150(ttot+1) 
	    	                    + (t.val - pm_ttot_val(ttot + 1)) * cm_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;  !! increase by cm_taxCO2_IncAfterPeakBudgYr per year 
          );
	      );
      
	    else   !! don't do anything if the peakBudgYr is already at the corner values (2040, 2100) or if the emissions in the peakBudgYr are close enough to 0 (within the range of +/- o45_change_totCO2emi_peakBudgYr)
             o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration)
       );
    );
    cm_peakBudgYr = o45_peakBudgYr_Itr(iteration+1);
    display cm_peakBudgYr;
  );       
  p45_taxCO2eq_anchor(ttot)$((ttot.val ge 2005) AND (ttot.val le cm_peakBudgYr)) = p45_taxCO2eq_anchor_until2150(ttot); !! until peakBudgYr, take the contiuous price trajectory
   
  !!-----B2: if there was a flip-floping of cm_peakBudgYr in the previous iterations, try to overome this by adjusting the CO2 price path after the peaking year	
  if (o45_delay_increase_peakBudgYear(iteration) = 1,   
    display "not shifting peakBudgYr right, instead adjusting CO2 price for following year";
    loop(ttot$(ttot.val eq cm_peakBudgYr),  !! set ttot to the current peakBudgYr 
      loop(t2$(t2.val eq pm_ttot_val(ttot+1)),  !! set t2 to the following time step
        o45_factorRescale_taxCO2_afterPeakBudgYr(iteration) = 1 + max(sum(regi2,vm_emiAll.l(ttot,regi2,"co2"))/sum(regi2,vm_emiAll.l("2015",regi2,"co2")),-0.75) ; 
	  !! this was inspired by Christoph's approach. This value is 1 if emissions in the peakBudgYr are 0; goes down to 0.25 if emissions are <0 and approaching the size of 2015 emissions, and > 1 if emissions > 0. 
         
	  !! in case the normal linear extension still is not enough to get emissions to 0 after the peakBudgYr, shift peakBudgYr right again:
        if( ( o45_reached_until2150pricepath(iteration-1) eq 1 ) AND ( o45_totCO2emi_peakBudgYr(iteration) > (0.1 + o45_change_totCO2emi_peakBudgYr(iteration)) ), 
          display "price in following year reached original path in previous iteration and is still not enough -> shift peakBudgYr to right";
          o45_delay_increase_peakBudgYear(iteration+1) = 0;  !! probably is not necessary
          o45_reached_until2150pricepath(iteration) = 0;
          o45_peakBudgYr_Itr(iteration+1) = t2.val;        !! shift PeakBudgYear to the following time step
          p45_taxCO2eq_anchor(t2) = p45_taxCO2eq_anchor_until2150(t2) ;  !! set CO2 price in t2 to value in the "continuous path"
    
	      elseif ( o45_totCO2emi_peakBudgYr(iteration) < (0.1 + o45_change_totCO2emi_peakBudgYr(iteration) ) ), 
          display "New intermediate price in timestep after cm_peakBudgYr is sufficient to stabilize peaking year - go back to normal loop";	
	        o45_delay_increase_peakBudgYear(iteration+1) = 0;  !! probably is not necessary
          o45_reached_until2150pricepath(iteration) = 0;
	        o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration);  

        else      !! either didn't reach the continued "until2150"-price path in last iteration, or the increase was high enough to get emissions to 0. 
	           !! in this case, keep PeakBudgYr, and adjust the price in the year after the peakBudgYr to get emissions close to 0,
	        o45_delay_increase_peakBudgYear(iteration+1) = 1; !! make sure next iteration peakBudgYr is not shifted right again
	        o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration);
          p45_taxCO2eq_anchor(t2) = max(p45_taxCO2eq_anchor(ttot), !! at least as high as the price in the peakBudgYr
                                     p45_taxCO2eq_anchor(t2) * (o45_factorRescale_taxCO2_afterPeakBudgYr(iteration) / p45_factorRescale_taxCO2_Funneled(iteration) ) !! the full path was already rescaled by p45_factorRescale_taxCO2_Funneled, so adjust the second rescaling
                                    );
          loop(regi,                   !! this loop is necessary to allow the <-comparison in the next if statement
            if( p45_taxCO2eq_anchor_until2150(t2) < p45_taxCO2eq_anchor(t2) ,   !! check if new price would be higher than the price if the peakBudgYr would be one timestep later 
              display "price increase reached price from path with cm_peakBudgYr one timestep later - downscale to 99%"; 
		        p45_taxCO2eq_anchor(t2) = 0.99 * p45_taxCO2eq_anchor_until2150(t2); !! reduce the new CO2 price to 99% of the price that it would be if the peaking year was one timestep later. The next iteration will show if this is enough, otherwise cm_peakBudgYr will be shifted right 
              o45_reached_until2150pricepath(iteration) = 1;             !! upward CO2 price correction reached the continued price path - check in next iteration if this is high enough.  
            );
          );
        );
       
        display o45_factorRescale_taxCO2_afterPeakBudgYr;
	      p45_taxCO2eq_anchor(t)$(t.val gt t2.val) = p45_taxCO2eq_anchor(t2) + (t.val - t2.val) * cm_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;  !! increase by cm_taxCO2_IncAfterPeakBudgYr per year
      ); !! loop t2$(t2.val eq pm_ttot_val(ttot+1)),  !! set t2 to the following time step
    );  !! loop ttot$(ttot.val eq cm_peakBudgYr),  !! set ttot to the current peakBudgYr 
    cm_peakBudgYr = o45_peakBudgYr_Itr(iteration+1);  !! this has to happen outside the loop, otherwise the loop condition might be true twice
  ); !! if o45_delay_increase_peakBudgYear(iteration) = 1,   !! if there was a flip-floping in the previous iterations, try to solve this
  display p45_taxCO2eq_anchor, p45_taxCO2eq_anchor_until2150, o45_delay_increase_peakBudgYear, o45_reached_until2150pricepath, o45_peakBudgYr_Itr, o45_pkBudgYr_flipflop, cm_peakBudgYr;
);   !! if cm_iterative_target_adj eq 9 OR (cm_iterative_target_adj eq 10) OR (cm_iterative_target_adj eq 11),

 if((cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 1), 
    p45_taxCO2_IncAfterPeakBudgYr_iter(iteration) = cm_taxCO2_IncAfterPeakBudgYr; !! Save slope from last iteration
    
    if(iteration.val gt 8, !! only change slope after iteration 12
      !! 1) if deviation within tolerance -> no rescaling needed, keep the 2100 Tax level from the previous iteration
       if (abs(p45_globalBudget2100_absDev_iter(iteration)) le cm_budget2100CO2_absDevTol,
         p45_new2100Value(iteration) =  sum(t2$(t2.val eq 2100), p45_taxCO2eq_anchor_iter(iteration,t2)); !! keep the 2100 value the same

      !! 2) if rescaling needed:
       else
          !! 2.a if there is information on the slope:
          if(s45_TaxBudget_ChangeSlopeBest <0 ,
          p45_new2100Value(iteration) =  - sm_globalBudget2100_absDev * s45_TaxBudget_ChangeSlopeBest  !! the needed change
                                              + sum(t2$(t2.val eq 2100), p45_taxCO2eq_anchor_iter(iteration,t2)); !! the tax from this iteration
          
          !! 2.b if there is no information on the slope:
          else 
            s45_postPeakRescalingFactor = s45_actualbudgetco2_2100 / cm_addbudgetCO2from2020to2100;
            p45_postPeakRescalingFactor_iter(iteration) = s45_postPeakRescalingFactor;
            p45_postPeakRescalingFactor_Funneled(iteration) = 
                  max(min( 2 * EXP( -0.15 * iteration.val ) + 1.01 , p45_postPeakRescalingFactor_iter(iteration)),
                  1/ ( 2 * EXP( -0.15 * iteration.val ) + 1.01)
                  ); 
            p45_new2100Value(iteration) =  p45_postPeakRescalingFactor_Funneled(iteration) * 
                                          sum(t2$(t2.val eq 2100), p45_taxCO2eq_anchor_iter(iteration,t2));
            );
        ); 
      cm_taxCO2_IncAfterPeakBudgYr = (p45_new2100Value(iteration) - sum(t2$(t2.val eq cm_peakBudgYr), p45_taxCO2eq_anchor(t2))) 
                                        / (2100 - cm_peakBudgYr) 
                                        / sm_DptCO2_2_TDpGtC;
      );
    
    p45_taxCO2eq_anchor(t)$(t.val gt cm_peakBudgYr) = sum(t2$(t2.val eq cm_peakBudgYr), p45_taxCO2eq_anchor_until2150(t2)) !! CO2 tax in peak budget year
                                                  + (t.val - cm_peakBudgYr) * cm_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;

    !! Always set carbon price constant after 2100 to prevent huge taxes after 2100 and the resulting convergence problems
    p45_taxCO2eq_anchor(t)$(t.val gt 2100) = p45_taxCO2eq_anchor("2100");

    !! Compute difference for debugging
    pm_taxCO2eq_anchor_iterationdiff(t) = p45_taxCO2eq_anchor(t) - p45_taxCO2eq_anchor_iter(iteration,t);
    o45_taxCO2eq_anchor_iterDiff_Itr(iteration) = pm_taxCO2eq_anchor_iterationdiff("2100");

    display p45_taxCO2eq_anchor, pm_taxCO2eq_anchor_iterationdiff, o45_taxCO2eq_anchor_iterDiff_Itr;
  ); !! if (cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 1)

if ((cm_iterative_target_adj eq 10) AND (cm_peakBudgYr lt 2100),
    !! Initialize CDR price anchor trajectory to CO2 tax anchor trajectory for
    !! all years.
    p45_taxCDR_anchor(t) = p45_taxCO2eq_anchor(t);

    !! Update post-peak CO2 and CDR prices
    loop(ttot$(ttot.val eq cm_peakBudgYr), !! set ttot to peak year
      loop(t2$(t2.val eq pm_ttot_val(ttot + 1)), !! set t2 to the next timestep after the peak year
        !! If the CDR price slope is negative, CDR and CO2 prices will differ
        !! startinging in peak-year+2: CO2 prices will be kept constant and CDR
        !! prices will decrease.
        if(s45_taxCDR_IncAfterPeakBudgYr le 0,
          loop(t$(t.val gt t2.val),
            !! Set CO2 price constant for all years starting in peak-year+2. In
            !! the time step directly after the peak year the CO2 price must
            !! not be changed, as it may have been increased to avoid
            !! flip-flopping of the peak year. Since
            !! `cm_taxCO2_IncAfterPeakBudgYr` was set to 0 in the
            !! `cm_iterative_target_adj eq 10` case, the CO2 price in
            !! peak-year+1 should be equal to the one in peak-year if it was
            !! not explicitly increased, so it is ok to start in peak-year+2 in
            !! any case.
            p45_taxCO2eq_anchor(t) = p45_taxCO2eq_anchor(t2); 

            !! Linearly decrease CDR prices starting in peak-year+2.
            p45_taxCDR_anchor(t)$(t.val gt t2.val) = max(
              cm_minimumCDRtaxAfterPeak * sm_DptCO2_2_TDpGtC, !! Make sure that the trajectory does not become too small.
              p45_taxCO2eq_anchor(t2) + (t.val - t2.val) * s45_taxCDR_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC 
            );
          );

        !! If the CDR price slope is positive, CDR and CO2 price trajectories
        !! will be equal.
        else
          !! If the CO2 price was increased in the time step directly after the
          !! peak year, only apply the updated slope starting in peak-year+2 to
          !! not overwrite the value in peak-year+1.
          if(p45_taxCO2eq_anchor(t2) gt p45_taxCO2eq_anchor(ttot),
            loop(t$(t.val gt t2.val),
              p45_taxCO2eq_anchor(t) = p45_taxCO2eq_anchor(t2) + (t.val - t2.val) * s45_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;
            );

          !! If the CO2 price was not increased in the time step directly after
          !! the peak year, directly apply the updated slope from peak-year+1
          !! on.
          else
            loop(t$(t.val gt ttot.val),
              p45_taxCO2eq_anchor(t) = p45_taxCO2eq_anchor(ttot) + (t.val - ttot.val) * s45_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;
            );
          );

          !! Set CDR price trajectory exactly to the CO2 price trajectory.
          p45_taxCDR_anchor(t) = p45_taxCO2eq_anchor(t);
        ); !! if(s45_taxCDR_IncAfterPeakBudgYr le 0,        
      ); !! loop(t2$(t2.val eq pm_ttot_val(ttot + 1)),
    ); !! loop(ttot$(ttot.val eq cm_peakBudgYr),

    !! Always set CO2 and CDR price constant after 2100 to prevent huge taxes
    !! after 2100 and the resulting convergence problems
    p45_taxCDR_anchor(t)$(t.val gt 2100)   = p45_taxCDR_anchor("2100");
    p45_taxCO2eq_anchor(t)$(t.val gt 2100) = p45_taxCO2eq_anchor("2100");

    display p45_taxCDR_anchor;
  ); !! if(cm_iterative_target_adj eq 10,

  display p45_taxCO2eq_anchor, p45_taxCO2eq_anchor_until2150, o45_delay_increase_peakBudgYear, o45_reached_until2150pricepath, o45_peakBudgYr_Itr, o45_pkBudgYr_flipflop, cm_peakBudgYr;

*** --------ALGORITHM for cm_iterative_target_adj eq 7 ----------------------------------------------------------------------------------------
*** Algorithm for ENGAGE peakBudg formulation that results in a peak budget with zero net CO2 emissions afterwards
if(cm_iterative_target_adj eq 7,
	  if(o_modelstat eq 2 AND ord(iteration)<cm_iteration_max AND s45_actualbudgetco2 > 0 AND abs(cm_budgetCO2from2020 - s45_actualbudgetco2) ge 0.5,   !!only for optimal iterations, and not after the last one, and only if budget still possitive, and only if target not yet reached
		display p45_taxCO2eq_anchor;		
*** make sure that iteration converges: 
*** use multiplicative for budgets higher than 1600 Gt; for lower budgets, use multiplicative adjustment only for first 3 iterations, 
			if(ord(iteration) lt 3 or cm_budgetCO2from2020 > 1600,
			    !! change in CO2 price through adjustment: new price - old price; needed for adjustment option 2
				pm_taxCO2eq_anchor_iterationdiff(t) = p45_taxCO2eq_anchor(t) * min(max((s45_actualbudgetco2/cm_budgetCO2from2020)** (25/(2 * iteration.val + 23)),0.5+iteration.val/208),2 - iteration.val/102)  - p45_taxCO2eq_anchor(t);
				p45_taxCO2eq_anchor(t)$(t.val le cm_peakBudgYr) = p45_taxCO2eq_anchor(t) + pm_taxCO2eq_anchor_iterationdiff(t) ;
				p45_taxCO2eq_anchor_until2150(t) = p45_taxCO2eq_anchor_until2150(t) + pm_taxCO2eq_anchor_iterationdiff(t) ;
*** then switch to triangle-approximation based on last two iteration data points			
			else
			    !! change in CO2 price through adjustment: new price - old price; the two instances of "p45_taxCO2eq_anchor" cancel out -> only the difference term
				!! until cm_peakBudgYr: expolinear price trajectory
				p45_taxCO2eq_anchor_iterationdiff_tmp(t) = 
				                      max(pm_taxCO2eq_anchor_iterationdiff(t) * min(max((cm_budgetCO2from2020 - s45_actualbudgetco2)/(s45_actualbudgetco2 - s45_actualbudgetco2_last),-2),2),-p45_taxCO2eq_anchor(t)/2);
				p45_taxCO2eq_anchor(t)$(t.val le cm_peakBudgYr) = p45_taxCO2eq_anchor(t) + 
				                      max(pm_taxCO2eq_anchor_iterationdiff(t) * min(max((cm_budgetCO2from2020 - s45_actualbudgetco2)/(s45_actualbudgetco2 - s45_actualbudgetco2_last),-2),2),-p45_taxCO2eq_anchor(t)/2);
			  p45_taxCO2eq_anchor_until2150(t) = p45_taxCO2eq_anchor_until2150(t) + 
				                      max(pm_taxCO2eq_anchor_iterationdiff(t) * min(max((cm_budgetCO2from2020 - s45_actualbudgetco2)/(s45_actualbudgetco2 - s45_actualbudgetco2_last),-2),2),-p45_taxCO2eq_anchor_until2150(t)/2);
				pm_taxCO2eq_anchor_iterationdiff(t) = p45_taxCO2eq_anchor_iterationdiff_tmp(t);
				!! after cm_peakBudgYr: adjustment so that emissions become zero: increase/decrease tax in each time step after cm_peakBudgYr by percentage of that year's total CO2 emissions of 2015 emissions
			);
      o45_taxCO2eq_anchor_iterDiff_Itr(iteration) = pm_taxCO2eq_anchor_iterationdiff("2100");
      display o45_taxCO2eq_anchor_iterDiff_Itr;
		else
			if(s45_actualbudgetco2 > 0 or abs(cm_budgetCO2from2020 - s45_actualbudgetco2) < 2, !! if model was not optimal, or if budget already reached, keep tax constant
			p45_taxCO2eq_anchor(t) = p45_taxCO2eq_anchor(t);
			else
*** if budget has turned negative, reduce CO2 price by 20%
			p45_taxCO2eq_anchor(t) = 0.8*p45_taxCO2eq_anchor(t);
			p45_taxCO2eq_anchor_until2150(t) = 0.8*p45_taxCO2eq_anchor_until2150(t);
			);	
		);
*** after cm_peakBudgYr: always adjust to bring emissions close to zero
		p45_taxCO2eq_anchor(t)$(t.val gt cm_peakBudgYr) = p45_taxCO2eq_anchor(t) + p45_taxCO2eq_anchor(t)*max(sum(regi2,vm_emiAll.l(t,regi2,"co2"))/sum(regi2,vm_emiAll.l("2015",regi2,"co2")),-0.75);

*** check if cm_peakBudgYr is correct: if global emissions already negative, move cm_peakBudgYr forward
*** similar code block as used in iterative-adjust 9 below (credit to RP)
    o45_diff_to_Budg(iteration) = (cm_budgetCO2from2020 - s45_actualbudgetco2);
    o45_totCO2emi_peakBudgYr(iteration) = sum(t$(t.val = cm_peakBudgYr), sum(regi2, vm_emiAll.l(t,regi2,"co2")) );
    o45_totCO2emi_allYrs(t,iteration) = sum(regi2, vm_emiAll.l(t,regi2,"co2") );
    o45_change_totCO2emi_peakBudgYr(iteration) = sum(ttot$(ttot.val = cm_peakBudgYr), (o45_totCO2emi_allYrs(ttot-1,iteration) - o45_totCO2emi_allYrs(ttot+1,iteration) )/4 );  !! Only gives a tolerance range, exact value not important. Division by 4 somewhat arbitrary - could be 3 or 5 as well. 

    display cm_peakBudgYr, o45_diff_to_Budg, o45_peakBudgYr_Itr, o45_totCO2emi_allYrs, o45_totCO2emi_peakBudgYr, o45_change_totCO2emi_peakBudgYr;

***if( sum(t,sum(regi2,vm_emiAll.l(t,regi2,"co2")$(t.val = cm_peakBudgYr))) < -0.1,
*** cm_peakBudgYr = tt.val(t - 1)$(t.val = cm_peakBudgYr);
***);		

    if( abs(o45_diff_to_Budg(iteration)) < 20,                      !! only think about shifting peakBudgYr if the budget is close enough to target budget
      display "close enough to target budget to check timing of peak year";
      loop(ttot$(ttot.val = cm_peakBudgYr),                               !! look at the peak timing
***        if(  ( (o45_totCO2emi_peakBudgYr(iteration) < -(0.1 + o45_change_totCO2emi_peakBudgYr(iteration)) ) AND (cm_peakBudgYr > 2040) ), !! no peaking time before 2040
        if(  ( (o45_totCO2emi_peakBudgYr(iteration) < -(0.1) ) AND (cm_peakBudgYr > 2040) ), !! no peaking time before 2040
        display "shift peakBudgYr left";
		  o45_peakBudgYr_Itr(iteration+1) =  pm_ttot_val(ttot - 1);                
***          p45_taxCO2eq_anchor(t)$(t.val gt pm_ttot_val(ttot - 1)) = p45_taxCO2eq_anchor_until2150(ttot-1) + (t.val - pm_ttot_val(ttot - 1)) * cm_taxCO2_IncAfterPeakBudgYr * sm_DptCO2_2_TDpGtC;  !! increase by cm_taxCO2_IncAfterPeakBudgYr per year after peakBudgYr
*** if tax after cm_peakBudgYr is higher than normal increase rate (exceeding a 20% tolerance): shift right
		elseif( ( sum(regi, sum(t2$(t2.val = pm_ttot_val(ttot+1)),p45_taxCO2eq_anchor(t2))) > sum(regi,sum(t2$(t2.val = pm_ttot_val(ttot+1)),p45_taxCO2eq_anchor_until2150(t2)))*1.2 ) AND (cm_peakBudgYr < 2100) ), !! if peaking time would be after 2100, keep 2100 budget year
          if(  (iteration.val > 2) AND ( o45_peakBudgYr_Itr(iteration - 1) > o45_peakBudgYr_Itr(iteration) ) AND ( o45_peakBudgYr_Itr(iteration - 2) = o45_peakBudgYr_Itr(iteration) ) , !! if the target year was just shifted left after being shifted right
            o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration); !! don't shift right again immediately
          else
		    display "shift peakBudgYr right";
            o45_peakBudgYr_Itr(iteration+1) =  pm_ttot_val(ttot + 1);  !! ttot+1 is the new peakBudgYr
			loop(t$(t.val ge pm_ttot_val(ttot + 1)),
              p45_taxCO2eq_anchor(t) = p45_taxCO2eq_anchor_until2150(t);
            );
		  );
        
		else   !! don't do anything if the peakBudgYr is already at the corner values (2040, 2100) or if the emissions in the peakBudgYr are close to 0
          o45_peakBudgYr_Itr(iteration+1) = o45_peakBudgYr_Itr(iteration)
        );
      );
      cm_peakBudgYr = o45_peakBudgYr_Itr(iteration+1);
      display cm_peakBudgYr;
    );
*** If functionalForm is linear, re-adjust global anchor trajectory to go through the point (cm_taxCO2_historicalYr, cm_taxCO2_historical) 
$ifThen.taxCO2functionalForm3 "%cm_taxCO2_functionalForm%" == "linear"
    p45_taxCO2eq_anchor_until2150(t) = s45_taxCO2_historical 
        + (sum(t2$(t2.val eq cm_peakBudgYr), p45_taxCO2eq_anchor_until2150(t2)) - s45_taxCO2_historical) / (cm_peakBudgYr - s45_taxCO2_historicalYr) !! Yearly increase of CO2 price that interpolates between cm_taxCO2_historical in cm_taxCO2_historicalYr and p45_taxCO2eq_anchor_until2150 in peak year
                                      * (t.val - s45_taxCO2_historicalYr) ;
    p45_taxCO2eq_anchor(t)$(t.val le cm_peakBudgYr) = p45_taxCO2eq_anchor_until2150(t);
    p45_taxCO2eq_anchor(t)$(t.val gt 2100) = p45_taxCO2eq_anchor("2100");
***TODO: CHECK IF ALGORITHM DOES WHAT IS EXPECTED. CURRENTLY NO RE-ADJUSTMENT OF GLOBAL ANCHOR TRAJECTORY BETWEEN PEAK YEAR AND 2100 AS NO SUCH ADJUSTMENT WAS CONTAINED IN ORIGINAL ALGORITHM
$endIf.taxCO2functionalForm3

display p45_taxCO2eq_anchor_until2150, p45_taxCO2eq_anchor;
); !! if(cm_iterative_target_adj eq 7,

*** Save s45_actualbudgetco2 for having it available in next iteration:
s45_actualbudgetco2_last = s45_actualbudgetco2;


***----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*** Part III (Regional differentiation): Re-compute p45_regiDiff_ratio, and 
***                                      re-create regional carbon price trajectories p45_taxCO2eq_regiDiff using p45_taxCO2eq_anchor (updated in parts I-II above) and p45_regiDiff_ratio
*** Re-create regional CO2 and potentially
*** CDR price trajectories p45_taxCO2eq_regiDiff using p45_taxCO2eq_anchor /
*** p45_taxCDR_anchor (updated in parts I-II above) and p45_regiDiff_convFactor
*** (computed in datainput)
***----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

*** Step III.1: Re-compute p45_regiDiff_initialRatio. This is necessary if p45_taxCO2eq_anchor was adjusted.

$ifThen.taxCO2regiDiff3 "%cm_taxCO2_regiDiff_startyearValue%" == "endogenous"
if( (cm_taxCO2_regiDiff = 0) or (cm_taxCO2_regiDiff = 3), !! none or gdpSpread
  !! Both parameters are not needed. No need to re-compute
elseif (cm_taxCO2_regiDiff = 1) or (cm_taxCO2_regiDiff = 2), !! initialSpread10 or initialSpread20
  !! Nothing to re-compute as p45_regiDiff_initialRatio does not depend on p45_taxCO2eq_anchor
elseif (cm_taxCO2_regiDiff = 5) or (cm_taxCO2_regiDiff = 6) or (cm_taxCO2_regiDiff = 7) or (cm_taxCO2_regiDiff = 8) or (cm_taxCO2_regiDiff = 10), !! ScenarioMIP2035, ScenarioMIP2050, ScenarioMIP2070, ScenarioMIP2100, or manual
  !! Re-compute p45_regiDiff_initialRatio based on regional carbon prices in p45_regiDiff_startYr
  p45_regiDiff_initialRatio(regi) = sum(ttot$(ttot.val eq p45_regiDiff_startYr(regi)), p45_taxCO2eq_path_gdx_ref(ttot,regi) / p45_taxCO2eq_anchor(ttot));
else
  abort "please choose a valid scenario via cm_taxCO2_regiDiff or set cm_taxCO2_regiDiff to manual"
);
$else.taxCO2regiDiff3
if( (cm_taxCO2_regiDiff = 0) or (cm_taxCO2_regiDiff = 1) or (cm_taxCO2_regiDiff = 2) or (cm_taxCO2_regiDiff = 3), !! none, initialSpread10, initialSpread20, or gdpSpread
  abort "Regional carbon prices can only be set manually via cm_taxCO2_regiDiff_startyearValue if cm_taxCO2_regiDiff equals (ScenarioMIP2035), (ScenarioMIP2050), (ScenarioMIP2070), (ScenarioMIP2100), or (manual)."
else
  !! Re-compute p45_regiDiff_initialRatio
  p45_regiDiff_initialRatio(regi) = sum(ttot$(ttot.val eq cm_startyear), p45_regiDiff_startyearValue(regi) / p45_taxCO2eq_anchor(ttot));
);
$endIf.taxCO2regiDiff3
display  p45_regiDiff_initialRatio;

*** Step III.3: Create ratio between regional carbon price and global anchor trajectory based on previously defined convergence

if( (cm_taxCO2_regiDiff = 0) or (cm_taxCO2_regiDiff = 3), !! none or gdpSpread
  !! Nothing to re-compute
else 
  !! Set convergence factor equal to p45_regiDiff_initialRatio before p45_regiDiff_startYr:
  p45_regiDiff_ratio(t,regi)$(t.val lt p45_regiDiff_startYr(regi)) = p45_regiDiff_initialRatio(regi);
  !! Set  convergence factor equal to 1 from p45_regiDiff_endYr:
  p45_regiDiff_ratio(t,regi)$(t.val ge p45_regiDiff_endYr(regi)) = 1;
  !! Create convergence between p45_regiDiff_startYr and p45_regiDiff_endYr:
  loop((t,regi)$((t.val ge p45_regiDiff_startYr(regi)) and (t.val lt p45_regiDiff_endYr(regi))),
    p45_regiDiff_ratio(t,regi) = p45_regiDiff_initialRatio(regi) 
                                + (1 - p45_regiDiff_initialRatio(regi)) * rPower( (t.val - p45_regiDiff_startYr(regi)) / (p45_regiDiff_endYr(regi) - p45_regiDiff_startYr(regi)), p45_regiDiff_exponent(regi));
  );
);
display p45_regiDiff_ratio;

*** Step III.4: Create regionally differentiated carbon price trajectories based on global anchor trajectory and p45_regiDiff_ratio

p45_taxCO2eq_regiDiff(t,regi) = p45_regiDiff_ratio(t,regi) * p45_taxCO2eq_anchor(t);
display p45_taxCO2eq_regiDiff;

if (cm_iterative_target_adj eq 10,
  p45_taxCDR_regiDiff(t,regi)   = p45_regiDiff_ratio(t,regi) * p45_taxCDR_anchor(t);
  display p45_taxCDR_regiDiff;
);

*** Step III.5: If regional carbon prices in cm_startyear where set manually via cm_taxCO2_regiDiff_startyearValue, ensure that convergence to global anchor trajectory does not lead to lower regional carbon prices in some timesteps (this could happen if regional carbon price in cm_startyear is much higher than global anchor price)

$ifThen.taxCO2regiDiffStartyearValue2 "%cm_taxCO2_regiDiff_startyearValue%" == "endogenous"
$else.taxCO2regiDiffStartyearValue2
  p45_taxCO2eq_regiDiff(t,regi) = max(p45_taxCO2eq_regiDiff(t,regi), p45_regiDiff_startyearValue(regi));
  display "Apply p45_regiDiff_startyearValue(regi) as lower bound for p45_taxCO2eq_regiDiff"
  display p45_taxCO2eq_regiDiff;
$endIf.taxCO2regiDiffStartyearValue2

***------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*** Part IV (Interpolation from path_gdx_ref): Re-create interpolation based on p45_taxCO2eq_regiDiff (updated in part III above), and s45_interpolation_startYr and s45_interpolation_endYr (computed in datainput)
***-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

*** Step IV.2: Re-create interpolation
pm_taxCO2eq(ttot,regi) = p45_taxCO2eq_path_gdx_ref(ttot,regi); !! Initialize pm_taxCO2eq with p45_taxCO2eq_path_gdx_ref. Then overwrite all time steps after cm_startyear
pm_taxCO2eq(t,regi)$(t.val le s45_interpolation_startYr) = p45_taxCO2eq_regiDiff(t,regi);
!! there is no adjustment if s45_interpolation_startYr = 2025 and s45_interpolation_endYr = 2030
pm_taxCO2eq(t,regi)$((t.val gt s45_interpolation_startYr) and (t.val lt s45_interpolation_endYr)) =
    sum(ttot2$(ttot2.val eq s45_interpolation_startYr), p45_taxCO2eq_path_gdx_ref(ttot2,regi)) !! value of p45_taxCO2eq_path_gdx_ref in s45_interpolation_startYr
    * (s45_interpolation_endYr - t.val) / (s45_interpolation_endYr - s45_interpolation_startYr)
  + sum(t2$(t2.val eq s45_interpolation_endYr), p45_taxCO2eq_regiDiff(t2,regi)) !! value of p45_taxCO2eq_regiDiff in s45_interpolation_endYr
    * (t.val - s45_interpolation_startYr) / (s45_interpolation_endYr - s45_interpolation_startYr);
pm_taxCO2eq(t,regi)$(t.val ge s45_interpolation_endYr) = p45_taxCO2eq_regiDiff(t,regi);

display pm_taxCO2eq;

*** Set the CDR tax exactly to the CO2 tax for all years in and before the peak
*** year, to use the same interpolation as for the CO2 tax and set it to the
*** derived regionally differentiated CDR tax levels afterwards.
if (cm_iterative_target_adj eq 10,
  pm_taxCDR(t,regi)$(t.val le cm_peakBudgYr) = pm_taxCO2eq(t,regi);
  pm_taxCDR(t,regi)$(t.val gt cm_peakBudgYr) = p45_taxCDR_regiDiff(t,regi);
  display pm_taxCDR;
);


**** derive the new CDR tax
if(cm_iterative_target_adj eq 11,
  if((iteration.val ge 12),
  !! Adjust CDR prices if additional 2100 target budget was set 
  !! 4) apply the rescaling factor to the CDR price --> exact approach depends on the CDR price shape
      !! 4A)_____ constant price_____
      if(cm_CDRpriceShape eq 1, 
          p45_taxCDR_anchor(ttot) = max(cm_minimumCDRtaxAfterPeak * sm_DptCO2_2_TDpGtC, !! lower bound on CDR tax.
                                     p45_taxCDR_anchor(ttot) * p45_factorRescale_CDRtax_Funneled(iteration));
          !! Todo: what if the tax is at the minimum but should further decrease?
          ); !! constant CDR price 

      !! 4B)_____ linearly falling CDR price_____
      if(cm_CDRpriceShape eq 2, 
        !! Adjust the 2100 value, respecting the minimum value!
         p45_taxCDR_anchor("2100") = max(cm_minimumCDRtaxAfterPeak * sm_DptCO2_2_TDpGtC, !! lower bound on CDR tax.
                                    p45_taxCDR_anchor("2100") * p45_factorRescale_CDRtax_Funneled(iteration)); !! 
        
        !! If the intended 2100 value was below the minimum, adjust the starting value by the refactoring value derived above  
                !! Also need to respect max and min values there though! 
                !! => TODO: introduce a check about what to do if the extremes have been reached (repeatedly -> failure)
                !! => TODO: introduce a check if the tax in cm_startYear was adjusted
        if ((cm_minimumCDRtaxAfterPeak * sm_DptCO2_2_TDpGtC) eq p45_taxCDR_anchor("2100"), !! equal because it was adjusted just before 
             p45_taxCDR_anchor(ttot)$(ttot.val eq cm_startyear) = max(
                    min( 1.5 * cm_CDRstartYearTax * sm_DptCO2_2_TDpGtC,   !! a) a maximum value for cm_startYear
                      p45_taxCDR_anchor(ttot) * p45_factorRescale_CDRtax_Funneled(iteration) ),          !! TODO: calculate a new adjustment value based on the price sensitivity at cm_startyear + 1 
                    0.5 * cm_CDRstartYearTax);   !! b) a minimum value for cm_startYear
            );             
                
        !! calculate the slope from starting year to end year (one of them is new)
        s45_taxCDR_slope = (p45_taxCDR_anchor("2100") -  
                            sum(ttot2$(ttot2.val eq cm_startyear), p45_taxCDR_anchor(ttot2))) / 
                            (2100 - cm_startyear );
        
        !! get the tax in between the two anchor points; (first value is actually cm_startyear)        
        loop(ttot$((ttot.val gt cm_startyear) AND (ttot.val lt 2100)),
          p45_taxCDR_anchor(ttot) = sum(ttot2$(ttot2.val eq cm_startyear), p45_taxCDR_anchor(ttot2)) + 
                              (ttot.val - cm_startyear) * s45_taxCDR_slope;
            ); !! fill tax between edge years
      ); !! linearly falling CDR price


      !! 4D)_____ linearly falling CDR price to 0 in 2100_____
      if(cm_CDRpriceShape eq 4, 
        !! Adjust the start Year value; no limits for now
         loop(ttot2$(ttot2.val eq cm_startYear),
         p45_taxCDR_anchor(ttot2) = p45_taxCDR_anchor(ttot2) * p45_factorRescale_CDRtax_Funneled(iteration); !! 
        );                       
        !! calculate the slope from starting year to end year, start year was adjusted
        s45_taxCDR_slope = (p45_taxCDR_anchor("2100") -  
                            sum(ttot2$(ttot2.val eq cm_startyear), p45_taxCDR_anchor(ttot2))) / 
                            (2100 - cm_startyear);
        
        !! set the tax in between the two anchor points       
        loop(ttot$((ttot.val gt cm_startyear) AND (ttot.val lt 2100)),
          p45_taxCDR_anchor(ttot) = sum(ttot2$(ttot2.val eq cm_startyear), p45_taxCDR_anchor(ttot2)) + 
                              (ttot.val - cm_startyear) * s45_taxCDR_slope;
            ); !! fill tax between edge years
      ); !! linearly falling CDR price
      
      
      !! 4C)_____ exponential increase to constant_____
      !! The joint price is kept and continues to increase if more NNE are needed
      !! update CDR tax to new carbon price trajectory 
      if(cm_CDRpriceShape eq 3, 
       s45_maxCDRtax = p45_taxCDR_anchor("2100") * p45_factorRescale_CDRtax_Funneled(iteration); !! take the 2100 value as it will always be the maximum value, incl. in iteration 1
       p45_maxCDRtax_iter(iteration) = s45_maxCDRtax;
       loop(regi,
        loop(ttot,
            !! will set the max. CDR price for a region once it is smaller than the CO2 price
           p45_taxCDR_anchor(ttot) =  min(p45_taxCO2eq_anchor(ttot),  
                                              s45_maxCDRtax) ;
        ); !! ttot
      ); !! regi
      ); !! price increases with CO2 tax until max reached
      
      !! ____________________________
      !! always set CDR price constant after 2100
      p45_taxCDR_anchor(ttot)$(ttot.val gt 2100)   = p45_taxCDR_anchor("2100");

  !! 5) assign the anchor CDR tax as new tax for next iteration)
  pm_taxCDR(t,regi)$(t.val ge cm_startYear) = p45_taxCDR_anchor(t);

else !! earlier iterations: keep what was set in the first 15 iterations
  pm_taxCDR(t,regi)$(t.val ge cm_startYear) =  pm_taxCDR(t,regi);
  );
); !! end of (cm_iterative_target_adj eq 11)

************************************************************************************************************
if((cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 2),
  pm_frac_NetNegEmi(iteration) = cm_frac_NetNegEmi; !! Save the fraction used in this iteration
  !! calculate the new cm_frac_NetNegEmi, starting after 12 iterations
  if((iteration.val ge 12),
  !! if the derived net negative price is > than the previous price: keep these 2 constant. !! In this case: would need an endogenously calculated increase of the CP ex post instead!
  if(p45_neweffectiveNNEtax(iteration) gt p45_taxCO2eq_anchor("2100"),
    cm_frac_NetNegEmi = 0;
  elseif(p45_neweffectiveNNEtax(iteration) le 0),
    cm_frac_NetNegEmi = 1;
  else 
    cm_frac_NetNegEmi = (p45_taxCO2eq_anchor("2100") - p45_neweffectiveNNEtax(iteration)) / p45_taxCO2eq_anchor("2100");
    );
  );
); !! (cm_iterative_target_adj eq 9) AND (cm_postPeakCpAdj eq 2)

************************************************************************************************************

*** Step IV.3: Re-introduce lower bound pm_taxCO2eq by p45_taxCO2eq_path_gdx_ref if switch cm_taxCO2_lowerBound_path_gdx_ref is on
if(cm_taxCO2_lowerBound_path_gdx_ref = 1,
  pm_taxCO2eq(t,regi) = max(pm_taxCO2eq(t,regi), p45_taxCO2eq_path_gdx_ref(t,regi));
  display pm_taxCO2eq;
  !! For CDR this should only be applied if the CDR tax does NOT have a
  !! negative slope. In case the CDR tax does not have a negative slope, it
  !! should be, by design, exactly the same as the CO2 tax. In the case it does
  !! have a negative slope, we explicitly want to allow for the CDR tax to go
  !! below the reference level. That should actually be quite a corner case,
  !! though.
  if (cm_iterative_target_adj eq 10 AND s45_taxCDR_IncAfterPeakBudgYr ge 0,
    pm_taxCDR(t,regi) = max(pm_taxCO2eq(t,regi), p45_taxCO2eq_path_gdx_ref(t,regi));
    display pm_taxCDR;
  );
  ); !! lower bound  
  
); !! if((cm_emiscen eq 9) AND ((cm_iterative_target_adj eq 5) OR (cm_iterative_target_adj eq 7) OR (cm_iterative_target_adj ge 9)),
*** EOF ./modules/45_carbonprice/functionalForm/postsolve.gms
