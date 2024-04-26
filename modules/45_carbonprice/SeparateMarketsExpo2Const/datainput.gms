*** |  (C) 2006-2022 Potsdam Institute for Climate Impact Research (PIK)
*** |  authors, and contributors see CITATION.cff file. This file is part
*** |  of REMIND and licensed under AGPL-3.0-or-later. Under Section 7 of
*** |  AGPL-3.0, you are granted additional permissions described in the
*** |  REMIND License Exception, version 1.0 (see LICENSE file).
*** |  Contact: remind@pik-potsdam.de
*** SOF ./modules/45_carbonprice/SeparateMarketsExpo2Const/datainput.gms
***----------------------------
*** CO2 Tax level
***----------------------------

*** CO2 tax level is calculated at an initial 5% exponential increase from the 2020 tax level exogenously defined

*GL: tax path in 10^12$/GtC = 1000 $/tC
*** according to Asian Modeling Excercise tax case setup, 30$/t CO2eq in 2020 = 0.110 k$/tC

if(cm_co2_tax_2020 lt 0,
abort "please choose a valid cm_co2_tax_2020"
elseif cm_co2_tax_2020 ge 0,
*** convert tax value from $/t CO2eq to T$/GtC
pm_taxCO2eq("2025",regi)= cm_co2_tax_2020 * sm_DptCO2_2_TDpGtC;
);

pm_taxCO2eq(ttot,regi)$(ttot.val ge 2025 AND ttot.val le c_peakBudgYr) = pm_taxCO2eq("2025",regi)*cm_co2_tax_growth**(ttot.val-2025);
pm_taxCO2eq(ttot,regi)$(ttot.val gt c_peakBudgYr) =sum(t$(t.val eq c_peakBudgYr),pm_taxCO2eq(t,regi)); !! keep taxes constant after cm_peakBudgYr

pm_taxCDR(ttot,regi) = pm_taxCO2eq(ttot,regi)


display pm_taxCDR, pm_taxCO2eq;

*** EOF ./modules/45_carbonprice/SeparateMarketsExpo2Const/datainput.gms
