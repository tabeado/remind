*** |  (C) 2006-2023 Potsdam Institute for Climate Impact Research (PIK)
*** |  authors, and contributors see CITATION.cff file. This file is part
*** |  of REMIND and licensed under AGPL-3.0-or-later. Under Section 7 of
*** |  AGPL-3.0, you are granted additional permissions described in the
*** |  REMIND License Exception, version 1.0 (see LICENSE file).
*** |  Contact: remind@pik-potsdam.de
*** SOF ./modules/01_macro/singleSectorGr/equations.gms


***---------------------------------------------------------------------------
*' Usable macroeconomic output - net of climate change damages - is calculated from the macroeconomic output, 
*' taking into account export and import of the final good, taking specific trade costs into account, 
*' which are assigned to the importer. The resulting output is used for consumption, 
*' for investments into the capital stock, and for the energy system cost components investments,
*' fuel costs and operation & maintenance. 
*' Other additional costs like non-energy related greenhouse gas abatement costs and
*' agricultural costs, which are delivered by the land use model MAgPIE, are deduced from disposable output. 
*' Net tax revenues and adjustment costs converge to zero in the optimal solution (equilibrium point).
***---------------------------------------------------------------------------
qm_budget(ttot,regi)$( ttot.val ge cm_startyear ) .. 
    vm_cesIO(ttot,regi,"inco") * vm_damageFactor(ttot,regi) 
  - vm_Xport(ttot,regi,"good") 
  + vm_Mport(ttot,regi,"good") * (1 - pm_tradecostgood(regi) - pm_risk_premium(regi))
  + vm_revenueFromSpecificGoods(ttot, regi)
  =g=
    vm_cons(ttot,regi)
  + sum(ppfKap(in), vm_invMacro(ttot,regi,in))
  + sum(ppfKap(in), v01_invMacroAdj(ttot,regi,in))  
  + sum(in, vm_invRD(ttot,regi,in))
  + sum(in, vm_invInno(ttot,regi,in))
  + sum(in, vm_invImi(ttot,regi,in))
  + sum(tradePe(enty)$(NOT tradeCap(enty)), pm_costsTradePeFinancial(regi,"Mport",enty) * vm_Mport(ttot,regi,enty))
  + sum(tradePe(enty)$(NOT tradeCap(enty)),
      (pm_costsTradePeFinancial(regi,"Xport",enty) * vm_Xport(ttot,regi,enty))
    * ( 1 
      + ( pm_costsTradePeFinancial(regi,"XportElasticity",enty)
        / sqr(pm_ttot_val(ttot)-pm_ttot_val(ttot-1))
        * ( vm_Xport(ttot,regi,enty) 
          / ( vm_Xport(ttot-1,regi,enty) + pm_costsTradePeFinancial(regi, "tradeFloor",enty) ) 
          - 1
          )
        )$( ttot.val ge max(2010, cm_startyear) )
      )
    )
  + sum(tradeSe, pm_MPortsPrice(ttot,regi,tradeSe) * vm_Mport(ttot,regi,tradeSe))
  - sum(tradeSe, pm_XPortsPrice(ttot,regi,tradeSe) * vm_Xport(ttot,regi,tradeSe))
  + sum(tradeCap, vm_costTradeCap(ttot,regi,tradeCap))
  + vm_taxrev(ttot,regi)$(ttot.val ge 2010)
  + vm_costAdjNash(ttot,regi)
  + sum(in_enerSerAdj(in), vm_enerSerAdj(ttot,regi,in))
  + sum(teEs, vm_esCapInv(ttot,regi,teEs))
  + vm_costpollution(ttot,regi)
  + pm_totLUcosts(ttot,regi)
*** agricultural MACs are part of pm_totLUcosts (see module 26_agCosts)
  + sum(enty$(emiMacSector(enty) AND (NOT emiMacMagpie(enty))), pm_macCost(ttot,regi,enty))  
  + vm_costEnergySys(ttot,regi)
;



***---------------------------------------------------------------------------
*' The labor available in every time step and every region comes from exogenous data. 
*' It is the population corrected by the population age structure,
*' which results in the labour force of people agged 15 to 65. 
*' The labor participation rate is not factored into the labour supply (as it would only imply a
*' rescaling of parameters without consequences for the model's dynamic). 
*' The labour market balance equation reads as follows:
***---------------------------------------------------------------------------
q01_balLab(t,regi)..
    vm_cesIO(t,regi,"lab") 
  =e= 
    pm_lab(t,regi)
;

***---------------------------------------------------------------------------
*' The production function is a nested CES (constant elasticity of substitution) production function. 
*' The macroeconomic output is generated by the inputs capital, labor, and total final energy (as a macro-ecoomic
*' aggregate in $US units). The generation of total final energy is described
*' by a CES production function as well, whose input factors are CES function outputs again. 
*' Hence, the outputs of CES nests are intermediates measured in $US units. 
*' According to the Euler-equation the value of the intermediate equals the sum of expenditures for the inputs. 
*' Sector-specific final energy types represent the bottom end of the `CES-tree'. These 'CES leaves' are
*' measured in physical units and have a price in $US per physical unit. 
*' The top of the tree is the total economic output measured in $US.
*' The following equation is the generic form of the production function. 
*' It treats the various CES nests separately and the nests are inter-connetected via mappings. 
*' This equation calculates the amount of intermediate output in a time-step and region 
*' from the associated factor input amounts according to:
*** Keep in mind to adjust the calculation of derivatives and shares 
*** in ./core/reswrite.inc if you change the structure of this function.
***---------------------------------------------------------------------------
q01_cesIO(t,regi,ipf(out))$( NOT ipf_putty(out) ) ..
  vm_cesIO(t,regi,out)
  =e=
  !! use exp(log(a) * b) = a ** b because the latter is not accurate in GAMS for
  !! very low values of a
  exp(
    log(
      sum(cesOut2cesIn(out,in),
        pm_cesdata(t,regi,in,"xi")
      * exp(
          log(
	    pm_cesdata(t,regi,in,"eff")
	  * vm_effGr(t,regi,in)
	  * vm_damageProdFactor(t,regi,in)
	  * vm_cesIO(t,regi,in)
	  )
	* pm_cesdata(t,regi,out,"rho")
	)
      )
    )
  * (1 / pm_cesdata(t,regi,out,"rho"))
  )
;

***---------------------------------------------------------------------------
*' Constraints for perfect complements in the CES tree
***---------------------------------------------------------------------------
q01_prodCompl(t,regi,in,in2) $ (complements_ref(in,in2) AND (( NOT in_putty(in2)) OR ppfIO_putty(in2))) ..
    vm_cesIO(t,regi,in) 
  =e= 
    pm_cesdata(t,regi,in2,"compl_coef")
  * vm_cesIO(t,regi,in2)
;



***---------------------------------------------------------------------------    
*' The capital stock is calculated recursively. Its amount in the previous time
*' step is devaluated by an annual depreciation factor and enlarged by investments. 
*' Both depreciation and investments are expressed as annual values,
*' so the time step length is taken into account.
***---------------------------------------------------------------------------
q01_kapMo(ttot,regi,ppfKap(in))$(
                             NOT in_putty(in)
                         AND ord(ttot) lt card(ttot)
                         AND pm_ttot_val(ttot+1) ge max(2010, cm_startyear)
                         AND pm_cesdata("2005",regi,in,"quantity") gt 0     ) ..
  vm_cesIO(ttot+1,regi,in)
  =e=
    vm_cesIO(ttot,regi,in)
  * (1 - pm_delta_kap(regi,in))
 ** (pm_ttot_val(ttot+1) - pm_ttot_val(ttot))
  + pm_cumDeprecFactor_old(ttot+1,regi,in) * vm_invMacro(ttot,regi,in)
  + pm_cumDeprecFactor_new(ttot+1,regi,in) * vm_invMacro(ttot+1,regi,in)
;

***---------------------------------------------------------------------------
*' Adjustment costs of macro economic investments:
***---------------------------------------------------------------------------
q01_invMacroAdj(ttot,regi,ppfKap(in))$( ttot.val ge max(2010, cm_startyear))..
    v01_invMacroAdj(ttot,regi,in)
  =e= 
    sqr( (vm_invMacro(ttot,regi,in)-vm_invMacro(ttot-1,regi,in)) 
      / (pm_ttot_val(ttot)-pm_ttot_val(ttot-1)) 
      / (vm_invMacro(ttot,regi,in)+0.0001)
    )
  * vm_cesIO(ttot,regi,in) / 11
*ML/RP* use "kap/11"  instead of "vm_invMacro" for the scaling to remove the "invest=0"-trap that sometimes appeared in delay scenarios; kap/11 corresponds to the global average ratio of investments to capital in 2005.
*** In some regions the ratio kap:invest is higher, in some it is lower.
;

***---------------------------------------------------------------------------
*' Initial conditions for capital:
***---------------------------------------------------------------------------
q01_kapMo0(t0(t),regi,ppfKap(in))$(pm_cesdata(t,regi,in,"quantity") gt 0)..
    vm_cesIO(t,regi,in) 
  =e= 
    pm_cesdata(t,regi,in,"quantity");

*' Limit the share of one ppfEn in total CES nest inputs:
q01_limitShPpfen(t,regi,out,in)$( pm_ppfen_shares(t,regi,out,in) ) ..
    vm_cesIO(t,regi,in) + pm_cesdata(t,regi,in,"offset_quantity")
  =l=
    pm_ppfen_shares(t,regi,out,in)
  * (sum(cesOut2cesIn(out,in2), vm_cesIO(t,regi,in2) + pm_cesdata(t,regi,in2,"offset_quantity")))
;

*' Limit the ratio of two ppfEn:
q01_limtRatioPpfen(t,regi,in,in2)$( p01_ppfen_ratios(t,regi,in,in2) ) ..
    vm_cesIO(t,regi,in) + pm_cesdata(t,regi,in,"offset_quantity")
  =l=
    p01_ppfen_ratios(t,regi,in,in2)
  * (vm_cesIO(t,regi,in2) + pm_cesdata(t,regi,in,"offset_quantity"))
;


***---------------------------------------------------------------------------                                
*** Start of Putty-Clay equations 
*' Putty-Clay production function:
***---------------------------------------------------------------------------
q01_cesIO_puttyclay(t,regi,ipf_putty(out)) ..
  vm_cesIOdelta(t,regi,out)
  =e=
    sum(cesOut2cesIn(out,in),
      pm_cesdata(t,regi,in,"xi")
    * ( 
        pm_cesdata(t,regi,in,"eff")
      * vm_effGr(t,regi,in)
      * vm_cesIOdelta(t,regi,in)
      )
   ** pm_cesdata(t,regi,out,"rho")
    )
 ** (1 / pm_cesdata(t,regi,out,"rho"))
;

*' Putty-Clay constraints for perfect complements in the CES tree:
q01_prodCompl_putty(t,regi,in,in2) $ (complements_ref(in,in2)
                                 AND ( in_putty(in2) AND  ( NOT ppfIO_putty(in2)))) ..
      vm_cesIOdelta(t,regi,in) =e=
                                pm_cesdata(t,regi,in2,"compl_coef")
                                * vm_cesIOdelta(t,regi,in2);

*' Correspondance between vm_cesIO and vm_cesIOdelta:
q01_puttyclay(ttot,regi,in_putty(in))$(ord(ttot) lt card(ttot)  AND (pm_ttot_val(ttot+1) ge max(2010, cm_startyear)))..
  vm_cesIO(ttot+1,regi,in)
  =e=
  vm_cesIO(ttot,regi,in)*(1- pm_delta_kap(regi,in))**(pm_ttot_val(ttot+1)-pm_ttot_val(ttot))
           +  pm_cumDeprecFactor_old(ttot+1,regi,in)* vm_cesIOdelta(ttot,regi,in)
           +  pm_cumDeprecFactor_new(ttot+1,regi,in)* vm_cesIOdelta(ttot+1,regi,in)
;

*' Capital motion equation for putty clay capital:
q01_kapMo_putty(ttot,regi,in_putty(in))$(ppfKap(in) AND (ord(ttot) le card(ttot)) AND (pm_ttot_val(ttot) ge max(2005, cm_startyear)) AND (pm_cesdata("2005",regi,in,"quantity") gt 0))..
    vm_cesIOdelta(ttot,regi,in)
    =e=
    vm_invMacro(ttot,regi,in)
;
***---------------------------------------------------------------------------
*** End of Putty-Clay equations
***---------------------------------------------------------------------------
*** EOF ./modules/01_macro/singleSectorGr/equations.gms
