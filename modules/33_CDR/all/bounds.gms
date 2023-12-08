*** |  (C) 2006-2023 Potsdam Institute for Climate Impact Research (PIK)
*** |  authors, and contributors see CITATION.cff file. This file is part
*** |  of REMIND and licensed under AGPL-3.0-or-later. Under Section 7 of
*** |  AGPL-3.0, you are granted additional permissions described in the
*** |  REMIND License Exception, version 1.0 (see LICENSE file).
*** |  Contact: remind@pik-potsdam.de
*** SOF ./modules/33_CDR/all/bounds.gms
vm_emiCdr.fx(t,regi,emi)$(not sameas(emi,"co2")) = 0.0;
v33_grindrock_onfield_tot.up(t,regi,rlf,rlf2) = s33_step;
v33_grindrock_onfield_tot.fx("2005",regi,rlf,rlf2) = 0.0;
v33_grindrock_onfield.fx(t,regi,rlf,rlf2)$(rlf2.val gt 10) = 0;
v33_grindrock_onfield_tot.fx(t,regi,rlf,rlf2)$(rlf2.val gt 10) = 0;
v33_emiDAC.up(t,regi) = 0; !! DAC has never positive emissions
v33_emiEW.up(t,regi) = 0;   !!EW has never positive emissions
v33_emiEW.lo(t,regi)$(t.val le 2030) = -0.03; !! [GtC/a]: EW until 2030 is limited to 100Mt CO2 per region
vm_emiCdr.up(t,regi,"co2")$(t.val gt 2015) = -0.0001;
vm_emiCdr.up(t,regi,"co2")$(t.val le 2015) = 0; !! no total positive emissions from DAC + EW before 2015
if (cm_emiscen ne 1,
    vm_cap.lo(t,regi,"dac",rlf)$(teNoTransform2rlf_dyn33("dac",rlf) AND (t.val ge max(2025,cm_startyear))) = 1e-7;  
);
*** EOF ./modules/33_CDR/all/bounds.gms
