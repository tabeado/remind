*** |  (C) 2006-2023 Potsdam Institute for Climate Impact Research (PIK)
*** |  authors, and contributors see CITATION.cff file. This file is part
*** |  of REMIND and licensed under AGPL-3.0-or-later. Under Section 7 of
*** |  AGPL-3.0, you are granted additional permissions described in the
*** |  REMIND License Exception, version 1.0 (see LICENSE file).
*** |  Contact: remind@pik-potsdam.de
*** SOF ./modules/33_CDR/DAC/datainput.gms

!! Beutler et al. 2019 (Climeworks)
!!fe demand electricity for ventilation
p33_dac_fedem_el("feels") = 5.28;
!!fe demand heat for material recovery
p33_dac_fedem_heat(regi,"fehes") = 21.12;
p33_dac_fedem_heat(regi,"fegas") = 21.12;
p33_dac_fedem_heat(regi,"feh2s") = 21.12;
p33_dac_fedem_heat(regi,"feels") = 21.12;
!!p33_dac_fedem_heat("MEA","feels") = 13.82;
!!p33_dac_fedem_heat("SSA","feels") = 15.01;
!!p33_dac_fedem_heat("IND","feels") = 15.86;
!!p33_dac_fedem_heat("LAM","feels") = 16.70;
!!p33_dac_fedem_heat("OAS","feels") = 16.91;
!!p33_dac_fedem_heat("CHA","feels") = 17.30;
!!p33_dac_fedem_heat("USA","feels") = 18.62;
!!p33_dac_fedem_heat("JPN","feels") = 18.63;
!!p33_dac_fedem_heat("EUR","feels") = 18.82;
!!p33_dac_fedem_heat("CAZ","feels") = 20.20;
!!p33_dac_fedem_heat("REF","feels") = 21.18;
!!p33_dac_fedem_heat("NEU","feels") = 23.08;

*** FS: sensitivity on DAC efficiency
$if not "%cm_DAC_eff%" == "off" parameter p33_dac_fedem_fac(entyFeStat) / %cm_DAC_eff% /;
$if not "%cm_DAC_eff%" == "off" p33_dac_fedem(entyFeStat) = p33_dac_fedem(entyFeStat) * p33_dac_fedem_fac(entyFeStat);

*** EOF ./modules/33_CDR/DAC/datainput.gms
