# Code for wave and current induced SPM resuspension in NEMO 4.2 - BAMHBI (Black Sea Ecosystem model). Deliverable for NECCTON

Author: M. Choblet, building on code built by A. Capet, L. Vandenbulcke and T. Lovato
Date: 19.03.2025 


Shared code are excerpts from entire BAMHBI MY_SRC structure (model description paper in preparation). It contains code for the benthic model, which can not be run without the entire MY_SRC structure:

Function call overview of the bamhbi_bentic.F90 functions


* trcsms_my_trc for BAMHBI (to be shared in a separate repository) calls:
	- trc_wbstress (wavestress computation)
	- bamhbi_benthic, which itself calls BOTTOMSTRESS for total bottom stress. Resuspension/Deposition calculation + parameterized benthic-pelagic fluxes
	- trc_sms_ben (source minus sink)
