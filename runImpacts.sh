#!/bin/bash
if [ -z $1 ]
then
    mainout='AMass_gauss'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi
#text2workspace.py datacard_full_${mainout}_mmmt.txt
#
#combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_${mainout}_mmmt.root -t -1 --expectSignal=1 --doInitialFit --robustFit 1 --X-rtd MINIMIZER_analytic
#
#combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_${mainout}_mmmt.root -t -1 --expectSignal=1 --doFits --robustFit 1 --X-rtd MINIMIZER_analytic
#
#combineTool.py -M Impacts -n  mass_a40_mmmt -d datacard_full_${mainout}_mmmt.root -m  40 -o testimpacts
#
#plotImpacts.py -i testimpacts -o testimpacts_12May2021

#text2workspace.py datacard_full_mass_a40_${mainout}_mmmt.txt

# combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_mass_a40_${mainout}_mmmt.root -t -1 --doInitialFit  --freezeParameters MH
#
# combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_mass_a40_${mainout}_mmmt.root -t -1 --doFits --freezeParameters MH
#
# combineTool.py -M Impacts -n  mass_a40_mmmt -d datacard_full_mass_a40_${mainout}_mmmt.root -m  40 -o testimpacts

text2workspace.py datacard_full_${mainout}.txt
combineTool.py -M Impacts -m 40 -n mass_a40_${mainout} -d datacard_full_${mainout}.root -t -1 --doInitialFit --expectSignal=1 --freezeParameters MH --X-rtd ADDNLL_RECURSIVE=0  --robustFit 1
#combineTool.py -M Impacts -m 40 -n mass_a40_${mainout} -d datacard_full_${mainout}.root -t -1 --doInitialFit --expectSignal=1 --freezeParameters MH --X-rtd ADDNLL_RECURSIVE=0  --robustFit 1
#combineTool.py -M Impacts -m 40 -n mass_a40_${mainout} -d datacard_full_${mainout}.root -t -1 --doInitialFit --expectSignal=1 --freezeParameters MH --X-rtd ADDNLL_RECURSIVE=0 --rMax 62  --robustFit 1 --autoMaxPOIs '"*"' --autoBoundsPOIs '"*"' --setParameterRanges '"c3_bkg_Nominal=0.0001"' --rMin 0.0001
#combineTool.py -M Impacts -m 40 -n mass_a40_${mainout} -d datacard_full_${mainout}.root -t -1 --doFits --expectSignal=1 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --freezeParameters MH
combineTool.py -M Impacts -m 40 -n mass_a40_${mainout} -d datacard_full_${mainout}.root -t -1 --doFits --expectSignal=1 --freezeParameters MH --X-rtd ADDNLL_RECURSIVE=0 --robustFit 1
#combineTool.py -M Impacts -m 40 -n mass_a40_${mainout} -d datacard_full_${mainout}.root -t -1 --doFits --expectSignal=1 --freezeParameters MH --X-rtd ADDNLL_RECURSIVE=0 --rMax 62 --robustFit 1 --autoMaxPOIs '"*"' --autoBoundsPOIs '"*"' --setParameterRanges '"c3_bkg_Nominal=0.0001"' --rMin 0.0001

combineTool.py -M Impacts -n  mass_a40_${mainout} -d datacard_full_${mainout}.root -m  40 -o testimpacts_${mainout}

plotImpacts.py -i testimpacts_${mainout} -o testimpacts_${mainout}


###other stuff
#combine -M MultiDimFit datacard_full_2017_mmmt.root --algo grid --cl=0.68 --expectSignal=1 -t 1 -m 40 --freezeParameters MH --points=1000 --rMax 62
