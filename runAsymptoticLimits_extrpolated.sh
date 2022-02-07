echo "attempting extrapolation"
datacardnum=20
if [ -z $1 ]
then
    mainout='test'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi
if [ -z $2 ]
then
    channel='mmmt'
    echo "User may want to provide a channel"
else
    channel=$2
    echo 'using '$channel' as channel'
fi
for i in {18..62}
do
if [ $((${i} % 5)) -eq 0 ]
then
    datacardnum=${i}
fi
echo "working on mass point "${i}
#combine -M AsymptoticLimits -m ${i} -n _mass_${i}_mmem --run blind --freezeParameters MH datacard_full_mass_${datacardnum}_${mainout}_mmem.txt
#combine -M AsymptoticLimits -m ${i} -n _mass_${i}_mmem --run blind datacard_full_mass_${datacardnum}_${mainout}_mmem.txt
#combine -M AsymptoticLimits -m ${i} -n _mass_${i}_${channel} --X-rtd ADDNLL_RECURSIVE=0 --freezeParameters MH --run blind datacard_full_${mainout}.txt
combine -M AsymptoticLimits -m ${i} -n _mass_${i}_${mainout} --X-rtd ADDNLL_RECURSIVE=0 --freezeParameters MH --run blind datacard_full_${mainout}.txt
done



echo "hadding the files"
rm higgsCombine_aa_${mainout}_best.root
hadd higgsCombine_aa_${mainout}_best.root higgsCombine_mass_*_${mainout}.AsymptoticLimits.mH*.root
root -l -b -q 'plotLimit.C+("aa", ${mainout} ,2)'
