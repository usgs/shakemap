#!/bin/sh

#
# This script runs a set of verification tests and then opens the 
# resulting plots for inspection. It should work on Mac OS X
# and linux. It probably won't work on anything else.
#
# This script must be run from the shakemap "bin" directory
# or it won't work.
#

CALLED_FROM_PYTEST='True'
export CALLED_FROM_PYTEST

basepath='../tests/data/eventdata'
filepath='current/products'

plotfiles=""

for i in 01 02 03 04 04b 05 06 07 08a 08b 08c 08d 08e 09 10 11
do
    vtest="verification_test_00$i"
    `shake $vtest assemble model`
    if [[ $i == 07 ]]; then
        `shake $vtest xtestplot_spectra`
    elif [[ $i == 11 ]]; then
        `shake $vtest xtestimage`
    elif [[ $i == 08e ]]; then
        `shake verification_test_0008 xtestplot_multi`
    elif [[ $i == 08a || $i == 08b || $i == 08c || $i == 08d ]]; then
        :
    else
        `shake $vtest xtestplot`
    fi
    if [[ $i == 07 ]]; then
        plotfiles="$plotfiles $basepath/$vtest/$filepath/${vtest}_spectra_plot.pdf"
    elif [[ $i == 08a ]]; then
        plotfiles="$plotfiles $basepath/$vtest/$filepath/verification_test_0008_PGA.pdf"
    elif [[ $i == 08b || $i == 08c || $i == 08d || $i == 08e ]]; then
        :
    elif [[ $i == 11 ]]; then
        plotfiles="$plotfiles $basepath/$vtest/$filepath/${vtest}_PSA3p0.pdf"
    else
        plotfiles="$plotfiles $basepath/$vtest/$filepath/${vtest}_PGA.pdf"
    fi
done

osname=`uname`

if [[ "$osname" == 'Darwin' ]]; then
    `open $plotfiles`
else
    `xdg-open $plotfiles`
fi
