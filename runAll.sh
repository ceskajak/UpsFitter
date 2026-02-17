#!/bin/bash

# Run for all bins in parallel
root -l -b -q 'fit_simul.cxx("BHT2", "int")' > /dev/null 2>&1 &
root -l -b -q 'fit_simul.cxx("BHT2", "0mult3")' > /dev/null 2>&1 &
root -l -b -q 'fit_simul.cxx("BHT2", "3mult6")' > /dev/null 2>&1 &
root -l -b -q 'fit_simul.cxx("BHT2", "6mult9")' > /dev/null 2>&1 &
root -l -b -q 'fit_simul.cxx("BHT2", "9mult12")' > /dev/null 2>&1 &
root -l -b -q 'fit_simul.cxx("BHT2", "12mult20")' > /dev/null 2>&1 &
root -l -b -q 'fit_simul.cxx("BHT2", "9mult20")' > /dev/null 2>&1 &

# Wait while showing job status
while jobs -r | grep -q .; do
    clear
    echo "Running jobs:"
    jobs
    sleep 2
done

echo "All bins finished."