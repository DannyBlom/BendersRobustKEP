#!/usr/bin/env bash

for method in {1..3}
do
    for cyclelen in 3 4
    do
        for chainlen in 2 3 4
	do
	    for lifting in TRUE FALSE
            do
	        for policy in 1 2
	        do
		    echo -e "kidney/method = $method\nkidney/maxcyclelength = $cyclelen\nkidney/maxchainlength = $chainlen\nkidney/liftbenderscuts = $lifting\nkidney/recoursepolicy = $policy" > method$method-cyc$cyclelen-chain$chainlen-lift$lifting-policy$policy.set
	        done
            done
        done
    done
done    
