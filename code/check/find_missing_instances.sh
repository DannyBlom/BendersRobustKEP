#!/usr/bin/bash

for attacks in 1 2 3 4
do
    for nodes in 20 50 100
    do
	for instance in {0..29}
        do
	    for cyclelen in 3 4
	    do
		for chainlen in 2 3 4
		do
		    for method in 1 2
		    do
		    	for lifting in TRUE FALSE
		    	do
			    for policy in 1 2
			    do
			
			    	#FILE=results/s145344.mcs.default.q.klimentova${nodes}.nattacks${attacks}.$((${instance}+1))_Klimentova_${nodes}_${instance}.txt.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method1-complexFALSE-cyc${cyclelen}-chain${chainlen}-lift${lifting}-policy${policy}.out
			    	
				if [[ ${method} -eq 1 ]]
				then
				    FILE=results/non_missing_results/s145344.mcs.default.q.klimentova${nodes}.nattacks${attacks}.$((${instance}+1))_Klimentova_${nodes}_${instance}.txt.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method1-complexFALSE-cyc${cyclelen}-chain${chainlen}-lift${lifting}-policy${policy}.out
			    	    if [ ! -f "$FILE" ]; then
				        echo "$FILE"
			    	    fi
				fi

				if [[ ${method} -eq 2 ]]
				then
				    FILE=results/non_missing_results/s145344.mcs.default.q.klimentova${nodes}.nattacks${attacks}.$((${instance}+1))_Klimentova_${nodes}_${instance}.txt.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method2-cyc${cyclelen}-chain${chainlen}-lift${lifting}-policy${policy}.out
				fi
			    done
			done
		    done
		
		    for method in 3
		    do
			#FILE=results/s145344.mcs.default.q.klimentova${nodes}.nattacks${attacks}.$((${instance}+1))_Klimentova_${nodes}_${instance}.txt.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method${method}-cyc${cyclelen}-chain${chainlen}-lift${lifting}-policy${policy}.out
				
			FILE=results/non_missing_results/s145344.mcs.default.q.klimentova${nodes}.nattacks${attacks}.$((${instance}+1))_Klimentova_${nodes}_${instance}.txt.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method3-cyc${cyclelen}-chain${chainlen}-liftFALSE-policy${policy}.out
			if [ ! -f "$FILE" ]; then
			    echo "$FILE"
			fi					    
		    done
		done
	    done
	done
    done
done
