respath = "../check/results"                   # path to the results directory

NNODES = [20, 50, 100]
ATTACKS = [1, 2, 3, 4]
CYCLELEN = [3, 4]
CHAINLEN = [2, 3, 4]
LIFT = ["TRUE","FALSE"]
METHOD = [1, 2]                 # don't compare with branch-and-bound
POLICY = [1, 2]

baseinstances = []
for a in NNODES:
    for b in ATTACKS:
        for c in CYCLELEN:
            for d in CHAINLEN:
                baseinstances.append((a,b,c,d))

results = {}

# get results for FR policy
for (a,b,c,d) in baseinstances:
    for e in LIFT:
        for f in METHOD:
            if f == 1:
                rf = open(respath+"/check.mcs.default.q.klimentova{}.nattacks{}.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method{}-complexFALSE-cyc{}-chain{}-lift{}-policy{}.out".format(a,b,f,c,d,e,1), 'r')
            else:
                rf = open(respath+"/check.mcs.default.q.klimentova{}.nattacks{}.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method{}-cyc{}-chain{}-lift{}-policy{}.out".format(a,b,f,c,d,e,1), 'r')

            name = ""
            obj = -1.0

            for line in rf:
                if line.startswith("@01"):
                    name = line.split()[1].split('/')[-1].split(".")[0].split('_')[-1]
                elif line.startswith("objective value"):
                    obj = float(line.split()[2])
                elif line.startswith("[optimal]"):

                    # store solution
                    if (a,b,c,d,name) not in results:
                        assert name != ""
                        assert obj != -1.0
                        results[(a,b,c,d,name)] = [obj]

                elif line.startswith("@04"):
                    name = ""
                    obj = -1.0

            rf.close()

# get results for FSE policy
for (a,b,c,d) in baseinstances:
    for e in LIFT:
        for f in METHOD:

            if f == 1:
                rf = open(respath+"/check.mcs.default.q.klimentova{}.nattacks{}.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method{}-complexFALSE-cyc{}-chain{}-lift{}-policy{}.out".format(a,b,f,c,d,e,2), 'r')
            else:
                rf = open(respath+"/check.mcs.default.q.klimentova{}.nattacks{}.fullrecoursekidney.linux.x86_64.gnu.opt.spx2.mcs-login001.method{}-cyc{}-chain{}-lift{}-policy{}.out".format(a,b,f,c,d,e,2), 'r')
            
            name = ""
            obj = -1.0

            for line in rf:
                if line.startswith("@01"):
                    name = line.split()[1].split('/')[-1].split(".")[0].split('_')[-1]
                elif line.startswith("objective value"):
                    obj = float(line.split()[2])
                elif line.startswith("[optimal]"):

                    # store solution
                    if len(results[(a,b,c,d,name)]) != 2:
                        results[(a,b,c,d,name)].append(obj)

                elif line.startswith("@04"):
                    name = ""
                    obj = -1.0

            rf.close()


# evaluate results

print(r"\begin{table}[htbp]")
print(r"  \begin{scriptsize}")
print(r"    \caption{Comparison of average objective values with FR and FSE for all instances that are solved by cut-CC and/or cut-PICEF for both recourse policies.}")
print(r"    \label{tab:fr-vs-fse}")
print(r"    \begin{tabular*}{\textwidth}{@{}l@{\;\;\extracolsep{\fill}}rrr@{}}\toprule")
print(r"      $|V|$ & $z^{*}_{FR} - z^{*}_{FSE}$ & $\#instances$ & Relative count \\")

for nnodes in [20, 50, 100]:
    print(r"      \midrule")
    print(r"      \multicolumn{4}{@{}l}{$|V|=%d$} \\" % nnodes)
    data = {}
    ninstances = 0
    for key in results.keys():
        if key[0] != nnodes:
            continue
        # skip problems which could not be solved by both policies
        if len(results[key]) < 2:
            continue
        nnodes = key[0]
        diff = round(results[key][0] - results[key][1]) if results[key][0] - results[key][1] > 0 else 0

        if (nnodes, diff) in data:
            data[nnodes, diff] = data[nnodes, diff] + 1
        else:
            data[nnodes, diff] = 1
        ninstances += 1


    # print("difference\tcount\trelative")
    for nnodes, diff in data:
        print(r"          & %2d & %2d & %5.6f \\" % (diff, data[nnodes, diff],data[nnodes, diff]/ninstances))
    print(r"      \midrule")
    print(r"        Total & & {} & \\".format(ninstances))
print(r"      \bottomrule")
print(r"    \end{tabular*}")
print(r"  \end{scriptsize}")
print(r"\end{table}")