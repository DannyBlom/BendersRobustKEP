import sys
import os

filebegin = sys.argv[1]
fileend = sys.argv[2]
maxnattacks = int(sys.argv[3])
avg = sys.argv[4]

nopt = 0
total_nscenarios = 0
avgtime = 0
avgsub = 0

# basic definitions
timelimit = 3600
timeshift = 10

for i in range(1, maxnattacks + 1):
    logfile = filebegin + str(i) + fileend
    # get results from test run
    os.system('grep "objective value:\\|@01\\|terminated after\\|\[scenario\]\|@06\|@95" ' + logfile + ' >> tmp_resfile.txt')

    f = open("tmp_resfile.txt", 'r')

    # lists and dictionaries to store information about instances in test set
    instances = []                  # instances that have been tested
    name = ""                       # name of current instance
    values = {}                     # optimal objective values (contains only information on optimally solved instances)
    time = {}                       # solution time per instance (timelimit also if instance hit memory limit)
    nscenarios = {}                 # number of explored scenarios for optimally solved instances
    nthirdstage = {}                # number of iterations in third stage (B&B nodes) for optimally solved instances

    nscen = 0
    for line in f:
        if line.startswith("@01"):
            name = line.split()[1].split('/')[-1].split(".txt")[0]
            instances.append(name)
            nscen = 0
        elif "terminated after" in line:
            content = line.split()
            totaltime = float(content[3])
            if "suboptimal" in line:
                time[name] = timelimit
            else:
                time[name] = min(float(totaltime), timelimit)

        elif "scenario" in line:
            nscen += 1
        elif line.startswith("@06"):
            nscenarios[name] = nscen
        elif line.startswith("objective value:"):
            val = float(line.split()[2])
            values[name] = val
        elif line.startswith("@95"):
            if not name in nthirdstage:
                nthirdstage[name] = 1
            else:
                nthirdstage[name] = nthirdstage[name] + 1

    f.close()

    os.system("rm tmp_resfile.txt")

    # get statistics
    arith = (sum(time[name] for name in values) + sum(timelimit for name in instances if not name in values)) / len(instances)

    geom = 1.0
    for name in instances:
        if name in values:
            geom *= time[name]**(1/len(instances))
        else:
            geom *= timelimit**(1/len(instances))

    shiftgeom = 1.0
    for name in instances:
        if name in values:
            shiftgeom *= (time[name] + timeshift)**(1/len(instances))
        else:
            shiftgeom *= (timelimit + timeshift)**(1/len(instances))
    shiftgeom -= timeshift

    cur_nopt = sum(1 for name in values)
    if cur_nopt > 0:
        avgscenario = sum(nscenarios[name] for name in values) / cur_nopt
        avgthirditer = sum(nthirdstage[name] for name in values) / cur_nopt
    else:
        avgscenario = 0
        avgthirditer = 0
    
    if avg == "g":
        avgtime += geom
    elif avg == "s":
        avgtime += shiftgeom
    else:
        avgtime += arith

    nopt += cur_nopt
    total_nscenarios += avgscenario
    avgsub += avgthirditer

avgtime = avgtime / maxnattacks
total_nscenarios = total_nscenarios / maxnattacks
avgsub = avgsub / maxnattacks

print("%7d & %7.1f & %7.1f & %7.1f" % (nopt, avgtime, total_nscenarios, avgsub))
