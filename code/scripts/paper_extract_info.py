import sys
import os

# basic definitions
timelimit = 3600
timeshift = 10

# get parameters on 1: log file, 2: average used for time measurement
logfile = sys.argv[1]
withstage3 = True
if int(sys.argv[2]) == 0:
    withstage3 = False
avg = 'a'
if len(sys.argv) > 3:
    avg = sys.argv[3]

# get results from test run
os.system('grep "objective value:\\|@01\\|terminated after\\|\[scenario\]\|@06\|@98\|@95" ' + logfile + ' >> tmp_resfile.txt')

f = open("tmp_resfile.txt", 'r')

# lists and dictionaries to store information about instances in test set
instances = []                  # instances that have been tested
name = ""                       # name of current instance
values = {}                     # optimal objective values (contains only information on optimally solved instances)
time = {}                       # solution time per instance (timelimit also if instance hit memory limit)
nscenarios = {}                 # number of explored scenarios for optimally solved instances
nfirststage = {}                # number of master iterations for optimally solved instances
nthirdstage = {}                # number of iterations in third stage (B&B nodes) for optimally solved instances
timestage1 = {}                 # percentage of time per instance spent for stage 1
timestage2 = {}                 # percentage of time per instance spent for stage 2
timestage3 = {}                 # percentage of time per instance spent for stage 3

nscen = 0
for line in f:
    if line.startswith("@01"):
        name = line.split()[1].split('/')[-1].split(".txt")[0]
        instances.append(name)
        nscen = 0
    elif "terminated after" in line:
        content = line.split()
        totaltime = float(content[3])
        time2 = float(content[7][:-1])
        time3 = float(content[10][:-2])
        if "suboptimal" in line:
            time[name] = timelimit
        else:
            time[name] = min(float(totaltime), timelimit)
        timestage2[name] = time2 / totaltime if totaltime > 0 else 0
        timestage3[name] = time3 / totaltime if totaltime > 0 else 0
                             
    elif "scenario" in line:
        nscen += 1
    elif line.startswith("@06"):
        nscenarios[name] = nscen
    elif line.startswith("objective value:"):
        val = float(line.split()[2])
        values[name] = val
    elif line.startswith("@98"):
        if not name in nfirststage:
            nfirststage[name] = 1
        else:
            nfirststage[name] = nfirststage[name] + 1
    elif line.startswith("@95"):
        if not name in nthirdstage:
            nthirdstage[name] = 1
        else:
            nthirdstage[name] = nthirdstage[name] + 1

f.close()

os.system("rm tmp_resfile.txt")

# print statistics
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

nopt = sum(1 for name in values)
if nopt > 0:
    avgscenario = sum(nscenarios[name] for name in values) / nopt
    avgmasteriter = sum(nfirststage[name] for name in values) / nopt
    avgthirditer = sum(nthirdstage[name] for name in values) / nopt

# If instances are missing "terminated after", we have hit a timeout of slurm job, meaning that the timelimit is reached indefinitely
avgtime2 = 100 * sum(timestage2[name] if name in timestage2 else 0 for name in instances) / len(instances)
avgtime3 = 100 * sum(timestage3[name] if name in timestage3 else 0 for name in instances) / len(instances)

avgtime = arith
if avg == "g":
    avgtime = geom
elif avg == "s":
    avgtime = shiftgeom

if withstage3:
    if nopt > 0:
        print(" \\num{%4d} & \\num{%7.2f} & \\SI{%4.2f}{\\percent} & \\SI{%4.2f}{\\percent} & \\num{%4.1f} & \\num{%7.1f} " % (len(values), avgtime, avgtime2, avgtime3, avgscenario, avgthirditer))
    else:
        print(" \\num{%4d} & \\num{%7.2f} & \\SI{%4.2f}{\\percent} & \\SI{%4.2f}{\\percent} &  %4s & %7s " % (len(values), avgtime, avgtime2, avgtime3, "---", "---"))
else:
    if nopt > 0:
        print(" \\num{%4d} & \\num{%7.2f} & \\SI{%4.2f}{\\percent} & \\num{%4.1f} & \\num{%7.1f} " % (len(values), avgtime, avgtime2, avgscenario, avgthirditer))
    else:
        print(" \\num{%4d} & \\num{%7.2f} & \\SI{%4.2f}{\\percent} &  %4s & %7s " % (len(values), avgtime, avgtime2, "---", "---"))
