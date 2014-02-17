""" This program cues up and executes multiple runs
of SEOIGE. Results of runs are parsed and placed in a
spreadsheet for easy visual analysis. 
  
Copyright (c) 2014 Michael Fenton
Hereby licensed under the GNU GPL v3."""

import os, evolver, analyser, shutil, datetime, xlwt
from operator import itemgetter
from xlutils.copy import copy
from xlrd import open_workbook, xldate_as_tuple

def execute(RUNS, FITNESS, GRAMMAR, LOAD, GENERATIONS, OPT, RE_GEN, BEST_STOP, OPT_ALL, FOLDER):
    """execute multiple runs"""
    first_progs = []
    first_fits = []
    time1 = datetime.datetime.now()
    print "Multi-Run Start:", time1
    for run in range(RUNS):
        # Settings for each evolutionary run
        evolver.FITNESSES = FITNESS
            # 0 = total_weight
            # 2 = max_displacement
            # 4 = cum_difference
            # 6 = total_cost
        evolver.GRAMMAR_FILE = "grammars/" + str(GRAMMAR)
        evolver.POPULATION_SIZE = 1000
        evolver.GENERATIONS = GENERATIONS
        evolver.MUTATION_PROBABILITY = 0.01
        evolver.CROSSOVER_PROBABILITY = 0.75
        evolver.SHOW_FINAL = False
        evolver.SAVE_DXF = False
        evolver.SAVE_BEST = False
        evolver.MAKE_GIF = False
        evolver.DEBUG = False
        evolver.MULTI_CORE = True
        evolver.POP_FULL_STOP = False
        evolver.REMOVE_DUPLICATES = False
        evolver.LOAD = LOAD
        evolver.GRADUAL_LOAD = False
        if len(FITNESS) > 1:
            evolver.NSGA = True
        else:
            evolver.NSGA = False
        evolver.OPTIMIZE = False
        analyser.OPT_ALL = OPT_ALL
        analyser.BEST_STOP = BEST_STOP
        analyser.OPTIMIZE = OPT
        analyser.STEPS = 5
        analyser.GENOME_REWRITE = RE_GEN
        answer = evolver.run_all()
        first_progs.append(answer)
    
    for ans, answer in enumerate(first_progs):
        # find the best individual from each run and add it to a list
        shutil.copyfile("EliteResults/" + str(answer) + ".py", str(answer) + ".py")
        alll = __import__(answer)
        alll.SHOW = False
        run = alll.run()
        total = [str(answer)]
        for i, thing in enumerate(run):
            total.append(run[i])
        first_fits.append(total)
        # Remove duplicate files
        if os.path.exists("/home/michael/Dropbox/Collij/Mike/truss/" + str(answer) + ".py"):
            os.remove("/home/michael/Dropbox/Collij/Mike/truss/" + str(answer) + ".py")        
        if os.path.exists("/home/michael/Dropbox/Collij/Mike/truss/" + str(answer) + ".pyc"):
            os.remove("/home/michael/Dropbox/Collij/Mike/truss/" + str(answer) + ".pyc")
    
    time2 = datetime.datetime.now()
    total_time = time2 - time1    
    if os.path.isdir("/home/michael/Dropbox/Collij/Mike/truss/EliteResults/" + str(RUN_FOLDER)):
        pass
    else:
        os.mkdir("/home/michael/Dropbox/Collij/Mike/truss/EliteResults/" + str(RUN_FOLDER))
    if os.path.isdir("/home/michael/Dropbox/Collij/Mike/truss/EliteResults/" + str(RUN_FOLDER) + "/" + str(FOLDER)):
        pass
    else:
        os.mkdir("/home/michael/Dropbox/Collij/Mike/truss/EliteResults/" + str(RUN_FOLDER) + "/" + str(FOLDER))
    
    # Write info about best indiv from each run to a file.
    filename = "EliteResults/" + str(RUN_FOLDER) + "/" + str(FOLDER) +".txt"
    savefile = open(filename, 'w')
    for ans, answer in enumerate(first_fits):
        savefile.write("Run " + str(ans) + "\tBest: " + str(answer) + "\n")
    first_fits.sort(key=itemgetter(1))
    savefile.write("\nBEST: " + str(first_fits[0]))
    savefile.write("\n\nTotal time taken for " + str(RUNS) + " runs: " + str(total_time))
    savefile.close()

    shutil.copyfile("EliteResults/" + str(first_fits[0][0]) + ".py", str(first_fits[0][0]) + ".py")
    print "\nTotal time taken for",RUNS,"runs:", total_time
    print "\nBEST:",first_fits[0]
    final = __import__(first_fits[0][0])
    final.SHOW = False
    final.run()
    
    # Write all info to a spreadsheet
    WB = xlwt.Workbook()
    sh = WB.add_sheet(str(FOLDER))
    for ans, answer in enumerate(first_progs):
        results = open("EliteResults/" + str(answer), 'r')
        lines = iter(results)
        n = 0
        for line in lines:
            result = line.split()
            if result == []:
                break
            while result[0] == "Gen:":
                sh.write(n, ans+2,int(float(result[7])))
                n+=1
                break
        WB.save("EliteResults/" + str(RUN_FOLDER) + "/" + str(FOLDER) + ".xls")
        shutil.copyfile("EliteResults/" + str(answer), "EliteResults/" + str(RUN_FOLDER) + "/" + str(FOLDER) + "/" + str(answer))
        shutil.copyfile("EliteResults/" + str(answer) + ".py", "EliteResults/" + str(RUN_FOLDER) + "/" + str(FOLDER) + "/" + str(answer) + ".py")
        os.remove("/home/michael/Dropbox/Collij/Mike/truss/EliteResults/" + str(answer) + ".py")
        os.remove("/home/michael/Dropbox/Collij/Mike/truss/EliteResults/" + str(answer))

if __name__ == '__main__':
    RUN_FOLDER = "Test 11 - Now properly analyzing trusses"
#    execute(30, [0], "Delaunay_cantilever_test.bnf", 444800, 100, False, False, False, False, "Percentage Regular")
#    execute(30, [0], "Delaunay_cantilever_test_2.bnf", 444800, 100, False, False, False, False, "Just Percentage")
    execute(30, [0], "Delaunay_cantilever_test_3.bnf", 444800, 100, False, False, False, False, "Just The Right Solution")
#    execute(30, [0], "Delaunay_cantilever.bnf", 444800, 100, False, False, False, False, "Random Task")
    
 #   execute(30, [0], "Delaunay_cantilever_test.bnf", 444800, 500, True, True, True, True, "Opt Re Best All")
 #   execute(30, [0], "Delaunay_cantilever.bnf", 444800, 500, True, True, True, False, "Opt Re Best Fit")
 #   execute(30, [0], "Delaunay_cantilever_test.bnf", 444800, 500, True, True, False, True, "Opt Re All All")
 #   execute(30, [0], "Delaunay_cantilever.bnf", 444800, 500, True, True, False, False, "Opt Re All Fit")
 #   execute(30, [0], "Delaunay_cantilever_test.bnf", 444800, 500, True, False, True, True, "Opt Best All")
 #   execute(30, [0], "Delaunay_cantilever.bnf", 444800, 500, True, False, True, False, "Opt Best Fit")
 #   execute(30, [0], "Delaunay_cantilever_test.bnf", 444800, 500, True, False, False, True, "Opt All All")
 #   execute(30, [0], "Delaunay_cantilever.bnf", 444800, 500, True, False, False, False, "Opt All Fit")
