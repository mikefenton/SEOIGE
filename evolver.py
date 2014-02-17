"""evolver: A standalone evolutionary algorithm
Copyright (c) 2013 Michael Fenton, Jonathan Byrne, Erik Hemberg and James McDermott
Hereby licensed under the GNU GPL v3."""

import sys, os, copy, random, math, datetime, shutil, subprocess, gc, xlrd
from operator import itemgetter
import analyser as AZR
import datetime
import multiprocessing
from images2gif import writeGif
from PIL import Image
#random.seed(1)
CROSSOVER = True
OPTIMIZE = False
time_list = []

class StructuralFitness():
    """Fitness function for testing generated mesh programs. """

    def __init__(self):
        self.maximise = False # false = smaller is better        

    def __call__(self, arse, unique_id, program, genome, used_codons_a, used_codons_b, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG = False, PRINT = False):
        import analyser as AZR
        analyser = AZR.Analyser(arse, unique_id, program, genome, used_codons_a, used_codons_b, ave, mats, generation, LOAD, MATERIALS_FILE)
        total_fitness = analyser.test_mesh(FITNESSES, LOAD, DEBUG, False, PRINT)
        fitness_results = []
        for i in total_fitness[0]:
            fitness_results.append(i)
        return analyser.used_codons_a, analyser.used_codons_b, fitness_results, total_fitness[1], analyser.node_list, analyser.edge_list

class Individual(object):
    """A GE individual"""

    def __init__(self, genome_a, genome_b, length=200):
        if genome_a == None and genome_b == None:
            self.genome_a = [random.randint(0, CODON_A_SIZE)
                           for _ in xrange(length)]
            genome_b = random.randint(0, CODON_B_SIZE)
            self.genome_b = [genome_b for _ in xrange(length)] # two loads chs [139, 109, 85, 110, 0, 109, 30, 0, 109, 38]# one load chs [133, 126, 86, 88, 2, 109, 68, 5, 90, 3]
        else:
            self.genome_a = copy.deepcopy(genome_a)
            self.genome_b = copy.deepcopy(genome_b)
        self.name = None
        self.bad = default_fitness(FITNESS_FUNCTION.maximise)
        self.uid = None
        self.phenotype = None
        self.rank = None
        self.distance = None
        self.good = False
        self.node_list = []
        self.edge_list = []
        self.fitness = []
        for i in range(len(FITNESSES)):
            self.fitness.append(int(self.bad))
        self.used_codons_a = len(self.genome_a)
        self.used_codons_b = len(self.genome_b)
        self.condon_list = []
        self.ave = False

    def __lt__(self, other):
        if FITNESS_FUNCTION.maximise:
            return self.fitness < other.fitness
        else:
            return other.fitness < self.fitness

    def __str__(self):
        return ("Individual: " + " uid: " + str(self.uid)
                + str(self.genome_a) + str(self.genome_b) + "; " + str(self.fitness))

    def evaluate(self, arse, fitness, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT):
        self.name = arse
        self.used_codons_a, self.used_codons_b, self.total_fitness, self.good, self.node_list, self.edge_list = fitness(self.name, self.uid, self.phenotype, self.genome_b, self.used_codons_a, self.used_codons_b, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT)
        self.fitness = []
        for i in self.total_fitness:
            b = round(i)
            self.fitness.append(b)
        return self.name, self.fitness, self.phenotype, self.used_codons_a, self.used_codons_b, self.node_list, self.edge_list, self.good

    def choose_genome(self):
        """returns either genome at random"""
        number = round(random.random())
        if number == 0:
            if self.used_codons_a != 0:
                return [self.genome_a, self.used_codons_a, CODON_A_SIZE]
            else:
                return [self.genome_b, self.used_codons_b, CODON_B_SIZE]
        else:
            if self.used_codons_b != 0:
                return [self.genome_b, self.used_codons_b, CODON_B_SIZE]
            else:
                return [self.genome_a, self.used_codons_a, CODON_A_SIZE]
    
    def set_values(self, values):
        self.phenotype = values['phenotype']
        self.genome_a = values['genome_a']
        self.genome_b = values['genome_b']
        self.used_codons_a = values['used_codons_a']
        self.used_codons_b = values['used_codons_b']

    def save_result(self, result):
        """assign values from result array to individual"""
        self.name, self.fitness, self.phenotype, self.used_codons_a, self.used_codons_b, self.node_list, self.edge_list, self.good = result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7]

def initialise_population(size=10):
    """Create a popultaion of size and return"""
    return [Individual(None, None) for _ in xrange(size)]


def print_stats(generation, individuals, best_ever, mats, best_nos, MATERIALS_FILE, DEBUG = False):
    """Print info about each run in real time"""
    time = datetime.datetime.now()
    time_list.append(time)
    def ave(values):
        return float(sum(values)) / len(values)
    def std(values, ave):
        return math.sqrt(float(sum((value - ave) ** 2
                                   for value in values)) / len(values))
    ave_fitnesses = []
    std_fitnesses = []
    fit_list = []
    for fit in range(len(FITNESSES)):
        ave_fitness = round(ave([i.fitness[fit]for i in individuals
                          if i.phenotype is not None]))
        ave_fitnesses.append(ave_fitness)
    for fit, n in enumerate(ave_fitnesses):    
        std_fitness = round(std([i.fitness[fit] for i in individuals
                          if i.phenotype is not None], n))
        std_fitnesses.append(std_fitness)
    ave_used_codons_a = ave([i.used_codons_a for i in individuals
                           if i.phenotype is not None])
    ave_used_codons_b = ave([i.used_codons_b for i in individuals
                           if i.phenotype is not None])
    max_used_codons_a = max([i.used_codons_a for i in individuals
                           if i.phenotype is not None])
    max_used_codons_b = max([i.used_codons_b for i in individuals
                           if i.phenotype is not None])
    min_used_codons_a = min([i.used_codons_a for i in individuals
                           if i.phenotype is not None])
    min_used_codons_b = min([i.used_codons_b for i in individuals
                           if i.phenotype is not None])
    std_used_codons_a = std([i.used_codons_a for i in individuals
                           if i.phenotype is not None], ave_used_codons_a)
    if len(time_list)>1:
        time_taken = time_list[-1]-time_list[-2]
    else:
        time_taken = time_list[0]
    print("Gen:%d\n  First: %s\n  Averages: %s\n  Used A Codons: [Min, Ave, Max]: [%.1f, %.1f, %.1f]\n  Used B Codons: [Min, Ave, Max]: [%.1f, %.1f, %.1f]\n  Time Taken:[%s]\n"
          % (generation, best_ever.fitness, ave_fitnesses, min_used_codons_a, ave_used_codons_a, max_used_codons_a, min_used_codons_b, ave_used_codons_b, max_used_codons_b, time_taken))
    if not DEBUG:
        filename = "./EliteResults/" + str(TIME_STAMP)
        savefile = open(filename, 'a')
        savefile.write("Gen:\t%d\tNo. Fit Indivs:\t%d\tFirst:\t%s\t\tAverage:\t%s\t\tUsed A Codons:\t%.1f\t%.1f\t%.1f\t\tUsed B Codons:\t%.1f\t%.1f\t%.1f\t\t%s"
          % (generation, best_nos, best_ever.fitness[0],  ave_fitnesses[0], min_used_codons_a, ave_used_codons_a, max_used_codons_a, min_used_codons_b, ave_used_codons_b, max_used_codons_b, time_taken)+ "\n")
        savefile.close()
    return best_ever.fitness[0]

def default_fitness(maximise=False):
    if maximise:
        return - DEFAULT_FIT
    else:
        return DEFAULT_FIT

def int_flip_mutation(individual, use_prob=True):
    """Mutate the individual by randomly chosing a new int with
    probability p_mut. Works per-codon, hence no need for
    "within_used" option."""
    if OPTIMIZE:
        genome = individual.genome_a
        codons = individual.used_codons_a
        CODON_SIZE = CODON_A_SIZE
    else:
        info = individual.choose_genome()
        genome = info[0]
        codons = info[1]
        CODON_SIZE = info[2]
    if use_prob:
        for i in xrange(len(genome)):
            if random.random() < MUTATION_PROBABILITY:
                genome[i] = random.randint(0, CODON_SIZE)
    else:
        idx = random.randint(0, codons - 1)
        genome[idx] = genome[idx] + 1
    return individual

def tournament_selection(population):
    """Given an entire population, draw <tournament_size> competitors
    randomly and return the best."""
    tournament_size = RETURN_PERCENT(3)
    winners = []
    while len(winners) < GENERATION_SIZE:
        competitors = random.sample(population, tournament_size)
        competitors.sort(reverse=True)
        winners.append(competitors[0])
    return winners

def onepoint_crossover(p1, p2, within_used=True):
    """Given two individuals, create two children using one-point
    crossover and return them."""
    # Get the chromosomes
    if CROSSOVER:
        GENOME_A = False
        GENOME_B = False
        p1a = [p1.genome_a, p1.used_codons_a, CODON_A_SIZE]
        p1b = [p1.genome_b, p1.used_codons_b, CODON_B_SIZE]
        p2a = [p2.genome_a, p2.used_codons_a, CODON_A_SIZE]  
        p2b = [p2.genome_b, p2.used_codons_b, CODON_B_SIZE]
        if OPTIMIZE:
            GENOME_A = True
            p1c,p2c = p1a,p2a
        elif random.random() > 0.5:
            GENOME_A = True
            p1c,p2c = p1a,p2a
        else:
            GENOME_B = True
            p1c,p2c = p1b,p2b
        # Uniformly generate crossover points.
        maxp1, maxp2 = p1c[1], p2c[1]
        use = min(maxp1, maxp2)
        pt_p1 = random.randint(1, use)
        # Make new chromosomes by crossover: these slices perform copies
        if random.random() < CROSSOVER_PROBABILITY:
            if GENOME_A:
                c = p1c[0][:pt_p1] + p2c[0][pt_p1:]
                d = p1b[0]
                e = p2c[0][:pt_p1] + p1c[0][pt_p1:]
                f = p2b[0]
            else:
                c = p1a[0]
                d = p1c[0][:pt_p1] + p2c[0][pt_p1:]
                e = p2a[0]
                f = p2c[0][:pt_p1] + p1c[0][pt_p1:]
        else:
            c, d = p1a[0], p1b[0]
            e, f = p2a[0], p2b[0]
    else:
        c, d = p1.genome_a, p1.genome_b
        e, f = p2.genome_a, p2.genome_b
    # Put the new chromosomes into new individuals
    return [Individual(c, d), Individual(e, f)]

def evaluate_fitness(individuals, grammar, fitness, ave, mats, generation, MATERIALS_FILE, LOAD, FITNESSES, DEBUG = False, MULTI_CORE = True, PRINT = False):
    """Perform the mapping and evaluate each individual across multiple available cores"""
    if MULTI_CORE:
        cores = multiprocessing.cpu_count() #   use all available cores
        pool = multiprocessing.Pool(processes=cores)
        for arse, ind in enumerate(individuals):
            bind = (arse, fitness, ind, grammar, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT)
            # Perform the mapping for each individual
            pool.apply_async(parallelize_indivs, args = (bind, ), callback = ind.save_result)
        pool.close()    
        pool.join()
    else:
        for arse, ind in enumerate(individuals):
            bind = (arse, fitness, ind, grammar, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT)
            parallelize_indivs(bind)
    counter = 0
    pounder = 0
    for ind in individuals:
        if ind.phenotype == None:
            bind = (1, fitness, ind, grammar, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT)
            parallelize_indivs(bind)
            if ind.phenotype == None:
                counter += 1
        if ind.good == True:
            pounder += 1
    if counter:
        print "Number of individuals with no phenotype:",counter

def parallelize_indivs(bind):
    """evaluates an individual using the Analyser class"""
    arse, fitness, ind, grammar, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT = bind[0], bind[1], bind[2], bind[3], bind[4], bind[5], bind[6], bind[7], bind[8], bind[9], bind[10], bind[11]
    values = grammar.generate(ind.genome_a, ind.genome_b)
    ind.set_values(values)
    if ind.phenotype != None:
        everything = ind.evaluate(arse, fitness, ave, mats, generation, LOAD, MATERIALS_FILE, FITNESSES, DEBUG, PRINT)
    else:
        print "BROKEN PHENOTYPE"
        AZR.log_error(ind.genome, "genotype could not be mapped")
    return everything

def step(percent, counter, previous_best, max_gens, parent_pop, grammar, selection, best_ever, fitness, ave, mats, generation, MATERIALS_FILE, LOAD, FITNESSES, DEBUG = False, MULTI_CORE = True):
    """perform single iteration and return next generation"""
    elites = []
    parent_pop.sort(reverse = True)
    while len(elites) < ELITE_SIZE:
        next = parent_pop[len(elites)]#.pop(0)
        elites.append(next)
    #Select parents using tournament selection, tournament size 3
    pop_size = len(parent_pop)
    parents = selection(parent_pop)
    #Crossover parents and add to the new population
    child_pop = []
    while len(child_pop) < GENERATION_SIZE:
        child_pop.extend(onepoint_crossover(*random.sample(parents, 2)))
    #Mutate the new child population
    child_pop = list(int_flip_mutation(child) for child in child_pop)
    #Evaluate the fitness of the new population
    total_pop = []
    papal_top = []
    papal_top.extend(child_pop)
    for i in papal_top:
        i.ave = ave
    evaluate_fitness(papal_top, grammar, fitness, ave, mats, generation, MATERIALS_FILE, LOAD, FITNESSES, DEBUG, MULTI_CORE)
    p = max(elites)
    for i, ind in enumerate(elites):
        ind.uid = i
    pre_elites = copy.deepcopy(elites)
    for i, elite in enumerate(elites):
        evaluate_fitness([elites[i]], grammar, fitness, elite.ave, mats, generation, MATERIALS_FILE, LOAD, FITNESSES, DEBUG, MULTI_CORE)
    post_elites = elites[0]
    papal_top.extend(elites)
    papal_top.sort(reverse = True)
    total_pop = papal_top
    new_pop = []
    total_pop.sort(reverse = True)
    for i, ind in enumerate(total_pop):
        if i < (GENERATION_SIZE):
            new_pop.append(ind)

    # Finds the number of fit individuals per generation
    power = 0
    for kid in total_pop:
        if kid.good == True:
            power += 1
    best_ever = max(total_pop)
    if MAKE_GIF:
        west_Mesh = AZR.Analyser(best_ever.name, "best", str(best_ever.phenotype), best_ever.genome_b, best_ever.used_codons_a, best_ever.used_codons_b, ave, MATS, generation, LOAD, MATERIALS_FILE)
        west_Mesh.create_graph(str(TIME_STAMP), str(best_ever.phenotype), best_ever.genome_b, LOAD, FITNESSES, DEBUG, FINAL=True)
        west_Mesh.create_mesh(str(generation))
    return new_pop, best_ever, power, percent

def assign_size(MATERIALS_FILE):
    """generates a list of materials with full section properties"""
    tables = open('./tables/' + str(MATERIALS_FILE) + '.txt', 'r')
    number = 0
    beams = []
    if MATERIALS_FILE == "CSSTables" or MATERIALS_FILE == "test1": # Units in lb (Imperial)
        for i, line in enumerate(tables):
            if line.startswith('#'):
                number = number + 1
                tables.readlines
            else:
                line = line.split()
                idx = i-number
                name = "cable"
                diameter = float(line[0])
                thickness = float(line[1])
                area = float(line[3])*100 # in millimeters squared
                iy = float(line[5]) # in millimeters to the four # = ix
                emod = 68950 #206843# in Newtons/millimeters or MegaPascals
                density = 2.76799 * 10 ** (-6)  #7.4182132 * 10 ** (-6)  #  in kg/millimeter cubed
                unitweight = area * density # in kg per millimeter
                one, one_point_five, two, two_point_five, three, three_point_five, four, five, six, seven, eight, nine, ten = 1,1.5,2,2.5,3,3.5,4,5,6,7,8,9,10
                material_properties = {'id':idx,'1':one,'1.5':one_point_five,'2':two,'2.5':two_point_five,
                                        '3':three,'3.5':three_point_five,'4':four,'5':five,'6':six,'7':seven,
                                        '8':eight,'9':nine,'10':ten, 'diameter':diameter,
                                        'thickness':thickness,'unitweight':unitweight,
                                        'area':area,'I':iy,'emod':emod,#'ix':ix,
                                        'density':density, 'name':name}
                if material_properties not in beams:
                    beams.append(material_properties)
        tables.close()
    else: # Units in kg (Metric)
        for i, line in enumerate(tables):
            if line.startswith('#'):
                number = number + 1
                tables.readlines
            else:
                line = line.split()
                idx = i-number
                diameter = float(line[0])
                thickness = float(line[1])
                unitweight = float(line[2])/1000 # in kg per meter
                area = float(line[3])* 100 # in millimeters squared
                iy = float(line[5])*10000 # in millimeters to the four = ix
                emod = 210000 # in Newtons/millimeters or MegaPascals
                density = 7.85 * 10 ** (-6)  # in kg/millimeter cubed
                if float(diameter) > 270:
                    two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, one, one_point_five, two_point_five, three_point_five = float(line[11]), float(line[12]),float(line[13]),float(line[14]),float(line[15]), float(line[16]),float(line[17]),float(line[18]),float(line[19]), float(line[20]),float(line[21]),float(line[22]),float(line[23]), 0, 0, 0, 0
                else:
                    one, one_point_five, two, two_point_five, three, three_point_five, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen = float(line[11]), float(line[12]),float(line[13]),float(line[14]),float(line[15]), float(line[16]),float(line[17]),float(line[18]),float(line[19]), float(line[20]),float(line[21]),float(line[22]),float(line[23]), 0, 0, 0, 0
                material_properties = {'id':idx,'1':one,'1.5':one_point_five,'2':two,'2.5':two_point_five,
                                        '3':three,'3.5':three_point_five,'4':four,'5':five,'6':six,'7':seven,
                                        '8':eight,'9':nine,'10':ten, '11':eleven, '12':twelve, '13':thirteen,
                                        '14':fourteen,'diameter':diameter,
                                        'thickness':thickness,'unitweight':unitweight,
                                        'area':area,'I':iy,'emod':emod,
                                        'density':density}
                if material_properties['thickness'] >= THICKNESS_LIMITER:
                    number = number + 1
                else:
                    beams.append(material_properties)
        tables.close()
    return beams

def search_loop(fitness_function, max_gens, individuals, grammar, selection, FITNESSES, DEFAULT_FIT, MATERIALS_FILE, LOAD, MATS, DEBUG = False, MULTI_CORE = True):
    """Loop over max generations"""
    #Evaluate initial population
    ave = False
    evaluate_fitness(individuals, grammar, fitness_function, ave, MATS, 1, MATERIALS_FILE, LOAD, FITNESSES, DEBUG, MULTI_CORE)
    best_ever = max(individuals)
    if MAKE_GIF:
        west_Mesh = AZR.Analyser(best_ever.name, "best", str(best_ever.phenotype), best_ever.genome_b, best_ever.used_codons_a, best_ever.used_codons_b, ave, MATS, 1, LOAD, MATERIALS_FILE)
        west_Mesh.create_graph(str(TIME_STAMP), str(best_ever.phenotype), best_ever.genome_b, LOAD, FITNESSES, DEBUG, FINAL=True)
        west_Mesh.create_mesh(str(1))
    last_ten = [12344587, 9876984732865, 432534458709, 987987958755, 3452089995, 1343125324, 13432424242, 9860980986, 465436665 ,best_ever.fitness[0]]
    individuals.sort(reverse=True)
    hound = 0
    counter = 0
    for ind in individuals:
        if ind.good == True:
            hound += 1
    previous_best = hound
    if hound == 0:
        counter += 1
    blenderation = 1 - counter
    percentage = float(blenderation)/float(max_gens)
    have = print_stats(1, individuals, best_ever, MATS, hound, MATERIALS_FILE, DEBUG)
    SWITCH = False
    for generation in xrange(2, (max_gens + 1)):
        if last_ten[0] == last_ten[-1]:
            ave = True
            SWITCH = True
        elif SWITCH:
            ave = True
        else:
            ave = False
        individuals, best_ever, best_nos, percent = step(
            percentage, counter, previous_best, max_gens, individuals, grammar, selection, best_ever, fitness_function, ave, MATS, generation, MATERIALS_FILE, LOAD, FITNESSES, DEBUG, MULTI_CORE)
        last_ten.pop(0)
        last_ten.append(best_ever.fitness[0])
        previous_best = best_nos
        percentage = percent
        if previous_best == 0:
            counter +=1
        have = print_stats(generation, individuals, best_ever, MATS, best_nos, MATERIALS_FILE, DEBUG)
    return best_ever, ave

def RETURN_PERCENT(num):
    percent = int(round(POPULATION_SIZE/100))
    if percent < num:
        return num
    else:
        return percent

# GE Properties
GRAMMAR_FILE = "grammars/Delaunay_cantilever_test_3.bnf"
MATERIALS_FILE = "CHSTables"
    # "CSSTables"
    # "CHSTables"
    # "CHSTablesSortedByArea"
    # "CHSTablesSortedByThickness"
    # "test1"
THICKNESS_LIMITER = 1000 # in millimeters
CODON_A_SIZE = 100000
MATS = assign_size(MATERIALS_FILE)
CODON_B_SIZE = len(MATS)-1
POPULATION_SIZE = 100
GENERATION_SIZE = POPULATION_SIZE
GENERATIONS = 10
ELITE_SIZE = RETURN_PERCENT(1)
FITNESSES = [0]
    # 0 = total_weight
    # 2 = max_displacement
    # 4 = cum_difference
    # 6 = total_cost
DEFAULT_FIT = 1000000000000
FITNESS_FUNCTION = StructuralFitness()
MUTATION_PROBABILITY = 0.01
CROSSOVER_PROBABILITY = 0.75
ACTUAL_GRAMMAR = GRAMMAR_FILE.split("/")[1].split(".")[0]
#create a timestamp
now = datetime.datetime.now()
hms = "%02d%02d%02d" % (now.hour, now.minute, now.second)
TIME_STAMP = (str(now.day) + "_" + str(now.month) + "_" + hms)
SHOW_FINAL = False
SAVE_DXF = False
SAVE_BEST = False
MAKE_GIF = False
DEBUG = True
MULTI_CORE = True
REMOVE_DUPLICATES = False
LOAD = 444800 # 8000000 # 12000000 # 17500000 # 667200 # 444800

def run_all():
    from optparse import OptionParser
    import grammar as GRAMMAR
    import datetime
    if os.path.isdir("/home/michael/Dropbox/Collij/Mike/truss/population"):
        shutil.rmtree("/home/michael/Dropbox/Collij/Mike/truss/population")
    os.mkdir("/home/michael/Dropbox/Collij/Mike/truss/population")
    if os.path.isdir("/home/michael/Dropbox/Collij/Mike/truss/slf"):
        shutil.rmtree("/home/michael/Dropbox/Collij/Mike/truss/slf")
    os.mkdir("/home/michael/Dropbox/Collij/Mike/truss/slf")
    if os.path.isdir("/home/michael/Dropbox/Collij/Mike/truss/Pics"):
        shutil.rmtree("/home/michael/Dropbox/Collij/Mike/truss/Pics")
    os.mkdir("/home/michael/Dropbox/Collij/Mike/truss/Pics")
    global MATS, MULTI_CORE, DEBUG, POPULATION_SIZE, FITNESS_FUNCTION, GENERATIONS, ELITE_SIZE, FITNESSES, DEFAULT_FIT, GENERATION_SIZE, MUTATION_PROBABILITY, CROSSOVER_PROBABILITY, ACTUAL_GRAMMAR, now, hms, TIME_STAMP, SHOW_FINAL, SAVE_DXF, SAVE_BEST, LOAD, GRAMMAR_FILE, MATERIALS_FILE, THICKNESS_LIMITER, CODON_A_SIZE, CODON_B_SIZE
    now = datetime.datetime.now()
    hms = "%02d%02d%02d" % (now.hour, now.minute, now.second)
    TIME_STAMP = (str(now.day) + "_" + str(now.month) + "_" + hms)
    time1 = datetime.datetime.now()
    print "start:", time1
    parser = OptionParser(usage="if nothing is specified, it uses the default values specified in evolver class")
    parser.set_defaults(pop_size=POPULATION_SIZE, generations=GENERATIONS,
                        elite_size=ELITE_SIZE, mutation=MUTATION_PROBABILITY,
                        bnf_grammar=GRAMMAR_FILE, crossover=CROSSOVER_PROBABILITY)
    parser.add_option("-p", "--population", dest="pop_size",
                      help=" Number of individuals in the population")
    parser.add_option("-g", "--generations", dest="generations",
                      help="Number of iterations of the algorithm")
    parser.add_option("-e", "--elite_size", dest="elite_size",
                      help=" How many get copied to next generation")
    parser.add_option("-m", "--mutation", dest="mutation",
                      help="probability of mutation on a per-codon basis")
    parser.add_option("-c", "--crossover", dest="crossover",
                      help="probability of crossover")
    parser.add_option("-b", "--bnf_grammar", dest="bnf_grammar",
                      help="bnf grammar for mapping")
    opts, args = parser.parse_args()
    POPULATION_SIZE = int(opts.pop_size)
    GENERATION_SIZE = int(opts.pop_size)
    GENERATIONS = int(opts.generations)
    ELITE_SIZE = int(opts.elite_size)
    MUTATION_PROBABILITY = float(opts.mutation)
    CROSSOVER_PROBABILITY = float(opts.crossover)
    GRAMMAR_FILE = opts.bnf_grammar

    # Read grammar
    BNF_GRAMMAR = GRAMMAR.Grammar(GRAMMAR_FILE)
    # Create Individual
    INDIVIDUALS = initialise_population(POPULATION_SIZE)
    # Loop
    LAST_POP, ave = search_loop(FITNESS_FUNCTION, GENERATIONS, INDIVIDUALS, BNF_GRAMMAR, tournament_selection, FITNESSES, DEFAULT_FIT, MATERIALS_FILE, LOAD, MATS, DEBUG, MULTI_CORE)
    time2 = datetime.datetime.now()
    print "end:", time2
    total_time = time2 - time1
    print "\ntime taken:", total_time
    print "\nBest fitness value:", LAST_POP.fitness,"\n"
    best_Mesh = AZR.Analyser(LAST_POP.name, "best", str(LAST_POP.phenotype), LAST_POP.genome_b, LAST_POP.used_codons_a, LAST_POP.used_codons_b, ave, MATS, GENERATION_SIZE, LOAD, MATERIALS_FILE)
    best_Mesh.test_mesh(FITNESSES, LOAD, DEBUG, FINAL = True)
    best_Mesh.create_mesh('best')
    if not DEBUG:
        newfilename = "./EliteResults/" + str(TIME_STAMP)
        savefile = open(newfilename, 'a')
        savefile.write("\nBest fitness values: " + str(LAST_POP.fitness) + "\n\n# Grammar = " + str(ACTUAL_GRAMMAR) + "\n# Population Size = " + str(POPULATION_SIZE) + "\n# Generation Size = " + str(GENERATION_SIZE) + "\n# Generations = " + str(GENERATIONS) + "\n# Mutation = " + str(MUTATION_PROBABILITY) + "\n# Crossover = " +str(CROSSOVER_PROBABILITY) + "\n# Codon A Size = " + str(CODON_A_SIZE) + "\n# Codon B Size = " + str(CODON_B_SIZE) + "\n")
        for i in range(len(FITNESSES)):
            savefile.write("\n# Fitness " + str(i+1) +": " + str(best_Mesh.fitness_selections[FITNESSES[i]+ 1]))
        best_Mesh.create_graph(str(TIME_STAMP), str(LAST_POP.phenotype), LAST_POP.genome_b, LOAD, FITNESSES, DEBUG, FINAL=True)
        savefile.write("\n\n# Span: " + str(best_Mesh.span) + " mm")
        if best_Mesh.UDL:
            savefile.write("\n# Load: " + str(best_Mesh.total_UDL) + " N UDL")
        elif best_Mesh.UDL_points:
            savefile.write("\n# Load: " + str(best_Mesh.total_UDL) + " N UDL (Point Load Approximation)")
        elif best_Mesh.point:
            savefile.write("\n# Load: " + str(best_Mesh.point_load) + " N Point Load")
        savefile.write("\n\n# Total time taken for run: " + str(total_time))
        filename = "./EliteResults/" + str(TIME_STAMP) + "_best"
        savefile.close()
    if SAVE_DXF:
        print "\nSaving best as DXF"
        best_Mesh.save_dxf(str(TIME_STAMP)) 
    if SHOW_FINAL:
        #using medit to show the graph
        meshName = 'population/best.mesh'
        cmd = 'ffmedit '+meshName
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        process.communicate()
        #using slffea to show the mesh
        best_Mesh.apply_stresses(best_Mesh.edge_list)
        best_Mesh.create_slf_file(best_Mesh.edge_list)
        best_Mesh.test_slf_file()
        best_Mesh.parse_results(best_Mesh.edge_list)
        best_Mesh.show_analysis()
    gc.enable()
    if MAKE_GIF:
        print "Making GIF image of best indivs"
        import glob
        pics = glob.glob('population/*.ppm')
        dicks = []
        for illest in pics:
            dick = illest.split("/")[1]
            wick = dick.split(".")[0]
            if str(wick) != "best":
                dicks.append(int(wick))
        bliss = sorted(dicks, key=int)
        images = [Image.open("population/"+str(fn)+".ppm") for fn in bliss]
        filename = "EliteResults/" + str(TIME_STAMP) + ".GIF"
        writeGif(filename, images, duration=0.1)
        print "GIF completed\n"   
        pics = []
    dirPath = "/home/michael/Dropbox/Collij/Mike/truss/population"
    fileList = os.listdir(dirPath)
    for fileName in fileList:
        os.remove(dirPath+"/"+fileName)
    shutil.rmtree("/home/michael/Dropbox/Collij/Mike/truss/slf")
    shutil.rmtree("/home/michael/Dropbox/Collij/Mike/truss/population")
    gc.collect()
    return TIME_STAMP

if __name__ == "__main__":
    run_all()
