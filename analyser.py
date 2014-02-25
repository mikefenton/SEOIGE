""" This class contains functions for parsing the programs generated
by the grammar, creating the resulting graphs, converting them to
slffea files, running slffea, collating the results and using them to
generate a fitness value. 

ALL UNITS IN MILLIMETERS, NEWTONS
  
Copyright (c) 2012
Michael Fenton, Jonathan Byrne, Erik Hemberg and James McDermott
Hereby licensed under the GNU GPL v3."""
import time, graph, random, matplotlib.delaunay, grammar
import evolver
#from grammar import generate
from subprocess import PIPE, Popen
from shutil import copyfile
from datetime import datetime
from geometry import interpolate, three_d_line_length
from math import atan, tan, cos, pi, sqrt
from operator import itemgetter
from os import remove

# global class variables
# Fitness variables
DEATH_PENALTY = False
DEFAULT_FIT = 1000000000000

# Analysis variables
BUCKLING_CHECK = True
REMOVE_UNSTRESSED_EDGES = False
SHOW_ANALYSIS = False

OPTIMIZE = True
# Turns on the pre-fitness function optimizer which will
# optimize the whole population
STEPS = 5
# Number of optimization steps
OPT_ALL = True
# Optimize all individuals in a generation
BEST_STOP = False
# Stop optimization once better fitness is discovered
GENOME_REWRITE = True
# Re-write Ch.B with new member sizes

def eval_or_exec(name, time, s, gen, ave, MAT_FILE, used_codons_a, 
                 LOAD, generation, Fitnesses, DEBUG = False, FINAL=False):
    """Handles different return vals from eval/exec"""
    s = python_filter(name, time, s, gen, ave, MAT_FILE, used_codons_a, LOAD, 
                      generation, Fitnesses, DEBUG, FINAL)
    dictionary = {"itemgetter": itemgetter, "genome": gen, "graph": graph,
                  "interpolate": interpolate, "sqrt":sqrt, "atan": atan,
                  "tan": tan, "cos":cos, "pi":pi, "random": random,
                  "ave": ave, "triang": matplotlib.delaunay}
    exec(s, dictionary)
    retval = dictionary['XXXeval_or_exec_outputXXX']
    return retval

def python_filter(name, time, txt, genome, ave, MAT_FILE, used_codons_a, LOAD,
                  generation, Fitnesses, DEBUG = False, FINAL=False):
    """Converts text into indented python code"""
    counter = 0
    if txt == None:
        log_error("None", "no program generated")
        return 0
    for char in txt:
        if char == "{":
            counter += 1
        elif char == "}":
            counter -= 1
        tabstr = "\n" + "    " * counter
        if char == "{" or char == "}":
            txt = txt.replace(char, tabstr, 1)
    txt = "\n".join([line for line in txt.split("\n")
                     if line.strip() != ""])
    if DEBUG or FINAL:
        if FINAL:
            blame = "EliteResults/" + str(time) + ".py"
        else:
            blame = "whatagrammar.py"
        temp = open("temp", "w")
        temp.write('SHOW = False\ndef run():\n')
        temp.write('    import analyser, evolver, subprocess, graph, operator')
        temp.write(', geometry, random, os\n')
        temp.write('    import matplotlib.delaunay as triang\n')
        temp.write('    from operator import itemgetter\n')
        temp.write('    from math import sqrt\n\n    ')
        temp.write('if os.path.isdir("/home/michael/Dropbox/Collij/Mike/')
        temp.write('truss/slf"):')
        temp.write('\n        pass\n    else:\n        ')
        temp.write('os.mkdir("/home/michael/Dropbox/Collij/Mike/truss/slf")\n')
        temp.write('    genome = ' + str(genome) + '\n')
        temp.write('    ave = ' + str(ave) + '\n')
        temp.write('    MAT_FILE = \"' + str(MAT_FILE) + '\"\n')
        temp.write('    LOAD = ' + str(LOAD) + '\n')
        temp.write(txt)
        temp.write("\n    testGraph = mutant()\n    ")
        temp.write("mats = evolver.assign_size(MAT_FILE)\n    ")
        temp.write("analyser = analyser.Analyser(" + str(name))
        temp.write(", 'test','testGraph', genome, " + str(used_codons_a))
        temp.write(", testGraph[1], " + str(ave) + ", mats, " + str(generation))
        temp.write(", LOAD, MAT_FILE)\n")
        temp.write("    analyser.my_graph=testGraph[0]")
        temp.write("\n    grammar_type = testGraph[2]\n    if grammar_type == ")
        temp.write("\"cant\":\n        analyser.cantilever = True\n    ")
        temp.write("if grammar_type == \"truss\":\n        ")
        temp.write("analyser.truss = True\n    fitness = analyser.run_graph(")
        temp.write(str(Fitnesses) + ", " + str(OPTIMIZE) + ", ")
        temp.write(str(GENOME_REWRITE) + ", SHOW)\n    return fitness[0]\n\n")
        temp.write("if __name__ == '__main__':\n    run()")
        temp.close()
        old = open("temp", "r")
        new = open(blame, "w")
        lines = old.readlines()
        for i in range(0, 15):
            new.write(lines[i])
        if evolver.GRAMMAR_FILE == "grammars/Delaunay_cantilever.bnf":
            for i in range(15, 128):
                new.write("    " + str(lines[i]))
            for i in range(128, len(lines)):
                new.write(lines[i])
        elif evolver.GRAMMAR_FILE == "grammars/Delaunay.bnf":
            for i in range(15, 128):
                new.write("    " + str(lines[i]))
            for i in range(128, len(lines)):
                new.write(lines[i])
        elif evolver.GRAMMAR_FILE == "grammars/Delaunay_cantilever_test.bnf":
            for i in range(15, 111):
                new.write("    " + str(lines[i]))
            for i in range(111, len(lines)):
                new.write(lines[i])
        else:
            for i in range(15, 111):
                new.write("    " + str(lines[i]))
            for i in range(111, len(lines)):
                new.write(lines[i])
        old.close()
        remove("temp")
    return txt

def log_error(phenotype, msg):
    """any problems dumped to err.log"""
    print "logging error(" + str(time.clock) + "):", msg
    errFile = open('err.log', 'a')
    errFile.write("error(" + str(time.clock) + "):" + msg + "\n")
    errFile.write(str(phenotype) + "\n")
    errFile.close()

class Analyser():
    """this class execs chromosomes and generates an slffea mesh.
    It is then analysed by slffea and the result is processed to
    generate a fitness value"""

    def __init__(self, arse, unique_id, program, genome_b, used_codons_a,
                 used_codons_b, ave, materials, generation, LOAD,
                 MATERIALS_FILE):
        """stores all the slffea values"""
        self.generation = generation
        self.name = arse
        self.unique_id = unique_id
        self.good = False
        self.used_codons_a = used_codons_a
        self.used_codons_b = used_codons_b
        self.genome_b = genome_b
        self.previous_gen_ave = ave
        self.my_graph = None
        self.program = program
        self.span = 0
        self.UDL_points = False
        self.UDL = False
        self.point = True
        self.truss = False
        self.cantilever = False
        self.total_UDL = 480000
        self.point_load = LOAD
        self.point_load_2 = 0 # 222400 #
        self.length_a = 0
        self.length_b = 0
        self.MATERIALS_FILE = MATERIALS_FILE
        self.beams = materials
        self.fitness_selections = []
        self.node_list = []
        self.edge_list = []
        self.opt_edge_list = []
        self.fixed_list = []
        self.nodeselfloads = []
        self.point_load_nodes = []
        self.corner_load_nodes = []
        self.load_elems = []
        self.new_node_list = []
        self.truss_type = None
        self.UDL_per_m = 0
        self.depth = 0
        self.range = 0

    def create_graph(self, time, program, genome, LOAD, Fitnesses,
                     DEBUG = False, FINAL=False):
        """execute program to create the graph """
        answer = eval_or_exec(self.name, time, program, genome,
                              self.previous_gen_ave, self.MATERIALS_FILE,
                              self.used_codons_a, LOAD, self.generation,
                              Fitnesses, DEBUG, FINAL)
        self.my_graph = answer[0]
        self.used_codons_b = answer[1]
        grammar_type = answer[2]
        if grammar_type == "cant":
            self.cantilever = True
        if grammar_type == "truss":
            self.truss = True
        self.parse_graph()

    def parse_graph(self):
        """gather nodes, edges and other information from output graph"""
        self.edge_list = []
        self.node_list = []
        answer = self.my_graph.return_graph_info()
        self.truss_type = answer[0]
        self.span = answer[1]
        self.UDL_per_m = self.total_UDL/self.span
        self.depth = answer[2]
        total_mass = 0
        for node in self.my_graph.nodes():
            xyz = self.my_graph.get_node_data(node)
            label = self.my_graph.node[node]['label']
            node = {'id': str(node), 'x': xyz[0], 'y': xyz[1],
                    'z': xyz[2], 'label': label}
            self.node_list.append(node)
        for idx, edge in enumerate(self.my_graph.edges_iter(data=True)):
            material = edge[2]['material']
            genome_id = edge[2]['genome_id']
            node = self.node_list[int(edge[0])]
            node1 = node
            label_1 = node['label']
            bode = self.node_list[int(edge[1])]
            node2 = bode
            label_2 = bode['label']
            if label_1 == label_2:
                name = label_1
            elif label_1 == "corner" and label_2 == "load" and node['z'] == bode['z']:
                name = 'load'
            elif label_2 == "corner" and label_1 == "load" and node['z'] == bode['z']:
                name = 'load'
            elif label_1 == "fixed" and label_2 == "load" and node['z'] == bode['z']:
                name = 'load'
            elif label_2 == "fixed" and label_1 == "load" and node['z'] == bode['z']:
                name = 'load'
            else:
                name = "crossbrace"
            length = three_d_line_length(node1, node2) # millimeters
            beam = self.beams[int(material)]
            area = beam['area']
            diameter = beam['diameter']
            thickness = beam['thickness']
            unitweight = beam['unitweight']
            I = beam['I']
            emod = beam['emod']
            density = beam['density']
            mass = length *  float(beam['unitweight']) # answer is in kg #
            total_mass = mass+total_mass
            max_c_s = self.get_max_c_s(diameter, area, length, beam)
            edge = {'id':idx, 'pt_a':int(edge[0]), 'pt_b':int(edge[1]),
                        'material':int(material), 'length':float(length),
                        'mass':float(mass), 'area':str(area),
                        'label':name, 'diameter':float(diameter),
                        'thickness':float(thickness),
                        'unitweight':float(unitweight),
                        'I':float(I), 'emod':float(emod),
                        'density':float(density),
                        'max_c_s': float(max_c_s),
                        'genome_id':int(genome_id)}
            self.edge_list.append(edge)
        return self.edge_list

    def get_max_c_s(self, diameter, area, length, beam):
        """returns the maximum allowable compressive stress
        for a given section"""
        if self.MATERIALS_FILE == "CHSTables":
            if diameter > 270:
                lengths = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
            else:
                lengths = [0, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10]
        else:
            lengths = [0, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10]
        while len(lengths) != 1:
            first = lengths.pop(0)
            second = lengths[0]
            if (1000*float(first)) < float(length) <= (1000*float(second)):
                max_c_s = 1000*float(beam[str(second)])/float(area)
                break
            else:
                max_c_s = 1000*float(beam[str(second)])/float(area)
        return max_c_s

    def create_mesh(self, name='indiv'):
        """produces mesh for medit"""
        if name == 'indiv':
            filename = "population/indiv." + str(self.name) + ".mesh"
        else:
            filename = "population/"+name + '.mesh'
        mesh = open(filename, 'w')
        mesh.write("MeshVersionFormatted 1\nDimension\n3 \n")
        mesh.write("Vertices\n" + str(len(self.node_list)) + " \n")
        for node in self.node_list:
            mesh.write(str(node['x']) + " " + str(node['y'])
                       + " " + str(node['z']) + " 0  \n")
        mesh.write("Edges\n" + str(len(self.edge_list)) + " \n")
        for edge in self.edge_list:
            pt_a, pt_b = int(edge['pt_a']), int(edge['pt_b'])
            mesh.write(str(pt_a + 1) + " " + str(pt_b + 1) + " 0 \n")
        mesh.write("End\n")
        mesh.close()
        if evolver.MAKE_GIF: # Makes an image of the individual
            self.mesh_to_pic(name) 
            filename = "Pics/"+str(name)+".ppm"
            copyfile("population/"+str(name)+".ppm", filename)

    def mesh_to_pic(self, name):
        """run the picture generating program linuxMedit, you must have
        compiled linuxMedit and added it to /usr/local/bin"""
        cmd = 'linuxMedit'
        process = Popen(cmd, shell=True, stdout=PIPE,
                                   stdin=PIPE)
        process.communicate("population/"+str(name)+".mesh")

    def save_dxf(self, name):
        """outputs nodes and edges in dxf format for other software"""
        filename = "CompletedRuns/" + str(name) + ".dxf"       
        DXF = file(filename, 'w')
        DXF.write('  0\n')
        DXF.write('SECTION\n')
        DXF.write('  2\n')
        DXF.write('ENTITIES\n')
        for edge in self.edge_list:
            node = self.node_list[int(edge['pt_a'])]
            X1, Y1, Z1 = node['x'], node['y'], node['z']
            bode = self.node_list[int(edge['pt_b'])]
            X2, Y2, Z2 = bode['x'], bode['y'], bode['z']
            DXF.write('  0\n')
            DXF.write('LINE\n')
            DXF.write('  8\n')
            DXF.write('Polygon\n')
            DXF.write(' 10\n')
            DXF.write(str(X1))
            DXF.write('\n 20\n')
            DXF.write(str(Y1))
            DXF.write('\n 30\n')
            DXF.write(str(Z1))
            DXF.write('\n 11\n')
            DXF.write(str(X2))
            DXF.write('\n 21\n')
            DXF.write(str(Y2))
            DXF.write('\n 31\n')
            DXF.write(str(Z2))
            DXF.write('\n')
        DXF.write('  0\n')
        DXF.write('ENDSEC\n')
        DXF.write('  0\n')
        DXF.write('EOF\n')
        DXF.close()

    def apply_stresses(self, edges):
        """build fixed points list and load list"""
        dnodeselfloads = []
        self.fixed_list = []
        self.nodeselfloads = []
        self.point_load_nodes = []
        self.corner_load_nodes = []
        self.load_elems = []
        self.range = 0
        for node in self.node_list:
            if  node['label'] == 'fixed':
                self.fixed_list.append([node, int(node['y'])])
            elif  node['label'] == 'corner/fixed':
                self.fixed_list.append([node, int(node['y'])])
                self.corner_load_nodes.append(int(node['id']))
                self.range = self.range + 1
            elif node['label'] == 'load':
                self.point_load_nodes.append(int(node['id']))
                self.range = self.range + 1
            elif node['label'] == 'corner':
                self.corner_load_nodes.append(int(node['id']))
                self.range = self.range + 1
        #SLFFEA applies load to edges, find edges connecting load_nodes
        for edge in edges:
            if edge['label'] == 'load':
                self.load_elems.append(edge['id'])
            pt_a, pt_b = int(edge['pt_a']), int(edge['pt_b'])
            self.add_self_loads(edge, pt_a, pt_b, dnodeselfloads)
        dnodeselfloads.sort(key=itemgetter(0))
        if dnodeselfloads:
            self.my_graph.collapse(dnodeselfloads, self.nodeselfloads)
        else:
            print "ERROR!!************NO DNODESELFLOADS THING!****************"
    
    def add_self_loads(self, edge, a, b, nodes):
        """SLFFEA doesn't consider the mass of the element;
        we have to compute this ourselves and add it as a point
        load to the nodes at each end of the element.
        Load per node is in newtons - to add self loading, 
        remove the zero. To remove self loading, comment out the 'load'."""
        load = 0 # float(edge['mass']) / 2
        loadA = [a, load]
        loadB = [b, load]
        nodes.append(loadA)
        nodes.append(loadB)
        return nodes

    def create_slf_file(self, edges):
        """outputs an slf file in truss format"""
        mesh = open("slf/"+str(self.name), 'w')
        mesh.write('numel numnp nmat nmode  (This is for a truss)\n')
        mesh.write(str(len(edges))+'\t'+str(len(self.node_list))
                    + '\t'+str(len(self.beams)) + '\t0\n')
        mesh.write('matl no., E modulus, density, and Area\n')
        for i in range(len(self.beams)):
            mesh.write(str(i)+'\t'+str(self.beams[i]['emod'])+'\t'
                    + str(self.beams[i]['density'])+'\t'
                    +str(self.beams[i]['area'])+ '\n')
        mesh.write('el no., connectivity, matl no\n')
        for i, edge in enumerate(edges):
            mesh.write(str(i)+'\t'+str(edge['pt_a'])+'\t'+str(edge['pt_b'])
                    + '\t'+str(edge['material'])+'\n')
        mesh.write('node no., coordinates\n')
        for node in self.node_list:
            mesh.write(node['id']+'\t'+str(node['x'])+'\t'+str(node['y'])+'\t'
                    +str(node['z'])+"\n")
        mesh.write('prescribed displacement x: node  disp value\n')
        for node in self.fixed_list:           
            mesh.write(node[0]['id']+"\t0.0\n")
        mesh.write('-10\nprescribed displacement y: node  disp value\n')
        for node in self.fixed_list:
           # if node[1] < 0: # comment when dealing with fixed-fixed structures
            mesh.write(node[0]['id']+"\t0.0\n")
        mesh.write('-10\nprescribed displacement z: node  disp value\n')
        for node in self.node_list:
            mesh.write(node['id']+"\t0.0\n")
        mesh.write('-10\nnode with point load and load vector in x, y, z\n')          
        if self.point:
            for node in self.nodeselfloads:
                if node[0] in self.point_load_nodes: 
                    node[1] = node[1] + self.point_load
                    mesh.write(str(node[0])+'\t0\t-'+str(round(node[1], 5))
                        +'\t0\n')
        mesh.write('-10\nelement with stress and tensile stress vector\n-10')
        mesh.close()

    def test_slf_file(self):
        """run the structural analysis software, you must have
        compiled slffea and added ts and tspost to /usr/local/bin"""
        cmd = 'ts'
        process = Popen(cmd, shell=True, stdout=PIPE, stdin=PIPE)
        process.communicate("slf/"+str(self.name))
            
    def show_analysis(self):
        """use tspost to show stresses"""
        cmd = "echo slf/" + str(self.name) + '.ots | tspost'
        process = Popen(cmd, shell=True, stdout=PIPE, stdin=PIPE)
        process.communicate("slf/" + str(self.name))

    def parse_results(self, edges):
        """ read the results of slffea output, check if it calculated
        the results correctly"""
        header1 = ("element no. with stress and tensile stress vector")
        header3 = ("node no., coordinates")
        self.length_a = len(edges)
        results = open("slf/" + str(self.name) + '.ots', 'r')
        # opens the results file
        lines = iter(results)
        self.new_node_list = []
        for line in lines:
            if line.startswith(header3):
                line = lines.next()
                while line[0] != ('prescribed'):
                    result = line.split()
                    if result[0] == ('prescribed'):
                        break
                    # find the new node locations
                    idx, x = int(result[0]), float(result[1])
                    y, z = float(result[2]), float(result[3])
                    node = {'id': idx, 'x': x, 'y': y, 'z': z}
                    self.new_node_list.append(node)
                    line = lines.next()
            if line.startswith(header1):
                line = lines.next()
                while  line.strip() != ('-10'):
                    result = line.split()
                    if 'nan' in result or '-nan' in result:
                        # analysis has failed :(
                        print "***********STRUCTURE FAILED*************"
                        stress = int(DEFAULT_FIT)
                        for edge in edges:
                            edge['stress'] = stress
                        break
                    else:
                        # analysis success, here are stress results
                        stressid = int(result[0])
                        stress = float(result[1])
                        edges[stressid]['stress'] = stress
                        if REMOVE_UNSTRESSED_EDGES:
                            if stress == 0:
                                edges.pop(stressid)
                        line = lines.next()
        results.close()
        self.length_b = len(edges) # how much of Ch.B is used

    def delete_all_files(self):
        """Delete all SLFFEA analysis files"""
        remove("slf/" + str(self.name))
        remove("slf/" + str(self.name) + '.ots')

    def run_optimization(self, steps, Fitnesses, previous, PRINT = False):
        """Optimizes the member sizes for a given structure, based on material
        stress. Reassassigns member sizings, then re-analyses the structure.
        Once re-sized, member stresses change and the structure will need
        to be optimized again. This loop continues for a specified number
        of steps, or until a better fitness is achieved."""
        edges = []
        edges = self.edge_list # create a copy of the edge list
        
        # First optimization pass
        self.reassign_materials_quickly(self.edge_list)
        self.apply_stresses(self.opt_edge_list)
        self.create_slf_file(self.opt_edge_list)
        self.test_slf_file()
        self.parse_results(self.opt_edge_list)
        answer = self.calculate_fitness(Fitnesses, self.opt_edge_list, PRINT)
        if PRINT:
            print "\nOptimization step 1 complete :", answer
            print "___________________________________________\n"
        if answer[0] < previous[0]:
            edges = self.opt_edge_list
            previous = answer
        n = 2
        if BEST_STOP: # Stop iterating once a better fitness is achieved
            while (answer[0] > previous[0]):
                self.reassign_size()
                answer = self.calculate_fitness(Fitnesses, self.opt_edge_list,
                                                PRINT)
                if PRINT:
                    print "\nOptimization step", n, "complete:", answer
                    print "___________________________________________\n"
                if answer[0] < previous[0]: # found a better solution
                    edges = self.opt_edge_list
                    previous = answer
                n += 1
                if n == STEPS:
                    break
        else: # keep going for a specified number of steps.
            while n < steps:
                self.reassign_size()
                answer = self.calculate_fitness(Fitnesses, self.opt_edge_list,
                                                PRINT)
                if PRINT:
                    print "\nOptimization step", n, "complete:", answer
                    print "___________________________________________\n"
                if answer[0] < previous[0]: # we've found a better solution
                    edges = self.opt_edge_list
                    previous = answer
                n += 1
        if GENOME_REWRITE: # re-writes the member sizing genome
            if PRINT:
                print "Before:", self.genome_b[:self.used_codons_b]
            for edge in edges:
                self.genome_b[edge['genome_id']] = edge['material']
            if PRINT:
                print "After:", self.genome_b[:self.used_codons_b]
        self.apply_stresses(edges)
        self.create_slf_file(edges)
        self.test_slf_file()
        self.parse_results(edges)
        return edges

    def reassign_size(self):
        """reassigns materials, then re-analyses the structure"""
        self.reassign_materials_quickly(self.opt_edge_list)
        self.apply_stresses(self.opt_edge_list)
        self.create_slf_file(self.opt_edge_list)
        self.test_slf_file()
        self.parse_results(self.opt_edge_list)

    def reassign_materials_quickly(self, edges):
        """takes a structure and calculates the optimum cross-sectional
        area for each member based on the actual stress in the member
        and its capacity. Reassigns that member to the nearest available
        section."""
        wedge_list = []
        for edge in edges:
            dredge = {}
            stress = edge['stress']
            if self.MATERIALS_FILE == "CSSTables" or self.MATERIALS_FILE == "test1":
                max_allowable_tensile_stress = 172.375
            else:
                max_allowable_tensile_stress = 215
            if stress > 0:
                allowable = max_allowable_tensile_stress
            else:
                allowable = edge['max_c_s']
            original_area = float(edge['area'])
            force = float(abs(stress))*original_area
            required_area = force/allowable
            if required_area > self.beams[156]['area']:
                beam = self.beams[156]
                dredge['id'] = edge['id']
                dredge['stress'] = edge['stress']
                dredge['genome_id'] = edge['genome_id']
                dredge['pt_a'] = edge['pt_a']
                dredge['pt_b'] = edge['pt_b']
                dredge['length'] = edge['length']
                dredge['label'] = edge['label']
                dredge['material'] = 156
                dredge['area'] = beam['area']
                dredge['unitweight'] = beam['unitweight']
                dredge['density'] = beam['density']
                dredge['diameter'] = beam['diameter']
                dredge['thickness'] = beam['thickness']
                dredge['emod'] = edge['emod']
                dredge['I'] = beam['I']
                dredge['max_c_s'] = self.get_max_c_s(dredge['diameter'],
                    dredge['area'], dredge['length'], beam)
                dredge['mass'] = float(dredge['length']) * float(dredge['unitweight'])
                # answer is in Newtons
                wedge_list.append(dredge)
            elif required_area < self.beams[0]['area']:
                beam = self.beams[0]
                dredge['id'] = edge['id']
                dredge['stress'] = edge['stress']
                dredge['genome_id'] = edge['genome_id']
                dredge['pt_a'] = edge['pt_a']
                dredge['pt_b'] = edge['pt_b']
                dredge['length'] = edge['length']
                dredge['label'] = edge['label']
                dredge['material'] = 0
                dredge['area'] = beam['area']
                dredge['unitweight'] = beam['unitweight']
                dredge['density'] = beam['density']
                dredge['diameter'] = beam['diameter']
                dredge['thickness'] = beam['thickness']
                dredge['emod'] = edge['emod']
                dredge['I'] = beam['I']
                dredge['max_c_s'] = self.get_max_c_s(dredge['diameter'],
                                        dredge['area'], dredge['length'], beam)
                dredge['mass'] = float(dredge['length']) * float(dredge['unitweight'])
                # answer is in Newtons
                wedge_list.append(dredge)
            else:
                for a, beam in enumerate(self.beams):
                    if a > 0:
                        if self.beams[a-1]['area'] < required_area < self.beams[a]['area']:
                            beam = self.beams[a]
                            dredge['id'] = edge['id']
                            dredge['stress'] = edge['stress']
                            dredge['genome_id'] = edge['genome_id']
                            dredge['pt_a'] = edge['pt_a']
                            dredge['pt_b'] = edge['pt_b']
                            dredge['length'] = edge['length']
                            dredge['label'] = edge['label']
                            dredge['material'] = a
                            dredge['area'] = beam['area']
                            dredge['unitweight'] = beam['unitweight']
                            dredge['density'] = beam['density']
                            dredge['diameter'] = beam['diameter']
                            dredge['thickness'] = beam['thickness']
                            dredge['emod'] = edge['emod']
                            dredge['I'] = beam['I']
                            dredge['max_c_s'] = self.get_max_c_s(dredge['diameter'], dredge['area'], dredge['length'], beam)
                            dredge['mass'] = float(dredge['length']) * float(dredge['unitweight'])
                            # answer is in Newtons
                            wedge_list.append(dredge)
                            break
        self.opt_edge_list = wedge_list

    def normalise(self, value, limit):
        # takes in a value and a limit and returns a percentage
        # of how far over that limit the value is
        normalised = ((float(value) - float(limit))/float(limit))*100.0
        return normalised
   
    def calculate_fitness(self, Fitnesses, edges, PRINT = False):
        """return values for length and average stress on the structure"""
        total_weight = 0
        if self.truss:
            max_allowable_displacement = self.span/250
        if self.cantilever:
            max_allowable_displacement = 50.8 # self.span/180 # 
        max_displacement = 0 # maximum actual displacement
        max_tensile_stress = 0 # maximum actual tensile stress
        max_comp_stress = 0 # maximum actual compressive stress
        max_norm_tens = 0 # maximum normalised tension
        max_norm_comp = 0 # maximum normalised compression
        total_cost = 0 # structure cost
        failure_counter = 0 # how many constraints have failed?
        cum_difference = 0 # the cumulative difference by which all 
                           # failed constraints have failed
        self.good = False # Need to set the fitness to bad before analysys
        
        # Set tensile stress limits
        if self.MATERIALS_FILE == "CSSTables" or self.MATERIALS_FILE == "test1":
            max_allowable_tensile_stress = 172.375 # 344.75 #
        else:
            max_allowable_tensile_stress = 215
        
        # Check 3 dimensional deflection of all nodes, find maximum
        for node in self.node_list:
            point = self.new_node_list[int(node['id'])]
            displacement = three_d_line_length(node, point)
            if displacement > max_displacement:
                max_displacement = displacement
            if displacement > max_allowable_displacement:
                failure_counter += 1
                diff = self.normalise(displacement, max_allowable_displacement)
                cum_difference = diff + cum_difference
                if PRINT:
                    print "Node", node['id'], " Fails in displacement by",
                    print displacement -max_allowable_displacement, "mm"
        
        # STEEL DESIGN TO BS 449 (TAKEN FROM STRUCTURAL ENGINEERS POCKETBOOK)
        # S355 STEEL!
        for edge in edges:
            total_weight = (total_weight + float(edge['mass']))
            
            # TENSILE LIMITS OF THE MATERIAL
            if edge['stress'] > 0: # edge in tension
                if edge['stress'] > max_tensile_stress:
                    max_tensile_stress = edge['stress']
                if edge['stress'] > max_allowable_tensile_stress:
                    failure_counter += 1
                    diff = self.normalise(edge['stress'],
                           max_allowable_tensile_stress)
                    if diff > max_norm_tens:
                        max_norm_tens = diff
                    cum_difference = diff + cum_difference
                    if PRINT:
                        print "Element", edge['id'], "(material",
                        print edge['material'], ") fails in tension by", diff,
                        print "%"
            
            else: # member is in compression
                # EULER BUCKLING CHECK
                press = abs(float(edge['stress']))
                if BUCKLING_CHECK:
                    limit = (pi**2)*(float(edge['emod']))*(float(edge['I']))/(float(edge['length'])**2)
                    if press > float(limit):
                        failure_counter += 1
                        diff = self.normalise(press, limit)
                        cum_difference = diff + cum_difference
                        if PRINT:
                            print "Element", edge['id'], "fails in buckling"
                
                # COMPRESSIVE LIMITS OF THE MATERIAL (UNITS IN KN)
                if self.MATERIALS_FILE == "CSSTables" or self.MATERIALS_FILE == "test1":
                    max_allowable_comp_stress = 172.375#344.75#
                else:
                    max_allowable_comp_stress = float(edge['max_c_s'])
                if press > max_comp_stress:
                    max_comp_stress = press
                if press > max_allowable_comp_stress:
                    failure_counter += 1
                    diff = self.normalise(press, max_allowable_comp_stress)
                    if diff > max_norm_comp:
                        max_norm_comp = diff
                    cum_difference = diff + cum_difference
                    if PRINT:
                        print "Element", edge['id'], "(material",
                        print edge['material'], ") fails in compression by",
                        print diff, "%"
        
        # Total structure weight in Imperial units
        total_weight = total_weight*2.20462
        max_norm_disp = 0
        
        # Compare all against limits
        if max_displacement > max_allowable_displacement:
            max_norm_disp = self.normalise(max_displacement,
                        max_allowable_displacement)
            if PRINT:
                print "Maximum Actual Displacement (normalised) =",
                print max_norm_disp, "% above limit"
        if PRINT:
            if max_norm_tens > 0:
                print "Maximum % Tensile Stress over the limit:",
                print max_norm_tens, "% above limit"
            if max_norm_comp > 0:
                print "Maximum % Compressive Stress over the limit:",
                print max_norm_comp, "% above limit"
        all_norm_fail = [["max_norm_disp", max_norm_disp],
            ["max_norm_tens", max_norm_tens], ["max_norm_comp", max_norm_comp]]
        all_norm_fail.sort(key=itemgetter(1))
        
        # Deal with unfit individuals
        if failure_counter >= 1:
            if PRINT:
                print "Total failed constraints:", failure_counter
                print "Cumulative % difference between limit and failed "
                print "constraints:", cum_difference, "% above limit"
                print "Worst failed constriaint:", all_norm_fail[-1]
            if DEATH_PENALTY:
                max_displacement = DEFAULT_FIT
                total_weight = DEFAULT_FIT
                total_cost = DEFAULT_FIT
            # Assigns a fitness which is a multiple of the greatest
            # normalised constraint violation (previously difference)
            else:
                max_displacement = 100 * (cum_difference +1) * max_displacement
                total_weight = 100 * (cum_difference +1) * total_weight
                total_cost = 100 * (cum_difference +1) * total_cost
        
        # what if we have a fit individual?
        elif failure_counter == 0:
            self.good = True
        
        # Returns all fitness values
        self.fitness_selections = [total_weight, "Self weight (kg)",
            max_displacement, "Maximum deflection (mm)", cum_difference,
            "Cumulative difference between limit and failed constraints:",
            total_cost,
            "the total estimated cost of construction of the structure"]
        final_fitnesses = []
        for i in Fitnesses:
            final_fitnesses.append(self.fitness_selections[i]) 
        return final_fitnesses, self.good

    def test_mesh(self, Fitnesses, LOAD, DEBUG = False,
                  FINAL = False, PRINT = False):
        """ calls all the functions in analyser"""
        self.create_graph("time", self.program,
                          self.genome_b, LOAD, Fitnesses, DEBUG)
        self.apply_stresses(self.edge_list)
        self.create_slf_file(self.edge_list)
        self.test_slf_file()
        if SHOW_ANALYSIS:
            self.show_analysis()
        self.parse_results(self.edge_list)
        if self.length_a != self.length_b:
        # can only occur if we're removing unstressed edges
            self.apply_stresses(self.edge_list)
            self.create_slf_file(self.edge_list)
            self.test_slf_file()
            self.parse_results(self.edge_list)
        answers = self.calculate_fitness(Fitnesses, self.edge_list, PRINT)
        if OPTIMIZE:
            if OPT_ALL:
                if FINAL or PRINT:
                    print "Optimizing...\nOriginal fitness:", answers[0], "\n"
                answers = self.calculate_fitness(Fitnesses,
                        self.run_optimization(STEPS, Fitnesses,
                        answers, PRINT), PRINT)
                if FINAL:
                    print "Optimized fitness:", answers
            elif answers[1]:
                if FINAL or PRINT:
                    print "Optimizing...\nOriginal fitness:", answers[0], "\n"
                answers = self.calculate_fitness(Fitnesses,
                        self.run_optimization(STEPS, Fitnesses,
                        answers, PRINT), PRINT)
                if FINAL or PRINT:
                    print "Optimized fitness:", answers
        self.delete_all_files()
        return answers
    
    def run_graph(self, Fitnesses, OPT, RE_GEN, show = False):
        """ calls all the functions in analyser"""
        self.parse_graph()
        self.apply_stresses(self.edge_list)
        self.create_slf_file(self.edge_list)
        self.test_slf_file()
        if show:
            self.show_analysis()
        self.parse_results(self.edge_list)
        answers = self.calculate_fitness(Fitnesses, self.edge_list)
        if (OPT):
            print "Optimizing...\nOriginal fitness:", answers[0]
            print "___________________________________________\n"
            answers = self.calculate_fitness(Fitnesses,
                    self.run_optimization(STEPS,
                    Fitnesses, answers, PRINT = True))
            print "Optimized fitness:", answers
            if show:
                self.show_analysis()
        self.delete_all_files()
        return answers
    
    def show_mesh(self, Fitnesses, LOAD, DEBUG = False):
        """generate and show mesh stresses using tspost"""
        self.name = str(datetime.now())[-15:]
        self.create_graph("time", self.program, self.genome_b,
                          LOAD, Fitnesses, DEBUG)
        self.apply_stresses(self.edge_list)
        self.create_slf_file(self.edge_list)
        self.test_slf_file()
        self.parse_results(self.edge_list)
        if self.length_a != self.length_b:
            self.apply_stresses(self.edge_list)
            self.create_slf_file(self.edge_list)
            self.test_slf_file()
            self.parse_results(self.edge_list)
        self.show_analysis()
        self.delete_all_files()
