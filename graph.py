""" creates graphs using methods from geometry class and turtle
graphics. built on the graph class from networkx
node methods:
pass in point[x,y,z], if new point then add to graph, return node_id
edge methods:
pass in node_ids, if no edge between nodes then add edge
Copyright (c) 2010 Jonathan Byrne, Erik Hemberg and James
McDermott Hereby licensed under the GNU GPL v3."""
import networkx as nx
import copy, math, time, os
from geometry import *


class GraphState(object):
    """This class is used for generating graphs using turtle graphics
    methods. It stores the current position, orientation, and node.
    Idea is: there's a "cursor", like a turtle, which moves around,
    rotates, etc. There's also current_node which isn't used for now.
    As in L-systems, we have save and restore functions for state,
    which is stored in a stack."""

    def __init__(self, pos=None, orientation=None):
        if pos is None:
            # x, y, z
            self.position = 0, 0, 0
        else:
            self.position = pos
        if orientation is None:
            # orientation is a vector
            self.orientation = 1, 0, 0
        else:
            self.orientation = orientation
        self.current_node = 0


class graph(nx.Graph):
    """ An extension of the networkx class that we use for
    representing our designs as a graph with xyz coordinates stored in
    the nodes. It was originally built for generating turtle-graphics"""

    def __init__(self, *args, **kwargs):
        super(graph, self).__init__(*args, **kwargs)
        self.states = [GraphState()]
        self.state = self.states[0]
        self.graph_info = None
        
        #You can debug your graphs by setting self.save=True. This
        #generates a graph in the population folder every time a node
        #is added, then movie.py will turn it into an avi
        self.unused_node = []
        self.save = False
        self.pop_folder = os.getcwd() + "/population/"
        self.frame_count = 0
        self.insulators = True
###########Turtle Methods############

    def initial_edge(self, dist):
        """This function ensures that every graph is non-empty."""
        self.add_node(0, xyz=self.state.position)
        self.project(dist)
        self.add_node(1, xyz=self.state.position)
        self.add_edge(0, 1)
        self.save_graph()

    def rotate(self, a, b, c):
        """"change the turtle's orientation. Not actual rotation"""
        self.state.orientation = (self.state.orientation[0] + a,
                                  self.state.orientation[1] + b,
                                  self.state.orientation[2] + c)

    def move(self, x_dist, y_dist, z_dist):
        """ Move the turtle, don't add any edge. """
        self.state.position = (self.state.position[0] + x_dist,
                               self.state.position[1] + y_dist,
                               self.state.position[2] + z_dist)

    def project(self, dist):
        """Move the turtle, taking account of orientation. Don't add
        an edge."""
        self.move(dist * self.state.orientation[0],
                  dist * self.state.orientation[1],
                  dist * self.state.orientation[2])

    def nearest_node_id(self):
        """ What node id is nearest turtle's current position? Unused in
        blender_graph.bnf for now."""
        pos = self.state.position
        dist = lambda k: euclidean_distance(self.node[k]['xyz'], pos)
        return min(self.node, key=dist)

    def incr_current_node(self):
        """New current_node."""
        self.state.current_node += 1
        self.state.current_node %= self.order()

    def add_node_connect2(self, id1, id2):
        """create a node and connect it to two existing nodes"""
        nodeid = self.get_unused_nodeid()
        self.add_node(nodeid, xyz=self.state.position)
        self.add_edge(nodeid, id1)
        self.add_edge(nodeid, id2)
        self.save_graph()

    def save_state(self):
        """save and restore methods for the stack of states."""
        self.states.append(copy.copy(self.state))

    def restore_state(self):
        """pop current state or create graph if None"""
        self.state = self.states.pop()
        if self.state is None:
            self.state = GraphState()
            
    def save_graph_info(self, info):
        """saves any information from a grammar/graph,
        to be returned at a later date"""
        self.graph_info = info
        
    def return_graph_info(self):
        """returns previously saved information"""
        return self.graph_info

######### GRAPH METHODS ########

    def get_unused_nodeid(self):
        """The order of a graph is the number of nodes. Nodes are
        numbered sequentially as we add them (see add_edge), so the
        order is the number of the next free id."""
        return self.order()

    def add_edge(self, node_a, node_b, attr_dict=None, **attr):
        """overriding graph add edge so that it can save graphs"""
        nx.Graph.add_edge(self, node_a, node_b, attr_dict, **attr)
        self.save_graph()

#    @profile
    def add_unique_node(self, coords, node_type):
        """Add a node and connect it to two others, ensuring that it's
        not lonesome"""
        #TODO Assert that all nodes are lists of ints and not empty, when? Catch exceptions        
        new = True
     #   print coords
        x, y, z = coords[0], coords[1], coords[2]
        a_coords = [int(x), int(y), int(z)]
        for node in self.node:
            b_coords = self.node[node]['xyz']
            if a_coords == b_coords:
                new = False
                node_id = node
                break
        if new:
            node_id = self.get_unused_nodeid()
            self.add_node(node_id, xyz=a_coords, label=node_type)
        return node_id

    def collapse(self, original, target):
        """Finds the nodes with the same node[0] (i.e. the same nodes)
        and adds up the individual loads (node[1]) to return the total
        load on that particular node."""
        for i, node in enumerate(original):
            pair = [node[0], node[1]]
            target.append(pair)
            if i > 0:
                penultimate = target.pop(-2)
                last = target.pop()
                if last[0] == penultimate[0]:
                    last[1] = last[1] + penultimate[1]
                    target.append(last)
                elif last[0] == (penultimate[0]+1):
                    target.append(penultimate)
                    target.append(last)

    def connect_by_height(self):
        """connect all nodes in the graph that share the
        same height (z-axis)"""
        for node_a in self.node:
            for node_b in self.node:
                if not self.has_edge(node_a, node_b):
                    pt1 = self.node[node_a]["xyz"]
                    pt2 = self.node[node_b]["xyz"]
                    if pt1[2] == pt2[2]:
                        self.add_edge(node_a, node_b)

    def connect_by_offset_height(self):
        """connect all nodes in the graph that share the
        same height (z-axis)but are offset"""
        for node_a in self.node:
            for node_b in self.node:
                if not self.has_edge(node_a, node_b):
                    pt1 = self.node[node_a]["xyz"]
                    pt2 = self.node[node_b]["xyz"]
                    if pt1[2] == pt2[2]:
                        if not pt1[1] == pt2[1]:
                            if not pt1[0] == pt2[0]:
                                self.add_edge(node_a, node_b)

    def connect_neighbours(self, node_list, length, range=False):
        """connect all nodes on a list  within a given range"""
        for node_a in node_list:
            for node_b in node_list:
                if not self.has_edge(node_a, node_b):
                    pt1 = self.node[node_a]["xyz"]
                    pt2 = self.node[node_b]["xyz"]
                    if range:
                        if int(distance(pt1, pt2)) < length:
                            self.add_edge(node_a, node_b)
                    else:
                        if int(distance(pt1, pt2)) == length:
                            self.add_edge(node_a, node_b)

    def connect_all_neighbours(self, length):
        """connect all nodes in the graph within a given range
        CAUTION: This is an expensive operation on big graphs"""
        for node_a in self.node:
            for node_b in self.node:
                if not self.has_edge(node_a, node_b):
                    pt1 = self.node[node_a]["xyz"]
                    pt2 = self.node[node_b]["xyz"]
                    if int(distance(pt1, pt2)) == length:
                        self.add_edge(node_a, node_b)

    def connect_nodes(self, node_list, x):
        """adds edges between nodes that are adjacent on a list"""
        for i in range(len(node_list) - 1):
            self.add_edge(node_list[i], node_list[i + 1], material=x)

    def connect_sorted_nodes(self, node_list):
        """adds edges between nodes that are sorted in order of creation"""
        node_list.sort()
        for i in range(len(node_list) - 1):
            self.add_edge(node_list[i], node_list[i + 1])

    def add_nodes(self, point_list, node_type, connect=False):
        """adds all coords on a list to the graph"""
        id_list = []
        for point in point_list:
            idx = self.add_unique_node(point, node_type)
            id_list.append(idx)
        if connect:
            self.connect_nodes(id_list, 0)
        return id_list

    def weave_nodes(self, a_list, b_list, offset):
        for i in range(len(a_list)):
            off = (i + offset) % len(b_list)
            self.add_edge(a_list[i], b_list[i])
            self.add_edge(a_list[i], b_list[off])
            self.add_edge(a_list[off], b_list[i])
    
    def sierpinski(self, tri, node_type, depth):      # fractals, baby!
        if depth > 1:
            node_ids = []
            triangles = []
            triangles.append((tri[0],midpoint(tri[0], tri[1]), 
                              midpoint(tri[0], tri[2])))
            triangles.append((tri[1],midpoint(tri[1], tri[0]), 
                              midpoint(tri[1], tri[2])))
            triangles.append((tri[2],midpoint(tri[2], tri[0]), 
                              midpoint(tri[2], tri[1])))
            triangles.append((midpoint(tri[0], tri[1]),
                              midpoint(tri[1], tri[2]), 
                              midpoint(tri[0], tri[2])))
            for idx, triangle in enumerate(triangles):
                ids = self.sierpinski(triangle, node_type, depth-1)
                node_ids.extend(ids)
            return node_ids
        else:
            triangle = [tri[0], tri[1], tri[2], tri[0]]
            ids = self.add_nodes(triangle, node_type, True)
            return ids

    def varinski(self, tri, node_type, depths): # variable sierpinski
        triangles = []
        node_ids = []
        triangles.append((tri[0],midpoint(tri[0], tri[1]), 
                          midpoint(tri[0], tri[2])))
        triangles.append((tri[1],midpoint(tri[1], tri[0]), 
                          midpoint(tri[1], tri[2])))
        triangles.append((tri[2],midpoint(tri[2], tri[0]), 
                          midpoint(tri[2], tri[1])))
        triangles.append((midpoint(tri[0], tri[1]),
                          midpoint(tri[1], tri[2]), 
                          midpoint(tri[0], tri[2])))

        for idx, triangle in enumerate(triangles):
            ids = self.sierpinski(triangle, node_type, depths[idx])
            node_ids.extend(ids)
        return node_ids

    def add_star_to_node(self, node, npts, node_type="none"):
        """Given a node id, and a radius, add a number of new nodes
        and edges from each to the central node."""
        node = self.get_node_idx_mod(node)
        xyz = self.get_node_data(node)
        id_list = []

        for i in range(npts):
            theta = 2 * math.pi * i / npts
            newxyz = (xyz[0] + math.cos(theta), xyz[1] + 3, xyz[2]
                      + math.sin(theta))
            new_nodeid = self.get_unused_nodeid()
            self.add_node(new_nodeid, xyz=newxyz, label=node_type)
            self.add_edge(node, new_nodeid)
            id_list.append(new_nodeid)
        return id_list

    def add_edge_between_existing_nodes(self, x_node, y_node):
        """this function assumes x and y exist."""
        self.add_edge(x_node, y_node)

    def get_node_idx_mod(self, x):
        """Given an integer, return a node id guaranteed to exist."""
        return x % self.order()

    def get_node_idx_mod_exclude_y(self, x, y):
        """Given an integer x, return a node id guaranteed to exist,
        and NOT node id y"""
        tmp = x % (self.order() - 1)
        if tmp >= y:
            tmp += 1
        return tmp

    def get_node_idx_float(self, x):
        """Given a number, return a node id guaranteed to exist."""
        return int(x * self.order())

    def get_node_idx_float_exclude_y(self, x, y):
        """Given a number, return an existing node id, not y"""
        tmp = int(x * (self.order() - 1))
        if tmp >= y:
            tmp += 1
        return tmp

    def get_nodes_with_n_edges(self, n):
        """What nodes have exactly n edges?"""
        return [node for node in self.adj
                if self.degree(node) == n]

    def get_nth_node_with_m_edges_mod(self, n, m):
        """Return a particular node of those with m edges"""
        nodes = self.get_nodes_with_n_edges(m)
        if len(nodes):
            return nodes[n % len(nodes)]
        else:
            return None

    def check_graph(self):
        """Scans graph for duplicates which would break slffea.
        This is an expensive method so its only used for debugging"""
        duplicate_count = 0
        for a_node in self.node:
            for b_node in self.node:
                if a_node != b_node:
                    a_coord = self.node[a_node]["xyz"]
                    b_coord = self.node[b_node]["xyz"]
                    if a_coord == b_coord:
                        duplicate_count += 1
        if duplicate_count > 0:
            print "duplicates!", duplicate_count
        else:
            print "no duplicates"

    def get_node_data(self, nodeid):
        """Get the cartesian coordinates of a node."""
        return self.node[nodeid]["xyz"]

    def get_label(self, nodeid):
        """return label info, used for identifying different parts of
        the structure"""
        return self.node[nodeid]["label"]

    def transpose_xz_coords(self):
        for node_id in self.node:
            xyz = self.node[node_id]['xyz']
            new_xyz = [xyz[2], xyz[1], xyz[0]]
            self.node[node_id]['xyz'] = new_xyz

    def rotate_points_around_xy_plane(self, points, sectors, materials):
        tangent = int(360 / sectors)
        for i in range(1, sectors):
            tangent_angle = tangent * i
            rotated_points = xy_rotate_points(points, tangent_angle)
            rotated_ids = self.add_nodes(rotated_points, "rotate")
            self.connect_nodes(rotated_ids, materials)

    def copy_and_rotate_around_xy_plane(self, original, tangent):
        #tangent = int(360 / sectors)
        orig_copy = original.copy()
        offset_copy = original.copy()

        for nodeid in offset_copy.node:
            xyz = offset_copy.node[nodeid]["xyz"]
            xyz = xy_rotate(xyz, tangent)
            xyz = [int(round(xyz[0], 0)), int(round(xyz[1], 0)),
                   int(round(xyz[2], 0))]
            offset_copy.node[nodeid]["xyz"] = xyz
        new_graph = nx.union(orig_copy, offset_copy, rename=("G-", "H-"))
        new_graph.frame_count = original.frame_count
        return new_graph
    
    def copy_and_rotate_around_yz_plane(self, original, tangent):
        #tangent = int(360 / sectors)
        orig_copy = original.copy()
        offset_copy = original.copy()

        for nodeid in offset_copy.node:
            xyz = offset_copy.node[nodeid]["xyz"]
            xyz = yz_rotate(xyz, tangent)
            xyz = [int(round(xyz[0], 0)), int(round(xyz[1], 0)),
                   int(round(xyz[2], 0))]
            offset_copy.node[nodeid]["xyz"] = xyz
        new_graph = nx.union(orig_copy, offset_copy, rename=("G-", "H-"))
        new_graph.frame_count = original.frame_count
        return new_graph

    def sanitise_pylon(self, original, width):
        new_graph = original.copy()
        new_graph.frame_count = original.frame_count
        excess_nodes = []
        insulator_nodes = []
        for nodeid in new_graph.node:
            xyz =  new_graph.node[nodeid]['xyz']
            if xyz[2] == 0:
                new_graph.node[nodeid]['label'] = 'base'

            if new_graph.node[nodeid]['label'] == ('line'):
                if abs(new_graph.node[nodeid]['xyz'][0]) > width:
                    found = self.traverse_node_type(new_graph,nodeid, 'arm')
                    excess_nodes.extend(found)
                elif self.insulators:
                    insulator_nodes.append(nodeid)
        #add insulators to the pylon
        for nodeid in insulator_nodes:
            xyz = new_graph.node[nodeid]['xyz']
            offset = pt_plus_pt(xyz,(0,0,-4000))
            new_id = new_graph.add_unique_node(offset, "insulator")
            new_graph.add_edge(nodeid, new_id)
        excess_nodes = list(set(excess_nodes))
        for node in excess_nodes:
            new_graph.remove_node(node)        
        arm_ids = []
        mirror_ids = []
        for nodeid in new_graph.node:
            if new_graph.node[nodeid]['label'] == ('arm'):
                arm_ids.append(nodeid)
                xyz = new_graph.node[nodeid]['xyz']
                mirror = [- xyz[0], xyz[1], xyz[2]]
                mirror_id = new_graph.add_unique_node(mirror, 'arm')
                mirror_ids.append(mirror_id)
                new_graph.add_edge(nodeid, mirror_id)
        new_graph = nx.convert_node_labels_to_integers(new_graph)
        return new_graph

    def traverse_node_type(self, new_graph, starting_node, node_type):
        unvisited = [starting_node]
        visited = []
        visited_edges = []
        while len(unvisited) > 0:
            next = unvisited.pop()
            visited.append(next)
            neighbours = new_graph.neighbors(next)
            for neighbour in neighbours:
                #check for self linking
                if not neighbour == next:
                    # only look at certain node type
                    if new_graph.node[neighbour]['label'] == node_type:
                        #check its not in the edgelist
                        if not set((neighbour, next)) in visited_edges:
                            visited_edges.append(set((neighbour, next)))
                            # nodes can have many incoming edges
                            if not neighbour in unvisited:
                                unvisited.append(neighbour)
        return visited

    def copy_and_offset_with_mirror(self, original, offset_val, reflect=False):
        """Add a copy of the graph, offsetting all nodes by a given
        vector. For nodes with the "rung" attribute, add an edge
        between existing node and its offset copy."""
        # make an unchanged copy and an offset/mirrored copy
        orig_copy = original.copy()
        offset_copy = original.copy()
        for nodeid in offset_copy.node:
            # perform an offset
            xyz = offset_copy.node[nodeid]["xyz"]
            xyz = pt_plus_pt(xyz, offset_val)
            if reflect:
                ## also perform a mirror in the y axis
                xyz = [xyz[0], - xyz[1], xyz[2]] # This should be xyz = [xyz[0], - xyz[1], - xyz[2]] 
            offset_copy.node[nodeid]["xyz"] = xyz

        # make a union of the original and copy, renaming nodes
        # note that this requires nx to be updated to svn 1520 or above
        # which fixes a bug where union discards node attributes
        new_graph = nx.union(orig_copy, offset_copy, rename=("G-", "H-"))
        new_graph.frame_count = original.frame_count
        return new_graph

    def replace_graph(self, new_graph):
        """replaces self.graph with sanitised new graph, removes any
        duplicates and only attaches only nodes that are connected to
        edges"""
        new_graph = nx.convert_node_labels_to_integers(new_graph)
        self.clear()

        #only add nodes that are connected to edges
     #   print "In graph"
        for edge in new_graph.edges_iter(data=True):
      #      print "edge:",edge
            if self.unused_node:
                if edge[0] == self.unused_node or edge[1] == self.unused_node:
                    print "GOT HIM!!!"
                    print edge
                    new_graph.remove_edge(edge[0], edge[1])
            elif edge[0] == edge[1]:
                new_graph.remove_edge(edge[0], edge[0])
            else:
                node_a = self.add_unique_node(new_graph.node[edge[0]]["xyz"],
                                              new_graph.node[edge[0]]["label"])
                node_b = self.add_unique_node(new_graph.node[edge[1]]["xyz"],
                                              new_graph.node[edge[1]]["label"])
                self.add_edge(node_a, node_b, edge[2])

    def sort_truss_support_nodes(self, new_graph, span, threed=False):
        """for the structural span and support grammars, makes sure there are reaction nodes"""
        for i, nodeid in enumerate(new_graph.node):
            if new_graph.node[nodeid] == {}:
                print "Deleting vacant node..."
                self.unused_node = i
                new_graph.node.pop(nodeid)
                print self.unused_node
                break   
            xyz =  new_graph.node[nodeid]['xyz']
            if xyz[2] == 0 and abs(xyz[1]) == span/2:
                new_graph.node[nodeid]['label'] = 'fixed'
            elif xyz[2] == 0 and xyz[1] == -span/2:
                new_graph.node[nodeid]['label'] = 'fixed'
        return new_graph

    def sort_cant_support_nodes(self, new_graph, span, depth, threed=False):
        """for the structural cantilever grammars, makes sure there are reaction nodes"""
        for i, nodeid in enumerate(new_graph.node):
            if new_graph.node[nodeid] == {}:
                print "\nDeleting vacant node...\n\n"
                self.unused_node = i
                new_graph.node.pop(nodeid)
                print self.unused_node
                break   
            xyz =  new_graph.node[nodeid]['xyz']
            if xyz[2] == 0 and abs(xyz[1]) == 0:
                new_graph.node[nodeid]['label'] = 'fixed'
            elif xyz[2] == depth and xyz[1] == 0:
                new_graph.node[nodeid]['label'] = 'fixed'
        return new_graph

    def sort_truss_load_nodes(self, new_graph, span, depth, threed=False):
        """for the structural span and support grammars, makes sure there are loaded nodes"""
        for nodeid in new_graph.node:
            if new_graph.node[nodeid] == {}:
                print "Deleting vacant node..."
                new_graph.node.pop(nodeid)
                break
            xyz =  new_graph.node[nodeid]['xyz']
            if xyz[2] == depth and xyz[1] == span/2:
                new_graph.node[nodeid]['label'] = 'corner'
            elif xyz[2] == depth and xyz[1] == -span/2:
                new_graph.node[nodeid]['label'] = 'corner'
         #   elif xyz[2] == 0 and new_graph.node[nodeid]['label'] != "fixed":
         #       new_graph.node[nodeid]['label'] = 'load'
        return new_graph
    
    def sort_cant_load_nodes(self, new_graph, span, depth, threed=False):
        """for the structural span and support grammars, makes sure there are loaded nodes"""
        for nodeid in new_graph.node:
            if new_graph.node[nodeid] == {}:
                print "Deleting vacant node..."
                break
            xyz =  new_graph.node[nodeid]['xyz']
            if xyz[2] == depth and xyz[1] == span:
                new_graph.node[nodeid]['label'] = 'corner'
            elif xyz[2] == depth and xyz[1] == 0:
                new_graph.node[nodeid]['label'] = 'corner/fixed'
            elif xyz[2] == depth:
                new_graph.node[nodeid]['label'] = 'load'
        return new_graph

######## UTILITY METHODS ##########
    def save_picture(self):
        """Save a 2d-visualisation of the graph (ignores 3D
        coordinates). Requires matplotlib"""
        nx.draw_graphviz(self)
        #plt.show()

    def create_mesh(self, name):
        """generates a .mesh file from the graph"""
        node_list = []
        edge_list = []
        for node_id in self.nodes():
            node = {}
            xyz = self.node[node_id]['xyz']
            x, y, z = xyz[0], xyz[1], xyz[2]
            label = self.node[node_id]['label']
            node = {'id': str(node), 'x': x, 'y': y, 'z': z, 'label': label}
            node_list.append(node)
        for idx, edge in enumerate(self.edges_iter()):
            edge = {'id': idx, 'pt_a': edge[0], 'pt_b': edge[1]}
            edge_list.append(edge)

        filename = name + '.mesh'
        mesh = open(filename, 'w')
        mesh.write("MeshVersionFormatted 1\nDimension\n3 \n")
        mesh.write("Vertices\n" + str(len(node_list)) + " \n")
        for node in node_list:
            mesh.write(str(node['x']) + " " + str(node['y'])
                       + " " + str(node['z']) + " 0  \n")
        mesh.write("Edges\n" + str(len(edge_list)) + " \n")
        for edge in edge_list:
            pt_a, pt_b = edge['pt_a'], edge['pt_b']
            if type(pt_a) == int and type(pt_b) == int:
                mesh.write(str(pt_a + 1) + " " + str(pt_b + 1) + " 0 \n")
        mesh.write("End\n")
        mesh.close()

    def save_graph(self, name=None):
        if not name == None:
            self.create_mesh(name)
        elif self.save:
            filename = self.pop_folder + "img%04d" % (self.frame_count)
            self.create_mesh(filename)
            self.frame_count += 1


########## PYLON GRAPH ############
class pylon_graph(graph):
    """Graph for a pylon that contains structures for line points and
    checks the validity of the graph when unique nodes are added."""

    LINE_LABEL = 'line'

    def __init__(self, *args, **kwargs):
        super(pylon_graph, self).__init__(*args, **kwargs)
        self.line_nodes = []
        self.valid = True

    def add_unique_node(self, coords, node_type, load):
        """Add a node and connect it to two others, ensuring that it's
        not lonesome. It is imperative for slffea that there are no
        duplicate nodes, this function has problems with floating
        points so try and set the xyz values as integers. Verifies
        that all the constraints are passed, otherwise an Exception is
        raised."""
        #TODO good import placement?
        import constraints
        new = True
        x, y, z = coords[0], coords[1], coords[2]
        a_coords = [int(round(x, 0)), int(round(y, 0)), int(round(z, 0))]
        #TODO use more memory and save the node_list structure in the
        #class (let the profiler decide)
        node_list = self.analyser_node_structure(self)
        for node in self.node:
            coords_n = self.node[node]['xyz']
            x, y, z = coords_n[0], coords_n[1], coords_n[2]
            b_coords = [int(round(x, 0)), int(round(y, 0)), int(round(z, 0))]
            if a_coords == b_coords:
                new = False
                node_id = node
                break
            #Check constraints
        if new:
            node_id = self.get_unused_nodeid()
            self.add_node(node_id, xyz=a_coords, label=node_type)
            new_node = {'id': str(node_id), 'x': a_coords[0], 'y': a_coords[1],'z': a_coords[2], 'label': node_type, 'load': load}
            if node_type == self.LINE_LABEL:
                self.line_nodes.append(new_node)
                self.valid = constraints.check_line_constraint(self.line_nodes)
                if not self.valid and __debug__:
#                    raise ValueError('Violating line constraint:',self.line_nodes, self.node)
                    print('Violating line constraint:',self.line_nodes, self.node)
                if self.valid and __debug__:
                        #Only need to verify the newly added node
                    self.valid = constraints.check_insulator_constraint([new_node], self)
                if not self.valid and __debug__:
#                    raise ValueError('Violating insulator constraint:',new_node)                
                    print('Violating insulator constraint:',new_node)                
                #Check if the line node violates the nodes
                self.valid = constraints.check_structure_constraint(node_list, [new_node])
                if not self.valid and __debug__:
#                    raise ValueError('Violating structure constraint:',new_node, node_list)
                    print('Violating structure constraint:',new_node, node_list)
            else:
                #Only check if the new node violates the line nodes
                self.valid = constraints.check_structure_constraint([new_node], self.line_nodes)
                if not self.valid and __debug__:
#                    raise ValueError('Violating structure constraint:',[new_node])   
                    print('Violating structure constraint:',[new_node])   
        return node_id

    def analyser_node_structure(self, my_graph):
        """Transfer my_graph nodes to a list of dictionaries, which
        contain the node information."""
        node_list = []
        for node in my_graph.nodes():
            xyz = my_graph.get_node_data(node)
            label = my_graph.node[node]['label']
            node = {'id': str(node), 'x': xyz[0], 'y': xyz[1],
                    'z': xyz[2], 'label': label}
            node_list.append(node)
        return node_list
