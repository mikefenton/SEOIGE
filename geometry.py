""" A module for calculating geometrical stuff, like distance from a
point to a line, whether a point is inside a polygon, and so on.
Methods in this class that a point [x,y,z] or pointList and return a
new point or pointList
Copyright (c) 2010 Jonathan Byrne, Michael Fenton, Erik Hemberg and James McDermott
Hereby licensed under the GNU GPL v3."""


from math import sin, cos, sqrt, pi, fabs, asin, acos, radians
import random


def mirror(pts, axis):
    """reflects points through x, y or z axis"""
    retval = list()
    for pt in pts:
        if axis == "x":
            inverse1, inverse2, inverse3 = - pt[0], pt[1], pt[2]
        if axis == "y":
            inverse1, inverse2, inverse3 = pt[0], - pt[1], pt[2]
        if axis == "z":
            inverse1, inverse2, inverse3 = pt[0], pt[1], - pt[2]
        inverse = (inverse1, inverse2, inverse3)
        retval.append(inverse)
    pts.reverse()
    for pt in pts:
        retval.append(pt)
    return retval

def check(number):
    if number%2==0:
        return "even"
    else:
        return "odd"

def check_point_on_xyline(a_xyz, b_xyz, xyz):
    if (a_xyz[0] <= xyz[0] and xyz[0] <= b_xyz[0]):
        if(a_xyz[1] <= xyz[1] and xyz[1] <= b_xyz[1]):
            return True
    else:
        return False

def midpoint(pt_a, pt_b):
    """calculate midpoint"""
    x = (pt_b[0] + pt_a[0]) / 2
    y = (pt_b[1] + pt_a[1]) / 2
    z = (pt_b[2] + pt_a[2]) / 2
    mid = (x, y, z)
    return mid

def offset(pts, off, axis):
    """moves points by by given off in x, y or z directions"""
    retval = list()
    for pt in pts:
        if axis == "x":
            inverse1, inverse2, inverse3 = pt[0] + off, pt[1], pt[2]
        if axis == "y":
            inverse1, inverse2, inverse3 = pt[0], pt[1] + off, pt[2]
        if axis == "z":
            inverse1, inverse2, inverse3 = pt[0], pt[1], pt[2] + off
        inverse = (inverse1, inverse2, inverse3)
        retval.append(inverse)
    return retval

def create_scaling(size_list):
    total = 0
    cumulative = 0
    ratio_list = []
    scaling = [0]

    for i in range(len(size_list)): size_list[i] = float(size_list[i])
    for size in size_list: total += size
    for size in size_list: ratio_list.append(float(1 / total) * size)
    for ratio in ratio_list:
        cumulative += ratio
        scaling.append(round(cumulative,3))
    return scaling

def scaled_list(min_range, max_range, size_list, reverse = False):
    diff = max_range - min_range
    scaled = []
    scaling = create_scaling(size_list)
    for scale in scaling:
        scaled.append(min_range + (diff * scale))
    if reverse:
        scaled.reverse()
    return scaled

def xy_rotate(point, degree_val):
    """rotate a point on the xy plane"""
    ang = radians(degree_val)
    rotate_x = (point[0] * cos(ang)) + (point[1] * sin(ang))
    rotate_y = (- point[0] * sin(ang)) + (point[1] * cos(ang))
    rotated = (rotate_x, rotate_y, point[2])
    return rotated
    
def xy_rotate_points(points, degree_val):
    "rotates list of points around xy axis"
    rotated_list = []
    for point in points:
        rotated = xy_rotate(point, degree_val)
        rotated_list.append(rotated)
    return rotated_list

def interpolate(p, xy):
    """given 2 points (in tuple format) and pval between 0 and 1
    return coords at p ratio"""
    p = 1 - p
    x, y = xy
    x0, y0, z0 = x
    x1, y1, z1 = y
    return [x0 * p + x1 * (1 - p),
            y0 * p + y1 * (1 - p),
            z0 * p + z1 * (1 - p)]

def subdivide_line(pt_a, pt_b, segments):
    point_list = []
    segs = segments + 1
    x0, y0, z0 = pt_a[0], pt_a[1], pt_a[2]
    x1, y1, z1 = pt_b[0], pt_b[1], pt_b[2]
    offset = ((x1 - x0) / segs, (y1 - y0) / segs,
              (z1 - z0) / segs)
    for i in range(0, segs+1):
        point = (x0 + (offset[0] * i),
                 y0 + (offset[1] * i),
                 z0 + (offset[2] * i))
        point_list.append(point)
    point_list.append(pt_b)
    return point_list

def bezier_form(t, p):
    p0, p1, p2, p3 = p
    p0x, p0y, p0z = p0
    p1x, p1y, p1z = p1
    p2x, p2y, p2z = p2
    p3x, p3y, p3z = p3
    return [((1 - t)**3 * p0x + 3 * (1 - t)**2 * t * p1x
             + 3 * (1 - t) * t ** 2 * p2x + t ** 3 * p3x),
            ((1 - t) ** 3 * p0y + 3 * (1 - t) ** 2 * t * p1y
             + 3 * (1 - t) * t ** 2 * p2y + t ** 3 * p3y),
            ((1 - t) ** 3 * p0z + 3 * (1 - t) ** 2 * t * p1z
             + 3 * (1 - t) * t ** 2 * p2z + t ** 3 * p3z)]


def pt_plus_pt(pt0, pt1):
    """get the dot-sum of pt0 and pt1"""
    retval = dot_operation(pt0, pt1, lambda x, y: x + y)
    #print "pt0", pt0, "pt1", pt1, "=", retval
    return retval

def sinusoid(pt_a, pt_b, height, length):
    """generates a sinusoid between two points, returns a list"""
    return_val = list()
    for i in range(pt_a[0], pt_b[0]):
        point = i, pt_a[1], sin(i / float(length)) * height
        return_val.append(point)
    return return_val


def square(origin, length):
    "creates a square on the xy axis given a starting point"
    point_list = []
    point_list.append(origin)
    point_list.extend(offset(point_list, length, "x"))
    point_list.extend(offset(point_list, length, "y"))
    return point_list


def rectangle(origin, width, length):
    "creates a rectangle on the xy axis given a starting point"
    point_list = []
    point_list.append(origin)
    point_list.extend(offset(point_list, width, "x"))
    point_list.extend(offset(point_list, length, "y"))
    return point_list


def cosine(pt_a, pt_b, height, length):
    """generates a sinusoid between two points, returns a list"""
    return_val = list()
    for i in range(pt_a[0], pt_b[0]):
        point = i, pt_a[1], cos(i / float(length)) * height
        return_val.append(point)
    return return_val


def shift_value(pts):
    """Shifts all z values up until theyare all positive, returns a list"""
    min_val = 0
    return_val = list()
    for pt in pts:
        if pt[2] < min_val:
            min_val = pt[2]
    if min_val < 0:
        for pt in pts:
            shift_point = pt[0], pt[1], pt[2] - min_val
            return_val.append(shift_point)
        return return_val
    else:
        return pts


def invert(pts, height):
    """flips a curve upside down for the bridge"""
    return_val = list()
    for pt in pts:
        inverted = pt[0], pt[1], height - pt[2]
        return_val.append(inverted)
    return return_val


def abs_value(pts):
    """takes a list and reflects any negative vals through the z axis"""
    return_val = list()
    for pt in pts:
        if pt[2] < 0:
            abs_point = pt[0], pt[1], 0 - pt[2]
            return_val.append(abs_point)
        else:
            return_val.append(pt)
    return return_val


def offset_list(pts, off):
    """This method takes a list of points (tuples) and an off tuple
    and returns a new list of points (tuples) off by the values
    given in the off tuple."""
    new_list = []
    for pt in pts:
        new_tup = (pt[0] + off[0], pt[1] + off[1], pt[2] + off[2])
        new_list.append(new_tup)
    return new_list


def distance(p, q):
    """returns euclidean distance between two points. The value is
    rounded to an integer"""
    dist = int(round(sqrt(sum([(p[i] - q[i]) ** 2
                             for i in range(len(p))])), 0))
    return dist


def euclidean_distance(p, q):
    """regular old euclidean distance"""
    return sqrt(sum([(p[i] - q[i]) ** 2 for i in range(len(p))]))


def spiral(t, radius, initial_phase, revs, fn):
    """fn is the "carrier function". fn(t) gives a point "in the centre of"
    the spiral. spiral() returns a single point on the outside of the
    spiral depending on the time parameter t. That is, call this
    multiple times with different values of t (same values for
    everything else) and you'll generate a spiral around the carrier
    curve."""
    epsilon = 0.001
    if t <= 1.0 - epsilon:
        pt0 = fn(t)
        pt1 = fn(t + epsilon)
    else:
        pt0 = fn(t - epsilon)
        pt1 = fn(t)
    phase = initial_phase + 2 * pi * t * revs
    x = disk_at_pt0_perp_to_line_pt0pt1(phase, pt0, pt1, radius)
    return x


def disk_at_pt0_perp_to_line_pt0pt1(theta, pt0, pt1, r):
    """calculating a disk in the plane perpendicular to the line
    between two points. from:
    http://local.wasp.uwa.edu.au/~pbourke/geometry/disk/"""
    rnorm, snorm = get_orthonormal_vectors(pt0, pt1)

    Qx = pt0[0] + r * cos(theta) * rnorm[0] + r * sin(theta) * snorm[0]
    Qy = pt0[1] + r * cos(theta) * rnorm[1] + r * sin(theta) * snorm[1]
    Qz = pt0[2] + r * cos(theta) * rnorm[2] + r * sin(theta) * snorm[2]

    return (Qx, Qy, Qz)


def square_at_pt0_perp_to_line_pt0pt1(pt0, pt1, side):
    """get the four corners of a square *centered* at pt0
    such that the square lies in the plane perpendicular
    to the line from pt0 to pt1"""
    return rect_at_pt0_perp_to_line_pt0pt1(pt0, pt1, side, side)


def rect_at_pt0_perp_to_line_pt0pt1(pt0, pt1, side1, side2):
    """get the four corners of a rectangle *centered* at pt0
    such that the rectangle lies in the plane perpendicular
    to the line from pt0 to pt1. Rectangle has sides of given size.
    """
    half1 = side1 / 2.0
    half2 = side2 / 2.0
    rnorm, snorm = get_orthonormal_vectors(pt0, pt1)

    # we return these points in a strange order, perhaps,
    # to satisfy makeBoard() in render.py
    return [
        pt_minus_pt(pt_minus_pt(pt0, scale(rnorm, half1)),
                    scale(snorm, half2)),
        pt_minus_pt(pt_plus_pt(pt0, scale(rnorm, half1)),
                    scale(snorm, half2)),
        pt_plus_pt(pt_minus_pt(pt0, scale(rnorm, half1)),
                   scale(snorm, half2)),
        pt_plus_pt(pt_plus_pt(pt0, scale(rnorm, half1)),
                   scale(snorm, half2)),
        ]


def get_orthonormal_vectors(pt0, pt1):
    """given two points, find two unit vectors which are orthogonal to
    the line between them. we use a couple of hacks to make sure
    they're nicely aligned to axes."""
    if sum(map(abs, pt_minus_pt(pt0, pt1))) < 0.000001:
        # Two points are coincident: any two orthogonal unit vectors
        # will do.
        return ((1, 0, 0), (0, 1, 0))

    # which way around should the results be? Want the vector which
    # will be used for short side of beam's cross-section to have the
    # zero z-component. Swap them around if necessary.
    swap_results = False
    if pt0[2] == pt1[2]:
        # in this condition, our heuristic below won't work. use a
        # special case. Arbitrary distant point.
        P = [pt1[0] + 100, pt1[1] + 57, pt0[2]]
        swap_results = False
    elif pt0[0] == pt1[0] and pt0[1] == pt1[1]:
        # same x and y values, differ only in z. another special case.
        return ((1, 0, 0), (0, 1, 0))
    else:
        # Set P to be pt1, dropped perpendicularly to the plane of pt0.
        # When we take R as being orthogonal to (pt1-pt0) and to (P-pt0),
        # the result is that R has a zero z-component. This means that the
        # beam's cross-section won't be canted.
        swap_results = True
        P = [pt1[0], pt1[1], pt0[2]]

    while True:
        #Calculate R and S as cross-products => they're orthogonal to the line
        R = cross_product(pt_minus_pt(P, pt0), pt_minus_pt(pt1, pt0))
        S = cross_product(R, pt_minus_pt(pt1, pt0))

        try:
            # Then normalise them. If R is zero this gives a
            # ZeroDivisionError. That happens if we're unlucky when
            # choosing P. So alter P. (Then we won't have the nice
            # axis-alignment)
            rnorm = normalised(R)
            snorm = normalised(S)
            break
        except ZeroDivisionError:
            print "ERR: Divide by zero. pt0, pt1, P, R, S:", pt0, pt1, P, R, S
            for i in range(3):
                P[i] += random.random()
            print "Now trying pt0, pt1, P, R, S:", pt0, pt1, P, R, S
            print "swap_results == " + str(swap_results)
    if swap_results:
        return snorm, rnorm
    else:
        return rnorm, snorm


def pt_minus_pt(pt0, pt1):
    """get the vector from pt0 to pt1"""
    return dot_operation(pt0, pt1, lambda x, y: x - y)


def dot_operation(pt0, pt1, fn):
    """use this to do (eg) (1, 1, 1) + (4, 4, 6)
    pass in lambda x, y: x + y as the operation."""
    return [fn(x0, x1) for x0, x1 in zip(pt0, pt1)]


def scale(pt, factor):
    """scale a vector by a scaling factor"""
    print "X", [x * factor for x in pt]
    return [x * factor for x in pt]

def scale_points(pts, factor):
    """scale a vector by a scaling factor"""
    for pt in pts:
        print "X", [x * factor for x in pt]
        return [x * factor for x in pt]


def cross_product(pt0, pt1):
    i = pt0[1] * pt1[2] - pt0[2] * pt1[1]
    j = pt0[2] * pt1[0] - pt0[0] * pt1[2]
    k = pt0[0] * pt1[1] - pt0[1] * pt1[0]
    return [i, j, k]


def vector_length(pt):
    return sqrt(sum([x ** 2 for x in pt]))

def three_d_line_length(pt1,pt2):
    """finds the three dimension length of a line
    from pt1 to pt2"""
    n1x,n1y,n1z = pt1['x'],pt1['y'],pt1['z']
    n2x,n2y,n2z = pt2['x'],pt2['y'],pt2['z']
    length = sqrt(((n1x-n2x)**2)+((n1y-n2y)**2)+((n1z-n2z)**2))
    return length

def normalised(pt):
    s = vector_length(pt)
    return [pt[0] / s, pt[1] / s, pt[2] / s]


def inside_polygon(pt, vertices):
    """Check whether a point is inside a polygon. I think the vertices
    have to be given in the anti-clockwise order, so we always want pt
    to be "on the left" of the line between the current vertex and the
    next vertex. This only works for convex polygons. They don't have
    to be regular though."""
    n = len(vertices)
    for i in range(n):
        if is_on_the_right(pt, (vertices[i], vertices[(i + 1) % n])):
            return False
    return True


def is_on_the_right(pt, vertices):
    """"On the right" means as we go from vertex 0 to vertex 1.  if
    the line is vertical, then pt is "on the right" in the v0->v1
    direction if 2nd vertex has greater y-value than the 1st and pt is
    actually "on the right" (ignoring direction). pt is "on the right"
    in the v0->v1 sense if it's actually on the left and v0 is above
    v1. Argh it's hard!"""
    if vertices[0][0] == vertices[1][0]:
        if vertices[1][1] > vertices[0][1]:
            return pt[0] > vertices[0][0]
        else:
            return pt[0] < vertices[0][0]
    # if the 2nd vertex has a smaller x-value than the 1st, then pt is
    # "on the right" *iff* it's *above* the line.
    above = above_line(pt, vertices)
    if (vertices[1][0] < vertices[0][0]):
        return above
    else:
        return not above


def above_line(pt, vertices):
    """Is a point (2-tuple) above a line (2-tuple of 2-tuples)?"""
    # if the line is vertical, the function is not defined
    if vertices[0][0] == vertices[1][0]:
        raise ValueError, "above_line() not defined for vertical lines"
    # if the line is horizontal, c won't be defined
    if vertices[0][1] == vertices[1][1]:
        return (pt[1] >= vertices[0][1])
    # otherwise use line equation: y - y1 = m(x - x1)
    # m is slope, c is x-intercept
    # m = (y2 - y1) / (x2 - x1)
    m = (vertices[1][1] - vertices[0][1]) / \
        float(vertices[1][0] - vertices[0][0])
    # deriving from the line equation: when y = 0 (x-intercept),
    # x = x1 - y1 / m
    y_val_at_x = m * (pt[0] - vertices[0][0]) + vertices[0][1]
    return (pt[1] >= y_val_at_x)


def point_line_dist(pt, line):
    """Distance from a point to a line. pt is a 2-tuple, (x0,
    y0). Line is a 2-tuple of two pts ((x1, y1), (x2, y2))."""
    x0, y0 = pt
    x1, y1, x2, y2 = line[0][0], line[0][1], line[1][0], line[1][1]
    if (x2 == x1):
        return fabs(x1 - x0)

    A = (y2 - y1) / (x2 - x1)
    B = -1
    C = y1 - x1 * (y2 - y1) / (x2 - x1)
    top = fabs(x0 * A + y0 * B + C)
    bottom = sqrt(A * A + B * B)
    return top / bottom


def point_line_intersection(pt, line):
    """Distance from a point to a line. pt is a 2-tuple, (x0,
    y0). Line is a 2-tuple of two pts ((x1, y1), (x2, y2))."""
    p1, p2 = line
    if p1 == p2:
        raise ValueError
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = pt

    d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))
    u = ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / (d * d)

    return (x1 + u * (x2 - x1), y1 + u * (y2 - y1))


def triangle_no(angle_val):
    """Given an angle in radians, says which of the six inscribed triangles
    in a hexagon (centred on the origin) the angle brings us to (going
    anticlockwise, counting from the rightmost point)."""
    return (int) (angle_val * 6 / (2 * pi))


def angle(pt):
    """pt is a 2-tuple. Returns an angle in radians."""
    x, y = pt
    try:
        if x >= 0.0 and y >= 0.0:
            return asin(y / sqrt(x * x + y * y))
        elif x < 0.0 and y >= 0.0:
            return pi - asin(y / sqrt(x * x + y * y))
        elif x < 0.0 and y < 0.0:
            return pi + asin(-y / sqrt(x * x + y * y))
        else:
            return 2 * pi - asin(-y / sqrt(x * x + y * y))
    except ZeroDivisionError:
        return 0.0


def get_hex_proportions(hex_pts, pt):
    """Given the vertices of a hexagon plus the centre point, plus a
    point in the hexagon, calculate the proportions (weights)
    corresponding to the vertices."""
    # which triangle is the point in?
    ang = angle(pt)
    tri = triangle_no(ang)

    # distances from pt to the three lines comprising that triangle
    dists = [0.0] * 7
    dists[tri] = point_line_dist(pt, (hex_pts[(tri + 1) % 6], hex_pts[6]))
    dists[(tri + 1) % 6] = point_line_dist(pt, (hex_pts[tri], hex_pts[6]))
    # centre proportion is dists[6]
    dists[6] = point_line_dist(pt, (hex_pts[tri], hex_pts[(tri + 1) % 6]))
    s = sum(dists)
    return [d / s for d in dists]


def get_tri_proportions(tri_pts, pt):
    """Given the vertices of a triangle, and a point within that
    triangle, calculate the proportions (weights) corresponding to the
    vertices."""
    # distance from pt to the three lines
    dists = [point_line_dist(pt, (tri_pts[i], tri_pts[(i + 1) % 3]))
             for i in range(3)]
    # scale the distances so that they add to 1.0
    s = sum(dists)
    return [d / s for d in dists]


def triangle_height(x):
    """Pass in the side of the equilateral triangle and this returns its
    height."""
    return sqrt(pow(x, 2) - pow(x / 2.0, 2))


def inside_unit_n_cube(pt):
    """Given a point in n-space, determine whether it's inside the unit
    n-cube."""
    for pti in pt:
        if pti < 0.0 or pti > 1.0:
            return False
    return True


def L1_distance_to_n_cube_boundary(pt):
    """The L1 distance just takes the minimum of the distances in each
    dimension."""
    return min([min([abs(pti), abs(1.0 - pti)]) for pti in pt])


def projection_to_boundary(cp, dp):
    """Given two 01-list presets (centre prs and a distance prs), find
    the point where the projection from one to the other meets the
    boundary.  There are better ways to do this, but it's not worth
    the effort here."""
    if cp == dp:
        raise ValueError, "Input points were identical, can't project"

    v = [dpi - cpi for cpi, dpi in zip(cp, dp)]

    # The point returned will be short of the boundary by up to this much:
    desired_accuracy = 0.001
    old_pt = dp[:]
    c_accuracy = 1.0
    scale_val = 0.5

    while c_accuracy > desired_accuracy and scale_val > 0.00001:
        print "in while, c_accuracy", c_accuracy, "scale_val", scale_val
        # Make a new point by adding a little bit (scale_val) to the last
        # pt inside the cube.

        new_pt = [opi + vi * scale_val for opi, vi in zip(old_pt, v)]

        # If it's still inside, continue from that point; else
        # decrease the amount to be added.
        if inside_unit_n_cube(new_pt):
            old_pt = new_pt
            c_accuracy = L1_distance_to_n_cube_boundary(old_pt)
        else:
            scale_val /= 2.0

    return old_pt


def generate_random_pt_inside_polygon(poly):
    """Make a new 2-tuple point inside the given polygon. The problem is
    to avoid generating a point outside the polygon. There are clever
    ways of doing it, using the right distribution, but simply looping
    is an ok solution."""

    maxx = max([pt[0] for pt in poly])
    minx = min([pt[0] for pt in poly])
    maxy = max([pt[1] for pt in poly])
    miny = min([pt[1] for pt in poly])

    outside = True
    while outside:
        pt = (random.uniform(minx, maxx), random.uniform(miny, maxy))
        if inside_polygon(pt, poly):
            outside = False
    return pt


def hex_vertices(size=1.0):
    """The vertices of a "unit" hexagon, starting with the right-most
    point and going anti-clockwise (in the direction of positive theta
    in the unit circle (cos(theta), sin(theta)))."""
    return [
        (size, 0.0),
        (0.5 * size, triangle_height(size)),
        (-0.5 * size, triangle_height(size)),
        (-size, 0.0),
        (-0.5 * size, -triangle_height(size)),
        (0.5 * size, -triangle_height(size)),
        ]


def tri_vertices(size=1.0):
    """The vertices of a unit equilateral triangle, anti-clockwise. Starts
    at origin, goes clockwise."""
    return [(0, 0),
            (1.0 * size, 0),
            (0.5 * size, triangle_height(1.0 * size))]


def dot_product(a, b):
    return sum([ai * bi for ai, bi in zip(a, b)])


def mag(a):
    return sqrt(sum([ai * ai for ai in a]))


def angle_between_two_vectors(a, b):
    try:
        # This can also raise a ZeroDivisionError: caller must catch
        # it and decide what to do with it.
        return acos(dot_product(a, b) / (mag(a) * mag(b)))
    except ValueError:
        # Floating-point rounding errors can make the X in
        # acos(X) slightly over 1.0 or slightly under -1.0: this
        # will raise a ValueError. When we catch this, check whether
        # the argument is positive or negative, and return the exact
        # value.
        if dot_product(a, b) > 0.0:
            return acos(1.0)
        else:
            return acos(-1.0)


def angle_between_two_vectors_given_base_pt(base, a, b):
    av = [ai - basei for ai, basei in zip(a, base)]
    bv = [bi - basei for bi, basei in zip(b, base)]
    return angle_between_two_vectors(av, bv)
