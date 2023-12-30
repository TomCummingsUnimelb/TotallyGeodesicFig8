from transform import *
from collections import defaultdict
import matplotlib.pyplot as plot
import argparse
import colorsys

from eisenstein import *

w = eisenstein(0,1)
genC = EisensteinMatrix(1,1,0,1)
genD = EisensteinMatrix(1,0,w,1)
genB = EisensteinMatrix(0,-w,1-w,-w-1)

verts = [eisenstein.convert(0), eisenstein.convert(1), w, w-1] #The vertices of the fundamental domain

'''
figureeight.py

Computes the elements of an orbit that intersect F.
Only run it on planes with Delta not an eisenstein norm, otherwise it will stall. Some examples

python figureeight.py 1 0 -2 -show
python figureeight.py 1 w -4 -show
python figureeight.py 1 0 -500 -show (this one might take a while, try zooming in on it though)

usage:
python figureeight.py a b c [-show] [-normal] [-drawsurface [-centre]]
arguments:
    a, b, c: The circle coefficients for the initial plane
    -show: Draws all of the circles in the orbit intersecting F
    -normal: Prints out the normal coordinates
    -drawsurface: draws a fundamental domain for the computed surface
    -centre: tries to centre the drawn fundamental domain by acting on tris and quads by elements of the stabiliser.
'''

class circle:
    def __init__(self, a, b, c, g) -> None:
        self.neighbours = defaultdict(lambda: None)
        if a.x < 0:
            self.coefficients = (-a,-b,-c)
        else:
            self.coefficients = (a,b,c)
        self.groupel = g

    #Returns a list of faces the circle cuts through. Faces are labelled by the fundamental group generator that passes through them.
    def boundaries(self):
        co = self.coefficients
        inset = set()
        for v in verts:
            if evaluate(*co, v) < 0:
                inset.add(v)
        outset = set(verts).difference(inset)
        outset.add("inf")
        faceset = set()
        for a in inset:
            for b in outset:
                pair = set([a,b])
                if pair == set([verts[0],verts[1]]):
                    faceset.add(genD)
                    faceset.add(genB)
                if pair == set([verts[0],verts[2]]):
                    faceset.add(genD.inverse())
                    faceset.add(genB)
                if pair == set([verts[0],verts[3]]):
                    faceset.add(genC.inverse())
                    faceset.add(genD.inverse())
                if pair == set([verts[1],verts[2]]):
                    faceset.add(genC)
                    faceset.add(genB)
                if pair == set([verts[2],verts[3]]):
                    faceset.add(genB.inverse())
                    faceset.add(genD.inverse())
                if pair == set([verts[0], "inf"]):
                    faceset.add(genD)
                    faceset.add(genC.inverse())
                if pair == set([verts[1], "inf"]):
                    faceset.add(genD)
                    faceset.add(genC)
                if pair == set([verts[2], "inf"]):
                    faceset.add(genC)
                    faceset.add(genB.inverse())
                if pair == set([verts[3], "inf"]):
                    faceset.add(genB.inverse())
                    faceset.add(genC.inverse())
        return faceset
    
    #Returns a list of tuples (circle, g), where circle is the neighbouring circle, and g is the generator leading to it
    #Note that my inverses might be mixed up
    def getneighbours(self):
        nset = self.boundaries()
        return [(circle(*transform(*self.coefficients, g.inverse()), g.inverse() * self.groupel), g) for g in nset]
    
    def addneighbour(self, g, n):
        self.neighbours[g] = n

    #Ideally, we shouldn't need to use this one, as the boundary detection should always provide valid neighbours
    def indomain(self):
        for vert in verts:
            if evaluate(*self.coefficients, vert) < 0:
                return True
        return False

#Returns a list of circles, along with a generating set of elements of the surface's fundamental group, i.e. the stabiliser subgroup.
def search(a, b, c):
    searchlist = [circle(a,b,c,EisensteinMatrix(1,0,0,1))]
    #Dictionaries are used to allow me to quickly determine if a circle is already in the list.
    out = defaultdict(lambda: None)
    fund = defaultdict(lambda: None)
    n=0
    while len(searchlist) != 0:
        n+=1
        # if n % 100 == 0:
        #     print("Searched {} circles".format(n))
        val = searchlist.pop(0)
        ns = val.getneighbours()
        for nb, el in ns:
            if out[nb.coefficients] == None:
                out[nb.coefficients] = nb
                searchlist.append(nb)
            else:
                g = nb.groupel.inverse()*out[nb.coefficients].groupel
                g = g.canonical()
                if fund[g] == None:
                    fund[g] = g
            neigh = out[nb.coefficients]
            val.neighbours[el] = neigh
            neigh.neighbours[el.inverse()] = val
    return (out, fund)


mina = 1
n = 0
transformlist = []
#Note: this is the old version. See above for the currently used version
def searchspace(circle):
    searchlist = [EisensteinMatrix(1,0,0,1)]
    map = defaultdict(lambda: None)
    successlist = []
    # transformlist = []
    n=0
    mina = circle[0].x
    while len(searchlist) != 0:
        n+=1
        next : EisensteinMatrix = searchlist.pop(0)
        a,b,c = transform(*circle, next)
        if a.x < 0:
            a = -a
            b = -b
            c = -c
        if map[(a,b,c)] != None:
            continue
        map[(a,b,c)] = "Added"
        flag = False
        for vert in verts:
            if evaluate(a,b,c,vert) < 0:
                flag = True
                break
        if a == eisenstein(0,0):
            flag = False
            for vert in verts:
                if c.x * evaluate(a,b,c,vert) < 0:
                    flag = True
                    break
        if not flag:
            continue
        if a.x < mina:
            mina = a.x
        successlist.append((a,b,c))
        transformlist.append(next)
        for gen in [genB, genC, genD]:
            mat1 : EisensteinMatrix = gen * next
            searchlist.append(mat1)
            mat2 : EisensteinMatrix = gen.inverse() * next
            searchlist.append(mat2)
    return successlist

def lerp(a,b, t):
    return (a[0]*(1-t) + b[0]*t, a[1](1-t), b[1]*t)

def toxy(a: eisenstein):
    return (a.x + a.y / 2, a.y * (math.sqrt(3) / 2))

def hsv2rgb(h,s,v):
    string = "#{:02x}{:02x}{:02x}".format(*tuple(round(abs(i) * 255) for i in colorsys.hsv_to_rgb(h,s,v)))
    return string if len(string) == 7 else "#000000"

def drawcircs(circlist, initial):
    circs = list()
    for a,b,c in circlist:
        if a.x != 0:
            circs.append(toprform(a,b,c))
        else: 
            print('line: {}'.format((b,c)))

    fig, ax = plot.subplots()
        
    #Drawing circles
    divs = 1000
    for circ in circlist:
        if circ[0] == eisenstein(0,0):
            continue
        p, r = toprform(*circ)
        x = p.x + (p.y / 2)
        y = p.y * (math.sqrt(3) / 2)
        xs = [x + r*math.cos(theta) for theta in [2 * math.pi * i / divs for i in range(divs + 1)]]
        ys = [y + r*math.sin(theta) for theta in [2 * math.pi * i / divs for i in range(divs + 1)]]
        # ax.plot(xs, ys, color=hsv2rgb((circ[0]).x / 200, 1.0, 1), linewidth=1)
        ax.plot(xs, ys, color="grey", linewidth=0.5)

    #Drawing straight lines
    for circ in circlist:
        if circ[0] != eisenstein(0,0):
            continue
        p0 = circ[1].conjugate() / (circ[1] * circ[1].conjugate())
        dp = (2*w - 1) * circ[1].conjugate() / 2
        x0, y0 = toxy(p0 + 2*dp)
        x1, y1 = toxy(p0 - 2*dp)
        xs = [x0, x1]
        ys = [y0, y1]
        ax.plot(xs, ys, 'c')
        
    #Drawing the fundamental domain
    for a,b in [(0,1), (0,2), (0,3), (1,2), (2,3)]:
        verta, vertb = toxy(verts[a]), toxy(verts[b])
        xs = [verta[0], vertb[0]]
        ys = [verta[1], vertb[1]]
        ax.plot(xs,ys,'b')
    
    # for i in range(4):
    #     v = verts[i]
    #     alignment = ['top', 'top', 'bottom', 'bottom'][i]
    #     ax.text(*toxy(v), str(v), horizontalalignment='center', verticalalignment=alignment)

    
    ax.set_title("Initial circle: {}".format(initial))

    ax.axis('equal')

    plot.savefig("plot.svg")
    plot.show()

def getnormalsurface(circlist):
    normalsurf = defaultdict(lambda: 0)
    for circ in circlist:
        tets = [[verts[0], verts[1], verts[2]], [verts[0], verts[2], verts[3]]]
        for i in range(2):
            tetrahedron = tets[i]
            inlist = list()
            for j in range(3):
                vert = tetrahedron[j]
                if evaluate(*circ, vert) < 0:
                    inlist.append(j)
            if len(inlist) == 0:
                continue
            if len(inlist) == 1:
                vert = inlist[0]
                normalsurf["{}: {}".format(i, vert)] += 1
                continue
            if len(inlist) == 2:
                complement = list()
                for k in range(4):
                    if k not in inlist:
                        complement.append(k)
                if inlist[0] == 0:
                    normalsurf["{}: {}{}/{}{}".format(i, inlist[0], inlist[1], complement[0], complement[1])] += 1
                else:
                    normalsurf["{}: {}{}/{}{}".format(i, complement[0], complement[1], inlist[0], inlist[1])] += 1
            if len(inlist) == 3:
                vert = inlist[0]
                normalsurf["{}: {}".format(i, 3)] += 1
                continue
    return normalsurf

def planetopoincare(p, r, x, y):
    px, py = toxy(p)
    val = max(r**2 - (px - x)**2 - (py - y)**2, 0)
    s = r / (r + math.sqrt(val))
    return (px + s*(x - px), py + s*(y - py))

def draw(initial, transformlist, generators):
    def lerp(a, b, t):
        return t*b + a*(1-t)
    fig, ax = plot.subplots()
    p, r = toprform(*initial)
    def addplots(tlist, colour):
        componentlist = list()
        for i in range(len(tlist)):
            current = tlist[i].inverse()
            tet1 = (0,1,w,'inf')
            tet2 = (0,w,w-1,'inf')
            componentlist.append(getintersection(tuple(map(lambda z: mobius(z, current), tet1)), initial))
            componentlist.append(getintersection(tuple(map(lambda z: mobius(z, current), tet2)), initial))
        for poly in componentlist:
            # if len(poly) == 3:
            #     continue
            for i in range(len(poly)):
                x1, y1 = poly[i].tocomplex()
                x2, y2 = poly[i-1].tocomplex()
                divs = 10
                points = [planetopoincare(p, r, lerp(x1, x2, t/divs), lerp(y1, y2, t/divs)) for t in range(divs + 1)]
                xs = [x for x, y in points]
                ys = [y for x, y in points]
                ax.plot(xs, ys, color=colour)
    n = 0
    for gen in generators:
        n += 1
        tlist = list(map(lambda x: x*gen, transformlist))
        addplots(tlist, hsv2rgb(n / 30, 0.7, 1))
        break # Removing this break allows the program to draw the neighboring copies of the fundamental domain. Much slower, only feasible for very small surfaces.
    x,y = toxy(p)
    xs = [x + r*math.cos(2 * math.pi * theta / 1000) for theta in range(1001)]
    ys = [y + r*math.sin(2 * math.pi * theta / 1000) for theta in range(1001)]
    ax.plot(xs,ys)
    ax.set_title("Initial circle: {}".format(initial))
    ax.axis('equal')
    plot.savefig("surface_fund.svg")
    plot.show()

def draw2(initial, transformlist, generators):
    def lerp(a, b, t):
        return t*b + a*(1-t)
    fig, ax = plot.subplots()
    p, r = toprform(*initial)
    def addplots(tlist, colour):
        componentlist = list()
        for i in range(len(tlist)):
            tet1 = (0,1,w,'inf')
            tet2 = (0,w,w-1,'inf')
            p1 = None
            p2 = None
            mindist1 = 'inf'
            mindist2 = 'inf'
            for g in generators:
                # print("here")
                current = (tlist[i]*g).inverse()
                tempp1 = getintersection(tuple(map(lambda z: mobius(z, current), tet1)), initial)
                dist1 = 0
                for a in tempp1:
                    dist1 += (p - a).norm()
                if mindist1 == 'inf' or dist1 < mindist1:
                    mindist1 = dist1
                    p1 = tempp1
                tempp2 = getintersection(tuple(map(lambda z: mobius(z, current), tet2)), initial)
                dist2 = 0
                for a in tempp2:
                    dist2 += (p - a).norm()
                if mindist2 == 'inf' or dist2 < mindist2:
                    mindist2 = dist2
                    p2 = tempp2
            componentlist.append(p1)
            componentlist.append(p2)
        for poly in componentlist:
            # if len(poly) == 3:
            #     continue
            for i in range(len(poly)):
                x1, y1 = poly[i].tocomplex()
                x2, y2 = poly[i-1].tocomplex()
                divs = 10
                points = [planetopoincare(p, r, lerp(x1, x2, t/divs), lerp(y1, y2, t/divs)) for t in range(divs + 1)]
                xs = [x for x, y in points]
                ys = [y for x, y in points]
                ax.plot(xs, ys, color=colour)
    addplots(transformlist, hsv2rgb(0, 0.7, 1))
    x,y = toxy(p)
    xs = [x + r*math.cos(2 * math.pi * theta / 1000) for theta in range(1001)]
    ys = [y + r*math.sin(2 * math.pi * theta / 1000) for theta in range(1001)]
    ax.plot(xs,ys)
    ax.set_title("Initial circle: {}".format(initial))
    ax.axis('equal')
    plot.savefig("surface_fund.svg")
    plot.show()


def getintersection(tet, circ):
    result = list()
    inset = set()
    for vert in tet:
        if vert == 'inf':
            continue
        if evaluate(*circ, vert) < 0:
            inset.add(vert)
    outset = set(tet).difference(inset)
    for invert in inset:
        for outvert in outset:
            if outvert == 'inf':
                result.append(invert)
                continue
            vec = ((evaluate(*circ, outvert) * invert) - (evaluate(*circ, invert) * outvert)) * (1 / (evaluate(*circ, outvert) - (evaluate(*circ, invert))))
            if vec == None:
                print("{}".format(((evaluate(*circ, outvert) * invert) - evaluate(*circ, invert) * outvert) / 0.25))
            result.append(vec)
    if len(result) == 4:
        result = [result[0], result[1], result[3], result[2]]
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "Totally Geodesic Surface Computer",
                                    description = "Computes totally geodesic surfaces in the figure-eight knot complement.")
    parser.add_argument('a', default=1)
    parser.add_argument('B', default=0)
    parser.add_argument('c', default=-2)
    parser.add_argument('-show', action='store_true')
    parser.add_argument('-normal', action='store_true')
    parser.add_argument('-drawsurface', action='store_true')
    parser.add_argument('-centre', action='store_true')
    args = parser.parse_args()

    a = eisenstein.convert(eval(args.a))
    b = eisenstein.convert(eval(args.B))
    c = eisenstein.convert(eval(args.c))
    circ = (a, b, c) #initial circle
    circdict, gdict = search(*circ)
    successlist = circdict.keys()

    print("Found {} images:".format(len(successlist)))
    print(successlist)
    print("Minimum value of a achieved: {}".format(mina))
    print("Generators of surface fundamental group:")
    print(gdict.values())
    
    # blist = list()
    # for val in successlist:
    #     if val[0] == 38:
    #         blist.append(-val[1].conjugate()/val[0])
    # fig, ax = plot.subplots()
    # for a,b in [(0,1), (0,2), (0,3), (1,2), (2,3)]:
    #     verta, vertb = toxy(verts[a]), toxy(verts[b])
    #     xs = [verta[0], vertb[0]]
    #     ys = [verta[1], vertb[1]]
    #     ax.plot(xs,ys,'b')
    # bpos = [toxy(b) for b in blist]
    # xs = [pos[0] for pos in bpos]
    # ys = [pos[1] for pos in bpos]
    # ax.scatter(xs, ys)
    # ax.axis('equal')
    # plot.show()

    if args.show:
        drawcircs(successlist, circ)

    if args.normal:
        print()
        #Triangulation Description. This is important if you want to compare with normal surfaces in Regina
        print("Gluing description")
        print("Tetrahedron \t Face (012) \t Face (013) \t Face (023) \t Face (123)")
        print("0 \t \t 1:(312) \t 1:(012) \t 1:(013) \t 1:(023)")
        print("1 \t \t 0:(013) \t 0:(023) \t 0:(123) \t 0:(120)")

        normalsurf = getnormalsurface(successlist)

        print()
        print("Normal surface description:")
        keys = list(normalsurf.keys())
        keys.sort(key=lambda x: "{}{}{}".format(x[0],len(x),x))
        string = ""
        for key in keys:
            string += "{} \t".format(key)
        print(string)
        string = ""
        for key in keys:
            string += "{} \t".format(normalsurf[key])
            if len(key) >= 8:
                string += "\t"
        print(string)

    if args.drawsurface:
        tlist = gdict.values()
        tlist2 = dict()
        for a in tlist:
            for b in tlist:
                tlist2[a*b]=True
        glist = list(tlist2.keys())
        if args.centre:
            draw2(circ, list(map(lambda x: x.groupel, circdict.values())), glist)
        else:
            draw(circ, list(map(lambda x: x.groupel, circdict.values())), glist)

    # mirror = list()
    # for val in successlist:
    #     a,b,c = val
    #     print(transform(a,b,c, EisensteinMatrix(0,w-1,1,0)))
    #     mirror.append(transform(a,b,c, EisensteinMatrix(0,w-1,1,0)))
    # drawcircs(mirror)
    # a,b,c = circ
    # if circ in mirror or (-a,-b,-c) in mirror:
    #     print("True")
            


# for circ in successlist():
#     for tri in [(eisenstein(0), eisenstein(1), w), (eisenstein(0), w, w-1)]:
#         pass


# for mat in transformlist:
#     for tet in [(eisenstein(0), eisenstein(1), w, 'inf'), (eisenstein(0), w, w-1, 'inf')]:
#         t1, t2, t3, t4 = mobius(tet[0], mat.inverse()), mobius(tet[1], mat.inverse()), mobius(tet[2], mat.inverse()), mobius(tet[3], mat.inverse())
#         inlist = []
#         outlist = []
#         for vert in [t1, t2, t3, t4]:
#             if vert != 'inf' and evaluate(circ[0], circ[1], circ[2], vert) < 0:
#                 inlist.append(vert)
#             else:
#                 outlist.append(vert)
#         shape = []
#         for vert1 in inlist:
#             for vert2 in outlist:
#                 if vert2 == 'inf':
#                     shape.append(vert1)
#                 else:
                    
