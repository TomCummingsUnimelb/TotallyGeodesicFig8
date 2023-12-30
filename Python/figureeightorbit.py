from eisenstein import *
from collections import defaultdict
from transform import *

w = eisenstein(0,1)
genC = EisensteinMatrix(1,1,0,1)
genD = EisensteinMatrix(1,0,w,1)
genB = EisensteinMatrix(0,1-w,1,1-2*w)

verts = [eisenstein.convert(0), eisenstein.convert(1), w, w-1]

class orbit:
    def __init__(self, circle, serialized = None) -> None:
        if serialized != None:
            self.deserialise(serialized)
            return
        self.initialcircle = circle
        self.isnormal = True
        searchlist = [EisensteinMatrix(1,0,0,1)]
        map = defaultdict(lambda: None)
        successlist = []
        transformlist = []
        n=0
        while len(searchlist) != 0:     #TODO: Figure out a proper way of simplifying fractions of quadratic integers
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
                self.isnormal = False
                flag = False
                for vert in verts:
                    if c.x * evaluate(a,b,c,vert) < 0:
                        flag = True
                        break
            if not flag:
                continue
            successlist.append((a,b,c))
            transformlist.append(next)
            for gen in [genB, genC, genD]:
                mat1 : EisensteinMatrix = gen * next
                searchlist.append(mat1)
                mat2 : EisensteinMatrix = gen.inverse() * next
                searchlist.append(mat2)
        self.orbitlist = successlist
    
    def mina(self) -> int:
        return min([eisenstein.convert(a[0]).x for a in self.orbitlist])

    def det(self) -> int:
        a,b,c = self.initialcircle
        return eisenstein.convert(b * b.conjugate() - a*c).x

    def normal(self) -> defaultdict():
        normalsurf = defaultdict(lambda: 0)
        circlist = self.orbitlist
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
    
    def triquadratio(self):
        normal = self.normal()
        tris = 0
        quads = 0
        for key in normal.keys():
            if '/' in key:
                quads += normal[key]
            else:
                tris += normal[key]
        return tris/quads
    
    def triquad(self):
        normal = self.normal()
        tris = 0
        quads = 0
        for key in normal.keys():
            if '/' in key:
                quads += normal[key]
            else:
                tris += normal[key]
        return (tris, quads)
    
    def euler(self) -> int:
        normal = self.normal()
        tris = 0
        quads = 0
        for key in normal.keys():
            if '/' in key:
                quads += normal[key]
            else:
                tris += normal[key]
        verts = (3*tris + 4*quads) // 6
        edges = (3*tris + 4*quads) // 2
        faces = tris + quads
        return verts - edges + faces
    
    def serialise(self) -> tuple:
        return (self.initialcircle, self.isnormal, self.orbitlist)
    
    def deserialise(self, string : str):
        init, isnormal, orbitlist = eval(string.replace("w", "*w"))
        self.initialcircle = init
        self.isnormal = isnormal
        self.orbitlist = orbitlist