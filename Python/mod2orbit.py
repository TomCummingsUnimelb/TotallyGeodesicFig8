from eisenstein import *
from transform import *

# This program was used to compute the orbits of hermitian matrices in M(2,F_4)

w = eisenstein(0,1)
def modcirc(a,b,c):
    return (a.x % 2, (b.x % 2) + (b.y % 2)*w, c.x % 2)

genA = EisensteinMatrix(1,1,0,1)
genB = EisensteinMatrix(1,0,w,1)
genC = EisensteinMatrix(0,-w,1-w,-w-1)
gens = [genA, genB, genC]
genname = {genA:"genA", genB:"genB", genC:"genC"}

o = [(0,0,0)]
even = [(1,0,0), (1,1,1), (1,w,1), (1,w+1,1), (0,0,1)]
odd1 = [(1,1,0), (0,w,0), (0,w,1), (1,0,1), (0,1,0)]
odd2 = [(0,1,1), (0,w+1,0), (1,w,0), (0,w+1,1), (1,w+1,0)]

#change this to one of even, odd1, odd2
parset = even

pairs = []
for i in range(len(parset)):
    circ = parset[i]
    for gen in gens:
        result = modcirc(*transform(*circ, gen))
        if result not in parset:
            print("Set is not closed")
        pairs.append((parset[i],genname[gen],result))
for pair in pairs:
    print("{1} sends {0} to {2}".format(*pair))