from eisenstein import *
from figureeight import searchspace, getnormalsurface
import  time
from figureeightorbit import *
from matplotlib import pyplot
import csv

'''
census.py

This was my first program for computing a full list of totally geodesic surfaces up to some bound.

It is very inefficient and was superceded by the version made with nim.
'''

def gcdlist(list):
    current = list[0]
    for entry in list:
        current = math.gcd(current, entry)
    return current

maxdet = 500
abound = maxdet
Bbound = math.ceil(math.sqrt(maxdet*2))
cbound = maxdet
arange = [eisenstein(1,0), eisenstein(3,0)] # [eisenstein(x,0) for x in range(1,abound + 1)] # Temporarily assuming mina = 1 or 3
#Centered hexagonal numbers
Brange = [eisenstein(x,y) if abs(x+y) <= Bbound else None for x in range(-Bbound, Bbound+1) for y in range(-Bbound, Bbound+1)]
Brange = list(filter(lambda x: x != None, Brange))
crange = [eisenstein(-x,0) for x in range(1,cbound + 1)]
totalcircs = len(arange) * len(Brange) * len(crange)
starttime = time.time()
orbitlist = list()
savetime = defaultdict(lambda: False)
for a in arange:
    for c in crange:
        shortlist = list()
        for b in Brange:
            if (b.norm() - a*c).x > maxdet:
                continue
            if b == None:
                continue
            if gcdlist([a.x, b.x, b.y, c.x]) > 1:
                continue
            if savetime[(eisenstein.convert(a),b,eisenstein.convert(c))]:
                continue
            current = orbit((a,b,c))
            for circ in current.orbitlist:
                savetime[circ] = True
            # for circ in current.orbitlist:
            #     if circ[0] == 0:
            #         continue
            orbitlist.append(current)

print("Searched {} starting circles and found {} distinct orbits".format(totalcircs, len(orbitlist)))

bset = set(map(lambda c: c.det(), orbitlist))
print("Number of distinct determinants: {}".format(len(bset)))
for orb in orbitlist:
    if orb.mina() % 3 == 0:
        continue
    if orb.det() % 2 == 0 and (1,0,-orb.det()) not in orb.orbitlist:
        print("exception")
    elif (1,w,-orb.det()+1) not in orb.orbitlist and (1,0,-orb.det()) not in orb.orbitlist:
        print("exception")


data1 = [circ.euler() for circ in orbitlist]
data2 = [circ.det() for circ in orbitlist]

fig, ax = pyplot.subplots()
ax.scatter(data2,data1)
pyplot.show()

f = open("outfiledet.txt", "w")
for circ in orbitlist:
    f.write(str(circ.serialise()) + "\n")
f.close()