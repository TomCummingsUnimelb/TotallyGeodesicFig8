include eisenstein
import std/sets
import std/times
import std/tables

const #Face pairing transformations and their inverses for the 4_1 knot complement
    genA : Matrix[Eisenstein] = [[i,i], [o,i]]
    genAinv : Matrix[Eisenstein] = [[i, -i], [o, i]]
    genD : Matrix[Eisenstein] = [[i,o], [w,i]]
    genDinv : Matrix[Eisenstein] = [[i,o],[-w,i]]
    genB : Matrix[Eisenstein] = [[o,-w],[i-w,-i-w]]
    genBinv : Matrix[Eisenstein] = genB.inverse()
    genarray = [genA, genAinv, genD, genDinv, genB, genBinv]

proc eval(m : Matrix[Eisenstein], z : Eisenstein) : int =
    let
        a = m[0,0]
        b = m[1,0]
        b2 = m[0,1]
        c = m[1,1]
    return a.real * z.norm() + (b*z + b2*z.conjugate()).real + c.real

const verts = [o,i,w,w-i]

proc intersects (m : Matrix[Eisenstein]) : bool =
    for vert in verts:
        if m.eval(vert) < 0:
            return true
    return false

# nim's ifdef sucks
when defined(compute_orientation):
    proc search (initcirc : Matrix[Eisenstein]): (seq[Matrix[Eisenstein]],bool) =
        assert initcirc.determinant() < 0
        var #Should this be let? Still learning nim...
            searchlist, donelist : seq[Matrix[Eisenstein]] #searchlist is circles left to be checked, donelist is the output list of circles
            searched = initHashSet[Matrix[Eisenstein]]() #Using a hashset allows us to save computational
                                                         #time on checking whether a circle has been encountered already
            orientation : Table[Matrix[Eisenstein],int]
            orientable = true
            current : Matrix[Eisenstein]
            newcirc: Matrix[Eisenstein]
        searchlist = @[initcirc]
        searched.incl(initcirc)
        orientation[initcirc] = 1
        while searchlist.len() > 0:
            current = searchlist.pop()
            donelist.add(current)
            for gen in genarray:
                newcirc = gen.adjoint() * current * gen #Transforms the current circle by a generator
                if newcirc[0,0].real == 0: #non-normal surface, this code isn't meant to handle this case
                    return (@[],true)
                if newcirc[0,0].real < 0: #Ensure this matrix is reduced
                    newcirc = -newcirc
                    if newcirc in searched:
                        if orientation[newcirc] == orientation[current]:
                            orientable = false
                    else:
                        orientation[newcirc] = -orientation[current]
                else:
                    orientation[newcirc] = orientation[current]
                if not (newcirc in searched):
                    searched.incl(newcirc)
                    if newcirc.intersects():
                        searchlist.add(newcirc)
        return (donelist,orientable)
else:
    proc search (initcirc : Matrix[Eisenstein]): seq[Matrix[Eisenstein]] =
        assert initcirc.determinant() < 0
        var #Should this be let? Still learning nim...
            searchlist, donelist : seq[Matrix[Eisenstein]] #searchlist is circles left to be checked, donelist is the output list of circles
            searched = initHashSet[Matrix[Eisenstein]]() #Using a hashset allows us to save computational
                                                         #time on checking whether a circle has been encountered already
            orientation : Table[Matrix[Eisenstein],int]
            orientable = true
            current : Matrix[Eisenstein]
            newcirc: Matrix[Eisenstein]
        searchlist = @[initcirc]
        searched.incl(initcirc)
        orientation[initcirc] = 1
        while searchlist.len() > 0:
            current = searchlist.pop()
            donelist.add(current)
            for gen in genarray:
                newcirc = gen.adjoint() * current * gen #Transforms the current circle by a generator
                if newcirc[0,0].real == 0: #non-normal surface, this code isn't meant to handle this case
                    return @[]
                if newcirc[0,0].real < 0: #Ensure this matrix is reduced
                    newcirc = -newcirc
                if not (newcirc in searched):
                    searched.incl(newcirc)
                    if newcirc.intersects():
                        searchlist.add(newcirc)
        return donelist

let
    quads = [4,5,6,11,12,13]
    tris = [0,1,2,3,7,8,9,10]
    vert0 = [o, i, w]
    vert1 = [o, w, w-i]
    ind = [[[3,4],[5,0]],[[6,1],[2,7]]]

#Produces the normal coordinates of a surface given a list of planes. Order is
#[0:0,0:1,0:2,0:3,0:01/23,0:02/13,0:03/12,1:0,1:1,1:2,1:3,1:01/23,1:02/13,1:03/12]#
proc normal(list : seq[Matrix[Eisenstein]]) : array[14,int] =
    var normalcoords : array[14,int]
    for mat in list:
        var ind1, ind2, ind3 : int
        ind1 = if mat.eval(vert0[0]) > 0: 1 else: 0
        ind2 = if mat.eval(vert0[1]) > 0: 1 else: 0
        ind3 = if mat.eval(vert0[2]) > 0: 1 else: 0
        var disk = ind[ind1][ind2][ind3]
        if disk != 7:
            normalcoords[disk] += 1
        ind1 = if mat.eval(vert1[0]) > 0: 1 else: 0
        ind2 = if mat.eval(vert1[1]) > 0: 1 else: 0
        ind3 = if mat.eval(vert1[2]) > 0: 1 else: 0
        disk = ind[ind1][ind2][ind3]
        if disk != 7:
            normalcoords[disk + 7] += 1
    return normalcoords

#Produces the Euler characteristic given the normal coordinates
proc euler(normal : array[14,int]) : int =
    var q = 0
    for i in quads:
        q += normal[i]
    return -q div 3

proc mina(list : seq[Matrix[Eisenstein]]) : int =
    assert list.len() > 0
    var min = list[0][0,0].real
    for mat in list:
        min = if mat[0,0].real < min : mat[0,0].real else : min
    return min

# let inittime = cpuTime()
# for num in 1..1000:
#     let mat : Matrix[Eisenstein] = [[i,o], [o,-num*i]]
#     var time = cpuTime()
#     discard search(mat)
#     echo num, " done in ", cpuTime() - time, " seconds"
# echo "Total time: ", cpuTime() - inittime