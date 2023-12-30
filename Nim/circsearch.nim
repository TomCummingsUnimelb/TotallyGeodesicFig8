include orbit
import std/sets
import std/math
import std/times
import std/strutils
import std/tables
import os

assert paramCount() > 0, "Remember to pass bound as an argument"

## circsearch.nim
## 
## This is a newer and faster implementation of the census program.
## 
## To build it, make sure you have nim installed, then run
## 
## nim c circsearch.nim
## 
## If you want to compute the orientability of a surface (a bit slower), run
## 
## nim c -d:compute_orientation
## 
## Then just run
## 
## circsearch.exe bound
## 
## where bound is the maximum magnitude of Delta to search.
## Recommended that you don't increase this bound above 1000 unless you have a lot of ram.
## This folder should come with precompiled binaries for both options, named "circsearch_orientation.exe" and "circsearch_no_orientation.exe"

proc min(a,b : int) : int =
    return if a < b : a else: b

proc max(a,b : int) : int =
    return if a < b : b else: a

proc gcd(a,b : int) : int =
    var
        low = min(a,b)
        high = max(a,b)
        temp : int
    while low != 0:
        temp = low
        low = high mod low
        high = temp
    return high

proc abs(a : int) : int =
    return if a < 0: -a else: a

# Computing up to 2000 used 13GB of memory, 8GB of which was virtual, so it was very slow
let bound = parseInt(paramStr(1))
let num = ceil(sqrt(bound.toFloat())).toInt()

# Precomputes a list of Eisenstein norms
# Note that we're talking about the field norm here
var eisnorms = initHashSet[int]()
for i in -num..num:
    for j in -num..num:
        eisnorms.incl(i*i + i*j + j*j)


let file = open("outfile.csv", fmwrite)
when defined(compute_orientation):
    file.writeline("Determinant,MinA,Euler,Orientable,A,B,C,0:0,0:1,0:2,0:3,0:01/23,0:02/13,0:03/12,1:0,1:1,1:2,1:3,1:01/23,1:02/13,1:03/12")
else:
    file.writeline("Determinant,MinA,Euler,A,B,C,0:0,0:1,0:2,0:3,0:01/23,0:02/13,0:03/12,1:0,1:1,1:2,1:3,1:01/23,1:02/13,1:03/12")

let t1 = cpuTime()
var count = 0
var circlist : seq[Matrix[Eisenstein]]
when defined(compute_orientation):
    var orientable : bool
var searched = initHashSet[Matrix[Eisenstein]]()
var normalset = initTable[array[14,int], (int, int)]()
for a in 1..bound: #ideally this would just be a in [1,3], but I haven't proven this to be the case
    let cbound = (bound div a) + 1
    for c in (-cbound)..(-1):
        if bound + a*c <= 0: #You wouldn't believe how big a difference this super obvious optimisation made
            continue
        let bbound = ceil(sqrt((bound + a*c).float)).toInt()
        for b1 in (-2*bbound)..(2*bbound):
            for b2 in (-2*bbound)..(2*bbound):
                let det = a*c - (b1*b1 + b1*b2 + b2*b2)
                if det < -bound or det >= 0 or -det in eisnorms: #only searching within a range of determinants, avoiding Eisenstein norms
                    continue
                if gcd(a, gcd(abs(b1), gcd(abs(b2), abs(c)))) != 1: #only searching cases where gcd(a,b1,b2,c)=1
                    continue
                let
                    A = a*i
                    B = b1*i+b2*w
                    C = c*i
                    mat : Matrix[Eisenstein] = [[A,B.conjugate()], [B,C]]
                if mat in searched:
                    continue
                count += 1
                if count mod 100 == 0:
                    echo count, " at ", cpuTime() - t1, "s"
                when defined(compute_orientation):
                    (circlist, orientable) = search(mat)
                else:
                    circlist = search(mat)
                for circ in circlist:
                    searched.incl(circ)
                let normalcoords = circlist.normal()
                let euler = normalcoords.euler()
                let mina = circlist.mina()
                
                ## Used to compute duplicate normal coordinates, currently disabled
                # if normalcoords in normalset:
                #     let val = normalset[normalcoords]
                #     echo "Duplicate for det1 ", det, ", mina1 ", mina, ", det2 ", val[0], ", mina2 ", val[1], ", coords ", normalcoords
                # else:
                #     normalset[normalcoords] = (det, mina)

                var str : string
                str = $(det) & ","
                str &= $(mina) & ","
                str &= $(euler) & ","
                when defined(compute_orientation):
                    str &= $(orientable) & ","
                var bstr = $(B)
                bstr = bstr.replace("*")
                str &= $(A) & "," & bstr & "," & $(C)
                for val in normalcoords:
                    str &= "," & $(val)
                file.writeline(str)

file.close()

echo "Time taken: ", cpuTime() - t1, "s"