from eisenstein import *
import math

def lensqrd(z: eisenstein):
    return z*z.conjugate()

def conjugate(z: eisenstein):
    return z.conjugate()

#Gives the equation coefficients associated with the transformed circle
#a,B,c are circle params, matrix is in PSL(2,C) NB: matrix convention is such that PSL(2,C) ⟳ C is a left action
def transform(a, B: eisenstein, c,  matrix: EisensteinMatrix):
    matrix = matrix.inverse()
    p, q, s, t = matrix.p, matrix.q, matrix.s, matrix.t
    newa = a*lensqrd(p) + B * p * conjugate(s) + conjugate(B * p) * s + c * lensqrd(s) # I can't explain these eqn's here, TODO: a write-up
    newB = a*p*conjugate(q) + B*p*conjugate(t) + s * conjugate(B*q) + c * s * conjugate(t)
    newc = a*lensqrd(q) + B * q * conjugate(t) + conjugate(B*q) * t + c * lensqrd(t)
    return (newa, newB, newc)

#Transforms z by an element matrix of PSL(2,C)
def mobius(z: eisenstein, matrix: EisensteinMatrix):
    if z == 'inf':
        if matrix.s != eisenstein(0,0):
            return matrix.p / matrix.s
        else:
            return 'inf'
    if (matrix.s * z + matrix.t) == eisenstein(0,0):
        return 'inf'
    return (matrix.p * z + matrix.q) / (matrix.s * z + matrix.t)

#Takes the coefficients to an equation a|z|^2 + 2ReBz + c = 0 and returns the pair (center, radius)
def toprform(a,B,c):
    p = -conjugate(B)/a
    r = math.sqrt((lensqrd(B) - a*c).x)/a.x
    return (p,r)

def evaluate(a, B, c, z):
    return (a * lensqrd(z) + B*z + conjugate(B*z) + c).x

def toxy(a: eisenstein):
    return (a.x + a.y / 2, a.y * (math.sqrt(3) / 2))