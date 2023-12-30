import math

YSCALE = math.sqrt(3)/2

class eisenstein:
    '''
    This class implements the eisenstein integers, along with overloading arithmetic operators.
    '''

    def convert(x):
        if isinstance(x, eisenstein):
            return x
        elif isinstance(x, int):
            return eisenstein(x,0)
        else:
            raise TypeError("unsupported conversion type: '{}'".format(type(x)))

    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def conjugate(self):
        return eisenstein(self.x + self.y,-self.y)

    def __neg__(self):
        return eisenstein(-self.x,-self.y)
    
    def __add__(self, other):
        if isinstance(other, eisenstein):
            return eisenstein(self.x + other.x, self.y + other.y)
        elif isinstance(other, int):
            return eisenstein(self.x + other, self.y)
        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))
    
    def __sub__(self, other):
        return self + (-other)
    
    def __radd__(self, other):
        return self + other
    
    def __rsub__(self, other):
        return (-self) + other
    
    def __mul__(self, other):
        if isinstance(other, eisenstein):
            return eisenstein(self.x * other.x - self.y * other.y, self.x * other.y + self.y * other.x + self.y * other.y)
        elif isinstance(other, int):
            return eisenstein(self.x * other, self.y * other)
        elif isinstance(other, float):
            return eisenstein(self.x * other, self.y * other)
        else:
            raise TypeError("unsupported operand type(s) for *: '{}' and '{}'".format(self.__class__, type(other)))
    
    def __rmul__(self, other):
        return self * other

    def __floordiv__(self, other): #other assumed to be an int
        return eisenstein(self.x // other, self.y // other)

    def __rtruediv__(self, other):
        selfconj = self.conjugate()
        lensqrd = (self * selfconj).x
        if isinstance(other, eisenstein) or isinstance(other, int):
            return (other * selfconj) / lensqrd
        else:
            raise TypeError("unsupported operand type(s) for /: '{}' and '{}'".format(type(other), self.__class__))
    
    def __truediv__(self, other):
        if isinstance(other, eisenstein):
            return other.__rtruediv__(self)
        if isinstance(other, int):
            if self.x % other == 0 and self.y % other == 0:
                return self // other
            else:
                return eisenstein(self.x / other, self.y / other) # this is really ugly, but I'm only using it in one spot, so fingers crossed it doesn't cause problems
    
    def __str__(self) -> str:
        if self.y == 0:
            return str(self.x)
        elif self.x == 0:
            return "{}w".format(self.y)
        if self.y > 0:
            return "{} + {}w".format(self.x, self.y)
        else:
            return "{} - {}w".format(self.x, -self.y)

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other):
        if isinstance(other, int):
            return self == eisenstein.convert(other)
        elif isinstance(other, eisenstein):
            return self.x == other.x and self.y == other.y
        return False

    def __hash__(self) -> int:
        return hash((self.x, self.y))
    
    #Note that this is the field norm, i.e. the magnitude squared
    def norm(self) -> int:
        return (self.x * self.x) + (self.x * self.y) + (self.y * self.y)
    
    def tocomplex(self) -> tuple:
        cx = self.x + (self.y / 2)
        cy = self.y * YSCALE
        return (cx, cy)


class EisensteinMatrix:
    '''
    Implements Matrices with eisenstein coefficients. 
    '''

    def __init__(self, p, q, s, t):
        self.p = eisenstein.convert(p)
        self.q = eisenstein.convert(q)
        self.s = eisenstein.convert(s)
        self.t = eisenstein.convert(t)

    def __mul__(self, other):
        if isinstance(other, eisenstein) or isinstance(other, int):
            return EisensteinMatrix(self.p * other, self.q * other, self.s * other, self.t * other)
        return EisensteinMatrix(self.p * other.p + self.q * other.s, self.p * other.q + self.q * other.t,
                                self.s * other.p + self.t * other.s, self.s * other.q + self.t * other.t)
    
    def __rmul__(self, other):
        if isinstance(other, eisenstein) or isinstance(other, int):
            return self * other
    
    def inverse(self):
        det = self.p*self.t - self.q*self.s
        return EisensteinMatrix(self.t / det, -self.q / det, -self.s / det, self.p / det)
    
    def det(self):
        return self.p*self.t - self.q*self.s
    
    def __repr__(self) -> str:
        return 'EisensteinMatrix({},{},{},{})'.format(self.p, self.q, self.s,self.t)
    
    def __hash__(self) -> int:
        return hash((self.p,self.q,self.s,self.t))
    
    def __eq__(self, o) -> bool:
        return isinstance(o, EisensteinMatrix) and o.p == self.p and o.q == self.q and o.s == self.s and o.t == self.t
    
    def canonical(self):
        factor = 1
        for val in [self.p, self.q, self.s, self.t]:
            if val != 0:
                if val.x < 0:
                    factor = -1
                    break
                elif val.x > 0:
                    factor = 1
                    break
                if val.x == 0 and val.y < 0:
                    factor = -1
                    break
                elif val.x == 0 and val.y > 0:
                    factor = 1
                    break
        return self * factor

