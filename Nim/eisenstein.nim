import std/strformat

type
    Eisenstein* = tuple
        real : int
        complex : int

proc `$` (x : Eisenstein) : string =
    if x.real != 0 and x.complex != 0:
        return fmt"{x.real}+{x.complex}*w"
    elif x.complex != 0:
        return fmt"{x.complex}*w"
    else:
        return fmt"{x.real}"

proc `+` (x : Eisenstein, y : Eisenstein) : Eisenstein =
    return (x[0] + y[0], x[1] + y[1])

proc `+` (x : int, y : Eisenstein) : Eisenstein = 
    return (x + y[0], y[1])

proc `+` (y : Eisenstein, x : int) : Eisenstein = 
    return (x + y[0], y[1])

proc `*` (x : Eisenstein, y : Eisenstein) : Eisenstein =
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] + x[1] * y[1])

proc `*` (x : int, y : Eisenstein) : Eisenstein =
    return (x * y[0], x * y[1])

proc `*` (y : Eisenstein, x : int) : Eisenstein =
    return (x * y[0], x * y[1])

proc `-` (x : Eisenstein, y : Eisenstein) : Eisenstein =
    return (x[0] - y[0], x[1] - y[1])

proc `-` (x : Eisenstein) : Eisenstein =
    return (-x[0], -x[1])

proc `-` (x : int, y : Eisenstein) : Eisenstein =
    return (x - y[0], - y[1])

proc `-` (x : Eisenstein, y : int) : Eisenstein =
    return (x[0] - y, x[1])

proc `==` (x : Eisenstein, y : Eisenstein) : bool =
    return x[0] == y[0] and x[1] == y[1]

proc `==` (x : Eisenstein, y : int) : bool =
    return x[0] == y and x[1] == 0

proc `==` (y : int, x : Eisenstein) : bool =
    return x[0] == y and x[1] == 0

proc conjugate(x: Eisenstein) : Eisenstein =
    return (x.real + x.complex, -x.complex)

proc norm(x : Eisenstein) : int =
    return x.real*x.real + x.real*x.complex + x.complex*x.complex



const i : Eisenstein = (1,0)
const w : Eisenstein = (0,1)
const o : Eisenstein = (0,0)

type
    Matrix[T] = array[2, array[2, T]]

const id: Matrix[Eisenstein] = [[i, o], [o, i]]

proc `$` (m : Matrix) : string =
    return &"{m[0,0]}, {m[0,1]}\n{m[1,0]}, {m[1,1]}"

proc `[]` (m : Matrix, a, b : static[int]) : m.T {.inline.} =
    m[a][b]

proc `*` (m : Matrix, n : Matrix ) : Matrix =
    return [[m[0,0]*n[0,0] + m[0,1]*n[1,0], m[0,0]*n[0,1] + m[0,1]*n[1,1]],[m[1,0]*n[0,0] + m[1,1]*n[1,0], m[1,0]*n[0,1] + m[1,1]*n[1,1]]]

proc `-` (m : Matrix) : Matrix =
    return [[-m[0,0], -m[0,1]],[-m[1,0], -m[1,1]]]

proc determinant(m : Matrix) : int =
    return (m[0,0]*m[1,1] - m[1,0]*m[0,1]).real

proc inverse (m : Matrix) : Matrix =
    assert m.determinant() == 1
    return [[m[1,1], -m[0,1]], [-m[1,0], m[0,0]]]

proc adjoint(m: Matrix[Eisenstein]) : Matrix[Eisenstein] =
    return [[m[0,0].conjugate(), m[1,0].conjugate()], [m[0,1].conjugate(), m[1,1].conjugate()]]
