MCALC2: A Matrix Calculator

Usage:
  mcalc2 ["expression"/"filename"]

Syntax:
  Use Matlab-style vectors & matrices, for example:
    [a2 a1 a0] or [a2, a1, a0] is a row vector
    [a11, a12;  a21, a22]      is a 2x2 matrix
  Use same square brackets to access matrix elements, for example:
    [1 2; 3 4][0][1] == 2
  or use parentheses (index starts at one):
    [1 2; 3 4](1,2) == 2
  Variable names can be up to 16 characters.
  Multiplication asterix '*' is always obligatory.

Operators by precedence (first to last):
  function call, parentheses, square brackets
  ', ^, member access
  + - (unary)
  * o / \  ./
  + -
  :
  < <= > >=
  == !=
  = += -= *= /= \= %=
Notes:
  '     is transpose (eg: M')
  o     is the tensor product
  \     is matrix division from left (first matrix should be square)
  .* ./ are element-wise operations

Keywords:
  help: Prints this info (works in cmd)
  clear: Clears all variables from memory
  vars: Shows all variables & answer count
  open: Opens a text file
  fractions: Toggle printing numbers as fractions
  tolerance: Set algorithm stopping tolerance
General functions:
  ans(n): The n-th answer
  cmd(w, h): Set console buffer size (in characters)
  printmode(n): n=0: Print decimals, n=1: Print fractions
Numbers:
  floor/ceil/round: Element-wise
  min/max/mean/sum: Reduce columns then rows (Matlab-style)
  frac(x, tolerance): Returns [floor(x), num, den]
  rand: Random number in [0, 1] (seed with RDTSC)
Vectors:
  cross(3D vec, 3D vec) -> 3D vec
  cross(2D vec, 2D vec) -> scalar
  For dot product use transpose (col'*col or row*row')
Matrices:
  ref: Row Echelon Form
  rref: Reduced Row Echelon Form (Gaussian Elimination)
  nullspace: Null space of a matrix
Square matrices:
  I(n): n-square identity matrix
  det: Determinant of square matrix
  inv: Inverse of square matrix
  tr: Trace of square matrix
  lu: LU Factorization of square matrix
  diag: Make a diagonal matrix
  egval: Eigenvalues of a matrix
  egvec(M, L): Eigenvectors, given matrix M and eigenvalues L
Polynomials:
  conv: Multiply polynomials
  polpow: Raise a polynomial power an integer
  roots: Find the roots of a polynomial
DSP:
  dft/idft: Discrete Fourier Transform
  dct/idct: Discrete Cosine Transforms II/III

