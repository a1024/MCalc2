MCALC2: A Matrix Calculator

Usage:
  mcalc2 ["expression"/"filename"]

Syntax:
  Use Matlab-style vectors & matrices.
  Elements are separated by commas or spaces, for example:
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
  \     is for square matrices (eg: A\B)
  .* ./ are element-wise operations

Keywords:
  help: Print this info (works in cmd)
  clear: Clears all variables from memory
  vars: Shows all variables & answer count
  open: Choose a text file to open
General functions:
  ans(n): The n-th answer
  cmd(w, h): set console buffer size (in characters)
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
  det: determinant of square matrix
  inv: inverse of square matrix
  tr: trace of square matrix
  lu: LU factorization of square matrix
  diag: make a diagonal matrix
  egval: Eigenvalues of a matrix
  egvec(M, L): Eigenvectors, given matrix M and eigenvalues L

