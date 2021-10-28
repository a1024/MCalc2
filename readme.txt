MCalc2: Matrix Calculator 2

Written in C++.

Usage:
  mcalc ["expression"/filename]

Syntax:
  Use Matlab-style vectors & matrices.
  Elements are separated by commas or spaces, for example:
    [a2 a1 a0] or [a2, a1, a0] is a row vector
    [a11, a12;  a21, a22]      is a 2x2 matrix
  For multiline statements:
    Open a parenthesis at the start of the first command, and close it at the
    end of the last command.
  Nested matrices and/or polynomials are not supported yet.
  Variable names can be up to 16 characters.
  Multiplication asterix '*' is always obligatory.

Operators by precedence (first to last):
  ^
  + - (unary)
  '
  / \ % * o - .- + .+
  = += -= *= /= \= %=
  ,
  ; [ ] ( )
Notes:
  '     is transpose (eg: M')
  o     is the tensor product
  \     is for square matrices (eg: A\B)
  .+ .- are element-wise operations

Keywords:
  help: Print this info (works in cmd)
  clear: Clears all variables from memory
  vars: Shows all variables & answer count
General functions:
  ans(n): The n-th latest answer
Vectors:
  cross(3D vec, 3D vec) -> 3D vec
  cross(2D vec, 2D vec) -> scalar
Matrices:
  ref: Row Echelon Form
  rref: Reduced Row Echelon Form (Gaussian Elimination)
Square matrices:
  I n: n-square identity matrix
  det: determinant of square matrix
  inv: inverse of square matrix
  diag: diagonal factorization of square matrix
  lu: LU factorization of square matrix
  tr: trace of square matrix
Fraction objects:
  sample(F): s to z domain
  ldiv(F, S): long division of fraction F by S steps
Utility:
  cmd(w, h): set console buffer size (in characters)

