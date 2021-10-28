#ifdef TOKEN
TOKEN(0, T_IGNORED)

//functions
TOKEN("ans", T_ANS) TOKEN("cmd", T_CMD)//general
TOKEN("conv", T_CONV) TOKEN("roots", T_ROOTS)//polynomials
TOKEN("ldiv", T_LDIV) TOKEN("sample", T_SAMPLE) TOKEN("invz", T_INVZ)//fractions
TOKEN("cross", T_CROSS)//row vectors
TOKEN("I", T_IDEN) TOKEN("sum", T_SUM) TOKEN("ref", T_REF) TOKEN("rref", T_RREF) TOKEN("det", T_DET) TOKEN("inv", T_INV) TOKEN("diag", T_DIAG) TOKEN("lu", T_LU) TOKEN("tr", T_TRACE)//matrices
TOKEN("dft", T_DFT) TOKEN("fft", T_FFT_UNUSED) TOKEN("idft", T_IDFT) TOKEN("ifft", T_IFFT_UNUSED)//matrices & polynomials

//operators
TOKEN("\'", T_TRANSPOSE)//unary postfix
TOKEN("^", T_POWER)
TOKEN("*", T_MUL) TOKEN("o", T_TENSOR) TOKEN("/", T_DIV) TOKEN("\\", T_DIV_BACK) TOKEN("%", T_MOD)
TOKEN(".*", T_MUL_EW) TOKEN("./", T_DIV_EW)//matrix-matrix element-wise operations
TOKEN("+", T_PLUS) TOKEN("-", T_MINUS)//can be unary		binary: if one is scalar then element-wise
//TOKEN(".+", T_ADD_EW) TOKEN(".-", T_SUB_EW)//scalar-matrix element-wise additions
TOKEN("<", T_LESS) TOKEN("<=", T_LESS_EQUAL) TOKEN(">", T_GREATER) TOKEN(">=", T_GREATER_EQUAL)
TOKEN("==", T_EQUAL) TOKEN("!=", T_NOT_EQUAL)
TOKEN("=", T_ASSIGN)

//control
TOKEN(",", T_COMMA) TOKEN(";", T_SEMICOLON)
TOKEN("[", T_LBRACKET) TOKEN("]", T_RBRACKET)
TOKEN("(", T_LPR) TOKEN(")", T_RPR)

//commands
TOKEN("help", T_HELP)
TOKEN("open", T_OPEN)
TOKEN("clear", T_CLEAR) TOKEN("vars", T_VARS)
TOKEN("gfset", T_GFSET)
TOKEN("exit", T_EXIT) TOKEN("quit", T_QUIT)

//constants
TOKEN("i", T_IMAG) TOKEN("j", T_IMAG_UNUSED) TOKEN("e", T_EULER) TOKEN("pi", T_PI) TOKEN("inf", T_INF) TOKEN("nan", T_NAN)

//lexer stuff
TOKEN(0, T_SPACE) TOKEN(0, T_NEWLINE) TOKEN(0, T_EOF)//consequent in this order
TOKEN(0, T_NUMBER)	//check lex_number
TOKEN(0, T_ID)		//check lex_id

//parser stuff
TOKEN(0, T_REAL) TOKEN(0, T_COMPLEX)

TOKEN(0, T_NTOKENS)
#endif