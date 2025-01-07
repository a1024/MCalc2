#ifdef TOKEN
TOKEN(0, T_IGNORED)

//functions
TOKEN("ans", T_ANS) TOKEN("cmd", T_CMD) //TOKEN("printmode", T_PRINTMODE)//general
TOKEN("conv", T_CONV) TOKEN("polpow", T_POLPOW) TOKEN("roots", T_ROOTS) TOKEN("poly", T_POLY)//polynomials
TOKEN("ldiv", T_LDIV) TOKEN("sample", T_SAMPLE) TOKEN("invz", T_INVZ)//fractions
TOKEN("cross", T_CROSS)//row vectors

//matrices
TOKEN("I", T_IDEN)
TOKEN("sum", T_SUM) TOKEN("mean", T_MEAN) TOKEN("min", T_MIN) TOKEN("max", T_MAX) TOKEN("product", T_PRODUCT)
TOKEN("ref", T_REF) TOKEN("rref", T_RREF) TOKEN("rref2", T_RREF2)
TOKEN("det", T_DET) TOKEN("tr", T_TRACE) TOKEN("inv", T_INV) TOKEN("lu", T_LU)
TOKEN("egval", T_EGVAL) TOKEN("egvec", T_EGVEC) TOKEN("nullspace", T_NULLSPACE) TOKEN("diag", T_DIAG)

//modular arithmetic
TOKEN("rdtsc", T_RDTSC) TOKEN("randmod", T_RANDMOD) TOKEN("randlog", T_RANDLOG)
TOKEN("gcd", T_GCD) TOKEN("lcm", T_LCM) TOKEN("invmod", T_INVMOD)
TOKEN("isprime", T_ISPRIME) TOKEN("mrtest", T_MRTEST) TOKEN("genprime", T_GENPRIME)
TOKEN("factorize", T_FACTORIZE) TOKEN("factorize2", T_FACTORIZE2) TOKEN("totient", T_TOTIENT) TOKEN("carmichael", T_CARMICHAEL) TOKEN("proots", T_PROOTS)
TOKEN("ntt", T_NTT) TOKEN("intt", T_INTT)
TOKEN("log2", T_LOG2) TOKEN("log10", T_LOG10) TOKEN("dlog", T_DLOG)

//generic functions
TOKEN("sqrt", T_SQRT) TOKEN("cbrt", T_CBRT)
//TOKEN("floor", T_FLOOR) TOKEN("ceil", T_CEIL) TOKEN("round", T_ROUND)
//TOKEN("exp", T_EXP) TOKEN("ln", T_LN) TOKEN("log", T_LOG)
//TOKEN("gamma", T_GAMMA) TOKEN("lngamma", T_LNGAMMA)
//TOKEN("cos", T_COS) TOKEN("acos", T_ACOS) TOKEN("cosd", T_COSD) TOKEN("acosd", T_ACOSD) TOKEN("cosh", T_COSH) TOKEN("acosh", T_ACOSH)
//TOKEN("sec", T_SEC) TOKEN("asec", T_ASEC) TOKEN("secd", T_SECD) TOKEN("asecd", T_ASECD) TOKEN("sech", T_SECH) TOKEN("asech", T_ASECH)
//TOKEN("sin", T_SIN) TOKEN("asin", T_ASIN) TOKEN("sind", T_SIND) TOKEN("asind", T_ASIND) TOKEN("sinh", T_SINH) TOKEN("asinh", T_ASINH)
//TOKEN("csc", T_CSC) TOKEN("acsc", T_ACSC) TOKEN("cscd", T_CSCD) TOKEN("acscd", T_ACSCD) TOKEN("csch", T_CSCH) TOKEN("acsch", T_ACSCH)
//TOKEN("tan", T_TAN) TOKEN("atan", T_ATAN) TOKEN("tand", T_TAND) TOKEN("atand", T_ATAND) TOKEN("tanh", T_TANH) TOKEN("atanh", T_ATANH)
//TOKEN("cot", T_COT) TOKEN("acot", T_ACOT) TOKEN("cotd", T_COTD) TOKEN("acotd", T_ACOTD) TOKEN("coth", T_COTH) TOKEN("acoth", T_ACOTH)

//crypto
TOKEN("hist", T_HIST)
TOKEN("sha224", T_SHA224) TOKEN("sha256", T_SHA256) TOKEN("sha384", T_SHA384) TOKEN("sha512", T_SHA512)
TOKEN("sha3-224", T_SHA3_224) TOKEN("sha3-256", T_SHA3_256) TOKEN("sha3-384", T_SHA3_384) TOKEN("sha3-512", T_SHA3_512)
TOKEN("shake128", T_SHAKE128) TOKEN("shake256", T_SHAKE256)
TOKEN("aes128", T_AES128) TOKEN("aes192", T_AES192) TOKEN("aes256", T_AES256)

//operators
TOKEN("\'", T_TRANSPOSE)//unary postfix
TOKEN("^", T_POWER) TOKEN(":", T_COLON)
TOKEN("*", T_MUL) TOKEN("o", T_TENSOR) TOKEN("/", T_DIV) TOKEN("\\", T_DIV_BACK) TOKEN("%", T_MOD)
TOKEN(".*", T_MUL_EW) TOKEN("./", T_DIV_EW)//matrix-matrix element-wise operations
TOKEN("+", T_PLUS) TOKEN("-", T_MINUS)//can be unary		binary: if one is scalar then element-wise
TOKEN("<", T_LESS) TOKEN("<=", T_LESS_EQUAL) TOKEN(">", T_GREATER) TOKEN(">=", T_GREATER_EQUAL)
TOKEN("==", T_EQUAL) TOKEN("!=", T_NOT_EQUAL) TOKEN(".==", T_EQUAL_EW) TOKEN(".!=", T_NOT_EQUAL_EW)
TOKEN("mod", T_MOD_LP)

TOKEN("+=", T_ASSIGN_ADD) TOKEN("-=", T_ASSIGN_SUB)
TOKEN("*=", T_ASSIGN_MUL) TOKEN("/=", T_ASSIGN_DIV) TOKEN("%=", T_ASSIGN_MOD)
TOKEN("=", T_ASSIGN)

//control
TOKEN(",", T_COMMA) TOKEN(";", T_SEMICOLON)
TOKEN("(", T_LPR) TOKEN(")", T_RPR)
TOKEN("[", T_LBRACKET) TOKEN("]", T_RBRACKET)	//matrix/vector
//TOKEN("[[", T_POLSTART) TOKEN("]]", T_POLEND)	//polynomial	//X  need nested structure support
//TOKEN("{", T_LBRACE) TOKEN("}", T_RBRACE)		//scope

//commands
TOKEN("help", T_HELP) TOKEN("keywords", T_KEYWORDS)
TOKEN("open", T_OPEN)
TOKEN("clear", T_CLEAR) TOKEN("vars", T_VARS) //TOKEN("fractions", T_FRACTIONS) TOKEN("tolerance", T_TOLERANCE)
TOKEN("gfset", T_GFSET)
TOKEN("benchmark", T_BENCHMARK)
TOKEN("bin", T_BIN) TOKEN("dec", T_DEC) TOKEN("hex", T_HEX)
TOKEN("ascii", T_ASCII)
TOKEN("exit", T_EXIT) TOKEN("quit", T_QUIT)

//constants
//TOKEN("i", T_IMAG) TOKEN("j", T_IMAG_UNUSED) TOKEN("e", T_EULER) TOKEN("pi", T_PI) TOKEN("inf", T_INF) TOKEN("nan", T_NAN)
TOKEN("rand", T_RAND)

//lexer stuff
TOKEN(0, T_SPACE) TOKEN(0, T_NEWLINE) TOKEN(0, T_EOF)//consequent in this order
TOKEN(0, T_NUMBER)	//check lex_number
TOKEN(0, T_ID)		//check lex_id

TOKEN(0, T_NTOKENS)
#endif