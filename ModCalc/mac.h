//mac.h - Main include file
//Copyright (C) 2021  Ayman Wagih Mohsen
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef MC2_H
#define MC2_H
#include	<string>
#include	<map>
#include	"mac_memory.h"
#ifdef __linux__
#define		sprintf_s	snprintf
#define		vsprintf_s	vsnprintf
#define		scanf_s		scanf
#define		__rdtsc		__builtin_ia32_rdtsc
#define		_HUGE		HUGE_VAL
char		_getch();
#endif
#define		SIZEOF(STATIC_ARRAY)	(sizeof(STATIC_ARRAY)/sizeof(*STATIC_ARRAY))

//mac_system
#ifdef __linux__
#define			set_console_buffer_size(...)
const char*		open_file_window();
#else
int				set_console_buffer_size(short w, short h);
const wchar_t*	open_file_window();
#endif
int				file_is_readablea(const char *filename);
int				file_is_readablew(const wchar_t *filename);


//mac_math
int		maximum(int a, int b);
int		minimum(int a, int b);
int		clamp(int lo, int x, int hi);
int		mod(int x, int n);
int		first_set_bit(unsigned long long n);//idx of LSB
int		first_set_bit16(unsigned short n);//idx of LSB
int		floor_log2(unsigned long long n);//idx of MSB
int		ceil_log2(unsigned long long n);
int		floor_log10(unsigned long long n);//idx of MSD
double	power(double x, int y);
double	_10pow(int n);

int		floor_log10(Int const &n);//index of MSD
int		floor_log2(Int const &n);//index of MSB
Int		impl_gcd(Int a, Int b);//by copy
Int		impl_lcm(Int *arr, int count);//array is destroyed
Int		impl_intEEA(Int const &x, Int const &n);
Int		impl_powmod(Int const &x, Int const &e, Int const &n);

void	dec2frac(double x, double error, int *i, int *num, int *den);
int		query_double(double x, int *point);
int		print_double(double x, int point_pos, int total);//point_pos==0: no leading spaces, total==0: no trailing spaces
int		print_double_frac(double x, int point_pos, int total);
void	print_matrix_debug(Int const *buf, int w, int h);

Int		isqrt(Int const &x);
Int		icbrt(Int const &x);
int		impl_isprime(Int &n);
bool	impl_mrtest(Int const &n, Int const &a);
Int		impl_genprime(Int start);
Int		impl_randmod(Int const &modulus);
Int		impl_randlog(Int const &L);
//void	impl_factorize(Int n, IntVec &factors, IntVec *powers=nullptr, IntVec *rfactors=nullptr);
Int		impl_totient(Int const &n, IntVec &unique_factors, IntVec &powers);
Int		impl_carmichael(Int const &n, IntVec &unique_factors, IntVec &powers);
bool	impl_proots(Int n, IntVec &roots);//primitive roots
Int		impl_dlog(Int const &x, Int const &b, Int const &n);//r=dlog(x, b, n)	b^r = x (mod n)
//Int	c_mod(Int const &x, Int const &n);
//Int	c_invmod(Int const &x, Int const &n);
//Int	c_mulmod(Int const &a, Int const &b, Int const &n);
//Int	c_divmod(Int const &a, Int const &b, Int const &n);

void	impl_addbuffers(Int *dst, Int const *a, Int const *b, int count);
void	impl_subbuffers(Int *dst, Int const *a, Int const *b, int count);
void	impl_negbuffer(Int *dst, Int const *a, int count);
void	impl_buf_plus_val(Int *dst, Int const *buf, Int const &val, int count);
void	impl_val_minus_buf(Int *dst, Int const &val, Int const *buf, int count);
void	impl_buf_mul_val(Int *dst, Int const *buf, Int const &val, int count);

void	impl_ref(Int *m, short dx, short dy, Int const &n);
void	impl_rref(Int *m, short dx, short dy, Int const &n);
void	impl_rref2(Int *m, short dx, short dy, Int const &n);
Int		impl_det(Int *m, int dx, Int const &n);//m is rref'ed
void	impl_matinv(Int *m, short dx, Int const &n);//resize m to (dy * 2dx) temporarily,		dx==dy always
void	impl_matmul(Int *dst, const Int *A, const Int *B, int h1, int w1h2, int w2);
void	impl_transpose(Int *dst, const Int *src, int src_dx, int src_dy);
void	impl_tensor(Int *dst, const Int *m1, const Int *m2, int dx1, int dy1, int dx2, int dy2);
void	impl_matdiv(Int *dst, const Int *num, Int *den, int num_dy, int dx, Int const &n);		//dst & num: num_dy*dx,  den: dx*(dx*2)		den is destroyed
void	impl_matdiv_back(Int *dst, Int *den, const Int *num, int dy, int num_dx, Int const &n);	//den: dy*(dy*2),  dst & num: dy*num_dx		den is destroyed
void	impl_matpow(Int *dst, Int *m1, int e, int dx, Int const &n);//dst: dx*dx,  m1: dx*(dx*2)		calculates m1^e,	m1 is destroyed

void	impl_lu(Int const *m, int n, Int *lower, Int *upper);
void	impl_egval2(Int const *M, Int *lambdas);
void	impl_egval3(Int const *M, Int *lambdas);
void	impl_egval(Int const *M, int n, Int *D, int it_limit);
int		impl_nullspace(Int *M, int dx, int dy, Int *solution, char *dep_flags, short *row_idx);
//int	impl_egvec(Int const *M, int n, Int const *lambdas, Int *S);

void	impl_polmul(Int *res, Int const *A, Int const *B, int asize, int bsize, int add);//res has correct size of (asize+bsize-1)


//mac.cpp
#define				G_BUF_SIZE	65536
extern const int	g_buf_size;
extern char			g_buf[G_BUF_SIZE];

extern char	gfmode;//later

enum		TokenType
{
#define		TOKEN(STR, LABEL)	LABEL,
#include	"mac_keywords.h"
#undef		TOKEN
	T_NKEYWORDS,
};
extern const char *keywords[T_NKEYWORDS];
struct		Token
{
	TokenType type;
	Int idata;
	char *sdata;
};

#define		STRING_LESS_THAN(LEFT, RIGHT)	(strcmp(LEFT, RIGHT)<0)
struct		StringLibrary
{
	std::vector<char*> v;
	bool binary_search(char *e, int &idx)const
	{
		int size=(int)v.size(), L=0, R=size-1;
		while(L<=R)
		{
			int middle=(L+R)>>1;
			if(STRING_LESS_THAN(v[middle], e))
				L=middle+1;
			else if(STRING_LESS_THAN(e, v[middle]))
				R=middle-1;
			else
			{
				idx=middle;
				return true;
			}
		}
		idx=L+(L<size&&STRING_LESS_THAN(v[L], e));
		return false;
	}
	char* add(const char *static_array, int len=0, bool *old=nullptr)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(!len)
			len=(int)strlen(static_array);
		char *p=new char[len+1];
		memcpy(p, static_array, len);
		p[len]=0;
		int idx=0;
		bool found=binary_search(p, idx);
		if(!found)
			v.insert(v.begin()+idx, p);
		if(old)
			*old=found;
		return v[idx];
	}
	void clear()
	{
		for(int k=0;k<(int)v.size();++k)
			delete v[k];
	}
};
extern StringLibrary strings;

extern Int	zero, ten;
enum		MatrixFlags
{
	M_UNSPECIFIED_RANGE	=0x00000001,
	M_POLYNOMIAL_MASK	=0x0000FFFF,
	M_SCALAR			=0x00010001,
};
struct		Matrix//8+4 bytes
{
	//TokenType type;
	union
	{
		//dx & dy == 0: unspecified range
		//dx==0: polynomial		X
		//otherwise: matrix
		unsigned flags;
		struct{unsigned short dx, dy;};
	};
	IntVec v;
	char *name;//freed by StringLibrary
	Matrix():dx(0), dy(0), name(nullptr){}
	//Matrix():dx(0), dy(0), data(nullptr), name(nullptr){}
	Matrix(Matrix const &other):dx(other.dx), dy(other.dy), v(other.v), name(other.name){}
	Matrix(Matrix &&other):dx(other.dx), dy(other.dy), v(other.v), name(other.name)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(this!=&other)
			MEMZERO(&other, 1);
	}
	~Matrix()
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		//if(data)
		//	FREE(data);
#ifdef DEBUG_MEMORY
		if(syscall_count==60)
			syscall_count=syscall_count;
#endif
	}
	Matrix& operator=(Matrix const &other)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(this!=&other)
		{
			//type=other.type;
			dx=other.dx;
			dy=other.dy;
			//if(data)
			//	FREE(data);
			int count=dx*dy;
			v=other.v;
			//ALLOC(data, dx*dy);
			//MEMCPY(data, other.data, count);
			name=other.name;
		}
		return *this;
	}
	Matrix& operator=(Matrix &&other)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(this!=&other)
		{
			//type=other.type;
			dx=other.dx;
			dy=other.dy;
			v=std::move(other.v);
			//if(data)
			//	FREE(data);
			//data=other.data;
			name=other.name;
			MEMZERO(&other, 1);
		}
		return *this;
	}
	int size()const
	{
		return dx*dy;//TODO: check for polynomial (not implemented yet)
	}
	void print()const;
	void resize()
	{
		v.resize(dx*dy);//TODO: check for polynomial (not implemented yet)
	}
	//void resize(short new_dx, short new_dy)
	//{
	//	if(new_dx!=-1)
	//		dx=new_dx;
	//	if(new_dy!=-1)
	//		dy=new_dy;
	//	v.resize(dx*dy);
	//}
	void move_in(IntVec &temp, short new_dx, short new_dy)
	{
		if(new_dx!=-1)
			dx=new_dx;
		if(new_dy!=-1)
			dy=new_dy;
		v=std::move(temp);
	}
	void move2temp(Matrix &m)//to be used on temp matrices (eg: 2*x, x looses its name)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(this!=&m)
		{
			dx=m.dx;
			dy=m.dy;
			v=std::move(m.v);
			//if(data)
			//	FREE(data);
			//data=m.data;
			name=nullptr;
			MEMZERO(&m, 1);
		}
	}
	Int*		data()		{return v.data();}
	Int const*	data()const	{return v.data();}
	Int&		get(int kx, int ky)		{return v[dx*ky+kx];}
	Int const&	get(int kx, int ky)const{return v[dx*ky+kx];}
	//Int*			end()		{return data+dx*dy;}
	//Int const*	end()const	{return data+dx*dy;}
	void alloc_ramp(short dx, short dy, Int const &start)
	{
		this->dy=dy, this->dx=dx;
		int size=dy*dx;
		v.resize(size);
		//ALLOC(data, size);
		for(int k=0;k<size;++k)
			v[k]=start+k;
	}
};
#define		GET(DATA, DX, KX, KY)	DATA[(DX)*(KY)+(KX)]
void		print_help();
bool		get_str_from_file(std::string &str);


//mac_lexer.cpp
extern Int			lex_number;
extern char			*lex_id;
extern int			text_size, idx;
extern const char	*text;
extern const char	*keywords[T_NTOKENS+1];
void		lex_init(const char *str, int len);
int			lex_skip_space();
TokenType	lex_get(bool space_sensitive);//space-sensitive: true: you need the leading space/newline token (before an actual token), false: space is ignored
//TokenType	lex_look_ahead(int k, bool space_sensitive);
Token*		lex_next(bool space_sensitive);
void		lex_save();
void		lex_restore();
void		print_sorted_keywords();


//mac_parser.cpp
extern bool	benchmark;
enum		SolveResultType
{
	SOLVE_OK,
	SOLVE_OK_NO_ANS,
	SOLVE_INCOMPLETE,
	SOLVE_PARSE_ERROR,
};
extern std::vector<Matrix> g_answers;
extern std::map<char*, Matrix> g_vars;
extern std::vector<std::string> errors;
int			solve(std::string &str, bool again);
bool		user_error(const char *format, ...);//full description and punctuation
bool		user_error2(int start, int end, const char *format, ...);//refer in text where the error happened, a period '.' is appended to each error
#endif