//mc2.h - Main include file
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
#include	<vector>
#include	<string>
#include	<map>
#include	"mc2_memory.h"
#ifdef __linux__
#define		sprintf_s	snprintf
#define		vsprintf_s	vsnprintf
#define		scanf_s		scanf
#define		__rdtsc		__builtin_ia32_rdtsc
#define		_HUGE		HUGE_VAL
extern "C" char _getch();
#endif

//mc2_system
#ifdef __linux__
#define			set_console_buffer_size(...)
const char*		open_file_window();
int				file_is_readablea(const char *filename);
#else
int				set_console_buffer_size(short w, short h);
const wchar_t*	open_file_window();
int				file_is_readablea(const char *filename);
int				file_is_readablew(const wchar_t *filename);
#endif


//mc2_math.c
extern "C"
{
	int		maximum(int a, int b);
	int		minimum(int a, int b);
	int		clamp(int lo, int x, int hi);
	int		mod(int x, int n);
	int		first_set_bit(unsigned long long n);//idx of LSB
	int		first_set_bit16(unsigned short n);//idx of LSB
	int		floor_log2(unsigned long long n);//idx of MSB
	int		ceil_log2(unsigned long long n);
	int		floor_log10(double x);
	double	power(double x, int y);
	double	_10pow(int n);
	
	void	dec2frac(double x, double error, int *i, int *num, int *den);
	int		query_double(double x, int *point);
	int		print_double(double x, int point_pos, int total);//point_pos==0: no leading spaces, total==0: no trailing spaces
	int		print_double_frac(double x, int point_pos, int total);
//	void	print_matrix(Comp const *data, int w, int h);
	void	print_matrix_debug(Comp const *buf, int w, int h);

	double	c_abs2(Comp *z);
	void	c_inv(Comp *dst, const Comp *z);
	void	c_mul(Comp *dst, const Comp *a, const Comp *b);
	void	c_mul_add(Comp *dst, const Comp *a, const Comp *b);
	void	c_div(Comp *dst, const Comp *a, const Comp *b);
	void	c_mod(Comp *dst, const Comp *a, const Comp *b);

	void	impl_addbuffers(Comp *dst, Comp const *a, Comp const *b, int count);//count of doubles
	void	impl_subbuffers(Comp *dst, Comp const *a, Comp const *b, int count);
	void	impl_negbuffer(Comp *dst, Comp const *a, int count);
	void	impl_buf_plus_val(Comp *dst, Comp const *buf, Comp val, int count);
	void	impl_val_minus_buf(Comp *dst, Comp val, Comp const *buf, int count);
	void	impl_buf_mul_val(Comp *dst, Comp const *buf, Comp val, int count);

	void	impl_ref(Comp *m, short dx, short dy);
	void	impl_rref(Comp *m, short dx, short dy);
	void	impl_rref2(Comp *m, short dx, short dy);
	Comp	impl_det(Comp *m, int dx);//m is destroyed
	void	impl_matinv(Comp *m, short dx);//resize m to (dy * 2dx) temporarily,		dx==dy always
	void	impl_matmul(Comp *dst, const Comp *A, const Comp *B, int h1, int w1h2, int w2);
	void	impl_transpose(Comp *dst, const Comp *src, int src_dx, int src_dy);
	void	impl_tensor(Comp *dst, const Comp *m1, const Comp *m2, int dx1, int dy1, int dx2, int dy2);
	void	impl_matdiv(Comp *dst, const Comp *num, Comp *den, int num_dy, int dx);		//dst & num: num_dy*dx,  den: dx*(dx*2)		den is destroyed
	void	impl_matdiv_back(Comp *dst, Comp *den, const Comp *num, int dy, int num_dx);	//den: dy*(dy*2),  dst & num: dy*num_dx		den is destroyed
	void	impl_matpow(Comp *dst, Comp *m1, int e, int dx);//dst: dx*dx,  m1: dx*(dx*2)		calculates m1^e,	m1 is destroyed

	void	impl_lu(Comp const *m, int n, Comp *lower, Comp *upper);
	void	impl_egval2(Comp const *M, Comp *lambdas);
	void	impl_egval3(Comp const *M, Comp *lambdas);
	void	impl_egval(Comp const *M, int n, Comp *D, int it_limit);
	int		impl_nullspace(Comp *M, int dx, int dy, Comp *solution, char *dep_flags, short *row_idx);
	int		impl_egvec(Comp const *M, int n, Comp const *lambdas, Comp *S);
//	int		impl_diag2(Comp const *M, Comp *invS, Comp *D, Comp *S);
	
	void	impl_polmul(Comp *res, Comp const *A, Comp const *B, int asize, int bsize, int add);//res has correct size of (asize+bsize-1)
}


//mc2.cpp
#define		G_BUF_SIZE	65536
extern const int	g_buf_size;
extern char			g_buf[G_BUF_SIZE];

extern char	gfmode;//later

enum		TokenType
{
#define		TOKEN(STR, LABEL)	LABEL,
#include	"mc2_keywords.h"
#undef		TOKEN
};

#define		STRING_LESS_THAN(LEFT, RIGHT)	(strcmp(LEFT, RIGHT)<0)
struct		StringLibrary
{
	std::vector<char*> v;
	bool binary_search(char *e, int &idx)const
	{
		int size=v.size(), L=0, R=size-1;
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
			len=strlen(static_array);
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

enum		MatrixFlags
{
	M_UNSPECIFIED_RANGE	=0x00000001,
	M_POLYNOMIAL_MASK	=0x0000FFFF,
	M_SCALAR			=0x00010001,
};
//#define	M_UNSPECIFIED_RANGE		0x0000000000000001
//#define	M_POLYNOMIAL_MASK		0xFFFF000000000000
//#define	M_SCALAR				0x0000000100010001
struct		Matrix//2*4+4+4=16 bytes
{
	union
	{
		unsigned flags;
		struct{unsigned short dx, dy;};
	};
	unsigned short dz, type;
	//unsigned short
	//	type,		//0: matrix, 1: fraction, 2: matrix of fractions
	//	dx, dy,
	//	dz;			//unused

	//TokenType type;
	//union
	//{
	//	//dx==1 & dy==0: unspecified range
	//	//dx>1, dy==0: polynomial		TODO
	//	//otherwise: matrix
	//	unsigned long long flags;
	//	struct{unsigned short type, dx, dy, dz;};
	//};
	union
	{
		Comp *data;//complex numbers are interleaved
		Matrix *poly;
	};
	char *name;//freed by StringLibrary
	Matrix():type(0), dx(0), dy(0), dz(0), data(nullptr), name(nullptr){}
	Matrix(Matrix const &other):type(other.type), dx(other.dx), dy(other.dy), dz(other.dz), data(other.data), name(other.name){}
	Matrix(Matrix &&other):type(other.type), dx(other.dx), dy(other.dy), dz(other.dz), data(other.data), name(other.name)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(this!=&other)
			MEMZERO(Matrix, &other, 1);
			//memset(&other, 0, sizeof(Matrix));
	}
	~Matrix()
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(data)
			CFREE(data);
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
			if(data)
				CFREE(data);
			int count=dx*dy;
			CALLOC(data, dx*dy);
			//data=(double*)malloc(bytesize);
			CMEMCPY(data, other.data, count);
			//memcpy(data, other.data, count*sizeof(double));
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
			if(data)
				CFREE(data);
			data=other.data;
			name=other.name;
			MEMZERO(Matrix, &other, 1);
			//memset(&other, 0, sizeof(Matrix));
		}
		return *this;
	}
	void print()const;
	void move2temp(Matrix &m)//to be used on temp matrices (eg: 2*x, x looses its name)
	{
#ifdef DEBUG_MEMORY
		const char file[]=__FILE__;
#endif
		if(this!=&m)
		{
			//type=m.type;
			dx=m.dx;
			dy=m.dy;
			if(data)
				CFREE(data);
			data=m.data;
			name=nullptr;
			MEMZERO(Matrix, &m, 1);
		}
	}
	Comp&		get(int kx, int ky)		{return data[dx*ky+kx];}
	Comp const& get(int kx, int ky)const{return data[dx*ky+kx];}
	Comp*		end()		{return data+dx*dy;}
	Comp const*	end()const	{return data+dx*dy;}
	void alloc_ramp(short dx, short dy, double start)
	{
		this->dy=dy, this->dx=dx;
		int size=dy*dx;
		CALLOC(data, size);
		for(int k=0;k<size;++k)
			data[k].r=start+k, data[k].i=0;
	}
	int size()const
	{
		return dx*dy;//TODO: check for polynomial (not implemented yet)
	}
	void resize()
	{
		CREALLOC(data, data, size());
	}
};
//struct		Fraction
//{
//	std::vector<Comp> num, den;
//};
#define		IN_RANGE(CONST_START, VAR, CONST_END_INCLUSIVE)	(unsigned((VAR)-(CONST_START))<unsigned((CONST_END_INCLUSIVE)+1-(CONST_START)))
#define		GET(DATA, DX, KX, KY)	DATA[(DX)*(KY)+(KX)]
void		print_help();
bool		get_str_from_file(std::string &str);


//mc2_lexer.cpp
extern double		lex_number;
extern char			*lex_id;
extern int			text_size, idx;
extern const char	*text;//
extern const char	*keywords[T_NTOKENS+1];
void		lex_init(const char *str, int len);
int			lex_skip_space();
TokenType	lex_get(bool space_sensitive);//space-sensitive: true: you need the leading space/newline token (before an actual token), false: space is ignored
TokenType	lex_look_ahead(int k, bool space_sensitive);
long long	acme_read_integer(const char *text, int size, int base, int start, int *ret_end, int *ret_logb);
char		acme_read_number(const char *text, int size, int k, int *advance, long long *ival, double *fval);//returns  0: integer, 1: float, 2: error


//mc2_parser.cpp
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