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
#include	"mc2_memory.h"
#include	<vector>
#include	<string>
#include	<map>

//mc2_system.c
extern "C"
{
	int				set_console_buffer_size(short w, short h);
	const wchar_t*	open_file_window();
}


//mc2_math.c
extern "C"
{
	int		maximum(int a, int b);
	int		minimum(int a, int b);
	int		mod(int x, int n);
	int		first_set_bit(unsigned long long n);//idx of LSB
	int		first_set_bit16(unsigned short n);//idx of LSB
	int		floor_log2(unsigned long long n);//idx of MSB
	int		ceil_log2(unsigned long long n);
	int		floor_log10(double x);
	double	power(double x, int y);
	double	_10pow(int n);
	
	void	print_matrix_debug(const double *buf, int w, int h);

	void	impl_addbuffers(double *dst, double const *a, double const *b, int size);
	void	impl_subbuffers(double *dst, double const *a, double const *b, int size);
	void	impl_negbuffer(double *dst, double const *a, int size);
	void	impl_buf_plus_val(double *dst, double const *buf, double val, int size);
	void	impl_val_minus_buf(double *dst, double val, double const *buf, int size);
	void	impl_buf_mul_val(double *dst, double const *buf, double val, int size);

	void	impl_ref(double *m, short dx, short dy);
	void	impl_rref(double *m, short dx, short dy);
	double	impl_det(double *m, int dx);//m is destroyed
	void	impl_matinv(double *m, short dx);//resize m to (dy * 2dx) temporarily,		dx==dy always
	void	impl_matmul(double *dst, const double *A, const double *B, int h1, int w1h2, int w2);
	void	impl_transpose(double *dst, const double *src, int src_dx, int src_dy);
	void	impl_tensor(double *dst, const double *m1, const double *m2, int dx1, int dy1, int dx2, int dy2);
	void	impl_matdiv(double *dst, const double *num, double *den, int num_dy, int dx);		//dst & num: num_dy*dx,  den: dx*(dx*2)		den is destroyed
	void	impl_matdiv_back(double *dst, double *den, const double *num, int dy, int num_dx);	//den: dy*(dy*2),  dst & num: dy*num_dx		den is destroyed
	void	impl_matpow(double *dst, double *m1, int e, int dx);//dst: dx*dx,  m1: dx*(dx*2)		calculates m1^e,	m1 is destroyed
	
	void	impl_polmul(double *res, double const *A, double const *B, int asize, int bsize, int add);//res has correct size of (asize+bsize-1)
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

struct		Matrix//12+4 bytes
{
	TokenType type;
	unsigned short dx, dy;
	double *data;//complex numbers are interleaved
	char *name;//freed by StringLibrary
	Matrix():type(T_IGNORED), dx(0), dy(0), data(nullptr), name(nullptr){}
	Matrix(Matrix const &other):type(other.type), dx(other.dx), dy(other.dy), data(other.data), name(other.name){}
	Matrix(Matrix &&other):type(other.type), dx(other.dx), dy(other.dy), data(other.data), name(other.name)
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
			free(data);
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
			type=other.type;
			dx=other.dx;
			dy=other.dy;
			int count=dx*dy;
			DALLOC(data, dx*dy);
			//data=(double*)malloc(bytesize);
			memcpy(data, other.data, count*sizeof(double));
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
			type=other.type;
			dx=other.dx;
			dy=other.dy;
			data=other.data;
			name=other.name;
			MEMZERO(Matrix, &other, 1);
			//memset(&other, 0, sizeof(Matrix));
		}
		return *this;
	}
	void print()
	{
		if(!data)
		{
			printf("Error: data == nullptr\n");
			return;
		}
		//if(name)
		//	printf("%s =\n", name);
		if(dx==1&&dy==1)//scalar
		{
			if(type==T_REAL)
				printf("%4g\n", *data);
			else//complex scalar
			{
				printf("%4g + j%4g\n", data[0], data[1]);
			}
		}
		else//matrix
		{
			printf("[\n");
			for(int ky=0;ky<dy;++ky)
			{
				for(int kx=0;kx<dx;++kx)
				{
					if(type==T_REAL)
						printf("%4g", data[dx*ky+kx]);
					else
						printf("%4g + j%4g", data[2*(dx*ky+kx)], data[2*(dx*ky+kx)+1]);
					if(kx<dx-1)
						printf(",");
				}
				if(ky<dy-1)
					printf(";\n");
			}
			printf("\n]\n");
		}
	}
	void reset()
	{
		type=T_IGNORED;
		dx=0, dy=0;
		data=nullptr;
		name=nullptr;
	}
};
extern std::vector<Matrix> g_answers;
extern std::map<char*, Matrix> g_vars;
void		print_help();
bool		get_str_from_file(std::string &str);


//mc2_lexer.cpp
extern double		lex_number;
extern char			*lex_id;
extern int			text_size, idx;
extern const char	*text;//
extern const char	*keywords[T_NTOKENS+1];
void		lex_init(const char *str, int len);
TokenType	lex_get();
TokenType	lex_look_ahead(int k);


//mc2_parser.cpp
enum		SolveResultType
{
	SOLVE_OK,
	SOLVE_OK_NO_ANS,
	SOLVE_INCOMPLETE,
	SOLVE_PARSE_ERROR,
};
int			solve(std::string &str, bool again);
#endif