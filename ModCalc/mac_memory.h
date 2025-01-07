//mac_memory.c - Memory operations
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

#ifndef MC2_CSTUFF
#define MC2_CSTUFF
#include		<string.h>
#include		"longer-int.h"
#include		<vector>
#include		<string>


//	#define	DEBUG_MEMORY

#ifdef __linux__
#include		<stdlib.h>
#include		<assert.h>
//#define		_aligned_malloc(SIZE, ALIGN)	aligned_alloc(ALIGN, SIZE)	//how to realloc?
//#define		_aligned_free(POINTER)			free(POINTER)
#define			GETPTR(PTR)		((void**)(PTR))[-1]
#define			GETSIZE(PTR)	((size_t*)(PTR))[-2]
static void*	_aligned_malloc(size_t size, size_t align)
{
	size_t addr=0, aa=0;
	assert(align>=sizeof(size_t));//
	addr=(size_t)malloc(size+align+2*sizeof(size_t));
	if(!addr)
		return 0;
	aa=addr+2*sizeof(size_t)+align;
	aa-=aa%align;
	GETPTR(aa)=(void*)addr;
	GETSIZE(aa)=size;
	return (void*)aa;
}
static void 	_aligned_free(void *p)
{
	free(GETPTR(p));
}
static void*	_aligned_realloc(void *oldp, size_t size, size_t align)//inherently slow
{
	int copysize=GETSIZE(oldp);
	void *p=_aligned_malloc(size, align);
	if(!p)
		return 0;
	if(copysize>size)
		copysize=size;
	memcpy(p, oldp, copysize);
	_aligned_free(oldp);
	return p;
}
#endif
#ifdef _MSC_VER
//#define		scanf	scanf_s
#endif
#ifdef DEBUG_MEMORY
#include	<conio.h>
#define		malloc(SIZE)							d_alloc(file, __LINE__, SIZE)
#define		realloc(POINTER, SIZE)					d_realloc(file, __LINE__, POINTER, SIZE)
#define		free(POINTER)							d_free(file, __LINE__, POINTER)

#define		_aligned_malloc(SIZE, ALIGN)			d_aligned_alloc(file, __LINE__, SIZE, ALIGN)
#define		_aligned_realloc(POINTER, SIZE, ALIGN)	d_aligned_realloc(file, __LINE__, POINTER, SIZE, ALIGN)
#define		_aligned_free(POINTER)					d_aligned_free(file, __LINE__, POINTER)

#define		memcpy(DST, SRC, SIZE)					d_memcpy(file, __LINE__, DST, SRC, SIZE)
#define		memmove(DST, SRC, SIZE)					d_memmove(file, __LINE__, DST, SRC, SIZE)
#define		memset(DST, VAL, SIZE)					d_memset(file, __LINE__, DST, VAL, SIZE)
//#ifdef __cplusplus
//extern "C"
//{
//#endif
	extern int	syscall_count, emergency_flag;
	void*	d_alloc				(const char *file, int line, size_t bytesize);
	void*	d_realloc			(const char *file, int line, void *p, size_t bytesize);
	int		d_free				(const char *file, int line, void *p);

	void*	d_aligned_alloc		(const char *file, int line, size_t bytesize, size_t alignment);
	void*	d_aligned_realloc	(const char *file, int line, void *p, size_t bytesize, size_t alignment);
	int		d_aligned_free		(const char *file, int line, void *p);

	void	d_memset			(const char *file, int line, void *dst, int val, size_t bytesize);
	void	d_memcpy			(const char *file, int line, void *dst, const void *src, size_t bytesize);
	void	d_memmove			(const char *file, int line, void *dst, const void *src, size_t bytesize);
//#ifdef __cplusplus
//}
//#endif
#endif
static void	memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
#ifdef DEBUG_MEMORY
	const char file[]=__FILE__;
#endif
	unsigned copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=(unsigned)srcbytes;
	memcpy(d, s, copied);
	while(((size_t)copied<<1)<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#define		ALLOC(NEWP, COUNT)						NEWP=(decltype(&*(NEWP)))malloc((COUNT)*sizeof(*NEWP))
#define		REALLOC(NEWP, OLDP, COUNT)				NEWP=(decltype(&*(NEWP)))realloc(OLDP, (COUNT)*sizeof(*NEWP))
#define		FREE(OLDP)								free(OLDP)
#define		MEMFILL(DST, SRC, DSTCOUNT, SRCCOUNT)	memfill(DST, SRC, (DSTCOUNT)*sizeof(*(DST)), (SRCCOUNT)*sizeof(*(SRC)))
#define		MEMZERO(DST, COUNT)						memset(DST, 0, (COUNT)*sizeof(*(DST)))
#define		MEMCPY(DST, SRC, COUNT)					memcpy(DST, SRC, (COUNT)*sizeof(*(DST)))
#define		MEMMOVE(DST, SRC, COUNT)				memmove(DST, SRC, (COUNT)*sizeof(*(DST)))
//template<typename T>inline void objzero(T *dst, size_t count)
//{
//	for(size_t k=0;k<count;++k)
//		dst[k]=T();
//}
template<typename T>inline T* objmove(T *dst, T *src, size_t count)
{
	if(dst<src)//move back
	{
		for(size_t k=0;k<count;++k)
			dst[k]=std::move(src[k]);
	}
	else if(dst>src)//move forward
	{
		for(int k=(int)(count-1);k>=0;--k)
			dst[k]=std::move(src[k]);
	}
	return dst;
}

struct		Int;
extern Int	zero, one, ten;
void		impl_ldiv(Int num, Int den, Int &Q, Int &R);
static Int	mulmod(Int const &a, Int const &b, Int const &n);
Int			impl_powmod(Int const &x, Int const &e, Int const &n);
struct		Int
{
	LINT *data;
	Int():data(new_lint_num_degree(0, 0)){}
	Int(Int const &other):data(clone_lint(other.data)){}
	Int(Int &&other):data(other.data)
	{
		other.data=nullptr;
	}
//	Int(LINT *data):data(data){}
	Int(const char *digits):data(new_lint_str(digits)){}
	Int(unsigned number, size_t exponent=0):data(new_lint_num_degree(number, exponent)){}
	~Int()
	{
		free_lint(data);
	}
	Int& operator=(Int const &other)
	{
		if(&other!=this)
			copy_lint(data, other.data);
		return *this;
	}
	Int& operator=(Int &&other)
	{
		if(&other!=this)
		{
			if(data)
				free_lint(data);
			data=other.data;
			other.data=nullptr;
		}
		return *this;
	}

	Int& operator+=(Int const &other)
	{
		add_lint(data, other.data);
		return *this;
	}
	Int& operator-=(Int const &other)
	{
		subtract_lint(data, other.data);
		return *this;
	}
	Int& operator*=(Int const &other)
	{
		mul_lint(data, other.data);
		return *this;
	}
	Int& operator/=(Int const &other)//use impl_ldiv
	{
		Int Q, R;
		impl_ldiv(*this, other, Q, R);
		*this=Q;
		//div_lint(data, other.data);
		return *this;
	}
	Int& operator%=(Int const &other)//use impl_ldiv
	{
		Int Q, R;
		impl_ldiv(*this, other, Q, R);
		*this=R;
		//mod_lint(data, other.data);
		return *this;
	}
	Int& operator^=(Int const &other)
	{
		pow_lint(data, other.data);
		return *this;
	}
	int mulmod(Int const &x, Int const &n);
	int powmod(Int const &x, Int const &n)
	{
		*this=impl_powmod(*this, x, n);
		return 1;//

		//return pow_lint_mod(data, x.data, n.data);
	}
	void negate()
	{
		if(data)
			data->sign=!data->sign;
	}
	void abs()
	{
		if(data)
			data->sign=0;
	}
	Int& operator>>=(size_t sh)
	{
		shr_lint(data, sh);
		return *this;
	}
	Int& operator<<=(size_t sh)
	{
		shl_lint(data, sh);
		return *this;
	}
	
	int construct(int number, size_t exponent=0)
	{
		data=new_lint_num_degree(number, exponent);
		return data!=0;
	}
	int construct(const char *digits)
	{
		data=new_lint_str(digits);
		return data!=0;
	}
	void setzero()
	{
		if(!data)
			data=new_lint_num_degree(0, 0);
		else
		{
			data->size=data->used_size=1;
			data->sign=0;
			REALLOC(data->x, data->x, 1);
			data->x[0]=0;
		}
		//if(data)
		//	free_lint(data);//need clear function
		//data=new_lint_num_degree(0, 0);
	}
	void setzero(int logn_32)
	{
		if(!data)
			data=new_lint_num_degree(0, 0);
		data->size=data->used_size=logn_32;
		data->sign=0;
		REALLOC(data->x, data->x, logn_32);
		MEMZERO(data->x, logn_32);
	}
	void set(long long n)
	{
		bool neg=n<0;
		if(neg)
			n=-n;
		if(!data)
			data=new_lint_num_degree(n>>32, 1);
		else
		{
			REALLOC(data->x, data->x, 2);
			data->x[1]=n>>32;
			data->used_size=data->size=2;
		}
		data->x[0]=(unsigned)n;
		data->sign=neg;
	}
	int set(const char *digits)
	{
		if(!data)
		{
			data=new_lint_str(digits);
			return data!=0;
		}
		return init_lint_str(data, digits);
	}
	int expand_size(){return ::expand_size(data);}
	int attempt_shrink(){return ::attempt_shrink(data);}

	int abs_cmp(Int const &other)const
	{
		return compare_lint_abs(data, other.data);
	}
	int getnibble(int idx)const
	{
		if(!data)
			return 0;
		return data->x[idx>>3]>>((idx&7)<<2)&15;
	}
	int getbit(int idx)const
	{
		if(!data)
			return 0;
		return data->x[idx>>5]>>(idx&31)&1;
	}
	void setbit(int idx)
	{
		if(!data)
			data=new_lint_num_degree(1<<(idx&31), idx>>5);
		else
		{
			int newsize=(idx>>5)+1;
			if((int)data->used_size<newsize)
			{
				REALLOC(data->x, data->x, newsize);
				MEMZERO(data->x+data->used_size, newsize-data->used_size);
				data->used_size=data->size=newsize;
			}
			data->x[idx>>5]|=1<<(idx&31);
		}
	}
	void clearbit(int idx)
	{
		if(data)
		{
			int newsize=(idx+31)>>5;
			if(newsize<(int)data->used_size)
				data->x[idx>>5]&=~(1<<(idx&31));
		}
	}
	void assignbit(int idx, int bit)
	{
		if(bit)
			setbit(idx);
		else
			clearbit(idx);
	}
	bool is_neg()const
	{
		if(!data)
			return false;
		return data->sign!=0;
	}
	bool is_odd()const
	{
		if(!data)
			return 0;
		return data->x[0]&1;
	}
	Int operator-()const
	{
		Int neg=*this;
		neg.data->sign=!neg.data->sign;
		return neg;
	}
	int to_int()const
	{
		if(!data||data->used_size>1)
			return 0;
		int ret=data->x[0];
		if(data->sign)
			ret=-ret;
		return ret;
	}
	long long to_longlong()const
	{
		if(!data||data->used_size>2)
			return 0;
		long long ret=data->x[0];
		if(data->used_size==2)
			ret|=(long long)data->x[1]<<32;
		if(data->sign)
			ret=-ret;
		return ret;
	}
	std::string to_string()const;
	int print()const
	{
		auto str=to_string();
		return printf("%s", str.c_str());

		//auto str=lint_itoa(data);
		//if(!str)
		//	return printf("nullptr");
		//int printed=printf("%s", str);
		//free(str);
		//return printed;

	//	print_lint(data);
	}
};
typedef std::vector<Int> IntVec;

inline void setzero(Int *dst, int count)
{
	for(int k=0;k<count;++k)
		dst[k].setzero();
}
inline Int& operator++(Int &n)
{
	add_lint(n.data, one.data);
	return n;
}
inline Int operator++(Int &n, int)//post-increment is bad
{
	Int temp=n;
	add_lint(n.data, one.data);
	return temp;
}
inline Int& operator--(Int &n)
{
	subtract_lint(n.data, one.data);
	return n;
}
inline Int operator--(Int &n, int)//post-decrement is bad
{
	Int temp=n;
	subtract_lint(n.data, one.data);
	return temp;
}
inline bool operator==(Int const &a, Int const &b)
{
	return !compare_lint(a.data, b.data);
}
inline bool operator!=(Int const &a, Int const &b)
{
	return compare_lint(a.data, b.data)!=0;
}
inline bool operator<(Int const &a, Int const &b)
{
	return compare_lint(a.data, b.data)<0;
}
inline bool operator>(Int const &a, Int const &b)
{
	return compare_lint(a.data, b.data)>0;
}
inline bool operator<=(Int const &a, Int const &b)
{
	int ret=compare_lint(a.data, b.data);
	return ret<=0;
}
inline bool operator>=(Int const &a, Int const &b)
{
	int ret=compare_lint(a.data, b.data);
	return ret>=0;
}
inline Int operator+(Int const &a, Int const &b)
{
	Int t=a;
	t+=b;
	return t;
}
inline Int operator-(Int const &a, Int const &b)
{
	Int t=a;
	t-=b;
	return t;
}
inline Int operator*(Int const &a, Int const &b)
{
	Int t=a;
	t*=b;
	return t;
}

int			floor_log2(Int const &n);//index of MSB
inline int	Int::mulmod(Int const &x, Int const &n)
{
	if(n==zero)
		return 0;
	*this*=x;
	Int Q, R;
	impl_ldiv(*this, n, Q, R);
	*this=R;
	return 1;
	//return mul_lint_mod(data, x.data, n.data);
}
static Int mulmod(Int const &a, Int const &b, Int const &n)
{
	Int t=a;
	t*=b;
	Int Q, R;
	impl_ldiv(t, n, Q, R);
	return R;
}
inline Int operator/(Int const &a, Int const &b)//use ldiv
{
	Int Q, R;
	impl_ldiv(a, b, Q, R);
	return Q;

	//Int t=a;
	//t/=b;
	//return t;
}
inline Int operator%(Int const &a, Int const &b)//use ldiv
{
	Int Q, R;
	impl_ldiv(a, b, Q, R);
	return R;

	//Int t=a;
	//t%=b;
	//return t;
}
inline Int operator^(Int const &a, Int const &b)
{
	Int t=a;
	t^=b;
	return t;
}
inline Int operator>>(Int const &a, size_t sh)
{
	Int t=a;
	shr_lint(t.data, sh);
	return t;
}
inline Int operator<<(Int const &a, size_t sh)
{
	Int t=a;
	shl_lint(t.data, sh);
	return t;
}

inline Int abs(Int const &x)
{
	Int t=x;
	t.data->sign=!t.data->sign;
	return t;
}
void		impl_factorize(Int n, IntVec &factors, IntVec *powers=nullptr, IntVec *rfactors=nullptr);//pass 'powers' vector to get unique factors
//enum		FactorizationType
//{
//	F_REPEATED,				//12: [2 2 3]
//	F_UNIQUE,				//12: [2 3]
//	F_UNIQUE_WITH_COUNT,	//12: [2 2;3 1]
//	F_POWERS,				//12: [4 3]
//	F_POWERS_WITH_COUNT,	//12: [4 2;3 1]
//};
//typedef LINT *Int;
//typedef long long Int;
/*union		Int2
{
	struct
	{
		Int lo, hi;
	};
	struct
	{
		unsigned u[4];
	};
	struct
	{
		int i[4];
	};
	Int2():lo(0), hi(0){}
	Int2(Int const &hi, Int const &lo):lo(lo), hi(hi){}
};//*/
#endif