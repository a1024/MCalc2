//mc2_math.c - MCalc2 math operations
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

#include	"mc2_memory.h"
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<tmmintrin.h>
int			maximum(int a, int b){return a>b?a:b;}
int			minimum(int a, int b){return a<b?a:b;}
int			clamp(int lo, int x, int hi)
{
	if(x<lo)
		x=lo;
	if(x>hi)
		x=hi;
	return x;
}
int			mod(int x, int n)
{
	x%=n;
	x+=n&-(x<0);
	return x;
}
int			first_set_bit(unsigned long long n)//idx of LSB
{
	int sh=((n&((1ULL<<32)-1))==0)<<5,	idx =sh;	n>>=sh;
		sh=((n&((1   <<16)-1))==0)<<4,	idx+=sh;	n>>=sh;
		sh=((n&((1   << 8)-1))==0)<<3,	idx+=sh;	n>>=sh;
		sh=((n&((1   << 4)-1))==0)<<2,	idx+=sh;	n>>=sh;
		sh=((n&((1   << 2)-1))==0)<<1,	idx+=sh;	n>>=sh;
		sh= (n&((1   << 1)-1))==0,		idx+=sh;
	return idx;
}
int			first_set_bit16(unsigned short n)//idx of LSB
{
	int sh=((n&((1<<8)-1))==0)<<3,	idx =sh;	n>>=sh;
		sh=((n&((1<<4)-1))==0)<<2,	idx+=sh;	n>>=sh;
		sh=((n&((1<<2)-1))==0)<<1,	idx+=sh;	n>>=sh;
		sh= (n&((1<<1)-1))==0,		idx+=sh;
	return idx;
}
int			floor_log2(unsigned long long n)//idx of MSB
{
	int sh=(n>=1ULL	<<32)<<5,	logn =sh; n>>=sh;
		sh=(n>=1	<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1	<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1	<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1	<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1	<< 1;		logn+=sh;
	return logn;
}
int			ceil_log2(unsigned long long n)
{
	int sh=(n>1ULL<<31)<<5,	logn =sh; n>>=sh;
		sh=(n>1   <<15)<<4;	logn+=sh, n>>=sh;
		sh=(n>1   << 7)<<3;	logn+=sh, n>>=sh;
		sh=(n>1   << 3)<<2;	logn+=sh, n>>=sh;
		sh=(n>1   << 1)<<1;	logn+=sh, n>>=sh;
		sh= n>1   << 0;		logn+=sh;
	return logn;

	//int fl=floor_log2(n);
	//for(int k=0;k<fl;++k)//redundant O(n) loop
	//{
	//	if(n>>k&1)
	//	{
	//		++fl;
	//		break;
	//	}
	//}
	//return fl;
}
int			floor_log10(double x)
{
	static const double pmask[]=//positive powers
	{
		1, 10,		//10^2^0
		1, 100,		//10^2^1
		1, 1e4,		//10^2^2
		1, 1e8,		//10^2^3
		1, 1e16,	//10^2^4
		1, 1e32,	//10^2^5
		1, 1e64,	//10^2^6
		1, 1e128,	//10^2^7
		1, 1e256	//10^2^8
	};
	static const double nmask[]=//negative powers
	{
		1, 0.1,		//1/10^2^0
		1, 0.01,	//1/10^2^1
		1, 1e-4,	//1/10^2^2
		1, 1e-8,	//1/10^2^3
		1, 1e-16,	//1/10^2^4
		1, 1e-32,	//1/10^2^5
		1, 1e-64,	//1/10^2^6
		1, 1e-128,	//1/10^2^7
		1, 1e-256	//1/10^2^8
	};
	int logn, sh;
	if(x<=0)
		return 0x80000000;
	if(x>=1)
	{
		logn=0;
		sh=(x>=pmask[17])<<8;	logn+=sh, x*=nmask[16+(sh!=0)];
		sh=(x>=pmask[15])<<7;	logn+=sh, x*=nmask[14+(sh!=0)];
		sh=(x>=pmask[13])<<6;	logn+=sh, x*=nmask[12+(sh!=0)];
		sh=(x>=pmask[11])<<5;	logn+=sh, x*=nmask[10+(sh!=0)];
		sh=(x>=pmask[9])<<4;	logn+=sh, x*=nmask[8+(sh!=0)];
		sh=(x>=pmask[7])<<3;	logn+=sh, x*=nmask[6+(sh!=0)];
		sh=(x>=pmask[5])<<2;	logn+=sh, x*=nmask[4+(sh!=0)];
		sh=(x>=pmask[3])<<1;	logn+=sh, x*=nmask[2+(sh!=0)];
		sh= x>=pmask[1];		logn+=sh;
		return logn;
	}
	logn=-1;
	sh=(x<nmask[17])<<8;	logn-=sh;	x*=pmask[16+(sh!=0)];
	sh=(x<nmask[15])<<7;	logn-=sh;	x*=pmask[14+(sh!=0)];
	sh=(x<nmask[13])<<6;	logn-=sh;	x*=pmask[12+(sh!=0)];
	sh=(x<nmask[11])<<5;	logn-=sh;	x*=pmask[10+(sh!=0)];
	sh=(x<nmask[9])<<4;		logn-=sh;	x*=pmask[8+(sh!=0)];
	sh=(x<nmask[7])<<3;		logn-=sh;	x*=pmask[6+(sh!=0)];
	sh=(x<nmask[5])<<2;		logn-=sh;	x*=pmask[4+(sh!=0)];
	sh=(x<nmask[3])<<1;		logn-=sh;	x*=pmask[2+(sh!=0)];
	sh= x<nmask[1];			logn-=sh;
	return logn;
}
double		power(double x, int y)
{
	double mask[]={1, 0}, product=1;
	if(y<0)
		mask[1]=1/x, y=-y;
	else
		mask[1]=x;
	for(;;)
	{
		product*=mask[y&1], y>>=1;	//67.7
		if(!y)
			return product;
		mask[1]*=mask[1];
	}
	return product;
}
double		_10pow(int n)
{
	static double *mask=0;
	int k;
//	const double _ln10=log(10.);
	if(!mask)
	{
		mask=(double*)malloc(616*sizeof(double));
		for(k=-308;k<308;++k)		//23.0
			mask[k+308]=power(10., k);
		//	mask[k+308]=exp(k*_ln10);//inaccurate
	}
	if(n<-308)
		return 0;
	if(n>307)
		return _HUGE;
	return mask[n+308];
}


//complex functions
void		print_matrix_debug(const double *data, int w, int h)
{
	int kx, ky;

	printf("\n");
	for(ky=0;ky<h;++ky)
	{
		for(kx=0;kx<w;++kx)
			printf("%4g + i%4g", data[(w*ky+kx)<<1], data[((w*ky+kx)<<1)+1]);
		printf("\n");
	}
	printf("\n");
}

double		c_abs2(Comp const *z){return z->r*z->r+z->i*z->i;}
void		c_inv(Comp *dst, const Comp *z)
{
	double invabs2=1/c_abs2(z);
	dst->r=z->r* invabs2;
	dst->i=z->i*-invabs2;
}
void		c_mul(Comp *dst, const Comp *a, const Comp *b)//(a0+a1i)(b0+b1i)=a0b0-a1b1+i(a0b1+a1b0)	//dst, a and b can point to the same address
{
	double
		r=a->r*b->r-a->i*b->i,
		i=a->r*b->i+a->i*b->r;
	dst->r=r, dst->i=i;
}
void		c_div(Comp *dst, const Comp *a, const Comp *b)
{
	double invabsb2=1/c_abs2(b);
	double
		r=(a->r*b->r+a->i*b->i)*invabsb2,
		i=(a->i*b->r-a->r*b->i)*invabsb2;
	dst->r=r, dst->i=i;
}
void		c_mod(Comp *dst, const Comp *a, const Comp *b)//dst=a-floor(a/b)*b
{
	Comp z;
	c_div(&z, a, b);
	z.r=floor(z.r);
	z.i=floor(z.i);
	c_mul(&z, &z, b);
	dst->r=a->r-z.r;
	dst->i=a->i-z.i;
}

void		impl_addbuffers(Comp *dst, Comp const *a, Comp const *b, int count)
{
	Comp *p=dst, *end=dst+count;
	__m128d va, vb;
	for(;p<end;++p, ++a, ++b)
	{
		va=_mm_load_pd((double*)a);
		vb=_mm_load_pd((double*)b);
		va=_mm_add_pd(va, vb);
		_mm_store_pd((double*)p, va);
	}
	//int k;
	//for(k=0;k<size;++k)
	//	dst[k]=a[k]+b[k];
}
void		impl_subbuffers(Comp *dst, Comp const *a, Comp const *b, int count)
{
	__m128d va, vb;
	Comp *p=dst, *end=dst+count;
	for(;p<end;p+=2, a+=2, b+=2)
	{
		va=_mm_load_pd((double*)a);
		vb=_mm_load_pd((double*)b);
		va=_mm_sub_pd(va, vb);
		_mm_store_pd((double*)p, va);
	}
	//int k;
	//for(k=0;k<size;++k)
	//	dst[k]=a[k]-b[k];
}
void		impl_negbuffer(Comp *dst, Comp const *a, int count)
{
	__m128d va, zero=_mm_setzero_pd();
	Comp *p=dst, *end=dst+count;
	for(;p<end;p+=2, a+=2)
	{
		va=_mm_load_pd((double*)a);
		va=_mm_sub_pd(zero, va);
		_mm_store_pd((double*)p, va);
	}
	//int k;
	//for(k=0;k<size;++k)
	//	dst[k]=-a[k];
}
void		impl_buf_plus_val(Comp *dst, Comp const *a, const Comp *val, int count)
{
	__m128d va, vb=_mm_load_pd((double*)val);
	Comp *p=dst, *end=dst+count;
	for(;p<end;p+=2, a+=2)
	{
		va=_mm_load_pd((double*)a);
		va=_mm_add_pd(va, vb);
		_mm_store_pd((double*)p, va);
	}
	//int k;
	//for(k=0;k<size;++k)
	//	dst[k]=buf[k]+val;
}
void		impl_val_minus_buf(Comp *dst, const Comp *val, Comp const *b, int count)
{
	__m128d va=_mm_load_pd((double*)val), vb;
	Comp *p=dst, *end=dst+count;
	for(;p<end;p+=2, b+=2)
	{
		vb=_mm_load_pd((double*)b);
		vb=_mm_sub_pd(va, vb);
		_mm_store_pd((double*)p, vb);
	}
	//int k;
	//for(k=0;k<size;++k)
	//	dst[k]=val-buf[k];
}
void		impl_buf_mul_val(Comp *dst, Comp const *a, const Comp *val, int count)
{
	Comp *p=dst, *end=dst+count;
	for(;p<end;p+=2, a+=2)
		c_mul(p, a, val);
	//int k;
	//for(k=0;k<size;++k)
	//	dst[k]=buf[k]*val;
}

void		impl_ref(Comp *m, short dx, short dy)
{
#ifdef _DEBUG
	Comp pivot;
#endif
	Comp coeff, temp;
	int mindim=dx<dy?dx:dy, it, ky, kx;
	for(it=0;it<mindim;++it)//iteration
	{
		for(ky=it;ky<dy;++ky)//find pivot
		{
			if(c_abs2(m+((dx*ky+it)<<1))>1e-10)
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				break;
			}
		}
		if(ky<dy)
		{
			if(ky!=it)
				for(kx=0;kx<dx;++kx)//swap rows
					coeff=m[dx*it+kx], m[dx*it+kx]=m[dx*ky+kx], m[dx*ky+kx]=coeff;
			for(++ky;ky<dy;++ky)//subtract pivot row
			{
				c_div(&coeff, m+dx*ky+it, m+dx*it+it);
				//coeff=m[dx*ky+it]/m[dx*it+it];
				for(kx=it;kx<dx;++kx)
				{
					c_mul(&temp, &coeff, m+dx*it+kx);
					m[dx*ky+kx].r-=temp.r;
					m[dx*ky+kx].i-=temp.i;
				}
					//m[dx*ky+kx]-=coeff*m[dx*it+kx];
			}
		}
	}
}
void		impl_rref(Comp *m, short dx, short dy)
{
#ifdef _DEBUG
	Comp pivot;
#endif
	Comp coeff, temp;
	int mindim=dx<dy?dx:dy, it, ky, kx;
	for(it=0;it<mindim;++it)//iteration
	{
		for(ky=it;ky<dy;++ky)//find pivot
		{
			if(c_abs2(m+((dx*ky+it)<<1))>1e-10)
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				break;
			}
		}
		if(ky<dy)
		{
			if(ky!=it)
			{
				for(kx=0;kx<dx;++kx)//swap rows
					coeff=m[dx*it+kx], m[dx*it+kx]=m[dx*ky+kx], m[dx*ky+kx]=coeff;
			}
			for(ky=0;ky<dy;++ky)
			{
				if(ky==it)//normalize pivot
				{
					c_inv(&coeff, m+dx*it+it);
					//coeff=1/m[dx*it+it];
					for(kx=it;kx<dx;++kx)
						c_mul(m+dx*it+kx, m+dx*it+kx, &coeff);
						//m[dx*it+kx]*=coeff;
				}
				else//subtract pivot row from all other rows
				{
					c_div(&coeff, m+dx*ky+it, m+dx*it+it);
					//coeff=m[dx*ky+it]/m[dx*it+it];
					for(kx=it;kx<dx;++kx)
					{
						c_mul(&temp, &coeff, m+dx*it+kx);
						m[dx*ky+kx].r-=temp.r;
						m[dx*ky+kx].i-=temp.i;
					}
						//m[dx*ky+kx]-=coeff*m[dx*it+kx];
				}
				//print_matrix_debug(m, dx, dy);
			}
		}
	}
}
Comp		impl_det(Comp *m, int dx)//m is destroyed
{
	int k, dxplus1=dx+1;
	Comp result;
	impl_ref(m, dx, dx);

	result=m[0];//accumulate diagonal
	for(k=1;k<dx;++k)
		c_mul(&result, &result, m+dxplus1*k);
	return result;
}
void		impl_matinv(Comp *m, short dx)//resize m to (dy * 2dx) temporarily,		dx==dy always
{
	int k, dy=dx, size=dx*dy;
			//print_matrix_debug(m, dx<<1, dy);//
	for(k=size-dx;k>=0;k-=dx)//expand M into [M, 0]
	{
		CMEMCPY(m+(k<<1), m+k, dx);
		//memcpy(m+(k<<1), m+k, dx*sizeof(double));
		CMEMZERO(m+(k<<1)+dx, dx);
		//memset(m+(k<<1)+dx, 0, dx*sizeof(double));
	}
			//print_matrix_debug(m, dx<<1, dy);//
	for(k=0;k<dx;++k)//add identity: [M, I]
		m[(dx<<1)*k+dx+k].r=1;
			//print_matrix_debug(m, dx<<1, dy);//
	impl_rref(m, dx<<1, dy);//[I, M^-1]
			//print_matrix_debug(m, dx<<1, dy);//
	for(k=0;k<size;k+=dx)//pack M^-1
		CMEMCPY(m+k, m+(k<<1)+dx, dx);
		//memcpy(m+k, m+(k<<1)+dx, dx*sizeof(double));
			//print_matrix_debug(m, dx<<1, dy);//
}
void		impl_matmul(Comp *dst, const Comp *m1, const Comp *m2, int h1, int w1h2, int w2)
{
	int kx, ky, kv;
	Comp *C, temp;
#ifdef DEBUG_MEMORY
	printf("impl_matmul():\n");
	print_matrix_debug(dst, w2, h1);
	printf("  =\n");
	print_matrix_debug(A, w1h2, h1);
	printf("  *\n");
	print_matrix_debug(B, w2, w1h2);
#endif
	for(ky=0;ky<h1;++ky)
	{
		for(kx=0;kx<w2;++kx)
		{
			C=dst+w2*ky+kx;
			C->r=0, C->i=0;
			for(kv=0;kv<w1h2;++kv)
			{
				c_mul(&temp, m1+w1h2*ky+kv, m2+w2*kv+kx);
				C->r+=temp.r, C->i+=temp.i;
			}
				//*C+=m1[w1h2*ky+kv]*m2[w2*kv+kx];
		}
	}
}
void		impl_transpose(Comp *dst, const Comp *src, int src_dx, int src_dy)
{
	int ky, kx;

	for(ky=0;ky<src_dy;++ky)
		for(kx=0;kx<src_dx;++kx)
			dst[src_dy*kx+ky]=src[src_dx*ky+kx];
}
void		impl_tensor(Comp *dst, const Comp *m1, const Comp *m2, int dx1, int dy1, int dx2, int dy2)
{
	int dx, kx, ky, kx2, ky2;
	Comp coeff;

	dx=dx1*dx2;
	for(ky=0;ky<dy1;++ky)
	{
		for(kx=0;kx<dx1;++kx)
		{
			coeff=m1[dx1*ky+kx];
			for(ky2=0;ky2<dy2;++ky2)
				for(kx2=0;kx2<dx2;++kx2)
					c_mul(dst+dx*(dy2*ky+ky2)+dx2*kx+kx2, &coeff, m2+dx2*ky2+kx2);
					//dst[dx*(dy2*ky+ky2)+dx2*kx+kx2]=coeff*m2[dx2*ky2+kx2];
		}
	}
}
void		impl_matdiv(Comp *dst, const Comp *num, Comp *den, int num_dy, int dx)//dst & num: num_dy*dx,  den: dx*(dx*2)		den is destroyed
{
	impl_matinv(den, dx);

	impl_matmul(dst, num, den, num_dy, dx, dx);
}
void		impl_matdiv_back(Comp *dst, Comp *den, const Comp *num, int dy, int num_dx)//den: dy*(dy*2),  dst & num: dy*num_dx		den is destroyed
{
	impl_matinv(den, dy);

	impl_matmul(dst, den, num, dy, dy, num_dx);
}
void		impl_matpow(Comp *dst, Comp *m1, int e, int dx)//dst: dx*dx,  m1: dx*(dx*2)		calculates m1^e,	m1 is destroyed
{
	int k, size=dx*dx;
	Comp *temp=m1+size;

	if(e<0)//negative exponent
	{
		e=-e;
		impl_matinv(m1, dx);
	}
	//memcpy(x, m1, size*sizeof(double));

	CMEMZERO(dst, size);
	//memset(dst, 0, size*sizeof(double));//identity matrix
	for(k=0;k<size;k+=dx+1)
		dst[k].r=1;

	for(;;)
	{
		if(e&1)
		{
			impl_matmul(temp, dst, m1, dx, dx, dx);
			CMEMCPY(dst, temp, size);
			//memcpy(dst, temp, size*sizeof(double));
		}
		e>>=1;
		if(!e)
			break;
		impl_matmul(temp, m1, m1, dx, dx, dx);
		CMEMCPY(m1, temp, size);
		//memcpy(m1, temp, size*sizeof(double));
	}
}

//polynomials are stored as they are read	p[0]*x^(n-1) + p[1]*x^(n-2) + ...p[n-1] of degree n-1
void		impl_polmul(Comp *res, Comp const *A, Comp const *B, int asize, int bsize, int add)//res has correct size of (asize+bsize-1)
{//int add:  -1: subtract from res;  0: overwrite res;  1: add to res
	int dst_size=asize+bsize-1, k, k2;
	Comp coeff, sign={1, 0}, temp;
	if(add==-1)
		sign.r=-1;
	if(!add)
		CMEMZERO(res, dst_size);
		//memset(res, 0, dst_size*sizeof(double));
	for(k=0;k<asize;++k)//schoolbook O(n2)
	{
		c_mul(&coeff, &sign, A+asize-1-k);
		//coeff=sign*A[asize-1-k];
		for(k2=0;k2<bsize;++k2)
		{
			c_mul(&temp, &coeff, B+bsize-1-k2);
			res[dst_size-1-k-k2].r+=temp.r;
			res[dst_size-1-k-k2].i+=temp.i;
		}
			//res[dst_size-1-k-k2]+=coeff*B[bsize-1-k2];
	}
}
#if 0
void		impl_poldiv(double const *num, double const *den, int nnum, int nden, double *q, double *r)//qsize=nnum-nden+1, rsize=nnum (actually min(num, nden-1))
{
	int qsize=nnum-nden+1, k, k2;
	CMEMCPY(r, num, nnum);
	//memcpy(r, num, nnum*sizeof(double));
	for(k=0;k<qsize;++k)//schoolbook O(n2)
	{
		q[k]=r[k]/den[0];
		for(k2=0;k2<nden;++k2)
			r[k+k2]-=q[k]*den[k2];
	}
}
int			impl_countleadzeros(double *pol, int oldsize)//returns newsize
{
	const double tolerance=1e-10;
	int k;
	for(k=0;k<oldsize;++k)
		if(fabs(pol[k])>tolerance)
			break;
	if(k==oldsize)
		return k-1;
	return k;
}
int			impl_polgcd(double *res, double const *A, double const *B, int asize, int bsize)//res size = min(asize, bsize), returns gcd size
{
	const double tolerance=1e-10;
	double *cvec, *q, *r[3];
	int qsize=maximum(asize, bsize), rsize[3]={asize, bsize, minimum(asize, maximum(bsize-1, 1))}, leadingzeros;
	V_CONSTRUCT(double, cvec, qsize+rsize[0]+rsize[1]+rsize[2], 0, 0);
	q=cvec;
	r[0]=cvec+qsize;
	r[1]=r[0]+rsize[0];
	r[2]=r[1]+rsize[1];
	CMEMCPY(r[0], A, asize);
	CMEMCPY(r[1], B, bsize);
	//memcpy(r[0], A, asize*sizeof(double));
	//memcpy(r[1], B, bsize*sizeof(double));
	for(;;)
	{
		impl_poldiv(r[0], r[1], rsize[0], rsize[1], q, r[2]);
		leadingzeros=impl_countleadzeros(r[2], rsize[2]);
		r[2]+=leadingzeros, rsize[2]-=leadingzeros;//advance pointer & decrease size, by leading zero count
		if(rsize[2]==1&&fabs(*r[2])<tolerance)
			break;
		q=r[2], r[2]=r[0], r[0]=r[1], r[1]=q, q=cvec;//cycle pointers & their sizes {r[0], <- r[1], <- r[2]}
		leadingzeros=rsize[2], rsize[2]=rsize[0], rsize[0]=rsize[1], rsize[1]=leadingzeros;
	}
	CMEMCPY(res, r[1], rsize[1]);
	//memcpy(res, r[1], rsize[1]*sizeof(double));
	v_destroy(&cvec);
	return rsize[1];
}
int			impl_fracremoveleadzeros(Object *frac)
{
	int numzeros=impl_countleadzeros(frac->r, frac->dx),
		denzeros=impl_countleadzeros(frac->r+frac->dx, frac->dy);
	if(frac->dx-numzeros<=0)
		--numzeros;
	if(numzeros>0)
	{
		mem_shiftback(frac->r, frac->r+numzeros, (frac->dx-numzeros)*sizeof(double));
		frac->dx-=numzeros;
	}
	if(numzeros+denzeros>0)
	{
		mem_shiftback(frac->r+frac->dx, frac->r+frac->dx+numzeros+denzeros, (frac->dy-denzeros)*sizeof(double));
		frac->dy-=denzeros;
		v_resize(&frac->r, frac->dx+frac->dy, 0);
		return 1;
	}
	return 0;
}
void		impl_fracsimplify(Object *frac)
{
	int gcdsize=minimum(frac->dx, frac->dy);
	double *cvec,
		*gcd, *q, *r;
	int k, size, newnumsize, newdensize, resized;
	double inv_a0;

	V_CONSTRUCT(double, cvec, gcdsize+frac->dx+frac->dy, 0, 0);
	gcd=cvec;

	gcdsize=impl_polgcd(gcd, frac->r, &v_at(frac->r, frac->dx), frac->dx, frac->dy);//find gcd(num, den)

	if(gcdsize>1)
	{
		q=cvec+gcdsize;
		newnumsize=frac->dx-gcdsize+1;
		r=q+newnumsize;
		impl_poldiv(frac->r, gcd, frac->dx, gcdsize, q, r);//num/gcd
		CMEMCPY(frac->r, q, newnumsize);
		//memcpy(frac->r, q, newnumsize*sizeof(double));
		
		newdensize=frac->dy-gcdsize+1;
		r=q+newdensize;
		impl_poldiv(frac->r+frac->dx, gcd, frac->dy, gcdsize, q, r);//den/gcd
		CMEMCPY(frac->r+newnumsize, q, newdensize);
		//memcpy(frac->r+newnumsize, q, newdensize*sizeof(double));

		frac->dx=newnumsize;
		frac->dy=newdensize;
	//	mem_shiftback(frac->r+qsize, frac->r+frac->dx, frac->dy*sizeof(double));//shift back den X
	}
	v_destroy(&cvec);

	resized=impl_fracremoveleadzeros(frac);//resize once
	if(!resized&&gcdsize>1)
		v_resize(&frac->r, frac->dx+frac->dy, 0);
	size=v_size(&frac->r);
	inv_a0=v_at(frac->r, frac->dx);
	if(inv_a0)
	{
		inv_a0=1/inv_a0;
		for(k=0;k<size;++k)
			v_at(frac->r, k)*=inv_a0;
	}
}
void		impl_fracadd(Object *dst, Object *A, Object *B, int subtract)
{
	int n1d2=A->dx+B->dy-1, n2d1=A->dy+B->dx-1;//n1/d1 + n2/d2 = (n1d2 + n2d1)/d1d2
	int rnum, rden;
	double *frac;
	int sign=!subtract-subtract;

	rnum=maximum(n1d2, n2d1);
	rden=A->dy+B->dy-1;
	V_CONSTRUCT(double, frac, rnum+rden, 0, 0);
	if(n1d2>n2d1)
	{
		impl_polmul(		frac,					A->r,		&v_at(	B->r, B->dx),	A->dx, B->dy, 0		);//n1d2 (larger)
		impl_polmul(&v_at(	frac, n1d2-n2d1), &v_at(A->r, A->dx),		B->r,			A->dy, B->dx, sign	);//n2d1
	}
	else
	{
		impl_polmul(&v_at(	frac, n2d1-n1d2),		A->r,		&v_at(	B->r, B->dx),	A->dx, B->dy, 0		);//n1d2
		impl_polmul(		frac,			&v_at(	A->r, A->dx),		B->r,			A->dy, B->dx, sign	);//n2d1 (larger)
	}
	impl_polmul(&v_at(frac, rnum), &v_at(A->r, A->dx), &v_at(B->r, B->dx), A->dy, B->dy, 0);

	dst->type=T_FRAC;
	dst->dx=rnum;
	dst->dy=rden;
	v_move(&dst->r, &frac);
//	impl_fracsimplify(dst);
}
void		impl_fracmul(Object *dst, Object *A, Object *B)
{
	int rnum, rden;
	double *frac;

	rnum=A->dx+B->dx-1;
	rden=A->dy+B->dy-1;
	V_CONSTRUCT(double, frac, rnum+rden, 0, 0);
	impl_polmul(		frac,				A->r,				B->r,			A->dx, B->dx, 0);
	impl_polmul(&v_at(	frac, rnum), &v_at(	A->r, A->dx), &v_at(B->r, B->dx),	A->dy, B->dy, 0);
	
	dst->type=T_FRAC;
	dst->dx=rnum;
	dst->dy=rden;
	v_move(&dst->r, &frac);
//	impl_fracsimplify(dst);
}
void		impl_fracdiv(Object *dst, Object *A, Object *B)
{
	int dx=A->dx+B->dy-1,
		dy=A->dy+B->dx-1;
	double *frac;
	V_CONSTRUCT(double, frac, dx+dy, 0, 0);
	impl_polmul(		frac,				A->r,		&v_at(	B->r, B->dx),	A->dx, B->dy, 0);
	impl_polmul(&v_at(	frac, dx), &v_at(	A->r, A->dx),		B->r,			A->dy, B->dx, 0);

	dst->type=T_FRAC;
	dst->dx=dx;
	dst->dy=dy;
	v_move(&dst->r, &frac);
//	impl_fracsimplify(dst);
}
void		impl_fracpow(Object *dst, Object *A, Object *B, Token *fn)
{
	int e=extract_int(B, fn);
	Object p={T_FRAC, 1, 1}, x;

	V_CONSTRUCT(double, p.r, 2, 0, 0);
	p.r[0]=p.r[1]=1;
	if(e>=0)
		obj_assign(&x, A);
	else//negative exponent
	{
		x.type=A->type;
		x.dx=A->dy;
		x.dy=A->dx;
		V_CONSTRUCT(double, x.r, x.dx+x.dy, 0, 0);
		CMEMCPY(x.r, &v_at(A->r, A->dx), A->dy);
		//memcpy(x.r, &v_at(A->r, A->dx), A->dy*sizeof(double));
		CMEMCPY(x.r+x.dx, A->r, A->dx);
		//memcpy(x.r+x.dx, A->r, A->dx*sizeof(double));
	}
	for(;;)
	{
		if(e&1)
			impl_fracmul(&p, &p, &x);
		e>>=1;
		if(!e)
			break;
		impl_fracmul(&x, &x, &x);
	}
	v_destroy(&x);
	obj_assign(dst, &p);
}
void		impl_ldiv(Object *dst, Object *frac, int nsteps)
{
	int nnum, nden, hsize, k, kh;
	double *num, *den, *q, *h, *h0;
	//int k2;//

	nnum=frac->dx, nden=frac->dy, hsize=maximum(nnum, nden);
	num=frac->r, den=frac->r+nnum;
	V_CONSTRUCT(double, q, nsteps+hsize, 0, 0);//qsize+hsize
	MEMZERO(q, nsteps+hsize);
	//memset(q, 0, (nsteps+hsize)*sizeof(double));
	h=q+nsteps;

	for(k=0;k<nsteps;++k)//standard programming
	{
		h0=h+mod(-k, hsize);
		*h0=k==0;//impulse
		for(kh=1;kh<nden;++kh)
			*h0-=den[kh]*h[mod(kh-k, hsize)];
		for(kh=0;kh<nnum;++kh)
			q[k]+=num[kh]*h[mod(kh-k, hsize)];
		//for(k2=0;k2<nsteps+hsize;++k2)//
		//	printf("  %g", q[k2]);
		//printf("\n");
	}
	v_resize(&q, nsteps, 0);
	dst->type=T_MATRIX;
	dst->dx=1;
	dst->dy=nsteps;
	v_move(&dst->r, &q);
}
void		impl_matpolsubs(double *res, const double *m, double *m_temp, const double *coeff, int ncoeff, int dx)
{
	int size=dx*dx, k, k2;

	//p(x) = (...(c[0]*x + c[1])*x + c[2])...)*x + c[n-1]
	memset(res, 0, size*sizeof(double));
	for(k=0;k<size;k+=dx+1)//res=c[0]*In;
		res[k]=coeff[0];
	for(k=1;k<ncoeff;++k)
	{
		CMEMCPY(m_temp, res, size);
		//memcpy(m_temp, res, size*sizeof(double));//res*=m
		impl_matmul(res, m_temp, m, dx, dx, dx);

		for(k2=0;k2<size;k2+=dx+1)//res+=c[k]*In
			res[k2]+=coeff[k];
	}
}
#endif