//mac_math.c - ModCalc math operations
//Copyright (C) 2021  Ayman Wagih Mohsen, unless source link provided
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

#define _CRT_SECURE_NO_WARNINGS
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<tmmintrin.h>
#include	<vector>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
#include	"mac_memory.h"
static const char file[]=__FILE__;

//	#define	DEBUG_NULLSPACE

#ifdef __linux__
#define		sprintf_s	snprintf
#define		vsprintf_s	vsnprintf
#define		scanf_s		scanf
#define		__rdtsc		__builtin_ia32_rdtsc
#define		_HUGE		HUGE_VAL
char		_getch();
#endif
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
int			floor_log10(unsigned long long n)//idx of MSD
{
	static const unsigned long long mask[]=
	{
		1, 10,
		1, 100,
		1, 10000,
		1, 100000000,
		1, 10000000000000000,
	};
	int logn=0, sh;
	sh=(n>=mask[9])<<4,	logn+=sh, n/=mask[8+(sh!=0)];
	sh=(n>=mask[7])<<3,	logn+=sh, n/=mask[6+(sh!=0)];
	sh=(n>=mask[5])<<2,	logn+=sh, n/=mask[4+(sh!=0)];
	sh=(n>=mask[3])<<1,	logn+=sh, n/=mask[2+(sh!=0)];
	sh= n>=mask[1],		logn+=sh;
	return logn;
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


//modular arithmetic
Int		zero, one(1), minusone(-1), ten(10);
std::vector<Int> powersof10;
int		floor_log10(Int const &n)
{
	if(n==zero)
		return 0x80000000;//
	if(!powersof10.size())
	{
		powersof10.push_back(ten);
		powersof10.push_back(Int(100));
	}
	int log2log10n=0;
	for(;;)//grow: keep squaring until surpasses n
	{
		if(n.abs_cmp(powersof10[log2log10n])<0)
			break;
		++log2log10n;
		while(log2log10n>=(int)powersof10.size())
		{
			powersof10.push_back(powersof10.back());
			powersof10.back()*=powersof10.back();
		}
	}
	if(!log2log10n)
		return 0;
	--log2log10n;
	int log10n=1<<log2log10n;
	Int ee[]={powersof10[log2log10n], powersof10[log2log10n]};
	for(--log2log10n;log2log10n>=0;--log2log10n)//binary search: converge to log10n
	{
		ee[1]*=powersof10[log2log10n];
		if(n.abs_cmp(ee[1])>0)
			ee[0]=ee[1], log10n|=1<<log2log10n;
	}
	return log10n;
}
int		floor_log2(Int const &n)
{
	if(!n.data||n==zero)
		return 0;
	int logn=(int)(n.data->used_size<<5)-1;
	for(;!n.getbit(logn);--logn);
	return logn;
	//int logn=n.data->used_size-1;
	//int MSS=n.data->x[logn];
	//logn<<=5;
	//logn|=31;//point at most significant possible bit
	//for(;!(MSS>>(logn&31)&1);--logn);//find MSB
	//return logn;
}
void	impl_ldiv(Int num, Int den, Int &Q, Int &R)//Q can't be num		don't do this at home: impl_ldiv(t, den, t, r)
{
	if(den==zero)
		return;
	char numneg=num.is_neg(), denneg=den.is_neg(), neg=numneg!=denneg;
	num.abs();
	den.abs();
	int lognum=floor_log2(num);
	Q.setzero();
	R.setzero();
	for(int i=lognum;i>=0;--i)//https://en.wikipedia.org/wiki/Division_algorithm
	{
		R<<=1;
		R.data->x[0]|=num.getbit(i);
		if(R>=den)
		{
			R-=den;
			Q.setbit(i);
		}
	}
	if(neg)
		Q.negate();
	if(neg)
		R-=den;
	if(numneg)
		R.negate();
}

void	rotate_back(Int *data, int size)
{
	Int temp=std::move(data[0]);
	for(int k=0;k<size-1;++k)
		data[k]=std::move(data[k+1]);
	data[size-1]=std::move(temp);
}
Int		impl_gcd(Int a, Int b)//by copy
{
	Int c=a%b;
	while(c!=zero)
		a=b, b=c, c=a%b;
	return b;
}
Int		impl_lcm(Int *arr, int count)//array is destroyed
{
	Int LCM=1, test=2, t;
	for(int it=0;;)
	{
		bool divided=false;
		bool all_ones=true;
		for(int k=0;k<count;++k)
		{
			auto &val=arr[k];
			if(val!=1)
			{
				all_ones=false;
				t=val/test;
				if(t*test==val)
				{
					divided=true;
					val=t;
				}
			}
		}
		if(all_ones)
			break;
		if(divided)
			LCM*=test;
		else
		{
			test+=2-!it;
			++it;
		}

	}
	return LCM;
}
Int		impl_intEEA(Int const &x, Int const &n)
{
//	Int q, r[3]={n, x}, t[3]={Int(), 1};

	Int q, r[3]={n, x%n}, t[3]={Int(), 1};
	if(x.is_neg()&&r[1]!=x)
		printf("Warning: %s was reduced mod %s before calculating the Extended Euclidean Algorithm.\n", x.to_string().c_str(), n.to_string().c_str());
	for(;;)
	{
		impl_ldiv(r[0], r[1], q, r[2]);
		if(r[2]==zero)
			break;
		t[2]=t[0]-q*t[1];

		rotate_back(r, 3);
		rotate_back(t, 3);
	}
	r[1].abs();
	if(r[1]!=1)
	{
		printf("Note: %s has no inverse mod %s\n", x.to_string().c_str(), n.to_string().c_str());
		return zero;
	}
	t[1]%=n;
	return t[1];
}
//void	impl_mulmodu(Int &dst, Int const &a, Int const &b, Int const &n)
//{
//	const Int m64=0xFFFFFFFFFFFFFFFF;
//	auto pa=(unsigned*)&a, pb=(unsigned*)&b;
//	Int t;
//	impl_mod(t, m64, n);
//	t+=1;
//
//	dst=((Int)pa[0]*pb[0]%n;
//}
Int		impl_powmod(Int const &x, Int const &e, Int const &n)
{
	Int mask, product(1);
	if(e.is_neg())
	{
		mask=impl_intEEA(x, n);
		if(mask==zero)//no inverse
			return zero;
	}
	else
		mask=x;
	int nbits=floor_log2(e)+1;
	Int Q, R;
	for(int k=0;k<nbits;++k)
	{
		//if(e.is_odd())
		if(e.getbit(k))
		{
			product*=mask;
			impl_ldiv(product, n, Q, R);
			product=R;
		}
		//	product.mulmod(mask, n);
		//shr_lint(e.data, 1);
		//if(e==zero)
		//	break;
		mask*=mask;
		impl_ldiv(mask, n, Q, R);
		mask=R;
		//mask.mulmod(mask, n);
	}
	return product;
}

Int		isqrt(Int const &x)//https://en.wikipedia.org/wiki/Integer_square_root
{
	Int L, R=x+1;
	int k=0;
	while(L!=R-1)
	{
		Int M=L+R;
		M>>=1;
		if(M*M<=x)
			L=std::move(M);
		else
			R=std::move(M);
	}
	return L;
}
Int		icbrt(Int const &x)
{
	Int L, R=x+1;
	while(L!=R-1)
	{
		Int M=L+R;
		M>>=1;
		if(M*M*M<=x)
			L=std::move(M);
		else
			R=std::move(M);
	}
	return L;
}
int		impl_isprime(Int &n)
{
	Int sqrtn=isqrt(n);
	for(Int test=2;test<sqrtn;++test)
		if(n%test==zero)
			return false;
	return true;
}
bool	impl_mrtest(Int const &n, Int const &a)
{
	int k=0;
	int nbits=(int)(n.data->used_size<<5);
	Int nm1=n-1, q=nm1;
	for(k=0;k<nbits&&!(q.data->x[0]&1);++k, q>>=1);
	Int aq=impl_powmod(a, q, n);
	if(aq==one)
		return true;
	for(int j=0;j<k;++j)
	{
		if(aq==nm1)
			return true;
		aq.mulmod(aq, n);
	}
	return false;
}
Int		impl_genprime(Int start)
{
	if(!start.is_odd())
		++start;
	Int witness=2+rand();
	for(;;start+=2)
	{
		if(impl_mrtest(start, witness))
			break;
	}
	return start;
}
Int		impl_randlog(Int const &L)
{
	srand((unsigned)__rdtsc());
	int nc=(L.to_int()+31)>>5;
	Int r;
	r.setzero(nc);
	for(int k=0;k<nc;++k)
		r.data->x[k]=rand()<<30|rand()<<15|rand();
	return r;
}
Int		impl_randmod(Int const &modulus)
{
	srand((unsigned)__rdtsc());
	int logn=floor_log2(modulus);
	Int r=impl_randlog(logn);
	Int Q, R;
	impl_ldiv(r, modulus, Q, R);
	return R;
}
void	impl_factorize(Int n, IntVec &factors, IntVec *powers, IntVec *rfactors)
{
	factors.clear();
	if(powers)
		powers->clear();
	if(rfactors)
		rfactors->clear();
	auto t1=__rdtsc();
	//for(Int test(2), latest=test;n>one;)
	for(Int test(2);n>one;)
	{
		auto t2=__rdtsc();
		if(t2-t1>0x100000000)
			printf("\rtesting: "), test.print(), printf("\t");
		if(n%test==zero)
		{
			if(powers)
			{
				if(!factors.size()||test!=factors.back())
				{
					factors.push_back(test);
					powers->push_back(1);
					if(rfactors)
						rfactors->push_back(test);
				}
				else
				{
					++powers->back();
					if(rfactors)
						rfactors->back()*=test;
				}
			}
			else
				factors.push_back(test);
		/*	switch(type)
			{
			case F_REPEATED:
				factors.push_back(test);
				break;
			case F_UNIQUE:
				if(!factors.size()||test!=factors.back())
					factors.push_back(test);
				break;
			case F_UNIQUE_WITH_COUNT:
				if(!factors.size()||test!=latest)
				{
					factors.push_back(test);
					factors.push_back(1);
					latest=test;
				}
				else
					++factors.back();
				break;
			case F_POWERS:
				if(!factors.size()||test!=latest)
				{
					factors.push_back(test);
					latest=test;
				}
				else
					factors.back()*=test;
				break;
			case F_POWERS_WITH_COUNT:
				if(!factors.size()||test!=latest)
				{
					factors.push_back(test);
					factors.push_back(1);
					latest=test;
				}
				else
				{
					factors[factors.size()-2]*=test;
					++factors.back();
				}
				break;
			}//*/
			n/=test;
		}
		else
			++test;
	}
	auto t3=__rdtsc();
	if(t3-t1>0x100000000)
		printf("\n");
}
Int		impl_totient(Int const &n, IntVec &unique_factors, IntVec &powers)
{
	Int phi=1, temp;
	for(int k=0;k<(int)unique_factors.size();++k)
	{
		if(powers[k]>1)
		{
			temp=unique_factors[k];
			temp^=powers[k]-1;
			temp*=unique_factors[k]-1;
			phi*=temp;
		}
		else
			phi*=unique_factors[k]-1;
	}
	//Int phi=unique_factors[0]-1;
	//for(int k=1;k<(int)unique_factors.size();++k)
	//	phi*=unique_factors[k]-1;
	return phi;
}
Int		impl_carmichael(Int const &n, IntVec &unique_factors, IntVec &powers)//unique_factors is destroyed
{
	Int lambda=1, phi;
	for(int k=0, nfactors=(int)unique_factors.size();k<nfactors;++k)
	{
		auto &prime=unique_factors[k], &power=powers[k];
		bool half=prime==2&&power>2;
		--power;
		if(power>zero)
		{
			phi=prime^power;
			--prime;
			phi*=prime;
		}
		else
		{
			--prime;
			phi=prime;
		}
		if(half)
			phi>>=1;
		unique_factors[k]=phi;
	}
	return impl_lcm(unique_factors.data(), (int)unique_factors.size());
}
bool	impl_proots(Int n, IntVec &roots)//https://en.wikipedia.org/wiki/Primitive_root_modulo_n#Finding_primitive_roots
{
	roots.clear();
	n.abs();
	if(n>0x10000)
	{
		auto str=n.to_string();
		printf(
			"Warning:\n"
			"  About to calculate the primitive roots of %s\n"
			"  This may take a long time and the result array may get large.\n"
			"Enter a ZERO to proceed: ", str.c_str());
		int c=0;
		scanf("%d", &c);
		if(c)
			return false;
	}
	IntVec factors, powers;
	impl_factorize(n, factors, &powers);
	Int phi=impl_totient(n, factors, powers),
		lambda=impl_carmichael(n, factors, powers);
	if(lambda<phi)//no primitive roots
	{
		auto sn=n.to_string(), sphi=phi.to_string(), slambda=lambda.to_string();
		printf(
			"Note:\n"
			"%s has no primitive roots\n"
			"  totient = %s\n"
			"  lambda  = %s\n", sn.c_str(), sphi.c_str(), slambda.c_str());
		return true;
	}
	impl_factorize(phi, factors, &powers);
	for(int k=0;k<(int)factors.size();++k)
	{
		auto &pfact=factors[k];
		pfact=phi/pfact;
	}
	for(Int g=2, prod;g<n;++g)
	{
		Int gcd=impl_gcd(g, n);
		if(gcd>1)
			continue;
		bool found=true;
		for(int k=0;k<(int)factors.size();++k)
		{
			prod=g;
			prod.powmod(factors[k], n);
			if(prod==1)
			{
				found=false;
				break;
			}
		}
		if(found)
			roots.push_back(g);
	}
	return true;
}
Int		impl_dlog(Int const &x, Int const &b, Int const &n)//r=dlog(x, b, n)	b^r = x (mod n)
{
	Int r, x2=1;
	for(;;)
	{
		if(x2==x)
			break;
		x2.mulmod(b, n);
		++r;
	}
	return r;
}


#define		G_BUF_SIZE	1024
static char	g_buf[G_BUF_SIZE]={0};
int			query_double(double x, int *point)
{
	if(point)
		*point=sprintf_s(g_buf, G_BUF_SIZE, "%lld", (long long)floor(x));
	if(fabs(x)<1e-10)
		return 1;
	return sprintf_s(g_buf, G_BUF_SIZE, "%g", x);
}
int			print_double(double x, int point_pos, int total)//point_pos==0: no leading spaces, total==0: no trailing spaces
{
	double a=abs(x);
	long long i=(long long)a;
	int nbefore=sprintf_s(g_buf, G_BUF_SIZE, "%lld", i);
	int nafter=sprintf_s(g_buf, G_BUF_SIZE, "%g", a-i)-2;
	int neg=x<0;
	int nspaces=point_pos-nbefore-neg;
	if(nspaces<0)
		nspaces=0;
	if(total&&nspaces+neg+nafter>total)
		return printf("%*s%.*lf", nspaces, "", total-(nspaces+neg), x);
	int printed=printf("%*s%g", nspaces, "", x);
	if(printed<total)
		printed+=printf("%*s", total-printed, "");
	return printed;
}
void		dec2frac(double x, double error, int *i, int *num, int *den)//https://stackoverflow.com/questions/5124743/algorithm-for-simplifying-decimal-to-fractions
{
	int lower_n, upper_n, middle_n,
		lower_d, upper_d, middle_d;
	int n=(int)floor(x);
	x-=n;
	if(x<error)
	{
		*i=n, *num=0, *den=1;
		return;
	}
	if(1-error<x)
	{
		*i=n+1, *num=0, *den=1;
		return;
	}
	lower_n=0, upper_n=1;
	lower_d=1, upper_d=1;
	for(;;)
	{
		middle_n=lower_n+upper_n;//The middle fraction is (lower_n + upper_n) / (lower_d + upper_d)
		middle_d=lower_d+upper_d;
		if(middle_d*(x+error)<middle_n)//If x + error < middle
			upper_n=middle_n, upper_d=middle_d;
		else if(middle_n<(x-error)*middle_d)//Else If middle < x - error
			lower_n=middle_n, lower_d=middle_d;
		else
			break;
	}
	*i=n, *num=middle_n, *den=middle_d;
}
int			print_double_frac(double x, int point_pos, int total)
{
	int i, num, den, printed=0;

	if(fabs(x)<1e-10)
		printed=printf("%*s0", (point_pos-1)&-(point_pos>0), "");
	//else if(x==floor(x))
	//	printed=printf("%*g", point_pos, x);
	else
	{
		dec2frac(x, 1e-10, &i, &num, &den);

		if(num)
		{
			if(i>-5&&i<5)
				printed=printf("%*d/%d", point_pos, num+i*den, den);
			else if(den<1000)
				printed=printf("%*d+%d/%d", point_pos, i, num, den);
			else
				printed=print_double(x, point_pos, total);
		}
		else
			printed=printf("%*d", point_pos, i);
	}
	if(printed<total)
		printf("%*s", total-printed, "");
	return printed;
}

void		print_matrix_debug(Int const *data, int w, int h)
{
	int kx, ky;

	printf("\n");
	for(ky=0;ky<h;++ky)
	{
		for(kx=0;kx<w;++kx)
			printf(" %7lld", data[w*ky+kx].to_longlong());
		printf("\n");
	}
	printf("\n");
}

void		impl_addbuffers(Int *dst, Int const *a, Int const *b, int count)
{
	//Comp *p=dst, *end=dst+count;
	//__m128i va, vb;
	//for(;p<end;++p, ++a, ++b)
	//{
	//	va=_mm_load_si128((__m128i*)a);
	//	vb=_mm_load_si128((__m128i*)b);
	//	va=_mm_add_epi64(va, vb);
	//	_mm_store_si218((__m128i*)p, va);
	//}
	int k;
	for(k=0;k<count;++k)
		dst[k]=a[k]+b[k];
}
void		impl_subbuffers(Int *dst, Int const *a, Int const *b, int count)
{
	//Comp *p=dst, *end=dst+count;
	//__m128i va, vb;
	//for(;p<end;++p, ++a, ++b)
	//{
	//	va=_mm_load_si128((__m128i*)a);
	//	vb=_mm_load_si128((__m128i*)b);
	//	va=_mm_sub_epi64(va, vb);
	//	_mm_store_si218((__m128i*)p, va);
	//}
	int k;
	for(k=0;k<count;++k)
		dst[k]=a[k]-b[k];
}
void		impl_negbuffer(Int *dst, Int const *a, int count)
{
	//__m128d va, zero=_mm_setzero_pd();
	//Comp *p=dst, *end=dst+count;
	//for(;p<end;p+=2, a+=2)
	//{
	//	va=_mm_load_pd((double*)a);
	//	va=_mm_sub_pd(zero, va);
	//	_mm_store_pd((double*)p, va);
	//}
	int k;
	for(k=0;k<count;++k)
		dst[k]=-a[k];
}
void		impl_buf_plus_val(Int *dst, Int const *buf, Int const &val, int count)
{
	//__m128d va, vb=_mm_load_pd((double*)val);
	//Comp *p=dst, *end=dst+count;
	//for(;p<end;p+=2, a+=2)
	//{
	//	va=_mm_load_pd((double*)a);
	//	va=_mm_add_pd(va, vb);
	//	_mm_store_pd((double*)p, va);
	//}
	int k;
	for(k=0;k<count;++k)
		dst[k]=buf[k]+val;
}
void		impl_val_minus_buf(Int *dst, Int const &val, Int const *buf, int count)
{
	//__m128d va=_mm_load_pd((double*)val), vb;
	//Comp *p=dst, *end=dst+count;
	//for(;p<end;p+=2, b+=2)
	//{
	//	vb=_mm_load_pd((double*)b);
	//	vb=_mm_sub_pd(va, vb);
	//	_mm_store_pd((double*)p, vb);
	//}
	int k;
	for(k=0;k<count;++k)
		dst[k]=val-buf[k];
}
void		impl_buf_mul_val(Int *dst, Int const *buf, const Int &val, int count)
{
	//Comp *p=dst, *end=dst+count;
	//for(;p<end;p+=2, a+=2)
	//	c_mul(p, a, val);
	int k;
	for(k=0;k<count;++k)
		dst[k]=buf[k]*val;
}

void		impl_ref(Int *m, short dx, short dy, Int const &n)
{
#ifdef _DEBUG
	Int pivot;
#endif
	Int coeff, temp;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it]!=zero)
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(ky<dy)
		{
			if(ky>it)
				for(kx=0;kx<dx;++kx)//swap rows
					coeff=m[dx*it+kx], m[dx*it+kx]=m[dx*ky+kx], m[dx*ky+kx]=coeff;
			for(++ky;ky<dy;++ky)//subtract pivot row
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
		}
	}
}
void		impl_rref(Int *m, short dx, short dy, Int const &n)
{
#ifdef _DEBUG
	Int pivot;
#endif
	Int coeff, temp;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		kpivot=-1;
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it]!=zero)
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(kpivot==-1)
			continue;
		if(kpivot>npivots-1)
		{
			for(kx=0;kx<dx;++kx)//swap rows
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)//normalize pivot row
			{
				coeff=1/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*kpivot+kx]*=coeff;
			}
			else//subtract pivot row from all other rows
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
			//print_matrix_debug(m, dx, dy);//
		}
	}
}
void		impl_rref2(Int *m, short dx, short dy, Int const &n)//may be buggy
{
	Int pivot;
	Int coeff, temp;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		kpivot=-1;
		//printf("Enter pivot ky: ");
		//scanf_s("%d", &ky);
		//pivot=m[dx*ky+it];
		//printf("pivot=");
		//print_value(&pivot);
		//printf("\n");
		for(ky=it;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it]!=zero)
			{
				pivot=m[dx*ky+it];
				printf("pivot=");
				pivot.print();
				printf("\n");
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(kpivot==-1)
			continue;
		if(kpivot>npivots-1)
		{
			for(kx=0;kx<dx;++kx)//swap rows
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)//normalize pivot
			{
				coeff=1/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*kpivot+kx]*=coeff;
			}
			else//subtract pivot row from all other rows
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
			print_matrix_debug(m, dx, dy);
		}
	}
}
Int			impl_det(Int *m, int dx, Int const &n)//m is destroyed
{
	int k, dxplus1=dx+1;
	Int result;

	//print_matrix_debug(m, dx, dx);//
	impl_ref(m, dx, dx, n);
	//print_matrix_debug(m, dx, dx);//

	result=m[0];//multiply diagonals
	for(k=1;k<dx;++k)
		result*=m[dxplus1*k];
	return result;
}
void		impl_matinv(Int *m, short dx, Int const &n)//resize m to (dy * 2dx) temporarily,		dx==dy always
{
	int k, dy=dx, size=dx*dy;
			//print_matrix_debug(m, dx<<1, dy);//
	for(k=size-dx;k>=0;k-=dx)//expand M into [M, 0]
	{
		objmove(m+((ptrdiff_t)k<<1), m+k, dx);
		//MEMCPY(m+(k<<1), m+k, dx);
		setzero(m+((ptrdiff_t)k<<1)+dx, dx);
		//MEMZERO(m+(k<<1)+dx, dx);
	}
			//print_matrix_debug(m, dx<<1, dy);//
	for(k=0;k<dx;++k)//add identity: [M, I]
		m[(dx<<1)*k+dx+k]=1;
			//print_matrix_debug(m, dx<<1, dy);//
	impl_rref(m, dx<<1, dy, n);//[I, M^-1]
			//print_matrix_debug(m, dx<<1, dy);//
	for(k=0;k<size;k+=dx)//pack M^-1
		objmove(m+k, m+((ptrdiff_t)k<<1)+dx, dx);
		//MEMCPY(m+k, m+(k<<1)+dx, dx);
			//print_matrix_debug(m, dx<<1, dy);//
}
void		impl_matmul(Int *dst, const Int *m1, const Int *m2, int h1, int w1h2, int w2)
{
	int kx, ky, kv;
	Int *C, temp;
//	MEMZERO(dst, h1*w2);
//#ifdef DEBUG_MEMORY
//	printf("impl_matmul():\n");
//	print_matrix_debug(dst, w2, h1);
//	printf("  =\n");
//	print_matrix_debug(m1, w1h2, h1);
//	printf("  *\n");
//	print_matrix_debug(m2, w2, w1h2);
//#endif
	for(ky=0;ky<h1;++ky)
	{
		for(kx=0;kx<w2;++kx)
		{
			C=dst+w2*ky+kx;
			C->setzero();
			for(kv=0;kv<w1h2;++kv)
				*C+=m1[w1h2*ky+kv]*m2[w2*kv+kx];
		}
	}
}
void		impl_transpose(Int *dst, const Int *src, int src_dx, int src_dy)
{
	int ky, kx;

	for(ky=0;ky<src_dy;++ky)
		for(kx=0;kx<src_dx;++kx)
			dst[src_dy*kx+ky]=src[src_dx*ky+kx];
}
void		impl_tensor(Int *dst, Int const *m1, const Int *m2, int dx1, int dy1, int dx2, int dy2)
{
	int dx, kx, ky, kx2, ky2;
	Int coeff;

	dx=dx1*dx2;
	for(ky=0;ky<dy1;++ky)
	{
		for(kx=0;kx<dx1;++kx)
		{
			coeff=m1[dx1*ky+kx];
			for(ky2=0;ky2<dy2;++ky2)
				for(kx2=0;kx2<dx2;++kx2)
					dst[dx*(dy2*ky+ky2)+dx2*kx+kx2]=coeff*m2[dx2*ky2+kx2];
		}
	}
}
void		impl_matdiv(Int *dst, Int const *num, Int *den, int num_dy, int dx, Int const &n)//dst & num: num_dy*dx,  den: dx*(dx*2)		den is destroyed
{
	impl_matinv(den, dx, n);

	impl_matmul(dst, num, den, num_dy, dx, dx);
}
void		impl_matdiv_back(Int *dst, Int *den, const Int *num, int dy, int num_dx, Int const &n)//den: dy*(dy*2),  dst & num: dy*num_dx		den is destroyed
{
	impl_matinv(den, dy, n);

	impl_matmul(dst, den, num, dy, dy, num_dx);
}
void		impl_matpow(Int *dst, Int *m1, int e, int dx, Int const &n)//dst: dx*dx,  m1: dx*(dx*2)		calculates m1^e,	m1 is destroyed
{
	int k, size=dx*dx;
	Int *temp=m1+size;

	if(e<0)//negative exponent
	{
		e=-e;
		impl_matinv(m1, dx, n);
	}
	//memcpy(x, m1, size*sizeof(double));

	setzero(dst, size);
	//MEMZERO(dst, size);//identity matrix
	for(k=0;k<size;k+=dx+1)
		dst[k]=1;

	for(;;)
	{
		if(e&1)
		{
			impl_matmul(temp, dst, m1, dx, dx, dx);
			objmove(dst, temp, size);
			//MEMCPY(dst, temp, size);
		}
		e>>=1;
		if(!e)
			break;
		impl_matmul(temp, m1, m1, dx, dx, dx);
		objmove(m1, temp, size);
		//MEMCPY(m1, temp, size);
	}
}

#if 0
void		impl_lu(Int const *m, int n, Int *lower, Int *upper)
{
	Int sum;
	int i, j, k, idx;
	int size=n*n;
	MEMZERO(lower, size), MEMZERO(upper, size);
 
	//Decomposing matrix into Upper and Lower triangular matrix
	for(i=0;i<n;++i)//https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/
	{
		for(k=i;k<n;++k)//Upper Triangular
		{
			//Summation of L(i, j) * U(j, k)
			sum=zero;
			for(j=0;j<i;++j)
				sum+=lower[n*i+j]*upper[n*j+k];
 
			//Evaluating U(i, k)
			idx=n*i+k;
			upper[idx]=m[idx]-sum;
			//upper[i][k] = mat[i][k] - sum;
		}
		for(k=i;k<n;++k)//Lower Triangular
		{
			if(i==k)
				lower[n*i+i]=1;//Diagonal as 1
			else
			{
				//Summation of L(k, j) * U(j, i)
				sum=zero;
				for(j=0;j<i;++j)
					sum+=lower[n*k+j]*upper[n*j+i];
 
				//Evaluating L(k, i)
				sum=m[n*k+i]-sum;
				lower[n*k+i]=(m[n*k+i]-sum)/upper[n*i+i];
			}
		}
	}
}
void		impl_qr(Int const *M, int n, Int *Q, Int *R, double *temp)//temp: a buffer of 2n doubles
{
	//factorizes a REAL square matrix M into Q and R
	//where Q is orthogonal and R is upper triangular
	//https://madrury.github.io/jekyll/update/statistics/2017/10/04/qr-algorithm.html
	//https://github.com/madrury/linalg
	int i, j, k;
	double *col=temp, *uvec=temp+n, d;

	MEMZERO(temp, n<<1);
//	ALLOC(col, n<<1);
//	MEMZERO(col, n<<1);
	uvec=col+n;
	for(i=0;i<n;++i)
	{
		for(k=0;k<n;++k)//copy column i from M
			col[k]=M[n*k+i].r;
		for(j=0;j<i;++j)
		{
			for(k=0;k<n;++k)//copy column j from Q
				uvec[k]=Q[n*k+j].r;
			d=0;
			for(k=0;k<n;++k)//projection of col on uvec
				d+=col[k]*uvec[k];
			for(k=0;k<n;++k)//col is now perpendicular to uvec
				col[k]-=d*uvec[k];
			R[n*i+j].r=d;
		}
		d=0;
		for(k=0;k<n;++k)
			d+=col[k]*col[k];
		d=sqrt(d);
		R[n*i+i].r=d;
		d=1/d;
		for(k=0;k<n;++k)
			Q[n*k+i].r=col[k]*d;
	}
//	FREE(col);
}
#endif

#if 0
Int			c_det22(Int const *M)
{
	return M[0]*M[3]-M[1]*M[2];
}
void		impl_egval2(Int const *M, Int *lambdas)//finds the eigenvalues of a complex 2x2 matrix
{
	Int C[2], temp;

	C[0]=c_det22(M);//C0=det(M)
	C[1]=-M[0]-M[3];//C1=-tr(M)

	//quadratic formula		lambda = C1*-0.5 +- sqrt(C1*C1*0.25-C0)
	C[1].r*=-0.5;
	C[1].i*=-0.5;
	c_mul(&temp, C+1, C+1);
	temp.r-=C->r;
	temp.i-=C->i;
	c_sqrt(&temp, &temp);
	lambdas[0].r=C[1].r-temp.r, lambdas[0].i=C[1].i-temp.i;
	lambdas[1].r=C[1].r+temp.r, lambdas[1].i=C[1].i+temp.i;
}
void		impl_solve_cubic(Int const *coeffs, Int *roots)//finds r[0], r[1], & r[2], the solutions of x^3 + c[2]x^2 + c[1]x + c[0] = 0
{
	//https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
	Comp p=coeffs[2], q=coeffs[1], r=coeffs[0],
		p2, p3, q2, prod, A, B;
	double _3sqrt3=3*sqrt(3.), inv3cbrt2=1./(3*cbrt(2.)), ninth=1./9;
	Comp cm={-0.5, -sqrt(3.)*0.5}, cp={-0.5, sqrt(3.)*0.5};

	c_mul(&p2, &p, &p);
	c_mul(&p3, &p2, &p);
	c_mul(&q2, &q, &q);

	c_mul(&A, &q2, &q);
	A.r*=4;
	A.i*=4;
	c_mul(&prod, &p2, &q2);
	A.r-=prod.r;
	A.i-=prod.i;
	c_mul(&prod, &p3, &r);
	A.r+=4*prod.r;
	A.i+=4*prod.i;
	c_mul(&prod, &p, &q);
	c_mul(&prod, &prod, &r);
	A.r-=18*prod.r;
	A.i-=18*prod.i;
	c_mul(&prod, &r, &r);
	A.r+=27*prod.r;
	A.i+=27*prod.i;
	c_sqrt(&A, &A);
	A.r*=_3sqrt3;
	A.i*=_3sqrt3;
	A.r-=27*r.r;
	A.i-=27*r.i;
	c_mul(&prod, &p, &q);
	A.r+=9*prod.r-2*p3.r;
	A.i+=9*prod.i-2*p3.i;
	c_cbrt(&A, &A);
	A.r*=inv3cbrt2;
	A.i*=inv3cbrt2;
	B.r=3*q.r-p2.r;
	B.i=3*q.i-p2.i;
	c_div(&B, &B, &A);
	B.r*=ninth;
	B.i*=ninth;
	roots[2].r=roots[1].r=roots[0].r=-p.r*(1./3);
	roots[2].i=roots[1].i=roots[0].i=-p.i*(1./3);
	roots[0].r+=A.r-B.r;
	roots[0].i+=A.i-B.i;
	c_mul_add(roots+1, &A, &cm);
	c_mul_sub(roots+1, &B, &cp);
	c_mul_add(roots+2, &A, &cp);
	c_mul_sub(roots+2, &B, &cm);
}
void		impl_egval3(Int const *M, Int *lambdas)//finds the eigenvalues of a complex 3x3 matrix
{
	Comp C[3], temp;
	
	C[2].r=-(M[0].r+M[4].r+M[8].r);//C[2] = -tr(M) = -(M[0]+M[4]+M[8]);
	C[2].i=-(M[0].i+M[4].i+M[8].i);

	//C[1] = cof0+cof4+cof8 = M[4]*M[8]-M[5]*M[7] + M[0]*M[8]-M[2]*M[6] + M[0]*M[4]-M[1]*M[3];
	c_mul(C+1, M+4, M+8);
	c_mul_sub(C+1, M+5, M+7), C[0]=C[1];
	c_mul_add(C+1, M+0, M+8);
	c_mul_sub(C+1, M+2, M+6);
	c_mul_add(C+1, M+0, M+4);
	c_mul_sub(C+1, M+1, M+3);

	//C[0] = -det(M) = -(M[0]*(M[4]*M[8]-M[5]*M[7]) - M[1]*(M[3]*M[8]-M[5]*M[6]) + M[2]*(M[3]*M[7]-M[4]*M[6]));
	c_mul(C, C, M);
	c_mul(&temp, M+3, M+8);
	c_mul_sub(&temp, M+5, M+6);
	c_mul_sub(C, &temp, M+1);
	c_mul(&temp, M+3, M+7);
	c_mul_sub(&temp, M+4, M+6);
	c_mul_add(C, &temp, M+2);
	C->r=-C->r, C->i=-C->i;

	impl_solve_cubic(C, lambdas);
}
int			is_upper_triangular(Int *M, int n)
{
	int k, k2;
	const double tolerance=1e-10;
	for(k=0;k<n;++k)
		for(k2=0;k2<k;++k)
			if(fabs(M[n*k+k2].r)>tolerance||fabs(M[n*k+k2].i)>tolerance)
				return 0;
	return 1;
}
void		impl_egval(Int const *M, int n, Int *D, int it_limit)
{
	int size=n*n, it=0;
	Comp *CALLOC(temp, size*2+n);
	Comp *Q=temp+n, *R=temp+n+size;
	CMEMCPY(D, M, size);
	do
	{
		impl_qr(D, n, Q, R, (double*)temp);
		impl_matmul(D, R, Q, n, n, n);
		++it;
	}while(it<it_limit&&!is_upper_triangular(D, n));
	CFREE(temp);
}
int			impl_nullspace(Int *M, int dx, int dy, Int *solution, char *dep_flags, short *row_idx)
{
	//M is rref'ed
	//solution allocated size dy*dy, pre-memset solution to zero
	//p_flags & row_idx both size dx
	//returns number of vectors in solution
	int kx, kxdst, ky, keq, kfree, idx, idx2, nvec;
#ifdef DEBUG_NULLSPACE
	printf("Before RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	impl_rref(M, dx, dy);
#ifdef DEBUG_NULLSPACE
	printf("After RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	memset(dep_flags, 0, dx);
	memset(row_idx, 0, dx*sizeof(short));
	for(ky=0;ky<dy;++ky)//find pivots (dependent variables)
	{
		for(kx=ky;kx<dx;++kx)
		{
			idx=dx*ky+kx;
			if(fabs(M[idx].r)>1e-10||fabs(M[idx].i)>1e-10)
				break;
		}
		if(kx<dx)
			dep_flags[kx]=1, row_idx[ky]=kx;
		else
			break;
	}
	nvec=dx-ky;
	for(ky=0, keq=0, kfree=0;ky<dx;++ky)
	{
		if(dep_flags[ky])//pivot, dependent variable
		{
			for(kx=0, kxdst=0;kx<dx;++kx)
			{
				if(dep_flags[kx])
					continue;
				idx=dx*ky+kxdst, idx2=dx*keq+kx;
				solution[idx].r=-M[idx2].r;
				solution[idx].i=-M[idx2].i;
				++kxdst;
			}
			++keq;
		}
		else//free variable
		{
			idx=dx*ky;
			CMEMZERO(solution+idx, nvec);
			solution[idx+kfree].r=1;
			++kfree;
		}
#ifdef DEBUG_NULLSPACE
		//printf("Nullspace row %d:\n", ky);
		print_matrix_debug(solution+dx*ky, nvec, 1);//
#endif
	}
#if 0
	for(kx=0, kfree=0;kx<dx;++kx)//find solution vectors
	{
		if(!p_flags[kx])//free variable
		{
			for(ky=0, ky2=0;ky<dx;++ky)
			{

			}
		/*	for(kxdst=0;kxdst<dx;++kxdst)
			{
				idx=dx*kxdst+kfree;
				if(kxdst==kfree)
				{
					solution[idx].r=1;
					solution[idx].i=0;
				}
				else
				{
					solution[idx].r=M[].r;
					solution[idx].i=M[].i;
				}
			}//*/
			++kfree;
		}
	}
#endif
	return nvec;
}
int			impl_egvec(Int const *M, int n, Int const *lambdas, Int *S)
{
	int kv, kx, nvec, size=n*n, again;
	Comp *CALLOC(temp, size);
	CMEMZERO(S, size);
	char *dep_flags=(char*)malloc(n);
	short *row_idx=(short*)malloc(n*sizeof(short));
	for(kv=0, nvec=0;kv<n&&nvec<n;++kv)//get eigenvectors
	{
		again=0;
		for(kx=0;kx<kv;++kx)//check for repeated eigenvalues
		{
			if(fabs(lambdas[kx].r-lambdas[kv].r)<1e-10&&fabs(lambdas[kx].i-lambdas[kv].i)<1e-10)
			{
				again=1;
				break;
			}
		}
		if(again)
			continue;
		CMEMCPY(temp, M, size);
		for(kx=0;kx<n;++kx)
		{
			temp[(n+1)*kx].r-=lambdas[kv].r;
			temp[(n+1)*kx].i-=lambdas[kv].i;
		}
		nvec+=impl_nullspace(temp, n, n, S+nvec, dep_flags, row_idx);
		//print_matrix_debug(S, n, n);//
		//for(ky=0;ky<n;++ky)
		//{
		//	for(kx=ky;kx<n;++kx)
		//		if(fabs(temp[n*ky+kx].r)>1e-10||fabs(temp[n*ky+kx].i)>1e-10)
		//			break;
		//}
	}
	free(dep_flags), free(row_idx);
	return nvec;
}
int			impl_diag2(Int const *M, Int *invS, Int *D, Int *S)//diagonalizes a complex 2x2 matrix
{
	int k, nvec;
	Comp C[2], temp, M2[4];

	c_det22(C, M);			//C0=det(M)
	C[1].r=-(M[0].r+M[3].r);//C1=tr(M)
	C[1].i=-(M[0].i+M[3].i);

	//quadratic formula		lambda = C1*-0.5 +- sqrt(C1*C1*0.25-C0)
	C[1].r*=-0.5;
	C[1].i*=-0.5;
	c_mul(&temp, C+1, C+1);
	temp.r-=C->r;
	temp.i-=C->i;
	c_sqrt(&temp, &temp);
	CMEMZERO(D, 4);
	D[0].r=C[1].r-temp.r, D[0].i=C[1].i-temp.i;
	D[3].r=C[1].r+temp.r, D[3].i=C[1].i+temp.i;

	nvec=0;
	for(k=0;k<2&&nvec<2;++k)//get eigenvectors
	{
		CMEMCPY(M2, M, 4);
		M2[0].r-=D[3*k].r, M2[0].i-=D[3*k].i;//kth diagonal
		M2[3].r-=D[3*k].r, M2[3].i-=D[3*k].i;
		impl_rref(M2, 2, 2);
		if(M2[0].r)
		{
			S[nvec].r=-M2[1].r;
			S[nvec].i=-M2[1].i;
			S[nvec+2].r=1;
			S[nvec+2].i=0;
			++nvec;
		}
		else if(M2[1].r||M2[1].i)
		{
			S[nvec].r=1;
			S[nvec].i=0;
			S[nvec+2].r=0;
			S[nvec+2].i=0;
			++nvec;
		}
		else
		{
			S[nvec  ].r=1, S[nvec  ].i=0;
			S[nvec+1].r=0, S[nvec+1].i=0;
			S[nvec+2].r=0, S[nvec+2].i=0;
			S[nvec+3].r=1, S[nvec+3].i=0;
			nvec+=2;
		}
	/*	printf("Lambda %d = ", k);
		print_value(D+3*k);
		printf("\n");
		print_matrix_debug(M2, 2, 2);
		printf("\n");//*/
	}

	for(k=0;k<2;++k)//normalize S
	{
		temp.r=0, temp.i=0;
		for(nvec=0;nvec<2;++nvec)
			c_mul_add(&temp, S+2*nvec+k, S+2*nvec+k);
		if(fabs(temp.r)<1e-10||fabs(temp.i)<1e-10)//eigenvector is zero
			continue;
		c_inv(&temp, &temp);
		for(nvec=0;nvec<2;++nvec)
			c_mul(S+2*nvec+k, S+2*nvec+k, &temp);
	}

	//invert S
	c_mul(C, S, S+3);
	c_mul(C+1, S+1, S+2);
	temp.r=C[0].r-C[1].r;
	temp.i=C[0].i-C[1].i;
	if(fabs(temp.r)<1e-10&&fabs(temp.i)<1e-10)//S is singular
	{
		CMEMZERO(invS, 4);
		return 0;
	}
	invS[0]=S[3], invS[3]=S[0];
	invS[1].r=-S[1].r;
	invS[1].i=-S[1].i;
	invS[2].r=-S[2].r;
	invS[2].i=-S[2].i;
	c_inv(&temp, &temp);
	for(k=0;k<4;++k)//divide by determinant
		c_mul(invS+k, invS+k, &temp);
	return 1;
}
#endif

//polynomials are stored as they are read	p[0]*x^(n-1) + p[1]*x^(n-2) + ...p[n-1] of degree n-1
void		impl_polmul(Int *res, Int const *A, Int const *B, int asize, int bsize, int add)//res has correct size of (asize+bsize-1)
{//int add:  -1: subtract from res;  0: overwrite res;  1: add to res
	int dst_size=asize+bsize-1, k, k2;
	Int coeff, sign={1, 0}, temp;
	if(add==-1)
		sign=-1;
	if(!add)
		setzero(res, dst_size);
		//MEMZERO(res, dst_size);
	for(k=0;k<asize;++k)//schoolbook O(n2)
	{
		coeff=sign*A[asize-1-k];
		for(k2=0;k2<bsize;++k2)
			res[dst_size-1-k-k2]+=coeff*B[bsize-1-k2];
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