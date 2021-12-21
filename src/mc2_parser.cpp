//mc2_parser.cpp - MCalc recursive parser
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

#include	<stdarg.h>
#include	<math.h>
#include	<map>
#include	<algorithm>
//#include	<conio.h>
#include	"mc2.h"
const char	file[]=__FILE__;

double		g_tolerance=1e-15;
bool		print_fractions=false;
void		Matrix::print()
{
	if(!data)
	{
		printf("nullptr\n");
		//printf("Error: matrix.data == nullptr\n");
		return;
	}
	auto prind_func=print_fractions?print_double_frac:print_double;
	if(flags==M_SCALAR)//no brackets
	{
		if(abs(data->i)>1e-10)//complex scalar
		{
			prind_func(data->r, 4, 0);//12
			if(data->i>0)
				printf("+");
			prind_func(data->i, 0, 0);
			printf("i\n");
		}
		//	printf("%4g+%4gi\n", data->r, data->i);
		else if(abs(data->r)>1e-10)
		{
			prind_func(data->r, 4, 0);
			printf("\n");
		}
		//	printf("%4g\n", data->r);
		else
			printf("   0\n");
	}
	else//matrix
	{
		auto colsizes=new char[dx*3];//{neg, point, total}
		memset(colsizes, 0, dx*3);
		for(int ky=0;ky<dy;++ky)
		{
			for(int kx=0;kx<dx;++kx)
			{
				auto &x=data[dx*ky+kx].r;
				int point=0;
				int total=query_double(x, &point);
				auto p=colsizes+3*kx;
				p[0]|=x<0;
				if(p[1]<p[0]-(x<0)+point)
					p[1]=p[0]-(x<0)+point;
				if(p[2]<p[0]-(x<0)+total)
					p[2]=p[0]-(x<0)+total;
			}
		}
		printf("[\n");
		for(int ky=0;ky<dy;++ky)
		{
			for(int kx=0;kx<dx;++kx)
			{
				auto &z=data[dx*ky+kx];
				auto p=colsizes+3*kx;
				printf(" ");
				if(abs(z.i)>1e-10)
				{
					prind_func(z.r, p[1], p[2]);
					if(z.i>0)
						printf("+");
					prind_func(z.i, 0, 0);
					printf("*i");
				}
				//	printf("%4g+%4gi", z.r, z.i);
				else if(abs(z.r)>1e-10)
					prind_func(z.r, p[1], p[2]);
				//	printf("%4g", z.r);
				else
					printf("%*s0", p[2]-1, "");
				//if(kx<dx-1)
				//	printf(",");
			}
			if(ky<dy-1)
				printf("\n");
				//printf(";\n");
		}
		printf("\n]\n");
		delete[] colsizes;
	}
}
std::vector<Matrix> g_answers;
std::map<char*, Matrix> g_vars;
std::vector<std::string> errors;
bool		user_error(const char *format, ...)
{
	va_list args;
	va_start(args, format);
	int printed=vsprintf_s(g_buf, g_buf_size, format, args);
	va_end(args);
	errors.push_back(std::string(g_buf, g_buf+printed));
	return false;
}
bool		user_error2(int start, int end, const char *format, ...)
{
	int printed=0;
	if(start>=end)
		end=start+4;
	if(end>text_size)
		end=text_size;
	if(start<end)
		printed=sprintf_s(g_buf, g_buf_size, "Error: \'%.*s\' ", end-start, text+start);
	else
		printed=sprintf_s(g_buf, g_buf_size, "Error: ");
	va_list args;
	va_start(args, format);
	printed+=vsprintf_s(g_buf+printed, g_buf_size-printed, format, args);
	va_end(args);
	printed+=sprintf_s(g_buf+printed, g_buf_size-printed, ".", end-start, text+start);
	errors.push_back(std::string(g_buf, g_buf+printed));
	return false;
}
bool		my_error(int start, int end)
{
	int printed=sprintf_s(g_buf, g_buf_size, "Error: \'%.*s\' Not implemented yet.", end-start, text+start);
	errors.push_back(std::string(g_buf, g_buf+printed));
	return false;
}

//utility
void		pack_rows(Comp *buffer, int dx_old, int dx_new, int dy)//dx_new < dx_old
{
	for(int k=1;k<dy;++k)
	{
		CMEMMOVE(buffer+dx_new*k, buffer+dx_old*k, dx_new);
	//	print_matrix_debug(buffer, dx_new, dy);//
	}
}

//inline math
const double _e=exp(1.), _pi=acos(-1.), _torad=_pi/180, _todeg=180/_pi;
inline double  cosd(double x){return cos(x*_torad);}
inline double acosd(double x){return acos(x)*_todeg;}
inline double  sec (double x){return 1/cos(x);}
inline double asec (double x){return acos(1/x);}
inline double  secd(double x){return 1/cos(x*_torad);}
inline double asecd(double x){return acos(1/x)*_todeg;}
inline double  sech(double x){return 1/cosh(x);}
inline double asech(double x){return acosh(1/x);}
inline double  sind(double x){return sin(x*_torad);}
inline double asind(double x){return asin(x)*_todeg;}
inline double  csc (double x){return 1/sin(x);}
inline double acsc (double x){return asin(1/x);}
inline double  cscd(double x){return 1/sin(x*_torad);}
inline double acscd(double x){return asin(1/x)*_todeg;}
inline double  csch(double x){return 1/sinh(x);}
inline double acsch(double x){return asinh(1/x);}
inline double  tand(double x){return tan(x*_torad);}
inline double atand(double x){return atan(x)*_todeg;}
inline double  cot (double x){return 1/tan(x);}
inline double acot (double x){return atan(1/x);}
inline double  cotd(double x){return 1/tan(x*_torad);}
inline double acotd(double x){return atan(1/x)*_todeg;}
inline double  coth(double x){return 1/tanh(x);}
inline double acoth(double x){return atanh(1/x);}

inline double c_distance(Comp const &a, Comp const &b)
{
	double r=a.r-b.r, i=a.i-b.i;
	return r*r+i*i;
}
inline void c_mul_i(Comp &z)
{
	double temp=z.r;
	z.r=-z.i, z.i=temp;
}
inline void c_mul_neg_i(Comp &z)
{
	double temp=z.r;
	z.r=z.i, z.i=-temp;
}
inline void c_sq(Comp &z)
{
	double
		r=z.r*z.r-z.i*z.i,
		i=z.r*z.i;
	z.r=r, z.i=i+i;
}
inline void	c_exp(Comp &z)
{
	double m=exp(z.r);
	z.r=m*cos(z.i);
	z.i=m*sin(z.i);
}
inline void	c_ln(Comp &z)
{
	double
		r=log(sqrt(c_abs2(&z))),
		i=atan2(z.i, z.r);
	z.r=r, z.i=i;
}
inline void c_sqrt(Comp &z)//sqrt(z)=exp(0.5lnz)
{
	if(z.r||z.i)
	{
		c_ln(z);
		z.r*=0.5, z.i*=0.5;
		c_exp(z);
	}
}
inline void c_cbrt(Comp &z)//sqrt(z)=exp(0.5lnz)
{
	if(z.r||z.i)
	{
		c_ln(z);
		z.r*=1./3, z.i*=1./3;
		c_exp(z);
	}
}
inline void	c_cos(Comp &z)
{
	Comp z2;
	c_mul_i(z);
	c_exp(z);
	c_inv(&z2, &z);
	z.r+=z2.r, z.i+=z2.i;
	z.r*=0.5, z.i*=0.5;

	//double
	//	r=cos(z.r)*cosh(z.i),
	//	i=-sin(z.r)*sinh(z.i);
	//z.r=r, z.i=i;
}
inline void	c_acos(Comp &z)
{
	Comp z2=z;
	c_sq(z);
	z.r-=1;
	c_sqrt(z);
	z.r+=z2.r, z.i+=z2.i;
	c_ln(z);
	c_mul_neg_i(z);
}
inline void	c_cosd(Comp &z)
{
	z.r*=_torad, z.i*=_torad;
	c_cos(z);
}
inline void	c_acosd(Comp &z)
{
	c_acos(z);
	z.r*=_todeg, z.i*=_todeg;
}
inline void c_cosh(Comp &z)
{
	Comp z2;
	c_exp(z);
	c_inv(&z2, &z);
	z.r+=z2.r, z.i+=z2.i;
	z.r*=0.5, z.i*=0.5;
}
inline void	c_acosh(Comp &z)
{
	Comp z2=z;
	c_sq(z);
	z.r-=1;
	c_sqrt(z);
	z.r+=z2.r, z.i+=z2.i;
	c_ln(z);
}
inline void c_sec(Comp &z)
{
	c_cos(z);
	c_inv(&z, &z);
}
inline void c_asec(Comp &z)
{
	c_inv(&z, &z);
	c_acos(z);
}
inline void c_secd(Comp &z)
{
	z.r*=_torad, z.i*=_torad;
	c_sec(z);
}
inline void c_asecd(Comp &z)
{
	c_asec(z);
	z.r*=_todeg, z.i*=_todeg;
}
inline void c_sech(Comp &z)
{
	c_cosh(z);
	c_inv(&z, &z);
}
inline void c_asech(Comp &z)
{
	c_inv(&z, &z);
	c_acosh(z);
}
inline void	c_sin(Comp &z)
{
	Comp z2;
	c_mul_i(z);
	c_exp(z);
	c_inv(&z2, &z);
	z.r-=z2.r, z.i-=z2.i;
	z.r*=0.5, z.i*=0.5;
}
inline void	c_asin(Comp &z)
{
	Comp z2=z;
	c_sq(z);
	z.r=1-z.r, z.i=-z.i;
	c_sqrt(z);
	c_mul_i(z2);
	z.r+=z2.r, z.i+=z2.i;
	c_ln(z);
	c_mul_neg_i(z);
}
inline void	c_sind(Comp &z)
{
	z.r*=_torad, z.i*=_torad;
	c_sin(z);
}
inline void	c_asind(Comp &z)
{
	c_asin(z);
	z.r*=_todeg, z.i*=_todeg;
}
inline void c_sinh(Comp &z)
{
	Comp z2;
	c_exp(z);
	c_inv(&z2, &z);
	z.r-=z2.r, z.i-=z2.i;
	z.r*=0.5, z.i*=0.5;
}
inline void	c_asinh(Comp &z)
{
	Comp z2=z;
	c_sq(z);
	z.r+=1;
	c_sqrt(z);
	z.r+=z2.r, z.i+=z2.i;
	c_ln(z);
}
inline void c_csc(Comp &z)
{
	c_sin(z);
	c_inv(&z, &z);
}
inline void c_acsc(Comp &z)
{
	c_inv(&z, &z);
	c_asin(z);
}
inline void c_cscd(Comp &z)
{
	z.r*=_torad, z.i*=_torad;
	c_csc(z);
}
inline void c_acscd(Comp &z)
{
	c_acsc(z);
	z.r*=_todeg, z.i*=_todeg;
}
inline void c_csch(Comp &z)
{
	c_sinh(z);
	c_inv(&z, &z);
}
inline void c_acsch(Comp &z)
{
	c_inv(&z, &z);
	c_asinh(z);
}
inline void c_tan(Comp &z)
{
	c_mul_i(z);
	z.r+=z.r, z.i+=z.i;
	c_exp(z);
	Comp z2=z;
	z.r-=1, z2.r+=1;
	c_div(&z, &z, &z2);
}
inline void c_atan(Comp &z)
{
	Comp z2=z;
	z.i+=1;
	z.r=-z.r, z.i=1-z.i;
	c_div(&z, &z, &z2);
	c_ln(z);
	c_mul_i(z);
	z.r*=0.5, z.i*=0.5;
}
inline void c_tand(Comp &z)
{
	z.r*=_torad, z.i*=_torad;
	c_tan(z);
}
inline void c_atand(Comp &z)
{
	c_atan(z);
	z.r*=_todeg, z.i*=_todeg;
}
inline void c_tanh(Comp &z)
{
	z.r+=z.r, z.i+=z.i;
	c_exp(z);
	Comp z2=z;
	z.r-=1;
	z2.r+=1;
	c_div(&z, &z, &z2);
}
inline void c_atanh(Comp &z)
{
	Comp z2=z;
	z.r+=1;
	z.r=1-z.r, z.i=-z.i;
	c_div(&z, &z, &z2);
	c_ln(z);
	z.r*=0.5, z.i*=0.5;
}
inline void c_cot(Comp &z)
{
	c_tan(z);
	c_inv(&z, &z);
}
inline void c_acot(Comp &z)
{
	c_inv(&z, &z);
	c_atan(z);
}
inline void c_cotd(Comp &z)
{
	z.r*=_torad, z.i*=_torad;
	c_cot(z);
}
inline void c_acotd(Comp &z)
{
	c_acot(z);
	z.r*=_todeg, z.i*=_todeg;
}
inline void c_coth(Comp &z)
{
	c_tan(z);
	c_inv(&z, &z);
}
inline void c_acoth(Comp &z)
{
	c_inv(&z, &z);
	c_atan(z);
}

void		dft_init(int size, Comp *&weights, Comp *&temp, bool inverse)
{
	CALLOC(weights, size*size);
	CALLOC(temp, size);
	double scale_re=1/sqrt(size), scale_im, freq=2*_pi/size;
	if(inverse)
		scale_im=-scale_re;
	else
		scale_im=scale_re;
	for(int ky=0;ky<size;++ky)
	{
		auto row=weights+size*ky;
		for(int kx=0;kx<size;++kx)
		{
			row[kx].r=scale_re*cos(freq*kx*ky);
			row[kx].i=scale_im*sin(freq*kx*ky);
		}
	}
}
void		dft_apply_1D(Comp *data, Comp const *weights, Comp *temp, int size, int stride)
{
	for(int k=0;k<size;++k)
	{
		auto row=weights+size*k;
		temp[k].r=temp[k].i=0;
		for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
			c_mul_add(temp+k, row+kd, data+ks);
	}
	for(int ks=0, kd=0;kd<size;++ks, kd+=stride)
		data[kd]=temp[ks];
}
void		dft_finish(Comp *&weights, Comp *&temp)
{
	CFREE(weights), weights=nullptr;
	CFREE(temp), temp=nullptr;
}

void		dct_init(int size, Comp *&weights, Comp *&temp)
{
	CALLOC(weights, size*size);
	CALLOC(temp, size);
	double scale=sqrt(2./size), freq=_pi/size;
	for(int ky=0;ky<size;++ky)
	{
		auto row=weights+size*ky;
		for(int kx=0;kx<size;++kx)
		{
			row[kx].r=scale*cos(freq*(kx+0.5)*ky);
			row[kx].i=0;
		}
	}
}
void		idct_init(int size, Comp *&weights, Comp *&temp)
{
	CALLOC(weights, size*size);
	CALLOC(temp, size);
	double scale=sqrt(2./size), freq=_pi/size;
	for(int ky=0;ky<size;++ky)
	{
		auto row=weights+size*ky;
		for(int kx=0;kx<size;++kx)
		{
			if(!kx)
				row[kx].r=scale*0.5;
			else
				row[kx].r=scale*cos(freq*kx*(ky+0.5));
			row[kx].i=0;
		}
	}
}

Comp		eval_poly(Comp const *coeff, int order, Comp &x)
{
	Comp result=coeff[0];
	for(int k=1;k<order;++k)
	{
		c_mul(&result, &result, &x);
		result.r+=coeff[k].r;
		result.i+=coeff[k].i;
	}
	return result;
}
bool		polmul(Matrix &dst, Matrix const &a, Matrix const &b, int idx0)
{
	if(a.dy!=1||b.dy!=1)
		return user_error2(idx0, idx, "Expected a row vector");
	auto CALLOC(temp, a.dx+b.dx-1);
	impl_polmul(temp, a.data, b.data, a.dx, b.dx, 0);
	if(dst.data)
		CFREE(dst.data);
	dst.data=temp;
	dst.dx=a.dx+b.dx-1;
}

inline bool	check_scalar(Matrix &m, int idx0)
{
	if(m.dx!=1||m.dy!=1)
		return user_error2(idx0, idx, "Expected a scalar");
	return true;
}
inline bool	check_real(Matrix &m, int idx0)
{
	if(!check_scalar(m, idx0))
		return false;
	if(abs(m.data->i)>1e-10)
		return user_error2(idx0, idx, "Expected a real value");
}
inline bool	get_int(Matrix &m, int idx0, int &i)
{
	if(m.dx!=1||m.dy!=1)
		return user_error2(idx0, idx, "Expected a scalar");
	if(abs(m.data->i)>1e-10)
		return user_error2(idx0, idx, "Expected a real value");
	i=(int)floor(m.data->r+0.5);
	return true;
}
inline void cmatrix_plus_cscalar(Matrix &m, Matrix &s)
{
	__m128d vb=_mm_load_pd((double*)s.data);
	for(auto p=m.data, end=m.end();p<end;++p)
	{
		__m128d va=_mm_load_pd((double*)p);
		va=_mm_add_pd(va, vb);
		_mm_store_pd((double*)p, va);
	}
}
inline bool	obj_add(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dy==1&&m2.dx==1)//addition of a scalar
	{
		cmatrix_plus_cscalar(m, m2);
		//__m128d vb=_mm_load_pd(m2.data);
		//for(auto p=m.data, end=m.end();p<end;p+=2)
		//{
		//	__m128d va=_mm_load_pd(p);
		//	va=_mm_add_pd(va, vb);
		//	_mm_store_pd(p, va);
		//}
		//for(int k=0, size=m.dy*m.dx;k<size;++k)
		//{
		//	__m128d va=_mm_load_pd(m.data+(k<<1));
		//	va=_mm_add_pd(va, vb);
		//	_mm_store_pd(m.data+(k<<1), va);
		//}
		//for(int k=0, size=m.dy*m.dx;k<size;++k)
		//{
		//	m.data[(k<<1)  ]+=m2.data[0];
		//	m.data[(k<<1)+1]+=m2.data[1];
		//}
	}
	else if(m.dy==1&&m.dx==1)//addition of a scalar
	{
		cmatrix_plus_cscalar(m2, m);
		//for(int k=0, size=m2.dy*m2.dx;k<size;++k)
		//	m2.data[k]+=m.data[0];
		m.move2temp(m2);
	}
	else//matrix addition
	{
		if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		impl_addbuffers(m.data, m.data, m2.data, m.dy*m.dx);
	}
	return true;
}
inline bool	obj_sub(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dy==1&&m2.dx==1)//subtraction of a scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
		{
			m.data[k].r-=m2.data->r;
			m.data[k].i-=m2.data->i;
		}
	}
	else if(m.dy==1&&m.dx==1)//scalar minus matrix
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
		{
			m2.data[k].r=m.data->r-m2.data[k].r;
			m2.data[k].i=m.data->i-m2.data[k].i;
		}
		m.move2temp(m2);
	}
	else
	{
		if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		impl_subbuffers(m.data, m.data, m2.data, m.dy*m.dx);
	}
	return true;
}
inline bool	obj_mul(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.flags==M_SCALAR)//multiplication by scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
			c_mul(m.data+k, m.data+k, m2.data);
	}
	else if(m.flags==M_SCALAR)//multiplication by scalar
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
			c_mul(m2.data+k, m2.data+k, m.data);
		m.move2temp(m2);
	}
	else//matrix multiplication
	{
		if(m.dx!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
		auto CALLOC(temp, m.dy*m2.dx);
		//double *temp=(double*)malloc(m.dy*m2.dx*sizeof(double));
		impl_matmul(temp, m.data, m2.data, m.dy, m.dx, m2.dx);
		CFREE(m.data);
		m.data=temp;
		m.dx=m2.dx;
	}
	return true;
}
inline bool	obj_div(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dx==1&&m2.dy==1)//division by scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
			c_div(m.data+k, m.data+k, m2.data);
	}
	else//matrix division
	{
		if(m2.dx!=m2.dy)//denominator matrix must be square
			return user_error2(idx0, idx, "Denominator matrix must be square: h2=%d != w2=%d", m2.dy, m2.dx);
		if(m.dx!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
		CREALLOC(m2.data, m2.data, m2.dx*m2.dy);
		//m2.data=(double*)realloc(m2.data, m2.dx*m2.dy*2*sizeof(double));
		auto CALLOC(temp, m.dx*m.dy);
		//double *temp=(double*)malloc(m.dx*m.dy*sizeof(double));
		impl_matdiv(temp, m.data, m2.data, m.dy, m2.dx);
		CFREE(m.data);
		m.data=temp;
	}
	return true;
}
inline bool	obj_mod(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dx==1&&m2.dy==1)//m2 can be scalar
	{
		for(int k=0, size=m.dx*m.dy;k<size;++k)
			c_mod(m.data+k, m.data+k, m2.data);
			//m.data[k]=m.data[k]-floor(m.data[k]/m2.data[0])*m2.data[0];
	}
	else
	{
		if(m2.dx!=m.dx||m2.dy!=m.dy)
			return user_error2(idx0, idx, "The modulus operator \'%%\' is element-wise: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		for(int k=0, size=m.dx*m.dy;k<size;++k)
			c_mod(m.data+k, m.data+k, m2.data+k);
			//m.data[k]=m.data[k]-floor(m.data[k]/m2.data[k])*m2.data[k];
	}
	return true;
}

int			parse_incomplete=false;

//precedence (first to last):
//	calls & brackets	f(a)	(a)		[a]		[[a]]
//	postfix & power		a'		a[b]	a^b
//	prefix				-a		+a
//	multicplicative		a*b		a o b	a/b		a\b		a%b		a.*b	a./b
//	additive			a+b		a-b
//	range				a:b
//	relational			a<b		a<=b	a>b		a>=b
//	equality			a==b	a!=b
//	assignment			a=b		a+=b	a-=b	a*=b	a/=b	a%=b

//parser declarations
bool		r_unary(Matrix &m, bool space_sensitive);
bool		r_equality(Matrix &m, bool space_sensitive);
bool		r_assign_expr(Matrix &m, bool space_sensitive);

//parser inspired by OpenC++
bool		r_row_vector(Matrix &m)
{
	lex_skip_space();
	if(!r_equality(m, true))
		return false;
	for(;;)
	{
		int idx0=idx;
		switch(auto t=lex_get(true))
		{
		case T_SPACE:
		case T_COMMA:
			{
				Matrix m2;
				int esize=errors.size();
				if(!r_equality(m2, true))
				{
					if(t==T_SPACE)//trailing whitespace
					{
						if(esize<(int)errors.size())
							errors.erase(errors.begin()+esize, errors.begin()+errors.size());//pop-back errors inserted by last r_equality()
						break;
					}
					return user_error2(idx0, idx, "Expected a matrix column");
				}

				//join matrices horizontally
				if(m.dy!=m2.dy)
					return user_error2(idx0, idx, "Matrices have different heights: h1=%d != h2=%d", m.dy, m2.dy);
				int newdx=m.dx+m2.dx;
				CREALLOC(m.data, m.data, newdx*m.dy);
				//m.data=(double*)realloc(m.data, newdx*m.dy*sizeof(double));
				for(int ky=m.dy-1;ky>=0;--ky)
				{
					CMEMCPY(m.data+newdx*ky, m.data+m.dx*ky, m.dx);
					//memcpy(m.data+newdx*ky, m.data+m.dx*ky, m.dx*sizeof(double));
					CMEMCPY(m.data+newdx*ky+m.dx, m2.data+m2.dx*ky, m2.dx);
					//memcpy(m.data+newdx*ky+m.dx, m2.data+m2.dx*ky, m2.dx*sizeof(double));
				}
				m.dx=newdx;
				m.name=nullptr;
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_postfix(Matrix &m, bool space_sensitive)//pre-allocated
{
	for(;;)
	{
		int idx0=idx;
		switch(lex_get(space_sensitive))
		{
		case T_TRANSPOSE:
			m.name=nullptr;
			if(m.dx!=1||m.dy!=1)
			{
				Comp *data=m.data;
				CALLOC(m.data, m.dx*m.dy);
				std::swap(m.dx, m.dy);
				for(int ky=0;ky<m.dy;++ky)
					for(int kx=0;kx<m.dx;++kx)
						m.get(kx, ky)=GET(data, m.dy, ky, kx);
				//for(int ky=0;ky<m.dy;++ky)
				//{
				//	for(int kx=0;kx<m.dx;++kx)
				//	{
				//		m.RDATA(m.dx, kx, ky)=RDATA(m.dy, ky, kx);
				//		m.IDATA(m.dx, kx, ky)=IDATA(m.dy, ky, kx);
				//	}
				//}
				//int ncomp=1+(m.type==T_COMPLEX);
				//for(int ky=0;ky<m.dy;++ky)
				//	for(int kx=0;kx<m.dx;++kx)
				//		for(int k=0;k<ncomp;++k)
				//			m.data[ncomp*(m.dx*ky+kx)+k]=data[ncomp*(m.dy*kx+ky)+k];
				CFREE(data);
			}
			continue;
		case T_LPR://Matlab-style subscript
			{
				m.name=nullptr;
				if(m.flags==M_UNSPECIFIED_RANGE)
					return user_error2(idx0, idx, "Expected a matrix, got an unspecified range");
				Matrix mky, mkx;
				if(!r_assign_expr(mky, false))
					return false;
				idx0=idx;
				if(lex_get(false)==T_COMMA)
				{
					if(!r_assign_expr(mkx, false))
						return false;
				}
				else
				{
					idx0=idx;
					if(m.dy==1)
					{
						mkx.move2temp(mky);
						mky.flags=M_UNSPECIFIED_RANGE;
					}
					else
						mkx.flags=M_UNSPECIFIED_RANGE;
				}
				if(lex_get(false)!=T_RPR)
					return user_error2(idx0, idx, "Expected a closing parenthesis \')\'");
				
				if(mky.flags==M_UNSPECIFIED_RANGE)
					mky.alloc_ramp(1, m.dy, 1);
				else
				{
					if(mky.dx>1&&mky.dy>1)
						return user_error2(idx0, idx, "Expected a vector");
					if(mky.dx>1)
						std::swap(mky.dx, mky.dy);
				}
				if(mkx.flags==M_UNSPECIFIED_RANGE)
					mkx.alloc_ramp(m.dx, 1, 1);
				else
				{
					if(mkx.dx>1&&mkx.dy>1)
						return user_error2(idx0, idx, "Expected a vector");
					if(mkx.dy>1)
						std::swap(mkx.dx, mkx.dy);
				}
				auto CALLOC(temp, mky.dy*mkx.dx);
				for(int ky=0;ky<mky.dy;++ky)
				{
					auto val=&mky.get(0, ky);
					if(abs(val->i)>1e-10)
					{
						CFREE(temp);
						return user_error2(idx0, idx, "Expected an index, got a complex number ky[%d] = %g + %gi", ky, val->r, val->i);
					}
					int ky2=clamp(0, (int)floor(val->r+0.5)-1, m.dy-1);
					for(int kx=0;kx<mkx.dx;++kx)
					{
						val=&mkx.get(kx, 0);
						if(abs(val->i)>1e-10)
						{
							CFREE(temp);
							return user_error2(idx0, idx, "Expected an index, got a complex number kx[%d] = %g + %gi", kx, val->r, val->i);
						}
						int kx2=clamp(0, (int)floor(val->r+0.5)-1, m.dx-1);
						temp[mkx.dx*ky+kx]=m.data[m.dx*ky2+kx2];
					}
				}
				m.dx=mkx.dx, m.dy=mky.dy;
				CFREE(m.data);
				m.data=temp;
			/*	Comp *temp=nullptr;
				if(m.dy>1)//select/permute rows		[r1; r2; r3]([2 1]) == [r2; r1]
				{
					CALLOC(temp, mky.dy*m.dx);
					for(int ky=0;ky<mky.dy;++ky)
					{
						for(int kx=0;kx<mky.dx;++kx)
						{
							auto &val=mky.get(kx, ky);
							if(abs(val.i)>1e-10)
							{
								CFREE(temp);
								return user_error2(idx0, idx, "Expected an index, got a complex number at [ky=%d][kx=%d]", ky, kx);
							}
							int idx=(int)floor(val.r+0.5);

							temp[mky.dx*ky+kx]=
						}
					}
				}
				else//select/permute columns		[c1 c2 c3]([2 1]) == [c2 c1]
				{
				}//*/
			}
			continue;
		case T_LBRACKET://C-style subscript
			{
				m.name=nullptr;
				if(m.flags==M_UNSPECIFIED_RANGE)
					return user_error2(idx0, idx, "Expected a matrix, got an unspecified range");
				Matrix midx;
				if(!r_assign_expr(midx, false))
					return false;
				if(lex_get(false)!=T_RBRACKET)
					return user_error2(idx0, idx, "Expected a closing bracket \']\'");

				if(midx.flags==M_SCALAR)
				{
					int a_idx=0;
					if(!get_int(midx, idx0, a_idx))
						return false;
					if(a_idx<0)
						a_idx=0;
					if(m.dy>1)//select row
					{
						if(a_idx>=m.dy)//TODO: out-of-bounds error?
							a_idx=m.dy-1;
						if(a_idx)
							CMEMCPY(m.data, m.data+m.dx*a_idx, m.dx);
						m.dy=1;
					}
					else//select column
					{
						if(a_idx>=m.dx)
							a_idx=m.dx-1;
						m.data[0]=m.data[a_idx];
						m.dx=1;
					}
					CREALLOC(m.data, m.data, m.dx*m.dy);
				}
			}
			continue;
		case T_POWER:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m.dx!=m.dy)
					return user_error2(idx0, idx, "Expected a square matrix");
				if(m2.dx==1&&m2.dy==1)
				{
					int e=0;
					if(!get_int(m2, idx0, e))
						return false;
					//int e=(int)floor(m2.data[0]+0.5);
					auto CALLOC(temp, m.dx*m.dy);
					CREALLOC(m.data, m.data, m.dx*m.dy*2);
					impl_matpow(temp, m.data, e, m.dx);
					CFREE(m.data);
					m.data=temp;
				}
				else//TODO: matrix power matrix
					return my_error(idx0, idx);
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_unary(Matrix &m, bool space_sensitive)
{
	for(;;)
	{
		int idx0=idx;
		auto t=lex_get(space_sensitive);
		switch(t)
		{
			//region - functions
#if 1
		case T_ANS:case T_PRINTMODE:
		case T_ROOTS:
		case T_SAMPLE:case T_INVZ:
		case T_IDEN:
		case T_SUM:case T_TRACE:
		case T_REF:case T_RREF:case T_RREF2:
		case T_DET:case T_INV:
		case T_LU:
		case T_EGVAL:case T_NULLSPACE:case T_DIAG:
		case T_SQRT:case T_CBRT:case T_EXP:case T_LN:case T_LOG:
	//	case T_GAMMA:case T_LNGAMMA:
		case T_COS:case T_ACOS:case T_COSD:case T_ACOSD:case T_COSH:case T_ACOSH:
		case T_SEC:case T_ASEC:case T_SECD:case T_ASECD:case T_SECH:case T_ASECH:
		case T_SIN:case T_ASIN:case T_SIND:case T_ASIND:case T_SINH:case T_ASINH:
		case T_CSC:case T_ACSC:case T_CSCD:case T_ACSCD:case T_CSCH:case T_ACSCH:
		case T_TAN:case T_ATAN:case T_TAND:case T_ATAND:case T_TANH:case T_ATANH:
		case T_COT:case T_ACOT:case T_COTD:case T_ACOTD:case T_COTH:case T_ACOTH:
		case T_DFT:case T_IDFT:
	//	case T_FFT:case T_IFFT:
		case T_DCT:case T_IDCT:
			m.name=nullptr;
			{
				int idx1=idx;
				if(lex_get(false)==T_LPR)
				{
					if(!r_assign_expr(m, false))
						return false;
					if(lex_get(false)!=T_RPR)
						return user_error2(idx0, idx, "Expected a closing parenthesis \')\' of function call");
				}
				else//parentheses are optional for unary functions
				{
					idx=idx1;
					if(!r_unary(m, space_sensitive))
						return false;
				}
					//return user_error2(idx0, idx, "Expected an opening parenthesis \'(\' of function call");
			}
			switch(t)
			{
			case T_ANS:
				{
					int a_idx=0;
					get_int(m, idx0, a_idx);
					//if(m.dx>1||m.dy>1)//expected a scalar
					//	return user_error2(idx0, idx, "Expected a scalar");
					//int a_idx=(int)floor(m.data[0]+0.5);
					if(a_idx<0)
						a_idx=0;
					if(a_idx>=(int)g_answers.size())
						a_idx=g_answers.size()-1;
					m=g_answers[a_idx];
				}
				break;
			case T_PRINTMODE:
				{
					int i=0;
					get_int(m, idx0, i);
					print_fractions=i!=0;
					if(print_fractions)
						printf("Fractions: ON\n");
					else
						printf("Fractions: OFF\n");
				}
				break;
			case T_ROOTS:
				{
					if(minimum(m.dx, m.dy)!=1)
						return user_error2(idx0, idx, "Expected a vector of polynomial coefficients");
					int size=maximum(m.dx, m.dy);
					int k=0;

					//remove leading zero coefficients
					for(;k<size;++k)
					{
						auto &z=m.data[k];
						if(abs(z.r*z.r+z.i*z.i)>1e-10)
							break;
					}
					--k;
					if(k>0)
					{
						size-=k;
						m.dx=size, m.dy=1;
						CREALLOC(m.data, m.data, size-k);
					}
					if(size<2)
						return user_error2(idx0, idx, "Degenerate equation");
					//if(!check_real(m2, idx0))
					//	return false;
					int nroots=size-1;
					auto CALLOC(der, nroots), CALLOC(roots, nroots), CALLOC(deltas, nroots);
					for(k=1;k<size;++k)
					{
						der[k-1].r=(size-k)*m.data[k].r;
						der[k-1].i=(size-k)*m.data[k].i;
					}
					double factor=2*_pi/nroots;
					srand((unsigned)__rdtsc());
					for(k=0;k<nroots;++k)
					{
						roots[k].r=(double)rand()/RAND_MAX;
						roots[k].i=(double)rand()/RAND_MAX;
					}
					//for(k=0;k<nroots;++k)//initialize roots on unit circle?
					//{
					//	roots[k].r=cos(factor*k);
					//	roots[k].i=sin(factor*k);
					//}
					for(k=0;k<1e6;++k)
					{
						for(int k2=0;k2<nroots;++k2)
						{
							Comp d=eval_poly(der, size-1, roots[k2]);
							Comp p=eval_poly(m.data, size, roots[k2]);
							Comp sum={};
							for(int k3=0;k3<nroots;++k3)
							{
								if(k3==k2)
									continue;
								Comp term={roots[k2].r-roots[k3].r, roots[k2].i-roots[k3].i};
								c_inv(&term, &term);
								sum.r+=term.r;
								sum.i+=term.i;
							}
							c_mul(&sum, &sum, &p);
							d.r-=sum.r;
							d.i-=sum.i;
							c_div(&p, &p, &d);
							deltas[k2]=p;
						}
						for(int k2=0;k2<nroots;++k2)
						{
							roots[k2].r-=deltas[k2].r;
							roots[k2].i-=deltas[k2].i;
						}
						//printf("%d:\n", k);//
						//print_matrix_debug(deltas, nroots, 1);//
						double max_delta=0;
						for(int k2=0;k2<nroots;++k2)
						{
							auto &dk=deltas[k2];
							double d=abs(dk.r*dk.r+dk.i*dk.i);
							if(max_delta<d)
								max_delta=d;
						}
						if(max_delta<g_tolerance)
							break;
					}
					if(k==1e6)
						user_error2(idx0, idx, "Algorithm is unstable");
					std::sort(roots, roots+nroots, [](Comp const &a, Comp const &b)
					{
						if(a.r==b.r)
							return a.i<b.i;
						return a.r<b.r;
					});
					//printf("%d:\n", k);//
					//print_matrix_debug(deltas, nroots, 1);//
					CFREE(der);
					CFREE(m.data);
					m.data=roots;
					m.dx=1, m.dy=nroots;
				}
				break;
			case T_SAMPLE://TODO
				return my_error(idx0, idx);
			case T_INVZ://TODO
				return my_error(idx0, idx);
			case T_IDEN:
				{
					int size=0;
					get_int(m, idx0, size);
					//if(m.dx>1||m.dy>1)//expected a scalar
					//	return user_error2(idx0, idx, "Expected a scalar");
					//int size=(int)floor(m.data[0]+0.5);
					if(size<1)//wrong size
						return user_error2(idx0, idx, "Wrong size of identity matrix");
					m.dx=m.dy=size;
					size*=size;
					CALLOC(m.data, size);
					CMEMZERO(m.data, size);
					//memset(m.data, 0, size*sizeof(double));
					for(int k=0;k<size;k+=m.dx+1)
						m.data[k].r=1;
				}
				break;
			case T_SUM://should work like in Matlab
				{
					if(m.dy==1)
					{
						for(int kx=1;kx<m.dx;++kx)
						{
							m.data[0].r+=m.data[kx].r;//SIMD
							m.data[0].i+=m.data[kx].i;
						}
						m.dx=1;
					}
					else
					{
						for(int ky=1;ky<m.dy;++ky)
							for(int kx=0;kx<m.dx;++kx)
							{
								m.data[kx].r+=m.data[m.dx*ky+kx].r;
								m.data[kx].i+=m.data[m.dx*ky+kx].i;
							}
						m.dy=1;
					}
					CREALLOC(m.data, m.data, m.dy*m.dx);
					//m.data=(double*)realloc(m.data, m.dy*m.dx*sizeof(double));
				}
				break;
			case T_TRACE:
				if(m.dy!=m.dx)
					return user_error2(idx0, idx, "Expected a square matrix");
				for(int k=1;k<m.dx;++k)
				{
					m.data[0].r+=m.data[(m.dx+1)*k].r;
					m.data[0].i+=m.data[(m.dx+1)*k].i;
				}
				m.dx=m.dy=1;
				CREALLOC(m.data, m.data, m.dx*m.dy);
				break;
			case T_REF:
				impl_ref(m.data, m.dx, m.dy);
				break;
			case T_RREF:
				impl_rref(m.data, m.dx, m.dy);
				break;
			case T_RREF2:
				impl_rref2(m.data, m.dx, m.dy);
				break;
			case T_DET:
				if(m.dx!=m.dy)//must be square
					return user_error2(idx0, idx, "Expected a square matrix");
				CREALLOC(m.data, m.data, m.dx*m.dy*2);
				//m.data=(double*)realloc(m.dx*m.dy*sizeof(double));
				*m.data=impl_det(m.data, m.dx);
				m.dx=m.dy=1;
				CREALLOC(m.data, m.data, m.dx*m.dy);
				break;
			case T_INV:
				if(m.dy!=m.dx)
					return user_error2(idx0, idx, "Expected a square matrix");
				CREALLOC(m.data, m.data, m.dx*m.dy*2);
				impl_matinv(m.data, m.dx);
				CREALLOC(m.data, m.data, m.dx*m.dy);
				break;
			case T_LU:
				{
					if(m.dy!=m.dx)
						return user_error2(idx0, idx, "Expected a square matrix");
					int size=m.dx*m.dy;
					CREALLOC(m.data, m.data, size*3);
					CMEMCPY(m.data+size*2, m.data, size);
					impl_lu(m.data+size*2, m.dx, m.data+size, m.data);
					m.dy<<=1;
					CREALLOC(m.data, m.data, m.dx*m.dy);
				}
				break;
			case T_EGVAL:
				{
					if(m.dy!=m.dx)
						return user_error2(idx0, idx, "Expected a square matrix");
					if(m.dx<=1)
						break;
					int size=m.dx*m.dx;
					auto CALLOC(temp, m.dx);
					if(m.dx==2)
						impl_egval2(m.data, temp);
					else if(m.dx==3)
						impl_egval3(m.data, temp);
					else
					{
						auto CALLOC(D, size);
						impl_egval(m.data, m.dx, D, 1000);
						for(int k=0, k2=0;k<size;k+=m.dx+1, ++k2)
							temp[k2]=D[k];
						CFREE(D);
					}
					CFREE(m.data);
					m.data=temp;
					m.dy=1;
				}
				break;
			case T_NULLSPACE:
				{
					int size=m.dx*m.dx;
					auto CALLOC(temp, size);
					CMEMZERO(temp, size);
					auto dep_flags=new char[m.dx];
					auto row_idx=new short[m.dx];
					int nvec=impl_nullspace(m.data, m.dx, m.dy, temp, dep_flags, row_idx);
					delete[] dep_flags, row_idx;
					pack_rows(temp, m.dx, nvec, m.dx);
					//for(int k=1;k<m.dx;++k)//pack rows
					//{
					//	CMEMMOVE(temp+nvec*k, temp+m.dx*k, nvec);
					////	print_matrix_debug(temp, nvec, m.dx);//
					//}
					CREALLOC(temp, temp, nvec*m.dx);
					CFREE(m.data);
					m.data=temp;
					m.dy=m.dx;
					m.dx=nvec;
				}
				break;
			case T_DIAG:
				if(m.dx*m.dy>1)
				{
					int dim=m.dx;
					if(dim<m.dy)
						dim=m.dy;
					int size=dim*dim;
					auto CALLOC(temp, size);
					CMEMZERO(temp, size);
					for(int ks=0, kd=0;kd<size;++ks, kd+=size+1)
						temp[kd]=m.data[ks];
					CFREE(m.data);
					m.data=temp;
					m.dx=m.dy=dim;
				}
				break;
			//	return my_error(idx0, idx);
		/*	case T_DIAG2:
				{
					if(m.dx!=2||m.dy!=2)
						return user_error2(idx0, idx, "Expected a 2x2 matrix, got %dx%d", m.dy, m.dx);
					int size=4;
					auto CALLOC(temp, size*3);
					int success=impl_diag22(m.data, temp+(size<<1), temp, temp+size);
					CFREE(m.data);
					m.data=temp;
					m.dy*=3;
					if(!success)
						printf("Matrix is not diagonalizable.\n");
					//	return user_error2(idx0, idx, "Matrix is not diagonalizable");
				}
				break;
			case T_DIAG3:
				break;//*/
#define		EW_FUNC(FUNC)	for(int k=0, size=m.dx*m.dy;k<size;++k)FUNC(m.data[k]);
			case T_SQRT:	EW_FUNC(c_sqrt)break;
			case T_CBRT:	EW_FUNC(c_cbrt)break;
			case T_EXP:		EW_FUNC(c_exp)break;
			case T_LN:		EW_FUNC(c_ln)break;
			//case T_GAMMA:		EW_FUNC(c_tgamma)break;
			//case T_LNGAMMA:	EW_FUNC(c_lgamma)break;
			case T_COS:		EW_FUNC(c_cos)break;
			case T_ACOS:	EW_FUNC(c_acos)break;
			case T_COSD:	EW_FUNC(c_cosd)break;
			case T_ACOSD:	EW_FUNC(c_acosd)break;
			case T_COSH:	EW_FUNC(c_cosh)break;
			case T_ACOSH:	EW_FUNC(c_acosh)break;
			case T_SEC:		EW_FUNC(c_sec)break;
			case T_ASEC:	EW_FUNC(c_asec)break;
			case T_SECD:	EW_FUNC(c_secd)break;
			case T_ASECD:	EW_FUNC(c_asecd)break;
			case T_SECH:	EW_FUNC(c_sech)break;
			case T_ASECH:	EW_FUNC(c_asech)break;
			case T_SIN:		EW_FUNC(c_sin)break;
			case T_ASIN:	EW_FUNC(c_asin)break;
			case T_SIND:	EW_FUNC(c_sind)break;
			case T_ASIND:	EW_FUNC(c_asind)break;
			case T_SINH:	EW_FUNC(c_sinh)break;
			case T_ASINH:	EW_FUNC(c_asinh)break;
			case T_CSC:		EW_FUNC(c_csc)break;
			case T_ACSC:	EW_FUNC(c_acsc)break;
			case T_CSCD:	EW_FUNC(c_cscd)break;
			case T_ACSCD:	EW_FUNC(c_acscd)break;
			case T_CSCH:	EW_FUNC(c_csch)break;
			case T_ACSCH:	EW_FUNC(c_acsch)break;
			case T_TAN:		EW_FUNC(c_tan)break;
			case T_ATAN:	EW_FUNC(c_atan)break;
			case T_TAND:	EW_FUNC(c_tand)break;
			case T_ATAND:	EW_FUNC(c_atand)break;
			case T_TANH:	EW_FUNC(c_tanh)break;
			case T_ATANH:	EW_FUNC(c_atanh)break;
			case T_COT:		EW_FUNC(c_cot)break;
			case T_ACOT:	EW_FUNC(c_acot)break;
			case T_COTD:	EW_FUNC(c_cotd)break;
			case T_ACOTD:	EW_FUNC(c_acotd)break;
			case T_COTH:	EW_FUNC(c_coth)break;
			case T_ACOTH:	EW_FUNC(c_acoth)break;
#undef		EW_FUNC
			case T_DFT:
			case T_IDFT:
				{
					int size=m.dx*m.dy;
					Comp *weights=nullptr, *temp=nullptr;
					dft_init(m.dx, weights, temp, t==T_IDFT);
					for(int ky=0;ky<m.dy;++ky)
						dft_apply_1D(m.data+m.dx*ky, weights, temp, m.dx, 1);
					if(m.dx!=m.dy)
					{
						dft_finish(weights, temp);
						dft_init(m.dy, weights, temp, t==T_IDFT);
					}
					for(int kx=0;kx<m.dx;++kx)
						dft_apply_1D(m.data+kx, weights, temp, m.dy, m.dx);
					dft_finish(weights, temp);
				}
				break;
			//case T_FFT:
			//	break;
			//case T_IFFT:
			//	break;
			case T_DCT:
			case T_IDCT:
				{
					int size=m.dx*m.dy;
					Comp *weights=nullptr, *temp=nullptr;
					if(t==T_IDCT)
						idct_init(m.dx, weights, temp);
					else
						dct_init(m.dx, weights, temp);
					for(int ky=0;ky<m.dy;++ky)
						dft_apply_1D(m.data+m.dx*ky, weights, temp, m.dx, 1);
					if(m.dx!=m.dy)
					{
						dft_finish(weights, temp);
						if(t==T_IDCT)
							idct_init(m.dy, weights, temp);
						else
							dct_init(m.dy, weights, temp);
					}
					for(int kx=0;kx<m.dx;++kx)
						dft_apply_1D(m.data+kx, weights, temp, m.dy, m.dx);
					dft_finish(weights, temp);
				}
				break;
			}
			return r_postfix(m, space_sensitive);
		case T_CMD:
		case T_FRAC:
		case T_POLPOW:
		case T_LDIV:
		case T_CROSS:
		case T_EGVEC:
			{
				m.name=nullptr;
				Matrix m2;
				if(lex_get(false)!=T_LPR)
					return user_error2(idx0, idx, "Expected an opening parenthesis \'(\' of function call");
				if(!r_assign_expr(m, false))
					return false;
				if(lex_get(false)!=T_COMMA)
					return user_error2(idx0, idx, "Expected a comma \',\' of binary function call");
				if(!r_assign_expr(m2, false))
					return false;
				if(lex_get(false)!=T_RPR)
					return user_error2(idx0, idx, "Expected a closing parenthesis \'(\' of function call");
				switch(t)
				{
				case T_CMD:
					{
						int w=0, h=0;
						if(!get_int(m, idx0, w)||!get_int(m2, idx0, h))
							return false;
						//if(m.dx!=1||m.dy!=1||m2.dx!=1||m2.dy!=1)//expected 2 scalars
						//	return user_error2(idx0, idx, "Expected two scalar arguments");
						//int w=(int)floor(*m.data+0.5), h=(int)floor(*m2.data+0.5);
						m.data->r=set_console_buffer_size(w, h), m.data->i=0;
						m.dx=m.dy=1;
						CREALLOC(m.data, m.data, 1);
					}
					break;
				case T_FRAC:
					{
						if(!check_scalar(m, idx0)||!check_scalar(m2, idx0))
							return false;
						if(abs(m.data->i)>1e-10||abs(m2.data->i)>1e-10)
							return user_error2(idx0, idx, "Expected two real arguments");
						int i=0, num=0, den=0;
						dec2frac(m.data->r, m2.data->r, &i, &num, &den);
						CREALLOC(m.data, m.data, 3);
						m.data[0].r=i, m.data[0].i=0;
						m.data[1].r=num, m.data[1].i=0;
						m.data[2].r=den, m.data[2].i=0;
						m.dx=3;
					}
					break;
				case T_POLPOW:
					{
						if(m.dy!=1)
							return user_error2(idx0, idx, "Expected a row vector");
						int p=0;
						if(!get_int(m2, idx0, p))
							return false;
						if(p<0)
							return user_error2(idx0, idx, "Expected a positive integer");
						Matrix product;
						product.name=nullptr;
						product.dx=product.dy=1;
						CALLOC(product.data, 1);
						product.data->r=1, product.data->i=0;
						for(;;)
						{
							if(p&1)
								polmul(product, product, m, idx0);
							p>>=1;
							if(!p)
								break;
							polmul(m, m, m, idx0);
						}
						m=std::move(product);
					}
					break;
				case T_EGVEC:
					{
						if(m.dy!=m.dx)
							return user_error2(idx0, idx, "Expected a square matrix");
						int size=m.dx*m.dy;
						auto CALLOC(temp, size);
						Comp *lambdas=nullptr;
						bool copy_lambdas=m2.dx==1||m2.dy==1;
						if(copy_lambdas)//vector
						{
							int vecsize=m2.dx*m2.dy;
							if(vecsize!=m.dx)
								return user_error2(idx0, idx, "Expected %d eigenvalues, got %d", m.dx, vecsize);
							lambdas=m2.data;
						}
						else//diagonal matrix
						{
							if(m2.dx!=m2.dy)
								return user_error2(idx0, idx, "Expected a square matrix");
							if(m2.dx!=m.dx)
								return user_error2(idx0, idx, "Expected a diagonal square matrix of size %d", m.dx);
							CALLOC(lambdas, m.dx);
							for(int k=0;k<m.dx;++k)
								lambdas[k]=m2.data[(m.dx+1)*k];
						}
						int nvec=impl_egvec(m.data, m.dx, lambdas, temp);
						if(nvec<m.dx)
						{
							pack_rows(temp, m.dx, nvec, m.dx);
							CREALLOC(temp, temp, nvec*m.dx);
						}
						CFREE(m.data);
						m.data=temp;
						m.dx=nvec;
					}
					break;
				case T_LDIV://TODO
					return my_error(idx0, idx);
				case T_CROSS:
					{
						int size1=m.dx*m.dy, size2=m2.dx*m2.dy;
						if(size1==2&&size2==2)
						{
							Comp t1, t2;
							c_mul(&t1, m.data, m2.data+1);
							c_mul(&t2, m.data+1, m2.data);
							t1.r-=t2.r, t1.i-=t2.i;
							m.data[0]=t1;
							m.dx=m.dy=1;
							//double temp=m.data[0]*m2.data[1]-m.data[1]*m2.data[0];
							CREALLOC(m.data, m.data, 1);
							//*m.data=temp;
						}
						else if(size1==3&&size2==3)
						{
							auto CALLOC(temp, 3);
							Comp t4;
							c_mul(temp  , m.data+1, m2.data+2), c_mul(&t4, m.data+2, m2.data+1), temp[0].r-=t4.r, temp[0].i-=t4.i;
							c_mul(temp+1, m.data+2, m2.data+0), c_mul(&t4, m.data+0, m2.data+2), temp[1].r-=t4.r, temp[1].i-=t4.i;
							c_mul(temp+2, m.data+0, m2.data+1), c_mul(&t4, m.data+1, m2.data+0), temp[2].r-=t4.r, temp[2].i-=t4.i;
							//temp[0]=m.data[1]*m2.data[2]-m.data[2]*m2.data[1];
							//temp[1]=m.data[2]*m2.data[0]-m.data[0]*m2.data[2];
							//temp[2]=m.data[0]*m2.data[1]-m.data[1]*m2.data[0];
							CFREE(m.data);
							m.data=temp;
						}
						else
							return user_error2(idx0, idx, "cross() expects two 2D or two 3D vectors");
					}
					break;
				}
			}
			return r_postfix(m, space_sensitive);
		case T_CONV://variadic function
			{
				m.name=nullptr;
				if(lex_get(false)!=T_LPR)
					return user_error2(idx0, idx, "Expected an opening parenthesis \'(\' of function call");
				if(!r_assign_expr(m, false))
					return false;
				if(m.dy!=1)
					return user_error2(idx0, idx, "Expected a row vector");
				for(t=lex_get(false);t==T_COMMA;)
				{
					Matrix m2;
					if(!r_assign_expr(m2, false))
						return false;

					polmul(m, m, m2, idx0);
					//if(m.dy!=1||m2.dy!=1)
					//	return user_error2(idx0, idx, "Expected a row vector");
					//auto CALLOC(temp, m.dx+m2.dx-1);
					//impl_polmul(temp, m.data, m2.data, m.dx, m2.dx, 0);
					//CFREE(m.data);
					//m.data=temp;
					//m.dx+=m2.dx-1;

					t=lex_get(false);
				}
				if(t!=T_RPR)
					return user_error2(idx0, idx, "Expected an closing parenthesis \')\' of function call");
			}
			break;
#endif
			//region - constants
#if 1
		case T_NUMBER:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 1);
			m.data->r=lex_number, m.data->i=0;
			return r_postfix(m, space_sensitive);
		case T_RAND:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 1);
			srand((unsigned)__rdtsc());
			m.data->r=(double)rand()/RAND_MAX, m.data->i=0;
			return r_postfix(m, space_sensitive);
		case T_IMAG:
		case T_IMAG_UNUSED:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 2);
			m.data->r=0, m.data->i=1;
			return r_postfix(m, space_sensitive);
		case T_EULER:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 1);
			m.data->r=_e, m.data->i=0;
			return r_postfix(m, space_sensitive);
		case T_PI:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 1);
			m.data->r=_pi, m.data->i=0;
			return r_postfix(m, space_sensitive);
		case T_INF:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 1);
			m.data->r=_HUGE, m.data->i=0;
			return r_postfix(m, space_sensitive);
		case T_NAN:
			m.name=nullptr;
			m.dx=m.dy=1;
			CALLOC(m.data, 1);
			m.data->r=_HUGE-_HUGE, m.data->i=0;
			return r_postfix(m, space_sensitive);
#endif
		case T_POLY:
			return my_error(idx0, idx);
		case T_MINUS:
			m.name=nullptr;
			if(!r_unary(m, space_sensitive))
				return false;
			if(!r_postfix(m, space_sensitive))
				return false;
			for(int k=0, size=m.dx*m.dy;k<size;++k)
			{
				m.data[k].r=-m.data[k].r;
				m.data[k].i=-m.data[k].i;
			}
			//for(int k=0, size=m.dx*m.dy*(1+(m.type==T_COMPLEX));k<size;++k)
			//	m.data[k]=-m.data[k];
			return true;
		case T_PLUS:
			m.name=nullptr;
			if(!r_unary(m, space_sensitive))
				return false;
			return r_postfix(m, space_sensitive);
		case T_LPR:
			m.name=nullptr;
			if(!r_assign_expr(m, false))
				return false;
			if(lex_get(false)!=T_RPR)
				return false;
			return r_postfix(m, space_sensitive);
		case T_LBRACKET://matrix
			{
				m.name=nullptr;
				if(!r_row_vector(m))
					return false;
				for(;;)
				{
					idx0=idx;
					switch(auto t=lex_get(true))
					{
					case T_SEMICOLON:
					case T_NEWLINE:
						{
							Matrix m2;
							int esize=errors.size();
							if(!r_row_vector(m2))
							{
								if(t==T_NEWLINE)//trailing whitespace
								{
									if(esize<(int)errors.size())
										errors.erase(errors.begin()+esize, errors.begin()+errors.size());//pop-back errors inserted by last r_row_vector()
									break;
								}
								return user_error2(idx0, idx, "Expected a matrix row");
							}

							//join matrices vertically
							if(m.dx!=m2.dx)
								return user_error2(idx0, idx, "Matices have different widths: w1=%d != w2=%d", m.dx, m2.dx);
							int newdy=m.dy+m2.dy;
							CREALLOC(m.data, m.data, m.dx*newdy);
							//m.data=(double*)realloc(m.data, m.dx*newdy*sizeof(double));
							CMEMCPY(m.data+m.dx*m.dy, m2.data, m2.dx*m2.dy);
							//memcpy(m.data+m.dx*m.dy, m2.data, m2.dx*m2.dy*sizeof(double));
							m.dy=newdy;
						}
						continue;
					}
					idx=idx0;
					break;
				}
				if(lex_get(false)!=T_RBRACKET)
					return user_error2(idx0, idx, "Expected a closing bracket \']\'");
			}
			return r_postfix(m, space_sensitive);
	/*	case T_POLSTART:
			{
				lex_skip_space();
				if(!r_equality(m, true))
					return false;
				for(;;)
				{
					int idx0=idx;
					switch(auto t=lex_get(true))
					{
					case T_SPACE:
					case T_COMMA:
						break;
					}
				}
			}
			return r_postfix(m, space_sensitive);//*/
		case T_ID://TODO: look up variable name
			{
				auto it=g_vars.find(lex_id);
				if(it==g_vars.end())
					return user_error2(idx0, idx, "Undefined");
				m=it->second;
			}
			return r_postfix(m, space_sensitive);
		default:
			user_error2(idx0, idx, "Expected an expression");
			idx=idx0;
			return false;
		}
		break;
	}
	return true;
}
bool		r_multiplicative(Matrix &m, bool space_sensitive)
{
	if(!r_unary(m, space_sensitive))
		return false;
	for(;;)
	{
		int idx0=idx;
		switch(lex_get(space_sensitive))
		{
		case T_MUL:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(!obj_mul(m, m2, idx0))
					return false;
			}
			continue;
		case T_TENSOR:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				auto CALLOC(temp, m.dx*m.dy*m2.dx*m2.dy);
				//double *temp=(double*)malloc(m.dx*m.dy*m2.dx*m2.dy*sizeof(double));
				impl_tensor(temp, m.data, m2.data, m.dx, m.dy, m2.dx, m2.dy);
				CFREE(m.data);
				m.data=temp;
				m.dx*=m2.dx;
				m.dy*=m2.dy;
			}
			continue;
		case T_DIV:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(!obj_div(m, m2, idx0))
					return false;
			}
			continue;
		case T_DIV_BACK:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m.flags==M_SCALAR)//multiplication by scalar
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						c_mul(m2.data+k, m2.data+k, m.data);
						//m2.data[k]*=m.data[0];
					m.move2temp(m2);
				}
				else//matrix division
				{
					if(m.dx!=m.dy)//denominator matrix must be square
						return user_error2(idx0, idx, "Denominator matrix must be square: h1=%d != w1=%d", m.dy, m.dx);
					if(m.dx!=m2.dy)//dimension mismatch
						return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
					CREALLOC(m.data, m.data, m2.dx*m2.dy);
					//m.data=(double*)realloc(m.data, m.dx*m.dy*2*sizeof(double));
					auto CALLOC(temp, m2.dx*m2.dy);
					//double *temp=(double*)malloc(m2.dx*m2.dy*sizeof(double));
					impl_matdiv_back(temp, m.data, m2.data, m.dy, m2.dx);
					CFREE(m.data);
					m.data=temp;
				}
			}
			continue;
		case T_MOD:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(!obj_mod(m, m2, idx0))
					return false;
			}
			continue;
		case T_MUL_EW:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m2.flags==M_SCALAR)//multiplication by scalar
				{
					for(int k=0, size=m.dy*m.dx;k<size;++k)
						c_mul(m.data+k, m.data+k, m2.data);
						//m.data[k]*=m2.data[0];
				}
				else if(m.flags==M_SCALAR)//multiplication by scalar
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						c_mul(m2.data+k, m2.data+k, m.data);
						//m2.data[k]*=m.data[0];
					m.move2temp(m2);
				}
				else
				{
					if(m2.dx!=m.dx||m2.dy!=m.dy)
						return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						c_mul(m.data+k, m.data+k, m2.data+k);
						//m.data[k]*=m2.data[k];
				}
			}
			continue;
		case T_DIV_EW:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m2.flags==M_SCALAR)//division by scalar
				{
					for(int k=0, size=m.dy*m.dx;k<size;++k)
						c_div(m.data+k, m.data+k, m2.data);
						//m.data[k]/=m2.data[0];
				}
				else if(m.flags==M_SCALAR)//element-wise division of scalar by matrix
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						c_div(m2.data+k, m.data, m2.data+k);
						//m2.data[k]=m.data[0]/m2.data[k];
					m.move2temp(m2);
				}
				else
				{
					if(m2.dx!=m.dx||m2.dy!=m.dy)
						return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						c_div(m.data+k, m.data+k, m2.data+k);
						//m.data[k]/=m2.data[k];
				}
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_additive(Matrix &m, bool space_sensitive)
{
	if(!r_multiplicative(m, space_sensitive))
		return false;
	for(;;)
	{
		int idx0=idx;
		switch(lex_get(space_sensitive))
		{
		case T_PLUS:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(!obj_add(m, m2, idx0))
					return false;
			}
			continue;
		case T_MINUS:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(!obj_sub(m, m2, idx0))
					return false;
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_range(Matrix &m, bool space_sensitive)
{
	int idx0=idx;
	if(lex_get(space_sensitive)==T_COLON)
	{
		m.flags=M_UNSPECIFIED_RANGE;
		return true;
	}
	idx=idx0;
	if(!r_additive(m, space_sensitive))
		return false;
	for(;;)
	{
		idx0=idx;
		switch(lex_get(space_sensitive))
		{
		case T_COLON://range
			{
				auto &mstart=m;
				Matrix mstep, mend;
				if(!check_scalar(mstart, idx0))
					return false;
				idx0=idx;
				if(!r_additive(mend, space_sensitive))
					return false;
				if(!check_scalar(mend, idx0))
					return false;
				idx0=idx;
				if(lex_get(space_sensitive)==T_COLON)
				{
					mstep.move2temp(mend);
					if(!r_additive(mend, space_sensitive))
						return false;
					if(!check_scalar(mend, idx0))
						return false;
				}
				else
				{
					idx=idx0;
					mstep.dx=1, mstep.dy=1;
					CALLOC(mstep.data, 1);
					mstep.data->r=1, mstep.data->i=0;
				}
				int guard=1000, size=0;
				auto &start=mstart.data[0], &step=mstep.data[0], &end=mend.data[0];
				Comp x=start;
				double dummy0=0, dummy1=1;
				double *rleft=&dummy0, *rright=&dummy1, *ileft=&dummy0, *iright=&dummy1;//in-loop <= comparisons
				//sign check
				if(start.r<end.r)
				{
					if(step.r<0)
						return user_error2(idx0, idx, "Wrong step sign: step.r < 0");
					rleft=&x.r, rright=&end.r;
				}
				else if(start.r>end.r)
				{
					if(step.r>0)
						return user_error2(idx0, idx, "Wrong step sign: step.r > 0");
					rleft=&end.r, rright=&x.r;
				}
				//else//start.r==end.r
				//{
				//}
				if(start.i<end.i)
				{
					if(step.i<0)
						return user_error2(idx0, idx, "Wrong step sign: step.i < 0");
					ileft=&x.i, iright=&end.i;
				}
				else if(start.i>end.i)
				{
					if(step.i>0)
						return user_error2(idx0, idx, "Wrong step sign: step.i > 0");
					ileft=&end.i, iright=&x.i;
				}
				//else//start.i==end.i
				//{
				//}
				if(start.r==end.r&&start.i==end.i)
				{
					if(!step.r&&!step.i)
						step.r=1;
					rleft=&x.r, rright=&end.r;
					ileft=&x.i, iright=&end.i;
				}
				if(!step.r&&!step.i)//TODO: better check for infinite loop
					return user_error2(idx0, idx, "Range step = 0");
				if(step.r)
					size=(int)floor((end.r-start.r)/step.r)+1;
				if(step.i)
				{
					int size2=(int)floor((end.i-start.i)/step.i)+1;
					if(size>size2)
						size=size2;
				}
				if(size>guard)
				{
					//display warning [y/n]
					printf("\nWarning: allocating a range of %d values (max is %d). Proceed? [Y/N] ", size, guard);
					char c=0;
					scanf("%c", &c);
					if((c&0xDF)!='Y')
						return false;
				}
				CREALLOC(m.data, m.data, size);
				m.dx=size, m.dy=1;
				for(int k=0;k<size&&*rleft<=*rright&&*ileft<=*iright;x.r+=step.r, x.i+=step.i, ++k)
					m.data[k]=x;
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_relational(Matrix &m, bool space_sensitive)
{
	if(!r_range(m, space_sensitive))
		return false;
	for(;;)
	{
		int idx0=idx;
		switch(lex_get(space_sensitive))
		{
		case T_LESS:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_range(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k].r=m.data[k].r<m2.data[k].r, m.data[k].i=0;
			}
			continue;
		case T_LESS_EQUAL:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_range(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k].r=m.data[k].r<=m2.data[k].r, m.data[k].i=0;
			}
			continue;
		case T_GREATER:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_range(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k].r=m.data[k].r>m2.data[k].r, m.data[k].i=0;
			}
			continue;
		case T_GREATER_EQUAL:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_range(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k].r=m.data[k].r>=m2.data[k].r, m.data[k].i=0;
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_equality(Matrix &m, bool space_sensitive)
{
	if(!r_relational(m, space_sensitive))
		return false;
	for(;;)
	{
		int idx0=idx;
		switch(auto t=lex_get(space_sensitive))
		{
		case T_EQUAL:
		case T_NOT_EQUAL:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_relational(m2, space_sensitive))
					return false;
				bool equal=m.dx==m2.dx&&m.dy==m2.dy;
				if(equal)
				{
					for(int k=0, size=m.dx*m.dy;k<size;++k)
					{
						if(c_distance(m.data[k], m2.data[k])>1e-10)
						{
							equal=false;
							break;
						}
					}
				}
				CREALLOC(m.data, m.data, 1);
				//m.data=(double*)realloc(m.data, sizeof(double));
				m.data[0].r=t==T_EQUAL==equal, m.data[0].i=0;
				m.dx=1;
				m.dy=1;
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_assign_expr(Matrix &m, bool space_sensitive)
{
	int idx0=idx;
	if(lex_get(space_sensitive)==T_ID)
	{
		char *name=lex_id;
		switch(auto t=lex_get(space_sensitive))
		{
		case T_ASSIGN:
			if(!r_assign_expr(m, space_sensitive))
				return false;
			if(m.flags==M_UNSPECIFIED_RANGE)
				return user_error2(idx0, idx, "Cannot assign an unspecified range");
			m.name=name;
			g_vars[name]=m;
			return true;
		case T_ASSIGN_ADD:
		case T_ASSIGN_SUB:
		case T_ASSIGN_MUL:
		case T_ASSIGN_DIV:
		case T_ASSIGN_MOD:
			{
				if(!r_assign_expr(m, space_sensitive))
					return false;
				if(m.flags==M_UNSPECIFIED_RANGE)
					return user_error2(idx0, idx, "Cannot assign an unspecified range");
				auto &mv=g_vars[name];
				switch(t)
				{
				case T_ASSIGN_ADD:
					if(!obj_add(mv, m, idx0))
						return false;
					break;
				case T_ASSIGN_SUB:
					if(!obj_sub(mv, m, idx0))
						return false;
					break;
				case T_ASSIGN_MUL:
					if(!obj_mul(mv, m, idx0))
						return false;
					break;
				case T_ASSIGN_DIV:
					if(!obj_div(mv, m, idx0))
						return false;
					break;
				case T_ASSIGN_MOD:
					if(!obj_mod(mv, m, idx0))
						return false;
					break;
				}
				m=mv;
				//m.name=name;
			}
			return true;
		}
		goto next;
	}
next:
	idx=idx0;
	return r_equality(m, space_sensitive);
}
int			solve(std::string &str, bool again)
{
	int nspaces=16;
	str.insert(str.size(), nspaces, ' ');
	//const char spaces[]="                ";
	//const int nspaces=sizeof(spaces)-1;//should be at least 16
	//str+=spaces;

	lex_init(str.c_str(), str.size()-nspaces);
	parse_incomplete=false;
	if(!g_answers.size()||!again&&g_answers.back().data)
		g_answers.push_back(Matrix());
	auto &m=g_answers.back();
	SolveResultType ret=SOLVE_OK;
	//try
	//{
	int idx0=idx;
	switch(lex_get(false))
	{
	case T_HELP:
		print_help();
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_OPEN:
		if(get_str_from_file(str))
		{
			str.insert(str.size(), nspaces, ' ');
			//str+=spaces;
			lex_init(str.c_str(), str.size()-nspaces);
			goto parse_file;
		}
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_CLEAR:
		g_vars.clear();
		g_answers.clear();
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_VARS:
		{
			bool printed=false;
			if(g_vars.size())
			{
				printf("Variables:\n");
				for(auto &var:g_vars)
				{
					printf("%s =\n", var.first);
					var.second.print();
				}
				printed=true;
			}
			if(g_answers.size()-1>0)
			{
				printf("Previous answers:\n");
				for(int k=0;k<(int)g_answers.size()-1;++k)
				{
					auto &ans=g_answers[k];
					printf("ans(%d) =\n", k);
					ans.print();
				}
				printed=true;
			}
			if(!printed)
				printf("Nothing to show\n");
		}
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_FRACTIONS:
		if(print_fractions=!print_fractions)
			printf("Fractions: ON\n");
		else
			printf("Fractions: OFF\n");
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_TOLERANCE:
		if(lex_get(false)!=T_NUMBER)
			printf("Expected a tolerance value\n");
		else
		{
			g_tolerance=lex_number;
			printf("algorithm tolerance = %g\n", g_tolerance);
		}
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_GFSET://TODO
		ret=SOLVE_PARSE_ERROR;
		break;
	case T_EXIT:
	case T_QUIT:
		exit(0);
		break;
	default:
		idx=idx0;
	parse_file:
		while(idx<text_size)
		{
			if(!r_assign_expr(m, false))
			{
				ret=SOLVE_PARSE_ERROR;
				break;
			}
			idx0=idx;
			switch(lex_get(false))
			{
			case T_SEMICOLON:
				continue;
			case T_EOF:
				break;
			default:
				user_error2(idx0, idx, "Unexpected text");
				ret=SOLVE_PARSE_ERROR;
				break;
			}
			break;
			//if(idx<text_size)
			//{
			//	auto t=lex_get(false);
			//	if(t!=T_SEMICOLON&&t!=T_EOF)
			//	{
			//		ret=SOLVE_PARSE_ERROR;
			//		break;
			//	}
			//}
		}
		if(m.flags==M_UNSPECIFIED_RANGE)
		{
			user_error2(idx0, idx, "Result is an unspecified range");
			ret=SOLVE_PARSE_ERROR;
		}
		break;
	}
	//}
	//catch(...)
	//{
	//	return SOLVE_INCOMPLETE;
	//}
	str.erase(str.size()-nspaces, nspaces);
	return ret;
}