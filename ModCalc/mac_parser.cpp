//mac_parser.cpp - ModCalc recursive parser
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

#define _CRT_SECURE_NO_WARNINGS
#include	<stdarg.h>
#include	<math.h>
#include	<map>
#include	<algorithm>
#include	<tmmintrin.h>
#include	"mac.h"
#include	"sha2.h"
#include	"sha3.h"
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
const char	file[]=__FILE__;

int			base=10;
bool		benchmark=false;
//double	g_tolerance=1e-15;
//bool		print_fractions=false;
std::string Int::to_string()const//don't forget to free buffer
{
	std::string str;
	if(!data||!data->x)
	{
		str+='0';
		return str;
	}
	switch(base)
	{
	case 2:
		{
			int nbits=floor_log2(*this)+1;
			str.reserve(nbits);
			if(is_neg())
				str+='-';
			str+="0b";
			for(int k=nbits-1;k>=0;--k)
				str+='0'+getbit(k);
		}
		break;
	case 16://print whole buffer
		{
			int ndigits=(int)(data->used_size<<3);
			str.reserve(ndigits);
			if(is_neg())
				str+='-';
			str+="0x";
			for(int k=ndigits-1;k>=0;--k)
			{
				int digit=getnibble(k);
				if(digit<10)
					str+='0'+digit;
				else
					str+='A'+digit-10;
			}
		}
		break;
	case 10:
		{
			Int val=*this;
			val.attempt_shrink();
			bool neg=val.data->sign!=0;
			val.abs();

			Int Q, R;
			//Int t2=val/ten;
			for(;;)
			{
				impl_ldiv(val, ten, Q, R);
				int digit=R.to_int();
				str.push_back('0'+digit);
				if(Q==zero)
					break;
				if(ten*Q+R!=val)
					impl_ldiv(val, ten, Q, R);
				val=Q;

				//int digit=(val-t2*ten).to_int();
				//str.push_back('0'+digit);
				//val=t2, t2/=ten;
			}
			//if(!str.size())
			//	str+='0';
			if(neg)
				str+='-';
			std::reverse(str.begin(), str.end());
		}
		break;
	}
	return str;

	//return lint_itoa(data);
}
void		Matrix::print()const
{
	if(!v.size())
	{
		printf("nullptr\n");
		//printf("Error: matrix.data == nullptr\n");
		return;
	}
	//auto print_func=print_fractions?print_double_frac:print_double;
	if(flags==M_SCALAR)//no brackets
	{
		printf(" ");
		v[0].print();
		printf("\n");
	}
	//	printf(" %lld\n", *data);
	else//matrix
	{
		int size=dy*dx;
		auto colsizes=new short[size+dx];
		memset(colsizes, 0, (size+dx)*sizeof(*colsizes));
		for(int ky=0;ky<dy;++ky)
		{
			for(int kx=0;kx<dx;++kx)
			{
				auto &x=get(kx, ky);
				auto &nd=colsizes[dx*ky+kx];
				nd=1;
				if(x!=zero)
					nd+=floor_log10(x)+(x<zero);
				if(colsizes[size+kx]<nd)
					colsizes[size+kx]=nd;
				//auto str=x.to_string();
				//int nd=strlen(str);
				//free(str);
			}
		}
		printf("[\n");
		for(int ky=0;ky<dy;++ky)
		{
			for(int kx=0;kx<dx;++kx)
			{
				int valsize=colsizes[dx*ky+kx], colsize=colsizes[size+kx];
				printf(" %*s", colsize-valsize, "");
				int printed=v[dx*ky+kx].print();
			}
			//	printf(" %*lld", colsizes[kx], get(kx, ky));
			if(ky<dy-1)
				printf("\n");
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
	printed+=sprintf_s(g_buf+printed, g_buf_size-printed, ".");
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
void		pack_rows(Int *buffer, int dx_old, int dx_new, int dy)//dx_new < dx_old
{
	for(int k=1;k<dy;++k)
	{
		objmove(buffer+dx_new*k, buffer+dx_old*k, dx_new);
	//	print_matrix_debug(buffer, dx_new, dy);//
	}
}

bool		ntt_init(int idx0, int size, Int const &weight, Int const &n, IntVec &weights, IntVec &phi, IntVec &temp, bool inverse, bool anti_cyclic)
{
	if(size<2)
		return user_error2(idx0, idx, "NTT size = %d < 2", size);
	weights.resize(size);
	temp.resize(size);
	Int root=weight;
	if(anti_cyclic)//weight should be sqrtw
	{
		phi.resize(size);
		phi[0]=1;
		if(inverse)//invphi = [1  -sqrtw^(size-1) ... -sqrtw^2  -sqrtw]
		{
			Int negsqrtw=n-root;
			for(int k=size-1;k>=1;--k)
			{
				phi[k]=phi[(k+1)%size];
				phi[k].mulmod(negsqrtw, n);
			}
			Int w_2size=phi[1]*negsqrtw%n;
			if(w_2size!=1)//check
			{
				auto s_root=root.to_string(), s_w_2size=w_2size.to_string();
				user_error2(idx0, idx, "Invalid anti-cyclic InvNTT parameters:\n\t(modulus-sqrtw)^size = %s^%d = %s != 1", s_root.c_str(), size, s_w_2size.c_str());
				return false;
			}
		}
		else//phi = [1  sqrtw  sqrtw^2 ... sqrtw^(size-1)]
		{
			for(int k=1;k<size;++k)
			{
				phi[k]=phi[k-1];
				phi[k].mulmod(root, n);
			}
			Int w_size=phi[size-1]*root%n, nm1=n-1;
			if(w_size!=nm1)//check
			{
				auto s_root=root.to_string(), s_w_size=w_size.to_string(), s_nm1=nm1.to_string();
				user_error2(idx0, idx, "Invalid anti-cyclic NTT parameters:\n\tsqrtw^size = %s^%d = %s\n\t!= modulus-1 = %s", s_root.c_str(), size, s_w_size.c_str(), s_nm1.c_str());
				return false;
			}
		}
		root.mulmod(weight, n);//w=sqrtw^2
	}
	weights[0]=1;
	if(inverse)//invweights = [1  w^(size-1) ... w^2  w] .* size^-1 mod n
	{
		for(int k=size-1;k>=1;--k)
		{
			weights[k]=weights[(k+1)%size];
			weights[k].mulmod(root, n);
		}
		if(!anti_cyclic)//check
		{
			Int wsize=weights[1]*root%n;
			wsize.mulmod(size, n);
			if(wsize!=1)
			{
				auto s_root=root.to_string(), s_wsize=wsize.to_string();
				user_error2(idx0, idx, "Invalid InvNTT parameters:\n\tw^size = %s^%d = %s != 1", s_root.c_str(), size, s_wsize.c_str());
				return false;
			}
		}
		Int invsize=impl_intEEA(Int(size), n);
		if(invsize==zero)
		{
			auto s_n=n.to_string();
			user_error2(idx0, idx, "Invalid InvNTT parameters:\n\tsize=%d has no inverse mod %s", size, s_n.c_str());
			return false;
		}
		for(int k=0;k<size;++k)
			weights[k].mulmod(invsize, n);
	}
	else//weights = [1  w  w^2 ... w^(size-1)]
	{
		for(int k=1;k<size;++k)
		{
			weights[k]=weights[k-1];
			weights[k].mulmod(root, n);
		}
		if(!anti_cyclic)//check
		{
			Int wsize=weights[size-1]*root%n;
			if(wsize!=1)
			{
				auto s_root=root.to_string(), s_wsize=wsize.to_string();
				user_error2(idx0, idx, "Invalid NTT parameters:\n\tw^size = %s^%d = %s != 1", s_root.c_str(), size, s_wsize.c_str());
				return false;
			}
		}
	}
/*	if(inverse)
	{
		Int wk=root;
		if(anti_cyclic)
		{
			for(int k=size-1;k>=0;--k)//inv_phi = {1, -phi^(n-1), ..., -phi}
			{
				phi[k]=wk;
				wk.mulmod(root, n);
			}
			root.mulmod(weight, n);
			wk=root;
		}
		Int invsize=impl_intEEA(Int(size), n);
		wk.mulmod(invsize, n);
		for(int k=size-1;k>=0;--k)
		{
			weights[k]=wk;
			wk.mulmod(root, n);
		}
		if(wk!=root)
		{
			auto str=wk.to_string();
			auto str2=weight.to_string();
			user_error2(idx0, idx, "Invalid root of unity: w^(%d+1)=%s != w=%s", size, str, str2);
			free(str);
			free(str2);
			return false;
		}
	}
	else
	{
		Int wk=1;
		if(anti_cyclic)
		{
			root.mulmod(weight, n);
		}
		for(int k=0;k<size;++k)
		{
			weights[k]=wk;
			wk.mulmod(root, n);
		}
		if(wk!=1)
		{
			auto str=wk.to_string();
			user_error2(idx0, idx, "Invalid root of unity: w^%d=%s != 1", size, str);
			free(str);
			return false;
		}
	}//*/
	return true;
}
void		ntt_apply_1D(Int *data, IntVec const &weights, IntVec const &phi, IntVec &temp, Int const &n, int size, int stride, bool anti_cyclic)
{
	for(int k=0;k<size;++k)
	{
		temp[k]=zero;
		for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
			temp[k]+=weights[kd*k%size]*data[ks];
		temp[k]%=n;
	}
	if(anti_cyclic)
	{
		for(int ks=0;ks<size;++ks)
			temp[ks].mulmod(phi[ks], n);
	}
	for(int ks=0, kd=0;ks<size;++ks, kd+=stride)
		data[kd]=temp[ks];
}

/*Int			eval_poly(Int const *coeff, int order, Int &x)
{
	Int result=coeff[0];
	for(int k=1;k<order;++k)
	{
		c_mul(&result, &result, &x);
		result.r+=coeff[k].r;
		result.i+=coeff[k].i;
	}
	return result;
}//*/
bool		polmul(Matrix &dst, Matrix const &a, Matrix const &b, int idx0)
{
	//if(a.dy!=1||b.dy!=1)
	//	return user_error2(idx0, idx, "Expected a row vector");
	//Int *ALLOC(temp, a.dx+b.dx-1);
	//impl_polmul(temp, a.data, b.data, a.dx, b.dx, 0);
	//if(dst.data)
	//	FREE(dst.data);
	//dst.data=temp;
	//dst.dx=a.dx+b.dx-1;
	return true;
}

char		find_quotes(int &start)
{
	start=idx;
	for(;start<text_size&&text[start]!='\"'&&text[start]!='\'';++start);
	if(start>=text_size)
		return 0;
	return text[start];
}
inline bool	check_scalar(Matrix &m, int idx0)
{
	if(m.flags!=M_SCALAR)
		return user_error2(idx0, idx, "Expected a scalar");
	return true;
}
inline bool	check_square(Matrix &m, int idx0)
{
	if(m.dy!=m.dx)
		return user_error2(idx0, idx, "Expected a square matrix");
	return true;
}
inline void cmatrix_plus_cscalar(Matrix &m, Matrix &s)
{
	int size=m.size();
	for(int k=0;k<size;++k)
		m.v[k]+=s.v[0];
	//__m128d vb=_mm_load_pd((double*)s.data.data());
	//for(auto p=m.data, end=m.end();p<end;++p)
	//{
	//	__m128d va=_mm_load_pd((double*)p);
	//	va=_mm_add_pd(va, vb);
	//	_mm_store_pd((double*)p, va);
	//}
}
inline bool	obj_add(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dy==1&&m2.dx==1)//addition of a scalar
	{
		cmatrix_plus_cscalar(m, m2);
	}
	else if(m.dy==1&&m.dx==1)//addition of a scalar
	{
		cmatrix_plus_cscalar(m2, m);
		m.move2temp(m2);
	}
	else//matrix addition
	{
		if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		impl_addbuffers(m.data(), m.data(), m2.data(), m.size());
	}
	return true;
}
inline bool	obj_sub(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dy==1&&m2.dx==1)//subtraction of a scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
		{
			m.v[k]-=m2.v[0];
			m.v[k]-=m2.v[0];
		}
	}
	else if(m.dy==1&&m.dx==1)//scalar minus matrix
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
		{
			m2.v[k]=m.v[0]-m2.v[k];
			m2.v[k]=m.v[0]-m2.v[k];
		}
		m.move2temp(m2);
	}
	else
	{
		if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		impl_subbuffers(m.data(), m.data(), m2.data(), m.size());
	}
	return true;
}
inline bool	obj_mul(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.flags==M_SCALAR)//multiplication by scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
			m.v[k]*=m2.v[0];
	}
	else if(m.flags==M_SCALAR)//multiplication by scalar
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
			m2.v[k]*=m.v[0];
		m.move2temp(m2);
	}
	else//matrix multiplication
	{
		if(m.dx!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
		IntVec temp(m.dy*m2.dx);
		//Int *ALLOC(temp, m.dy*m2.dx);
		impl_matmul(temp.data(), m.data(), m2.data(), m.dy, m.dx, m2.dx);
		m.move_in(temp, m2.dx, -1);
		//m.v=std::move(temp);
		//FREE(m.data);
		//m.data=temp;
		//m.dx=m2.dx;
	}
	return true;
}
#if 0
inline bool	obj_div(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dx==1&&m2.dy==1)//division by scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
			m.data[k]/=m2.data[0];
			//c_div(m.data+k, m.data+k, m2.data);
	}
	else//matrix division
	{
		if(m2.dx!=m2.dy)//denominator matrix must be square
			return user_error2(idx0, idx, "Denominator matrix must be square: h2=%d != w2=%d", m2.dy, m2.dx);
		if(m.dx!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
		REALLOC(m2.data, m2.data, m2.dx*m2.dy);
		Int *ALLOC(temp, m.dx*m.dy);
		impl_matdiv(temp, m.data, m2.data, m.dy, m2.dx);
		FREE(m.data);
		m.data=temp;
	}
	return true;
}
#endif
inline bool	obj_mod(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dx==1&&m2.dy==1)//m2 can be scalar
	{
		for(int k=0, size=m.dx*m.dy;k<size;++k)
			m.v[k]%=m2.v[0];
			//c_mod(m.data+k, m.data+k, m2.data);
			//m.data[k]=m.data[k]-floor(m.data[k]/m2.data[0])*m2.data[0];
	}
	else
	{
		if(m2.dx!=m.dx||m2.dy!=m.dy)
			return user_error2(idx0, idx, "The modulus operator \'%%\' is element-wise: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		for(int k=0, size=m.dx*m.dy;k<size;++k)
			m.v[k]%=m2.v[k];
			//c_mod(m.data+k, m.data+k, m2.data+k);
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
				int esize=(int)errors.size();
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
				int dx0=m.dx;
				m.dx+=m2.dx;
				m.resize();
				for(int ky=m.dy-1;ky>=0;--ky)
				{
					objmove(m.data()+m.dx*ky, m.data()+dx0*ky, dx0);
					objmove(m.data()+m.dx*ky+dx0, m2.data()+m2.dx*ky, m2.dx);
				}
			//	int newdx=m.dx+m2.dx;
			//	m.v.resize(newdx*m.dy);
			//	//REALLOC(m.data, m.data, newdx*m.dy);
			//	for(int ky=m.dy-1;ky>=0;--ky)
			//	{
			//		MEMMOVE(m.data()+newdx*ky, m.data()+m.dx*ky, m.dx);
			//		MEMMOVE(m.data()+newdx*ky+m.dx, m2.data()+m2.dx*ky, m2.dx);
			//	}
			//	m.dx=newdx;
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
				auto temp=std::move(m.v);
				m.v.resize(m.dx*m.dy);
			//	Int *temp=m.data;
			//	ALLOC(m.data, m.dx*m.dy);
				std::swap(m.dx, m.dy);
				for(int ky=0;ky<m.dy;++ky)
					for(int kx=0;kx<m.dx;++kx)
						m.get(kx, ky)=GET(temp, m.dy, ky, kx);
			//	FREE(temp);
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
				IntVec temp(mky.size());
				//Int *ALLOC(temp, mky.dy*mkx.dx);
				for(int ky=0;ky<mky.dy;++ky)
				{
					auto val=&mky.get(0, ky);
					int ky2=clamp(0, val->to_int()-1, m.dy-1);
					for(int kx=0;kx<mkx.dx;++kx)
					{
						val=&mkx.get(kx, 0);
						int kx2=clamp(0, val->to_int()-1, m.dx-1);
						temp[mkx.dx*ky+kx]=m.v[m.dx*ky2+kx2];
					}
				}
				m.move_in(temp, mkx.dx, mky.dy);
				//m.dx=mkx.dx, m.dy=mky.dy;
				//FREE(m.data);
				//m.data=temp;
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
					int a_idx=midx.v[0].to_int();
					if(a_idx<0)
						a_idx=0;
					if(m.dy>1)//select row
					{
						if(a_idx>=m.dy)//TODO: out-of-bounds error?
							a_idx=m.dy-1;
						if(a_idx)
							objmove(m.data(), m.data()+m.dx*a_idx, m.dx);
							//MEMCPY(m.data(), m.data()+m.dx*a_idx, m.dx);
						m.dy=1;
					}
					else//select column
					{
						if(a_idx>=m.dx)
							a_idx=m.dx-1;
						m.v[0]=m.v[a_idx];
						m.dx=1;
					}
					m.resize();
				//	REALLOC(m.data, m.data, m.dx*m.dy);
				}
			}
			continue;
		case T_POWER:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;

				if(!check_scalar(m, idx0)||!check_scalar(m2, idx0))
					return false;
				m.v[0]^=m2.v[0];
				//if(m.dx!=m.dy)
				//	return user_error2(idx0, idx, "Expected a square matrix");
				//if(m2.flags==M_SCALAR)
				//{
				//	int e=m2.data->to_int();
				//	Int *ALLOC(temp, m.dx*m.dy);
				//	REALLOC(m.data, m.data, m.dx*m.dy*2);
				//	impl_matpow(temp, m.data, e, m.dx);
				//	FREE(m.data);
				//	m.data=temp;
				//}
				//else//TODO: matrix power matrix
				//	return my_error(idx0, idx);
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
		case T_ANS:
		case T_ROOTS:
		case T_SAMPLE:case T_INVZ:
		case T_IDEN:
		case T_SUM:case T_MEAN:case T_MIN:case T_MAX:case T_PRODUCT:
		case T_REF:case T_RREF:case T_RREF2:
		case T_DET:case T_TRACE:case T_INV:
		case T_LU:
		case T_EGVAL:case T_NULLSPACE:case T_DIAG:
		case T_GCD:case T_LCM:
		case T_RANDMOD:case T_RANDLOG:
		case T_ISPRIME:case T_FACTORIZE:case T_TOTIENT:case T_CARMICHAEL:case T_PROOTS:
		case T_LOG2:case T_LOG10:
		case T_SQRT:case T_CBRT:
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
					if(m.dx>1||m.dy>1)//expected a scalar
						return user_error2(idx0, idx, "Expected a scalar");
					int a_idx=m.v[0].to_int();
					if(a_idx<0)
						a_idx=0;
					if(a_idx>=(int)g_answers.size())
						a_idx=(int)(g_answers.size()-1);
					m=g_answers[a_idx];
				}
				break;
			//case T_PRINTMODE:
			//	{
			//		if(m.dx>1||m.dy>1)//expected a scalar
			//			return user_error2(idx0, idx, "Expected a scalar");
			//		int i=m.data->to_int();
			//		print_fractions=i!=0;
			//		if(print_fractions)
			//			printf("Fractions: ON\n");
			//		else
			//			printf("Fractions: OFF\n");
			//	}
			//	break;
			//case T_ROOTS:
			//	break;
			case T_SAMPLE://TODO
				return my_error(idx0, idx);
			case T_INVZ://TODO
				return my_error(idx0, idx);
			case T_IDEN:
				{
					if(m.dx>1||m.dy>1)//expected a scalar
						return user_error2(idx0, idx, "Expected a scalar");
					int size=m.v[0].to_int();
					if(size<1)//wrong size
						return user_error2(idx0, idx, "Wrong size of identity matrix");
					m.dx=m.dy=size;
					size*=size;
					m.resize();
					//ALLOC(m.data, size);
					//MEMZERO(m.data(), size);
					//memset(m.data, 0, size*sizeof(double));
					for(int k=0;k<size;k+=m.dx+1)
						m.v[k]=1;
				}
				break;
			case T_SUM://should work like in Matlab
				if(m.dy==1)
				{
					for(int kx=1;kx<m.dx;++kx)
						m.v[0]+=m.v[kx];
					m.dx=1;
				}
				else
				{
					for(int ky=1;ky<m.dy;++ky)
					{
						for(int kx=0;kx<m.dx;++kx)
							m.v[kx]+=m.v[m.dx*ky+kx];
					}
					m.dy=1;
				}
				m.resize();
			//	REALLOC(m.data, m.data, m.dy*m.dx);
				break;
			case T_MEAN:
				if(m.dy==1)
				{
					for(int kx=1;kx<m.dx;++kx)
						m.v[0]+=m.v[kx];
					m.v[0]/=m.dx;
					m.dx=1;
				}
				else
				{
					for(int ky=1;ky<m.dy;++ky)
					{
						for(int kx=0;kx<m.dx;++kx)
							m.v[kx]+=m.v[m.dx*ky+kx];
					}
					for(int kx=0;kx<m.dx;++kx)
						m.v[kx]/=m.dy;
					m.dy=1;
				}
				m.resize();
			//	REALLOC(m.data, m.data, m.dy*m.dx);
				break;
			case T_MIN:
				{
					if(m.dy==1)
					{
						for(int kx=1;kx<m.dx;++kx)
							if(m.v[0]>m.v[kx])
								m.v[0]=m.v[kx];
						m.dx=1;
					}
					else
					{
						for(int ky=1;ky<m.dy;++ky)
						{
							for(int kx=0;kx<m.dx;++kx)
								if(m.v[kx]>m.get(kx, ky))
									m.v[kx]=m.get(kx, ky);
						}
						m.dy=1;
					}
					m.resize();
				//	REALLOC(m.data, m.data, m.dy*m.dx);
				}
				break;
			case T_MAX:
				{
					if(m.dy==1)
					{
						for(int kx=1;kx<m.dx;++kx)
							if(m.v[0]<m.v[kx])
								m.v[0]=m.v[kx];
						m.dx=1;
					}
					else
					{
						for(int ky=1;ky<m.dy;++ky)
						{
							for(int kx=0;kx<m.dx;++kx)
								if(m.v[kx]<m.get(kx, ky))
									m.v[kx]=m.get(kx, ky);
						}
						m.dy=1;
					}
					m.resize();
				//	REALLOC(m.data, m.data, m.dy*m.dx);
				}
				break;
			case T_PRODUCT:
				if(m.dy==1)
				{
					for(int kx=1;kx<m.dx;++kx)
						m.v[0]*=m.v[kx];
					m.dx=1;
				}
				else
				{
					for(int ky=1;ky<m.dy;++ky)
					{
						for(int kx=0;kx<m.dx;++kx)
							m.v[kx]*=m.v[m.dx*ky+kx];
					}
					m.dy=1;
				}
				m.resize();
				break;
			case T_TRACE:
				if(!check_square(m, idx0))
					return false;
				for(int k=1;k<m.dx;++k)
					m.v[0]+=m.v[(m.dx+1)*k];
				m.dx=m.dy=1;
				m.resize();
			//	REALLOC(m.data, m.data, m.dx*m.dy);
				break;
			//case T_REF:
			//	impl_ref(m.data, m.dx, m.dy);
			//	break;
			//case T_RREF:
			//	impl_rref(m.data, m.dx, m.dy);
			//	break;
			//case T_RREF2:
			//	impl_rref2(m.data, m.dx, m.dy);
			//	break;
			//case T_DET:
			//	if(m.dx!=m.dy)//must be square
			//		return user_error2(idx0, idx, "Expected a square matrix");
			//	REALLOC(m.data, m.data, m.dx*m.dy*2);
			//	*m.data=impl_det(m.data, m.dx);
			//	m.dx=m.dy=1;
			//	REALLOC(m.data, m.data, m.dx*m.dy);
			//	break;
			//case T_INV:
			//	if(m.dy!=m.dx)
			//		return user_error2(idx0, idx, "Expected a square matrix");
			//	REALLOC(m.data, m.data, m.dx*m.dy*2);
			//	impl_matinv(m.data, m.dx);
			//	REALLOC(m.data, m.data, m.dx*m.dy);
			//	break;
			//case T_LU:
			//	{
			//		if(m.dy!=m.dx)
			//			return user_error2(idx0, idx, "Expected a square matrix");
			//		int size=m.dx*m.dy;
			//		REALLOC(m.data, m.data, size*3);
			//		MEMCPY(m.data+size*2, m.data, size);
			//		impl_lu(m.data+size*2, m.dx, m.data+size, m.data);
			//		m.dy<<=1;
			//		REALLOC(m.data, m.data, m.dx*m.dy);
			//	}
			//	break;
			//case T_EGVAL:
			//	{
			//		if(m.dy!=m.dx)
			//			return user_error2(idx0, idx, "Expected a square matrix");
			//		if(m.dx<=1)
			//			break;
			//		int size=m.dx*m.dx;
			//		Int *ALLOC(temp, m.dx);
			//		if(m.dx==2)
			//			impl_egval2(m.data, temp);
			//		else if(m.dx==3)
			//			impl_egval3(m.data, temp);
			//		else
			//		{
			//			Int *ALLOC(D, size);
			//			impl_egval(m.data, m.dx, D, 1000);
			//			for(int k=0, k2=0;k<size;k+=m.dx+1, ++k2)
			//				temp[k2]=D[k];
			//			FREE(D);
			//		}
			//		FREE(m.data);
			//		m.data=temp;
			//		m.dy=1;
			//	}
			//	break;
			//case T_NULLSPACE:
			//	{
			//		int size=m.dx*m.dx;
			//		Int *ALLOC(temp, size);
			//		MEMZERO(temp, size);
			//		auto dep_flags=new char[m.dx];
			//		auto row_idx=new short[m.dx];
			//		int nvec=impl_nullspace(m.data, m.dx, m.dy, temp, dep_flags, row_idx);
			//		delete[] dep_flags, row_idx;
			//		pack_rows(temp, m.dx, nvec, m.dx);
			//		REALLOC(temp, temp, nvec*m.dx);
			//		FREE(m.data);
			//		m.data=temp;
			//		m.dy=m.dx;
			//		m.dx=nvec;
			//	}
			//	break;
			case T_DIAG:
				if(m.dx*m.dy>1)
				{
					int dim=m.dx;
					if(dim<m.dy)
						dim=m.dy;
					int size=dim*dim;
					IntVec temp(size);
				//	Int *ALLOC(temp, size);
				//	MEMZERO(temp.data(), size);
					for(int ks=0, kd=0;kd<size;++ks, kd+=size+1)
						temp[kd]=m.v[ks];
					m.move_in(temp, dim, dim);
				//	FREE(m.data);
				//	m.data=temp;
				//	m.dx=m.dy=dim;
				}
				break;
			case T_GCD:
				for(int k=1, size=m.size();k<size;++k)
					m.v[0]=impl_gcd(m.v[0], m.v[k]);
				m.dx=m.dy=1;
				m.resize();
				break;
			case T_LCM:
				m.v[0]=impl_lcm(m.data(), m.size());
				m.dx=m.dy=1;
				m.resize();
				break;
			case T_ISPRIME:
				for(int k=0, size=m.size();k<size;++k)
				{
					auto &val=m.v[k];
					val.data->x[0]=impl_isprime(val);
					val.data->used_size=1;
					val.attempt_shrink();
				}
				break;
			case T_FACTORIZE:
				{
					if(!check_scalar(m, idx0))
						return false;
					IntVec temp;
					impl_factorize(m.v[0], temp);
					m.move_in(temp, (short)temp.size(), 1);
				}
				break;
			case T_TOTIENT:
				{
					IntVec factors, powers;
					for(int k=0, size=m.size();k<size;++k)
					{
						auto &val=m.v[k];
						impl_factorize(val, factors, &powers);
						val=impl_totient(val, factors, powers);
					}
				}
				break;
			case T_CARMICHAEL:
				{
					IntVec factors, powers;
					for(int k=0, size=m.size();k<size;++k)
					{
						auto &val=m.v[k];
						impl_factorize(val, factors, &powers);
						val=impl_carmichael(val, factors, powers);
					}
				}
				break;
			case T_PROOTS:
				{
					if(!check_scalar(m, idx0))
						return false;
					IntVec temp;
					impl_proots(m.v[0], temp);
					m.move_in(temp, (short)temp.size(), 1);
				}
				break;
#define		EW_FUNC(FUNC)	for(int k=0, size=m.size();k<size;++k)m.v[k]=FUNC(m.v[k]);
			case T_RANDMOD:	EW_FUNC(impl_randmod)break;
			case T_RANDLOG:	EW_FUNC(impl_randlog)break;
			case T_GENPRIME:EW_FUNC(impl_genprime)break;
			case T_LOG2:	EW_FUNC(floor_log2)break;
			case T_LOG10:	EW_FUNC(floor_log10)break;
			case T_SQRT:	EW_FUNC(isqrt)break;
			case T_CBRT:	EW_FUNC(icbrt)break;
#undef		EW_FUNC
			default:
				return my_error(idx0, idx);
			}
			return r_postfix(m, space_sensitive);
		case T_MRTEST:
		case T_INVMOD:
		case T_CMD:
	//	case T_FRAC:
		case T_POLPOW:
		case T_LDIV:
		case T_CROSS:
		case T_EGVEC:
		case T_NTT:case T_INTT:
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
				case T_MRTEST:
					{
						int step=0;
						if(m2.flags!=M_SCALAR)
						{
							if(m.dx!=m2.dx||m.dy!=m2.dy)
								return user_error2(idx0, idx, "Dimension mismatch: m1:%dx%d != m2:%dx%d. Second matrix should be scalar of match the size of first one", m.dx, m.dy, m2.dx, m2.dy);
							step=1;
						}
						for(int k=0, size=m.size(), k2=0;k<size;++k, k2+=step)
						{
							auto &val=m.v[k];
							val.data->x[0]=impl_mrtest(val, m2.v[k2]);
							val.data->used_size=1;
							val.attempt_shrink();
						}
					}
					break;
				case T_INVMOD:
					{
						if(!check_scalar(m2, idx0))
							return false;
						for(int k=0, size=m.size();k<size;++k)
							m.v[k]=impl_intEEA(m.v[k], m2.v[0]);
					}
					break;
#ifndef __linux__
				case T_CMD:
					{
						if(!check_scalar(m, idx0)||!check_scalar(m2, idx0))
							return false;
						int w=m.v[0].to_int(), h=m2.v[0].to_int();
						m.v[0]=set_console_buffer_size(w, h);
						m.dx=m.dy=1;
						m.resize();
					//	REALLOC(m.data, m.data, 1);
					}
					break;
#endif
				case T_POLPOW:
					{
						if(m.dy!=1)
							return user_error2(idx0, idx, "Expected a row vector");
						if(!check_scalar(m, idx0))
							return false;
						int p=m.v[0].to_int();
						if(p<0)
							return user_error2(idx0, idx, "Expected a positive integer");
						Matrix product;
						product.name=nullptr;
						product.dx=product.dy=1;
						product.resize();
					//	ALLOC(product.data, 1);
						product.v[0]=1;
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
			/*	case T_EGVEC:
					{
						if(m.dy!=m.dx)
							return user_error2(idx0, idx, "Expected a square matrix");
						int size=m.dx*m.dy;
						Int *ALLOC(temp, size);
						Int *lambdas=nullptr;
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
							ALLOC(lambdas, m.dx);
							for(int k=0;k<m.dx;++k)
								lambdas[k]=m2.data[(m.dx+1)*k];
						}
						int nvec=impl_egvec(m.data, m.dx, lambdas, temp);
						if(nvec<m.dx)
						{
							pack_rows(temp, m.dx, nvec, m.dx);
							REALLOC(temp, temp, nvec*m.dx);
						}
						FREE(m.data);
						m.data=temp;
						m.dx=nvec;
					}
					break;//*/
				case T_LDIV://TODO
					return my_error(idx0, idx);
				case T_CROSS:
					{
						int size1=m.dx*m.dy, size2=m2.dx*m2.dy;
						if(size1==2&&size2==2)
						{
							m.v[0]=m.v[0]*m2.v[1]-m.v[1]*m2.v[0];
							m.dx=m.dy=1;
							m.resize();
						}
						else if(size1==3&&size2==3)
						{
							IntVec temp(3);
							temp[0]=m.v[1]*m2.v[2]-m.v[2]*m2.v[1];
							temp[1]=m.v[2]*m2.v[0]-m.v[0]*m2.v[2];
							temp[2]=m.v[0]*m2.v[1]-m.v[1]*m2.v[0];
							m.move_in(temp, -1, -1);
						}
						else
							return user_error2(idx0, idx, "cross() expects two 2D or two 3D vectors");
					}
					break;
				case T_NTT://ntt(M, [N wx wy anti_cyclic])
				case T_INTT:
					{
						m.name=nullptr;
						if(m2.dx!=4||m2.dy!=1)
							return user_error2(idx0, idx, "Second argument should be [modulus wx wy anti_cyclic]");
						auto &modulus=m2.v[0], &wx=m2.v[1], &wy=m2.v[2];
						bool anti_cyclic=m2.v[3]!=zero;
						if(modulus==zero)
							return user_error2(idx0, idx, "modulus is zero");
						wx%=modulus;
						wy%=modulus;
						int size=m.dx*m.dy;
						IntVec weights, temp, phi;
						bool initialized=false;
						if(m.dx>1)
						{
							initialized=true;
							if(!ntt_init(idx0, m.dx, wx, modulus, weights, phi, temp, t==T_INTT, anti_cyclic))
								return false;
							for(int ky=0;ky<m.dy;++ky)
								ntt_apply_1D(m.data()+m.dx*ky, weights, phi, temp, modulus, m.dx, 1, anti_cyclic);
						}
						if(m.dy>1)
						{
							if(!initialized||m.dx!=m.dy||wx!=wy)
								if(!ntt_init(idx0, m.dy, wy, modulus, weights, phi, temp, t==T_INTT, anti_cyclic))
									return false;
							for(int kx=0;kx<m.dx;++kx)
								ntt_apply_1D(m.data()+kx, weights, phi, temp, modulus, m.dy, m.dx, anti_cyclic);
						}
					}
					break;
				}
			}
			return r_postfix(m, space_sensitive);
		case T_DLOG://dlog(x, base, modulus)
			{
				Matrix base, modulus;
				if(lex_get(false)!=T_LPR)
					return user_error2(idx0, idx, "Expected an opening parenthesis \'(\' of function call");
				if(!r_assign_expr(m, false))
					return false;
				if(lex_get(false)!=T_COMMA)
					return user_error2(idx0, idx, "Expected a comma \',\'");
				if(!r_assign_expr(base, false))
					return false;
				if(lex_get(false)!=T_COMMA)
					return user_error2(idx0, idx, "Expected a comma \',\'");
				if(!r_assign_expr(modulus, false))
					return false;
				if(lex_get(false)!=T_RPR)
					return user_error2(idx0, idx, "Expected a closing parenthesis \'(\' of function call");
				m.v[0]=impl_dlog(m.v[0], base.v[0], modulus.v[0]);
			}
			break;
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

					t=lex_get(false);
				}
				if(t!=T_RPR)
					return user_error2(idx0, idx, "Expected an closing parenthesis \')\' of function call");
			}
			break;
		case T_HIST:
			{
				int start=0;
				char quotes=find_quotes(start);
				if(!quotes)
					return user_error2(idx0, idx, "Expected quoted text");
				char quote=text[start];
				++start;
				int end=start;
				for(;end<text_size&&text[end]!=quote;++end);
				if(text[end]!=quote)
				{
					if(quote=='\'')
						return user_error2(idx0, idx, "Expected a closing quote");
					return user_error2(idx0, idx, "Expected closing double quotes");
				}

				auto histogram=new int[256];
				memset(histogram, 0, 256*sizeof(int));
				for(int k=start;k<end;++k)
					++histogram[text[k]];
				int histmax=0;
				for(int k=0;k<256;++k)
					if(histmax<histogram[k])
						histmax=histogram[k];
				if(histmax)
				{
					printf("Symbol  Freq\n");
					for(int k=0;k<256;++k)
					{
						if(histogram[k])
						{
							if(k>32&&k<127)
								printf("%3d \'%c\': %3d ", k, (char)k, histogram[k]);
							else
								printf("%3d     : %3d ", k, histogram[k]);
							for(int k2=0, nstars=histogram[k]*64/histmax;k2<nstars;++k2)
								printf("*");
							printf("\n");
						}
					}
				}
				delete[] histogram;
				idx=end+1;
			}
			return r_postfix(m, space_sensitive);
		case T_SHA224:
		case T_SHA256:
		case T_SHA384:
		case T_SHA512:
		case T_SHA3_224:
		case T_SHA3_256:
		case T_SHA3_384:
		case T_SHA3_512:
			{
				int start=0;
				char quotes=find_quotes(start);
				//int start=idx;
				//for(;start<text_size&&text[start]!='\"'&&text[start]!='\'';++start);
				//if(text[start]!='\"'&&text[start]!='\'')
				if(quotes)
				{
					char quote=text[start];
					++start;
					int end=start;
					for(;end<text_size&&text[end]!=quote;++end);
					if(text[end]!=quote)
					{
						if(quote=='\'')
							return user_error2(idx0, idx, "Expected a closing quote");
						return user_error2(idx0, idx, "Expected closing double quotes");
					}
					unsigned hash[16]={};
					m.dx=m.dy=1;
					m.resize();
					auto &val=m.v[0];
					switch(t)
					{
					case T_SHA224:
						sha256((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(7);
						memcpy(val.data->x, hash, 7*sizeof(int));
						break;
					case T_SHA256:
						sha256((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(8);
						memcpy(val.data->x, hash, 8*sizeof(int));
						break;
					case T_SHA384:
						sha384((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(12);
						memcpy(val.data->x, hash, 12*sizeof(int));
						break;
					case T_SHA512:
						sha512((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(16);
						memcpy(val.data->x, hash, 16*sizeof(int));
						break;
					case T_SHA3_224:
						FIPS202_SHA3_224((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(7);
						memcpy(val.data->x, hash, 7*sizeof(int));
						break;
					case T_SHA3_256:
						FIPS202_SHA3_256((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(8);
						memcpy(val.data->x, hash, 8*sizeof(int));
						break;
					case T_SHA3_384:
						FIPS202_SHA3_384((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(12);
						memcpy(val.data->x, hash, 12*sizeof(int));
						break;
					case T_SHA3_512:
						FIPS202_SHA3_512((unsigned char*)(text+start), end-start, (unsigned char*)hash);
						val.setzero(16);
						memcpy(val.data->x, hash, 16*sizeof(int));
						break;
					}
					std::reverse((unsigned char*)val.data->x, (unsigned char*)(val.data->x+val.data->used_size));
					idx=end+1;
				}
				else//hash elements of a matrix		TODO: single hash for whole matrix
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
					unsigned hash[16]={};
					for(int k=0;k<m.size();++k)
					{
						auto &val=m.v[k];
						if(val.data)
						{
							std::reverse((unsigned char*)val.data->x, (unsigned char*)(val.data->x+val.data->used_size));
							switch(t)
							{
							case T_SHA224:
								sha224((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(7);
								memcpy(val.data->x, hash, 7*sizeof(int));
								break;
							case T_SHA256:
								sha256((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(8);
								memcpy(val.data->x, hash, 8*sizeof(int));
								break;
							case T_SHA384:
								sha384((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(12);
								memcpy(val.data->x, hash, 12*sizeof(int));
								break;
							case T_SHA512:
								sha512((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(16);
								memcpy(val.data->x, hash, 16*sizeof(int));
								break;
							case T_SHA3_224:
								FIPS202_SHA3_224((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(7);
								memcpy(val.data->x, hash, 7*sizeof(int));
								break;
							case T_SHA3_256:
								FIPS202_SHA3_256((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(8);
								memcpy(val.data->x, hash, 8*sizeof(int));
								break;
							case T_SHA3_384:
								FIPS202_SHA3_384((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(12);
								memcpy(val.data->x, hash, 12*sizeof(int));
								break;
							case T_SHA3_512:
								FIPS202_SHA3_512((unsigned char*)val.data->x, (unsigned)val.data->used_size, (unsigned char*)hash);
								val.setzero(16);
								memcpy(val.data->x, hash, 16*sizeof(int));
								break;
							}
							std::reverse((unsigned char*)val.data->x, (unsigned char*)(val.data->x+val.data->used_size));
						}
					}
				}
			}
			return r_postfix(m, space_sensitive);
		case T_SHAKE128://f(m, m2/quoted_text)
		case T_SHAKE256:
		case T_AES128:
		case T_AES192:
		case T_AES256:
			{
				if(lex_get(false)!=T_LPR)
					return user_error2(idx0, idx, "Expected an opening parenthesis \'(\'");
				if(!r_unary(m, space_sensitive))
					return false;
				if(lex_get(false)!=T_COMMA)
					return user_error2(idx0, idx, "Expected a comma \',\'");
				int start=0;
				char quotes=find_quotes(start);
				if(quotes)
				{
				}
				else
				{
				}
				if(lex_get(false)!=T_LPR)
					return user_error2(idx0, idx, "Expected a closing parenthesis \')\'");
			}
			break;
		//case T_POLY:
		//	return my_error(idx0, idx);
#endif
			//region - constants
#if 1
		case T_NUMBER:
			m.name=nullptr;
			m.dx=m.dy=1;
			m.resize();
		//	ALLOC(m.data, 1);
			m.v[0]=std::move(lex_number);
			return r_postfix(m, space_sensitive);
		case T_RDTSC:
			m.name=nullptr;
			m.dx=m.dy=1;
			m.resize();
			m.v[0].set(__rdtsc());
			return r_postfix(m, space_sensitive);
		case T_RAND:
			m.name=nullptr;
			m.dx=m.dy=1;
			m.resize();
		//	ALLOC(m.data, 1);
			srand((unsigned)__rdtsc());
			m.v[0]=rand()<<15|rand();
			return r_postfix(m, space_sensitive);
#endif
		case T_MINUS:
			m.name=nullptr;
			if(!r_unary(m, space_sensitive))
				return false;
			if(!r_postfix(m, space_sensitive))
				return false;
			for(int k=0, size=m.dx*m.dy;k<size;++k)
				m.v[k].negate();
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
							int esize=(int)errors.size();
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
								return user_error2(idx0, idx, "Matrices have different widths: w1=%d != w2=%d", m.dx, m2.dx);
							int dy0=m.dy;
							m.dy+=m2.dy;
							m.resize();
							objmove(m.data()+m.dx*dy0, m2.data(), m2.dx*m2.dy);
							//MEMCPY(m.data()+m.dx*dy0, m2.data(), m2.dx*m2.dy);
							
							//int newdy=m.dy+m2.dy;
							//REALLOC(m.data, m.data, m.dx*newdy);
							//MEMCPY(m.data+m.dx*m.dy, m2.data, m2.dx*m2.dy);
							//m.dy=newdy;
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
				IntVec temp(m.size()*m2.size());
			//	Int *ALLOC(temp, m.dx*m.dy*m2.dx*m2.dy);
				impl_tensor(temp.data(), m.data(), m2.data(), m.dx, m.dy, m2.dx, m2.dy);
				m.move_in(temp, m.dx*m2.dx, m.dy*m2.dy);
				//FREE(m.data);
				//m.data=temp;
				//m.dx*=m2.dx;
				//m.dy*=m2.dy;
			}
			continue;
		case T_DIV:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				return user_error2(idx0, idx, "Division isn't supported yet");//TODO
				//if(!obj_div(m, m2, idx0))
				//	return false;
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
						m2.v[k]*=m.v[0];
					m.move2temp(m2);
				}
				else//matrix division
				{
					if(m.dx!=m.dy)//denominator matrix must be square
						return user_error2(idx0, idx, "Denominator matrix must be square: h1=%d != w1=%d", m.dy, m.dx);
					if(m.dx!=m2.dy)//dimension mismatch
						return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
					return user_error2(idx0, idx, "Division isn't supported yet");//TODO
					//REALLOC(m.data, m.data, m2.dx*m2.dy);
					//Int *ALLOC(temp, m2.dx*m2.dy);
					//impl_matdiv_back(temp, m.data, m2.data, m.dy, m2.dx);
					//FREE(m.data);
					//m.data=temp;
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
						m.v[k]*=m2.v[0];
				}
				else if(m.flags==M_SCALAR)//multiplication by scalar
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						m2.v[k]*=m.v[0];
					m.move2temp(m2);
				}
				else
				{
					if(m2.dx!=m.dx||m2.dy!=m.dy)
						return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.v[k]*=m2.v[k];
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
						m.v[k]/=m2.v[0];
				}
				else if(m.flags==M_SCALAR)//element-wise division of scalar by matrix
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						m2.v[k]=m.v[0]/m2.v[k];
					m.move2temp(m2);
				}
				else
				{
					if(m2.dx!=m.dx||m2.dy!=m.dy)
						return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.v[k]/=m2.v[k];
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
					mstep.resize();
				//	ALLOC(mstep.data, 1);
					mstep.v[0]=1;
				//	mstep.data->construct(1);
				}
				int guard=1000, size=0;
				int start=mstart.v[0].to_int(), step=mstep.v[0].to_int(), end=mend.v[0].to_int();
				int x=start;
				int rleft=0, rright=1;//in-loop <= comparisons
				//sign check
				if(start<end)
				{
					if(step<zero)
						return user_error2(idx0, idx, "Wrong step sign: step.r < 0");
					rleft=x, rright=end;
				}
				else if(start>end)
				{
					if(step>zero)
						return user_error2(idx0, idx, "Wrong step sign: step.r > 0");
					rleft=end, rright=x;
				}
				if(step==zero)//TODO: better check for infinite loop
					return user_error2(idx0, idx, "Range step = 0");
				if(step!=zero)
					size=(end-start)/step+1;
				if(size>guard)
				{
					//display warning [y/n]
					printf("\nWarning: allocating a range of %d values (max is %d). Proceed? [Y/N] ", size, guard);
					char c=0;
					scanf("%c", &c);
					if((c&0xDF)!='Y')
						return false;
				}
			//	REALLOC(m.data, m.data, size);
				m.dx=size, m.dy=1;
				m.resize();
				for(int k=0;k<size&&rleft<=rright;x+=step, ++k)
					m.v[k]=x;
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
					m.v[k]=m.v[k]<m2.v[k];
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
					m.v[k]=m.v[k]<=m2.v[k];
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
					m.v[k]=m.v[k]>m2.v[k];
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
					m.v[k]=m.v[k]>=m2.v[k];
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
						if(m.v[k]!=m2.v[k])
						{
							equal=false;
							break;
						}
					}
				}
				m.dx=m.dy=1;
				m.resize();
			//	REALLOC(m.data, m.data, 1);
				m.v[0]=t==T_EQUAL==equal;
			}
			continue;
		case T_EQUAL_EW:
		case T_NOT_EQUAL_EW:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_relational(m2, space_sensitive))
					return false;
				if(t==T_EQUAL)
				{
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.v[k]=m.v[k]==m2.v[k];
				}
				else
				{
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.v[k]=m.v[k]!=m2.v[k];
				}
			}
			continue;
		}
		idx=idx0;
		break;
	}
	return true;
}
bool		r_mod_lp(Matrix &m, bool space_sensitive)
{
	if(!r_equality(m, space_sensitive))
		return false;
	for(;;)
	{
		int idx0=idx;
		switch(auto t=lex_get(space_sensitive))
		{
		case T_MOD_LP:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_equality(m2, space_sensitive))
					return false;
				if(!obj_mod(m, m2, idx0))
					return false;
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
					return user_error2(idx0, idx, "Division isn't supported yet");//TODO
					//if(!obj_div(mv, m, idx0))
					//	return false;
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
	return r_mod_lp(m, space_sensitive);
}
int			solve(std::string &str, bool again)
{
	auto t1=__rdtsc();
	int nspaces=16;
	str.insert(str.size(), nspaces, ' ');
	//const char spaces[]="                ";
	//const int nspaces=sizeof(spaces)-1;//should be at least 16
	//str+=spaces;

	lex_init(str.c_str(), (int)(str.size()-nspaces));
	parse_incomplete=false;
	if(!g_answers.size()||!again&&g_answers.back().v.size()&&g_answers.back().v[0].data)
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
	case T_KEYWORDS:
		print_sorted_keywords();
		//for(int k=0;k<T_NKEYWORDS;++k)
		//{
		//	if(keywords[k])
		//	{
		//		printf("%s%c", keywords[k], !((k+1)&3)?'\n':'\t');
		//		//if((k+1)&3)
		//		//	printf("%s\t", keywords[k]);
		//		//else
		//		//	printf("%s\n", keywords[k]);
		//	}
		//}
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_OPEN:
		if(get_str_from_file(str))
		{
			str.insert(str.size(), nspaces, ' ');
			//str+=spaces;
			lex_init(str.c_str(), (int)(str.size()-nspaces));
			goto parse_file;
		}
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_BIN:
		base=2;
		printf("Binary\n");
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_DEC:
		base=10;
		printf("Decimal\n");
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_HEX:
		base=16;
		printf("Hexadecimal\n");
		ret=SOLVE_OK_NO_ANS;
		break;
	case T_CLEAR:
		g_vars.clear();
		g_answers.clear();
		printf("Previous operations cleared.\n");
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
	//case T_SHA256:
	//	ret=SOLVE_PARSE_ERROR;
	//	break;
	//case T_SHA512:
	//	ret=SOLVE_PARSE_ERROR;
	//	break;
	//case T_AES128:
	//	ret=SOLVE_PARSE_ERROR;
	//	break;
	//case T_AES256:
	//	ret=SOLVE_PARSE_ERROR;
	//	break;
	//case T_FRACTIONS:
	//	if(print_fractions=!print_fractions)
	//		printf("Fractions: ON\n");
	//	else
	//		printf("Fractions: OFF\n");
	//	ret=SOLVE_OK_NO_ANS;
	//	break;
	//case T_TOLERANCE:
	//	if(lex_get(false)!=T_NUMBER)
	//		printf("Expected a tolerance value\n");
	//	else
	//	{
	//		g_tolerance=lex_number.to_int();
	//		printf("algorithm tolerance = %g\n", g_tolerance);
	//	}
	//	ret=SOLVE_OK_NO_ANS;
	//	break;
	case T_GFSET://TODO
		ret=SOLVE_PARSE_ERROR;
		break;
	case T_ASCII:
		{
			printf("ASCII table:\n");
			const char *nonprintables[]=
			{
				"Null terminator",
				"Start of Heading",
				"Start of Text",
				"End of Text",
				"End of Transmission",
				"Enquiry",
				"Acknowledgement",
				"\\a Bell",
				"\\b Backspace",
				"\\t Horizontal Tab",
				"\\n Line Feed",
				"\\v Vertical Tab",
				"\\f Form Feed",
				"\\r Carriage Return",
				"Shift Out",
				"Shift In",
				"Data Link Escape",
				"Device Control 1 (often XON)",
				"Device Control 2",
				"Device Control 3 (often XOFF)",
				"Device Control 4",
				"Negative Acknowledgement",
				"Synchronous Idle",
				"End of Transmission Block",
				"Cancel",
				"End of Medium",
				"Substitute",
				"\\e Escape",
				"File Separator",
				"Group Separator",
				"Record Separator",
				"Unit Separator",
			};
			for(int k=0;k<32;++k)
				printf("%3d\t0x%02X\t%s\n", k, k, nonprintables[k]);
			//printf("  0\t0x00\t\'\\0\'\n");
			//printf("  7\t0x07\t\'\\a\'\n");
			//printf("  8\t0x08\t\'\\b\'\n");
			//printf("  9\t0x09\t\'\\t\'\n");
			//printf(" 10\t0x0A\t\'\\n\'\n");
			//printf(" 11\t0x0B\t\'\\v\'\n");
			//printf(" 13\t0x0D\t\'\\r\'\n");
			printf("\n");
			for(int k=32;k<127;++k)
				printf("%3d\t0x%X\t\'%c\'\n", k, k, (char)k);
			printf("127\t0x7F\tDelete\n");
			ret=SOLVE_OK_NO_ANS;
		}
		break;
	case T_BENCHMARK:
		if(benchmark=!benchmark)
			printf("Benchmark ON\n");
		else
			printf("Benchmark OFF\n");
		ret=SOLVE_OK_NO_ANS;
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
	auto t2=__rdtsc();
	if(benchmark)
		printf("%lld cycles\n", t2-t1);
	return ret;
}