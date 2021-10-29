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

#include	"mc2.h"
#include	<stdarg.h>
#include	<map>
#include	<math.h>
static const char file[]=__FILE__;

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
inline bool	obj_add(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dy==1&&m2.dx==1)//addition of a scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
			m.data[k]+=m2.data[0];
	}
	else if(m.dy==1&&m.dx==1)//addition of a scalar
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
			m2.data[k]+=m.data[0];
		//if(m2_is_temp)
			m.move2temp(m2);
		//else
		//	m=m2;
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
			m.data[k]-=m2.data[0];
	}
	else if(m.dy==1&&m.dx==1)//scalar minus matrix
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
			m2.data[k]=m.data[0]-m2.data[k];
		//if(m2_is_temp)
			m.move2temp(m2);
		//else
		//	m=m2;
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
	if(m2.dx==1&&m2.dy==1)//multiplication by scalar
	{
		for(int k=0, size=m.dy*m.dx;k<size;++k)
			m.data[k]*=m2.data[0];
	}
	else if(m.dx==1&&m.dy==1)//multiplication by scalar
	{
		for(int k=0, size=m2.dy*m2.dx;k<size;++k)
			m2.data[k]*=m.data[0];
		//if(m2_is_temp)
			m.move2temp(m2);
		//else
		//	m=m2;

		//m=std::move(m2);

		//free(m.data);
		//m.data=m2.data;
		//m.dx=m2.dx, m.dy=m2.dy;
	}
	else//matrix multiplication
	{
		if(m.dx!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
		auto DALLOC(temp, m.dy*m2.dx);
		//double *temp=(double*)malloc(m.dy*m2.dx*sizeof(double));
		impl_matmul(temp, m.data, m2.data, m.dy, m.dx, m2.dx);
		free(m.data);
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
			m.data[k]/=m2.data[0];
	}
	else//matrix division
	{
		if(m2.dx!=m2.dy)//denominator matrix must be square
			return user_error2(idx0, idx, "Denominator matrix must be square: h2=%d != w2=%d", m2.dy, m2.dx);
		if(m.dx!=m2.dy)//dimension mismatch
			return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
		REALLOC(double, m2.data, m2.data, m2.dx*m2.dy*2);
		//m2.data=(double*)realloc(m2.data, m2.dx*m2.dy*2*sizeof(double));
		auto DALLOC(temp, m.dx*m.dy);
		//double *temp=(double*)malloc(m.dx*m.dy*sizeof(double));
		impl_matdiv(temp, m.data, m2.data, m.dy, m2.dx);
		free(m.data);
		m.data=temp;
	}
	return true;
}
inline bool	obj_mod(Matrix &m, Matrix &m2, int idx0)
{
	if(m2.dx==1&&m2.dy==1)//m2 can be scalar
	{
		for(int k=0, size=m.dx*m.dy;k<size;++k)
			m.data[k]=m.data[k]-floor(m.data[k]/m2.data[0])*m2.data[0];
	}
	else
	{
		if(m2.dx!=m.dx||m2.dy!=m.dy)
			return user_error2(idx0, idx, "The modulus operator \'%%\' is element-wise: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
		for(int k=0, size=m.dx*m.dy;k<size;++k)
			m.data[k]=m.data[k]-floor(m.data[k]/m2.data[k])*m2.data[k];
	}
	return true;
}

int			parse_incomplete=false;

//precedence (first to last):
//	function call		f(a)
//	transpose			a'			< these two are swapped, but it doesn't matter, or does it?
//	signs				-a	+a		<
//	multicplicative		a*b		a o b	a/b		a\b		a%b		a.*b	a./b
//	additive			a+b		a-b
//	relational			a<b		a<=b	a>b		a>=b
//	equality			a==b	a!=b
//	assignment			a=b

//parser declarations
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
				DREALLOC(m.data, m.data, newdx*m.dy);
				//m.data=(double*)realloc(m.data, newdx*m.dy*sizeof(double));
				for(int ky=m.dy-1;ky>=0;--ky)
				{
					DMEMCPY(m.data+newdx*ky, m.data+m.dx*ky, m.dx);
					//memcpy(m.data+newdx*ky, m.data+m.dx*ky, m.dx*sizeof(double));
					DMEMCPY(m.data+newdx*ky+m.dx, m2.data+m2.dx*ky, m2.dx);
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
			if(m.dx!=1||m.dy!=1)
			{
				double *data=m.data;
				DALLOC(m.data, m.dx*m.dy);
				//m.data=(double*)malloc(m.dx*m.dy*sizeof(double));
				std::swap(m.dx, m.dy);
				int ncomp=1+(m.type==T_COMPLEX);
				for(int ky=0;ky<m.dy;++ky)
					for(int kx=0;kx<m.dx;++kx)
						for(int k=0;k<ncomp;++k)
							m.data[ncomp*(m.dx*ky+kx)+k]=data[ncomp*(m.dy*kx+ky)+k];
				free(data);
				m.name=nullptr;
			}
			continue;
		case T_LBRACKET:
			{
				Matrix m2;
				if(!r_assign_expr(m2, false))
					return false;
				if(lex_get(space_sensitive)!=T_RBRACKET)
					return user_error2(idx0, idx, "Expected a closing bracket \']\'");
				if(m2.dx!=1||m2.dy!=1)
					return user_error2(idx0, idx, "Expected a scalar");
				int a_idx=(int)floor(m2.data[0]+0.5);
				if(a_idx<0)
					a_idx=0;
				if(m.dy>1)//select row
				{
					if(a_idx>=m.dy)//TODO: out-of-bounds error?
						a_idx=m.dy-1;
					if(a_idx)
						DMEMCPY(m.data, m.data+m.dx*a_idx, m.dx);
					m.dy=1;
				}
				else//select column
				{
					if(a_idx>=m.dx)
						a_idx=m.dx-1;
					m.data[0]=m.data[a_idx];
					m.dx=1;
				}
				DREALLOC(m.data, m.data, m.dx*m.dy);
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
		switch(auto t=lex_get(space_sensitive))
		{
			//region - functions
#if 1
		case T_ANS:
		case T_ROOTS:
		case T_SAMPLE:
		case T_INVZ:
		case T_IDEN:
		case T_SUM:
		case T_REF:
		case T_RREF:
		case T_DET:
		case T_INV:
		case T_DIAG:
		case T_LU:
		case T_TRACE:
		case T_COS:case T_ACOS:case T_COSD:case T_ACOSD:case T_COSH:case T_ACOSH:
		case T_SEC:case T_ASEC:case T_SECD:case T_ASECD:case T_SECH:case T_ASECH:
		case T_SIN:case T_ASIN:case T_SIND:case T_ASIND:case T_SINH:case T_ASINH:
		case T_CSC:case T_ACSC:case T_CSCD:case T_ACSCD:case T_CSCH:case T_ACSCH:
		case T_TAN:case T_ATAN:case T_TAND:case T_ATAND:case T_TANH:case T_ATANH:
		case T_COT:case T_ACOT:case T_COTD:case T_ACOTD:case T_COTH:case T_ACOTH:
		case T_GAMMA:
			m.name=nullptr;
			if(lex_get(false)!=T_LPR)
				return user_error2(idx0, idx, "Expected an opening parenthesis \'(\' of function call");
			if(!r_assign_expr(m, false))
				return false;
			if(lex_get(false)!=T_RPR)
				return user_error2(idx0, idx, "Expected a closing parenthesis \')\' of function call");
			switch(t)
			{
			case T_ANS:
				{
					if(m.dx>1||m.dy>1)//expected a scalar
						return user_error2(idx0, idx, "Expected a scalar");
					int a_idx=(int)floor(m.data[0]+0.5);
					if(a_idx<0)
						a_idx=0;
					if(a_idx>=(int)g_answers.size())
						a_idx=g_answers.size()-1;
					m=g_answers[a_idx];
				}
				break;
			case T_ROOTS://TODO
				return my_error(idx0, idx);
			case T_SAMPLE://TODO
				return my_error(idx0, idx);
			case T_INVZ://TODO
				return my_error(idx0, idx);
			case T_IDEN:
				{
					if(m.dx>1||m.dy>1)//expected a scalar
						return user_error2(idx0, idx, "Expected a scalar");
					int a_idx=(int)floor(m.data[0]+0.5);
					if(a_idx<1)//wrong size
						return user_error2(idx0, idx, "Wrong size of identity matrix");
					m.dx=m.dy=a_idx;
					a_idx*=a_idx;
					DALLOC(m.data, a_idx);
					//m.data=(double*)malloc(a_idx*sizeof(double));
					memset(m.data, 0, a_idx*sizeof(double));
					for(int k=0;k<a_idx;k+=m.dx+1)
						m.data[k]=1;
				}
				break;
			case T_SUM://should work like in Matlab
				{
					if(m.dy==1)
					{
						for(int kx=1;kx<m.dx;++kx)
							m.data[0]+=m.data[kx];
						m.dx=1;
					}
					else
					{
						for(int ky=1;ky<m.dy;++ky)
							for(int kx=0;kx<m.dx;++kx)
								m.data[kx]+=m.data[m.dx*ky+kx];
						m.dy=1;
					}
					DREALLOC(m.data, m.data, m.dy*m.dx);
					//m.data=(double*)realloc(m.data, m.dy*m.dx*sizeof(double));
				}
				break;
			case T_REF:
				impl_ref(m.data, m.dx, m.dy);
				break;
			case T_RREF:
				impl_rref(m.data, m.dx, m.dy);
				break;
			case T_DET:
				if(m.dx!=m.dy)//must be square
					return user_error2(idx0, idx, "Expected a square matrix");
				DREALLOC(m.data, m.data, m.dx*m.dy*2);
				//m.data=(double*)realloc(m.dx*m.dy*sizeof(double));
				impl_det(m.data, m.dx);
				m.dx=m.dy=1;
				DREALLOC(m.data, m.data, m.dx*m.dy);
				break;
			case T_INV:
				if(m.dy!=m.dx)
					return user_error2(idx0, idx, "Expected a square matrix");
				DREALLOC(m.data, m.data, m.dx*m.dy*2);
				impl_matinv(m.data, m.dx);
				DREALLOC(m.data, m.data, m.dx*m.dy);
				break;
			case T_DIAG://TODO
			case T_LU:
				return my_error(idx0, idx);
			case T_TRACE:
				if(m.dy!=m.dx)
					return user_error2(idx0, idx, "Expected a square matrix");
				for(int k=1;k<m.dx;++k)
					m.data[0]+=m.data[(m.dx+1)*k];
				m.dx=m.dy=1;
				DREALLOC(m.data, m.data, m.dx*m.dy);
				break;
#define		EW_FUNC(FUNC)	for(int k=0, size=m.dx*m.dy;k<size;++k)m.data[k]=FUNC(m.data[k]);
			case T_COS:		EW_FUNC(cos)break;
			case T_ACOS:	EW_FUNC(acos)break;
			case T_COSD:	EW_FUNC(cosd)break;
			case T_ACOSD:	EW_FUNC(acosd)break;
			case T_COSH:	EW_FUNC(cosh)break;
			case T_ACOSH:	EW_FUNC(acosh)break;
			case T_SEC:		EW_FUNC(sec)break;
			case T_ASEC:	EW_FUNC(asec)break;
			case T_SECD:	EW_FUNC(secd)break;
			case T_ASECD:	EW_FUNC(asecd)break;
			case T_SECH:	EW_FUNC(sech)break;
			case T_ASECH:	EW_FUNC(asech)break;
			case T_SIN:		EW_FUNC(sin)break;
			case T_ASIN:	EW_FUNC(asin)break;
			case T_SIND:	EW_FUNC(sind)break;
			case T_ASIND:	EW_FUNC(asind)break;
			case T_SINH:	EW_FUNC(sinh)break;
			case T_ASINH:	EW_FUNC(asinh)break;
			case T_CSC:		EW_FUNC(csc)break;
			case T_ACSC:	EW_FUNC(acsc)break;
			case T_CSCD:	EW_FUNC(cscd)break;
			case T_ACSCD:	EW_FUNC(acscd)break;
			case T_CSCH:	EW_FUNC(csch)break;
			case T_ACSCH:	EW_FUNC(acsch)break;
			case T_TAN:		EW_FUNC(tan)break;
			case T_ATAN:	EW_FUNC(atan)break;
			case T_TAND:	EW_FUNC(tand)break;
			case T_ATAND:	EW_FUNC(atand)break;
			case T_TANH:	EW_FUNC(tanh)break;
			case T_ATANH:	EW_FUNC(atanh)break;
			case T_COT:		EW_FUNC(cot)break;
			case T_ACOT:	EW_FUNC(acot)break;
			case T_COTD:	EW_FUNC(cotd)break;
			case T_ACOTD:	EW_FUNC(acotd)break;
			case T_COTH:	EW_FUNC(coth)break;
			case T_ACOTH:	EW_FUNC(acoth)break;
			case T_GAMMA:	EW_FUNC(tgamma)break;
			case T_LNGAMMA:	EW_FUNC(lgamma)break;
#undef		EW_FUNC
				break;
			}
			return r_postfix(m, space_sensitive);
		case T_CMD:
		case T_CONV:
		case T_LDIV:
		case T_CROSS:
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
						if(m.dx!=1||m.dy!=1||m2.dx!=1||m2.dy!=1)//expected 2 scalars
							return user_error2(idx0, idx, "Expected two scalar arguments");
						int w=(int)floor(*m.data+0.5), h=(int)floor(*m2.data+0.5);
						*m.data=set_console_buffer_size(w, h);
					}
					break;
				case T_CONV:
					{
						if(m.dy!=1||m2.dy!=1)//expected 2 row vectors
							return user_error2(idx0, idx, "Expected two row vectors");
						auto DALLOC(temp, m.dx+m2.dx-1);
						impl_polmul(temp, m.data, m2.data, m.dx, m2.dx, 0);
						free(m.data);
						m.data=temp;
						m.dx+=m2.dx-1;
					}
					break;
				case T_LDIV://TODO
					return my_error(idx0, idx);
				case T_CROSS:
					{
						int size1=m.dx*m.dy, size2=m2.dx*m2.dy;
						if(size1==2&&size2==2)
						{
							double temp=m.data[0]*m2.data[1]-m.data[1]*m2.data[0];
							DREALLOC(m.data, m.data, 1);
							*m.data=temp;
						}
						else if(size1==3&&size2==3)
						{
							auto DALLOC(temp, 3);
							temp[0]=m.data[1]*m2.data[2]-m.data[2]*m2.data[1];
							temp[1]=m.data[2]*m2.data[0]-m.data[0]*m2.data[2];
							temp[2]=m.data[0]*m2.data[1]-m.data[1]*m2.data[0];
							free(m.data);
							m.data=temp;
						}
						else
							return user_error2(idx0, idx, "cross() expects two 2D or two 3D vectors");
					}
					break;
				}
			}
			return r_postfix(m, space_sensitive);
#endif
			//region - constants
#if 1
		case T_NUMBER:
			m.name=nullptr;
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=lex_number;
			return r_postfix(m, space_sensitive);
		case T_IMAG:
		case T_IMAG_UNUSED:
			m.name=nullptr;
			m.type=T_COMPLEX;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 2);
			//m.data=(double*)malloc(2*sizeof(double));
			m.data[0]=0;
			m.data[1]=1;
			return r_postfix(m, space_sensitive);
		case T_EULER:
			m.name=nullptr;
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_e;
			return r_postfix(m, space_sensitive);
		case T_PI:
			m.name=nullptr;
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_pi;
			return r_postfix(m, space_sensitive);
		case T_INF:
			m.name=nullptr;
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_HUGE;
			return r_postfix(m, space_sensitive);
		case T_NAN:
			m.name=nullptr;
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_HUGE-_HUGE;
			return r_postfix(m, space_sensitive);
#endif
		case T_MINUS:
			m.name=nullptr;
			if(!r_unary(m, space_sensitive))
				return false;
			for(int k=0, size=m.dx*m.dy*(1+(m.type==T_COMPLEX));k<size;++k)
				m.data[k]=-m.data[k];
			return r_postfix(m, space_sensitive);
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
		case T_LBRACKET:
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
							DREALLOC(m.data, m.data, m.dx*newdy);
							//m.data=(double*)realloc(m.data, m.dx*newdy*sizeof(double));
							DMEMCPY(m.data+m.dx*m.dy, m2.data, m2.dx*m2.dy);
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
				auto DALLOC(temp, m.dx*m.dy*m2.dx*m2.dy);
				//double *temp=(double*)malloc(m.dx*m.dy*m2.dx*m2.dy*sizeof(double));
				impl_tensor(temp, m.data, m2.data, m.dx, m.dy, m2.dx, m2.dy);
				free(m.data);
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
				if(m.dx==1&&m.dy==1)//multiplication by scalar
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						m2.data[k]*=m.data[0];
					m.move2temp(m2);
				}
				else//matrix division
				{
					if(m.dx!=m.dy)//denominator matrix must be square
						return user_error2(idx0, idx, "Denominator matrix must be square: h1=%d != w1=%d", m.dy, m.dx);
					if(m.dx!=m2.dy)//dimension mismatch
						return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
					REALLOC(double, m.data, m.data, m2.dx*m2.dy);
					//m.data=(double*)realloc(m.data, m.dx*m.dy*2*sizeof(double));
					auto DALLOC(temp, m2.dx*m2.dy);
					//double *temp=(double*)malloc(m2.dx*m2.dy*sizeof(double));
					impl_matdiv_back(temp, m.data, m2.data, m.dy, m2.dx);
					free(m.data);
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
				if(m2.dx==1&&m2.dy==1)//multiplication by scalar
				{
					for(int k=0, size=m.dy*m.dx;k<size;++k)
						m.data[k]*=m2.data[0];
				}
				else if(m.dx==1&&m.dy==1)//multiplication by scalar
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						m2.data[k]*=m.data[0];
					m.move2temp(m2);
				}
				else
				{
					if(m2.dx!=m.dx||m2.dy!=m.dy)
						return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.data[k]*=m2.data[k];
				}
			}
			continue;
		case T_DIV_EW:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m2.dx==1&&m2.dy==1)//division by scalar
				{
					for(int k=0, size=m.dy*m.dx;k<size;++k)
						m.data[k]/=m2.data[0];
				}
				else if(m.dx==1&&m.dy==1)//element-wise division of scalar by matrix
				{
					for(int k=0, size=m2.dy*m2.dx;k<size;++k)
						m2.data[k]=m.data[0]/m2.data[k];
					m.move2temp(m2);
				}
				else
				{
					if(m2.dx!=m.dx||m2.dy!=m.dy)
						return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.data[k]/=m2.data[k];
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
bool		r_relational(Matrix &m, bool space_sensitive)
{
	if(!r_additive(m, space_sensitive))
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
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k]=m.data[k]<m2.data[k];
			}
			continue;
		case T_LESS_EQUAL:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k]=m.data[k]<=m2.data[k];
			}
			continue;
		case T_GREATER:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k]=m.data[k]>m2.data[k];
			}
			continue;
		case T_GREATER_EQUAL:
			{
				m.name=nullptr;
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k]=m.data[k]>=m2.data[k];
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
						if(abs(m.data[k]-m2.data[k])>1e-10)
						{
							equal=false;
							break;
						}
					}
				}
				DREALLOC(m.data, m.data, 1);
				//m.data=(double*)realloc(m.data, sizeof(double));
				m.data[0]=t==T_EQUAL==equal;
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
	if(!g_answers.size()||!again&&g_answers.back().type!=T_IGNORED)
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
			switch(lex_get(false))
			{
			case T_SEMICOLON:
				continue;
			case T_EOF:
				break;
			default:
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