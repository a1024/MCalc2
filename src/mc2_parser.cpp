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
static const char file[]=__FILE__;
const double _e=exp(1.), _pi=acos(-1.);

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
			}
			continue;
		default:
			idx=idx0;
			break;
		}
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
						a_idx=g_answers.size();
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
			}
			break;
		case T_CMD:
		case T_CONV:
		case T_LDIV:
		case T_CROSS:
			{
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
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=lex_number;
			return r_postfix(m, space_sensitive);
		case T_IMAG:
		case T_IMAG_UNUSED:
			m.type=T_COMPLEX;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 2);
			//m.data=(double*)malloc(2*sizeof(double));
			m.data[0]=0;
			m.data[1]=1;
			return r_postfix(m, space_sensitive);
		case T_EULER:
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_e;
			return r_postfix(m, space_sensitive);
		case T_PI:
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_pi;
			return r_postfix(m, space_sensitive);
		case T_INF:
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_HUGE;
			return r_postfix(m, space_sensitive);
		case T_NAN:
			m.type=T_REAL;
			m.dx=1, m.dy=1;
			DALLOC(m.data, 1);
			//m.data=(double*)malloc(sizeof(double));
			m.data[0]=_HUGE-_HUGE;
			return r_postfix(m, space_sensitive);
#endif
		case T_MINUS:
			if(!r_unary(m, space_sensitive))
				return false;
			for(int k=0, size=m.dx*m.dy*(1+(m.type==T_COMPLEX));k<size;++k)
				m.data[k]=-m.data[k];
			return r_postfix(m, space_sensitive);
		case T_PLUS:
			if(!r_unary(m, space_sensitive))
				return false;
			return r_postfix(m, space_sensitive);
		case T_LPR:
			if(!r_assign_expr(m, false))
				return false;
			if(lex_get(false)!=T_RPR)
				return false;
			return r_postfix(m, space_sensitive);
		case T_LBRACKET:
			{
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
			return false;
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
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Dimension mismatch in matrix multiplication: w1=%d != h2=%d", m.dx, m2.dy);
				auto DALLOC(temp, m.dy*m2.dx);
				//double *temp=(double*)malloc(m.dy*m2.dx*sizeof(double));
				impl_matmul(temp, m.data, m2.data, m.dy, m.dx, m2.dx);
				free(m.data);
				m.data=temp;
				m.dx=m2.dx;
			}
			continue;
		case T_TENSOR:
			{
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
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
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
			continue;
		case T_DIV_BACK:
			{
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
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
			continue;
		case T_MOD:
			{
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
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
			}
			continue;
		case T_MUL_EW:
			{
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m2.dx==1&&m2.dy==1)//m2 can be scalar
				{
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.data[k]*=m2.data[0];
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
				Matrix m2;
				if(!r_unary(m2, space_sensitive))
					return false;
				if(m2.dx==1&&m2.dy==1)//m2 can be scalar
				{
					for(int k=0, size=m.dx*m.dy;k<size;++k)
						m.data[k]/=m2.data[0];
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
		default:
			idx=idx0;
			break;
		}
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
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				impl_addbuffers(m.data, m.data, m2.data, m.dy*m.dx);
			}
			continue;
		case T_MINUS:
			{
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				impl_subbuffers(m.data, m.data, m2.data, m.dy*m.dx);
			}
			continue;
		default:
			idx=idx0;
			break;
		}
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
				Matrix m2;
				if(!r_multiplicative(m2, space_sensitive))
					return false;
				if(m.dx!=m2.dx||m.dy!=m2.dy)//dimension mismatch
					return user_error2(idx0, idx, "Element-wise operation: %dx%d != %dx%d", m.dy, m.dx, m2.dy, m2.dx);
				for(int k=0, size=m.dx*m.dy;k<size;++k)
					m.data[k]=m.data[k]>=m2.data[k];
			}
			continue;
		default:
			idx=idx0;
			break;
		}
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
		default:
			idx=idx0;
			break;
		}
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
		if(lex_get(space_sensitive)!=T_ASSIGN)
			goto next;
		if(!r_equality(m, space_sensitive))
			return false;
		m.name=name;
		g_vars[name]=m;
		return true;
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
	if(!again)
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
		printf("Variables:\n");
		for(auto &var:g_vars)
		{
			printf("%s =\n", var.first);
			var.second.print();
		}
		printf("Previous answers:\n");
		for(int k=0;k<(int)g_answers.size()-1;++k)
		{
			auto &ans=g_answers[k];
			printf("ans(%d) =\n", k);
			ans.print();
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