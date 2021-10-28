//mc2_lexer.cpp - MCalc lexer
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
#include	"mc2_aligned_vector.h"
#include	<string.h>
#include	<tmmintrin.h>
#include	<assert.h>
#include	<vector>
#include	<algorithm>
const char	file[]=__FILE__;
const char	*keywords[]=
{
#define		TOKEN(STR, LABEL)	STR,
#include	"mc2_keywords.h"
#undef		TOKEN
};

StringLibrary strings;

double		lex_number=0;
char		*lex_id=nullptr;

static bool	lexer_initialized=false;
const char	*text=nullptr;
int			text_size=0, idx=0;

typedef unsigned char byte;
#define			CASE_MASK		0xDF
#define			IN_RANGE(CONST_START, VAR, CONST_END_INCLUSIVE)	(unsigned((VAR)-(CONST_START))<unsigned((CONST_END_INCLUSIVE)+1-(CONST_START)))
inline char		is_whitespace(byte c)
{
	return c>=' '||c<='\t'||c>='\r'||c<='\n';
}
inline char		is_idstart(byte c)
{
	return IN_RANGE('A', c, 'Z')||IN_RANGE('a', c, 'z')||c=='_';
}
inline char		is_alphanumeric(byte c)
{
	return IN_RANGE('0', c, '9')||IN_RANGE('A', c, 'Z')||IN_RANGE('a', c, 'z')||c=='_';
}

struct			LexOption
{
	__m128i option, mask;
	int token;
	short len;
	union
	{
		short flags;//endswithnd<<1|endswithnan
		struct{short endswithnan:1, endswithnd:1;};//followed by non-alphanumeric, followed by a non-digit
	};
	const char *name;
};
typedef std::vector<LexOption, AlignedAllocator<LexOption, 16>> LexSlot;
LexSlot		slots[256];

inline int	strcmp_idx(const char *s1, const char *s2)//returns -1 if identical, otherwise the index where character start to differ
{
	int k=0;
	for(;*s1&&*s2&&*s1==*s2;++k, ++s1, ++s2);
	if(!*s1&&!*s2)
		return -1;
	return k;
}
bool		operator<(LexOption const &a, LexOption const &b)
{
	int k=strcmp_idx(a.name, b.name);
	if(k==-1)
		return 0;
	char ac=a.name[k], bc=b.name[k];
	if(!ac)
		ac=127;
	if(!bc)
		bc=127;
	return ac<bc;
}
void		lex_init(const char *str, int len)
{
	if(!lexer_initialized)
	{
		lexer_initialized=true;

		__m128i ones=_mm_set1_epi32(-1);
		for(int kt=0;kt<T_NTOKENS;++kt)
		{
			auto keyword=keywords[kt];
			if(keyword)
			{
				short len=(short)strlen(keyword);
				char endswithnd=len==1&&keyword[0]=='.';//only '.' should end with a non-digit
				LexOption opt=
				{
					_mm_setzero_si128(),
					_mm_setzero_si128(),
					kt,
					len,
					endswithnd<<1|is_alphanumeric(keyword[0]),//only alphanumeric keywords should be followed by non-alphanumeric
					keyword,
				};
				memcpy(&opt.option, keyword, opt.len);
				opt.mask=_mm_cmpeq_epi8(opt.option, _mm_setzero_si128());
				opt.mask=_mm_xor_si128(opt.mask, ones);
				slots[keyword[0]].push_back(opt);
			}//if keyword
		}
		for(int k=0;k<128;++k)
		{
			auto &slot=slots[k];
			if(slot.size()>1)
				std::sort(slot.begin(), slot.end());
		}
	}
	text=str;
	text_size=len;
	idx=0;
}

//number reader
long long		acme_read_integer(const char *text, int size, int base, int start, int *ret_end)
{
	int k;
	byte temp;
	long long ival=0;
	int hex=-(base==16);

	for(k=start;k<size;++k)
	{
		temp=text[k];
		byte c=temp-'0';
		if(c<10)
		{
			ival*=base;
			ival+=c;
			continue;
		}
		else if(hex&&(c=(temp&0xDF)-'A', c<6))
		{
			ival*=base;
			ival+=10+c;
			//ival+=10+temp-'A';
			continue;
		}
		//else if(hex&(temp-'a'<6))
		//{
		//	ival*=base;
		//	ival+=10+temp-'a';
		//	continue;
		//}
		else if(temp=='\'')//allow quote as separator
			continue;
		break;
	}
	*ret_end=k;
	return ival;
}
static double	acme_read_tail(const char *text, int size, double inv_base, int start, int *ret_end)
{
	int k;
	byte temp;
	double fval=0;

	for(*ret_end=start;*ret_end<size&&(text[*ret_end]>='0'&&text[*ret_end]<='9'||text[*ret_end]>='A'&&text[*ret_end]<='F'||text[*ret_end]>='a'&&text[*ret_end]<='f'||text[*ret_end]=='\'');++*ret_end);
	for(k=*ret_end-1;k>=start;--k)
	{
		temp=text[k];
		byte c=temp-'0';
		if(c<10)
			fval+=c;
		else if(c=(temp&0xDF)-'A', c<6)
			fval+=c+10;
		//else if(temp>='A'&&temp<='F')
		//	fval+=temp-'A'+10;
		//else if(temp>='a'&&temp<='f')
		//	fval+=temp-'a'+10;
		else if(temp=='\'')//allow quote as separator
			continue;
		fval*=inv_base;
	}
	return fval;
}
static char		acme_read_number_pt2(const char *text, int size, int base, double invbase, int start, int *advance, long long *ival, double *fval)
{
	char isfloat=0, neg_exponent;
	int end;
	int exponent;
	char temp;

	*ival=acme_read_integer(text, size, base, start, &end);
	if(end<size&&text[end]=='.')
	{
		int start2=end+1;
		*fval=acme_read_tail(text, size, invbase, start2, &end);
		*fval+=(double)*ival;
		isfloat=1;
	}
	if(end+1<size)
	{
		temp=text[end]&0xDF;
		if(base==10&&temp=='E'||temp=='P')
		{
			start=end+1;
			neg_exponent=text[start]=='-';
			start+=text[start]=='-'||text[start]=='+';

			exponent=(int)acme_read_integer(text, size, base, start, &end);
			if(!isfloat)
				*fval=(double)*ival;
			if(neg_exponent)
				exponent=-exponent;
			*fval*=power((double)base, exponent);
			isfloat=1;
		}
	}
	*advance=end-start;
	return isfloat;
}
char			acme_read_number(const char *text, int size, int k, int *advance, long long *ival, double *fval)//returns  0: integer, 1: float, 2: error
{
	char temp, isfloat;

	assert(ival&&fval);

	isfloat=0;
	*ival=0;
	if(text[k]=='0'&&(k+1==size||is_alphanumeric(text[k+1])))//identify base prefix
	{
		temp=text[k+1]&0xDF;
		if(temp=='X')//hex
			isfloat=acme_read_number_pt2(text, size, 16, 1./16, k+2, advance, ival, fval);
		else if(temp=='B')//bin
			isfloat=acme_read_number_pt2(text, size, 2, 1./2, k+2, advance, ival, fval);
		else if(text[k+1]>='0'&&text[k+1]<='7')//oct	//TODO: "09" should generate an error
			isfloat=acme_read_number_pt2(text, size, 8, 1./8, k+1, advance, ival, fval);
		else if(text[k+1]=='.')//dec	//eg: "0.1"
			isfloat=acme_read_number_pt2(text, size, 10, 1./10, k, advance, ival, fval);
		else
			isfloat=2;
	}
	else//dec
		isfloat=acme_read_number_pt2(text, size, 10, 1./10, k, advance, ival, fval);
	return isfloat;
}

TokenType	lex_get()
{
//again:
	if(idx>=text_size)
		return T_EOF;
	auto &slot=slots[(byte)text[idx]];
	int nopts=slot.size();
	if(nopts)
	{
		__m128i val=_mm_loadu_si128((__m128i*)(text+idx));//text is over-allocated
		for(int ko=0;ko<nopts;++ko)
		{
			auto &opt=slot[ko];
			__m128i v2=_mm_and_si128(val, opt.mask);
			v2=_mm_cmpeq_epi8(v2, opt.option);
			int result=_mm_movemask_epi8(v2);
			if(result==0xFFFF)
			{
				char next=text[idx+opt.len];
				if(!opt.endswithnan||!(IN_RANGE('0', next, '9')|IN_RANGE('A', next, 'Z')|IN_RANGE('a', next, 'z')|(next=='_')))
				{
					idx+=opt.len;
					return (TokenType)opt.token;
				}
			}
		}
	}
	switch(text[idx])
	{
	case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':case '8':case '9':case '.':
		{
			long long ival;
			double fval;
			int advance;
			char ret=acme_read_number(text, text_size, idx, &advance, &ival, &fval);
			if(ret==2||!advance)
				goto unrecognized;
			if(ret)
				lex_number=fval;
			else
				lex_number=(double)ival;
			idx+=advance;
			return T_NUMBER;
		}
		break;
	case ' ':case '\t'://whitespace
	case '\r':case '\n':
		{
			char newline=0;
			for(;idx<text_size&&text[idx]&&(text[idx]==' '||text[idx]=='\t'||text[idx]=='\r'||text[idx]=='\n');newline|=text[idx]=='\r'||text[idx]=='\n', ++idx);
			return TokenType(T_SPACE+newline);
		}
	default:
		if(!is_idstart(text[idx]))
			goto unrecognized;
		else
		{
			int start=idx;
			for(;idx<text_size&&text[idx]&&is_alphanumeric(text[idx]);++idx);
			lex_id=strings.add(text+start, idx-start);
			return T_ID;
		}
	}
unrecognized:
	printf("Error: unrecognized text at %d: \'%.*s\'\n", idx, idx+50>text_size?50:text_size-idx, text+idx);
	return T_IGNORED;
}
TokenType	lex_look_ahead(int k)
{
	int idx0=idx;
	TokenType type=T_IGNORED;
	for(int k2=-1;k2<k;++k2)
		type=lex_get();
	idx=idx0;
	return type;
}