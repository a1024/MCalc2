#include	<immintrin.h>
#include	"aes.h"
bool		aes_ni_supported=false;
void		(*AES::enc)(unsigned char *text)=nullptr, (*AES::dec)(unsigned char *text)=nullptr;
unsigned char	AES::s_box[256], AES::s_box_1[256], AES::mult_by_2[256], AES::x_pow_i_4_1[11];
unsigned int	AES::DK0[256], AES::DK1[256], AES::DK2[256], AES::DK3[256],
	AES::E0[256], AES::E1[256], AES::E2[256], AES::E3[256], AES::SBS8[256], AES::SBS16[256], AES::SBS24[256],
	AES::D0[256], AES::D1[256], AES::D2[256], AES::D3[256], AES::SB1S8[256], AES::SB1S16[256], AES::SB1S24[256];
unsigned char	(*AES::key)[4][4], AES::Dkey[11][4][4];
unsigned char	AES::mult_gf_2_8(unsigned char a, unsigned char b)
{
	int result=0;
	for(int k=0;k<8;++k)
	{
		if(b&1<<k)
			result^=a;
		a=a&0x80?a<<1^0x1B:a<<1;
	}
	return result;
}
int			AES::leftmost_up_bit_pos(int x)
{
	int k=31;
	for(;k>=0;--k)
		if(x&1<<k)
			break;
	return k;
}
int			AES::mult_gf_2(int a, int b)
{
	int result=0;
	for(int k=0;k<32;++k)
		if(b&1<<k)
			result^=a<<k;
	return result;
}
int			AES::divide_gf_2(int a, int b, int *r=0)
{
	int q=0;
	for(int xb=leftmost_up_bit_pos(b);a>=b;)
	{
		int xa_b=leftmost_up_bit_pos(a)-xb;
		q^=1<<xa_b, a^=b<<xa_b;
	}
	if(r)
		*r=a;
	return q;
}
bool		AES::mult_inv_gf_2_8(unsigned char x, unsigned char &x_1)
{
	int Q, A[3]={1, 0, 0x11B}, B[3]={0, 1, x}, T[3];
	for(;B[2]!=1&&B[2]!=0;)
	{
		Q=divide_gf_2(A[2], B[2]);
		T[0]=A[0], T[1]=A[1], T[2]=A[2];
		A[0]=B[0], A[1]=B[1], A[2]=B[2];
		B[0]=T[0]^mult_gf_2(Q, B[0]), B[1]=T[1]^mult_gf_2(Q, B[1]), B[2]=T[2]^mult_gf_2(Q, B[2]);
	}
	if(B[2])
	{
		x_1=B[1];
		return true;
	}
	return false;
}
void		AES::s_box_step_4(unsigned char &x)
{
	int result=0, c=0x63;
	for(int k=0;k<8;++k)
		result^=(x>>k&1^x>>(k+4)%8&1^x>>(k+5)%8&1^x>>(k+6)%8&1^x>>(k+7)%8&1^c>>k&1)<<k;
	x=result;
}
void		AES::s_box_1_step_3(unsigned char &x)
{
	int result=0, d=0x05;
	for(int k=0;k<8;++k)
		result^=(x>>(k+2)%8&1^x>>(k+5)%8&1^x>>(k+7)%8&1^d>>k&1)<<k;
	x=result;
}
void		AES::initiate()
{
	select_AES_NI(aes_ni_supported);
	for(int k=0;k<256;++k)
	{
		s_box[k]=k;
		mult_inv_gf_2_8(s_box[k], s_box[k]);
		s_box_step_4(s_box[k]);

		s_box_1[k]=k;
		s_box_1_step_3(s_box_1[k]);
		mult_inv_gf_2_8(s_box_1[k], s_box_1[k]);
	/*	printf("%02X ", (unsigned char)s_box_1[k]);
		if(!((k+1)%16))
			printf("\n");*/

		mult_by_2[k]=mult_gf_2_8(k, 2);
	}
	for(int k=0;k<256;++k)//2113	3211	1321	1132
	{
		unsigned char c=s_box[k], c2=mult_by_2[c], c3=c2^c;
		E0[k]=c2|c<<8|c<<16|c3<<24, E1[k]=c3|c2<<8|c<<16|c<<24, E2[k]=c|c3<<8|c2<<16|c<<24, E3[k]=c|c<<8|c3<<16|c2<<24, SBS8[k]=c<<8, SBS16[k]=c<<16, SBS24[k]=c<<24;
	}
	for(int k=0;k<256;++k)//E9DB	BE9D	DBE9	9DBE
	{
		unsigned char c=k, c2=mult_by_2[c], c3=c2^c, c4=mult_by_2[c2], c8=mult_by_2[c4], c9=c8^c, cB=c8^c3, cC=c8^c4, cD=cC^c, cE=cC^c2;
		DK0[k]=cE|c9<<8|cD<<16|cB<<24, DK1[k]=cB|cE<<8|c9<<16|cD<<24, DK2[k]=cD|cB<<8|cE<<16|c9<<24, DK3[k]=c9|cD<<8|cB<<16|cE<<24;
	}
	for(int k=0;k<256;++k)//E9DB	BE9D	DBE9	9DBE
	{
		unsigned char c=s_box_1[k], c2=mult_by_2[c], c3=c2^c, c4=mult_by_2[c2], c8=mult_by_2[c4], c9=c8^c, cB=c8^c3, cC=c8^c4, cD=cC^c, cE=cC^c2;
		D0[k]=cE|c9<<8|cD<<16|cB<<24, D1[k]=cB|cE<<8|c9<<16|cD<<24, D2[k]=cD|cB<<8|cE<<16|c9<<24, D3[k]=c9|cD<<8|cB<<16|cE<<24, SB1S8[k]=c<<8, SB1S16[k]=c<<16, SB1S24[k]=c<<24;
	}
	{
		unsigned char smiley=0x8D;
		for(int k=0;k<11;++k)
			smiley=mult_by_2[x_pow_i_4_1[k]=smiley];
	}
}
void		AES::expand_key(unsigned char key[11*16])
{
	*(unsigned char**)&AES::key=key;
	for(int k=16;k<176;k+=4)
	{
		if(!(k%16))
			key[k  ]=s_box[key[k-3]]^x_pow_i_4_1[k/16]^key[k-16],
			key[k+1]=s_box[key[k-2]]^key[k-15],
			key[k+2]=s_box[key[k-1]]^key[k-14],
			key[k+3]=s_box[key[k-4]]^key[k-13];
		else
			key[k  ]=key[k-4]^key[k-16],
			key[k+1]=key[k-3]^key[k-15],
			key[k+2]=key[k-2]^key[k-14],
			key[k+3]=key[k-1]^key[k-13];
	}
	for(int k=0;k<16;k+=4)
		*(int*)((char*)Dkey+k)=*(int*)(key+160+k);
	for(int k=16;k<160;k+=16)
		for(int k2=0;k2<16;k2+=4)
			*(int*)((char*)Dkey+k+k2)=DK0[key[160-k+k2]]^DK1[key[160-k+k2+1]]^DK2[key[160-k+k2+2]]^DK3[key[160-k+k2+3]];
	for(int k=0;k<16;k+=4)
		*(int*)((char*)Dkey+160+k)=*(int*)(key+k);

/*	for(int k=0;k<16;k+=4)
		*(int*)((char*)Dkey+k)=*(int*)(key+k);
	for(int k=16;k<160;k+=4)
		*(int*)((char*)Dkey+k)=DK0[key[k]]^DK1[key[k+1]]^DK2[key[k+2]]^DK3[key[k+3]];
	for(int k=160;k<176;k+=4)
		*(int*)((char*)Dkey+k)=*(int*)(key+k);*/
}
void		AES::encrypt(unsigned char *text, unsigned char *key){expand_key(key), enc(text);}
void		AES::encrypt_ni(unsigned char text[16])
{
	__m128i m_text=_mm_loadu_si128((__m128i*)text);
	__m128i m_key=_mm_loadu_si128((__m128i*)key);
	//std::cout<<"key:\t", print_buffer(key, 16);//
	//std::cout<<"text:\t", print_buffer(text, 16);//
	m_text=_mm_xor_si128(m_text, m_key);
	//print_buffer(&m_text, 16);//
	m_key=_mm_loadu_si128((__m128i*)key+1);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+2);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+3);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+4);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+5);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+6);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+7);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+8);
	m_text=_mm_aesenc_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)key+9);
	m_text=_mm_aesenc_si128(m_text, m_key);
	//for(int r=1;r<10;++r)
	//{
	//	m_key=_mm_loadu_si128((__m128i*)key+r);
	//	m_text=_mm_aesenc_si128(m_text, m_key);
	//	print_buffer(&m_text, 16);//
	//}
	m_key=_mm_loadu_si128((__m128i*)key+10);
	m_text=_mm_aesenclast_si128(m_text, m_key);
	//print_buffer(&m_text, 16);//
	_mm_storeu_si128((__m128i*)text, m_text);
}
void		AES::encrypt_ia32(unsigned char text[16])//https://web.archive.org/web/20120126130601/http://software.intel.com/en-us/articles/optimizing-performance-of-the-aes-algorithm-for-the-intel-pentiumr-4-processor/
{
	//unsigned char text2[16];
	//for(int k=0;k<16;++k)
	//	text2[k]=text[k];
	unsigned char temp0[4][4], temp1[4][4];

	*(int*)temp0[0]=*(int*)key[0][0]^*(int*) text    ;
	*(int*)temp0[1]=*(int*)key[0][1]^*(int*)(text+ 4);
	*(int*)temp0[2]=*(int*)key[0][2]^*(int*)(text+ 8);
	*(int*)temp0[3]=*(int*)key[0][3]^*(int*)(text+12);
	//std::cout<<"key:\t", print_buffer(key, 16);//
	//std::cout<<"text:\t", print_buffer(text, 16);//
	//std::cout<<"xor k:\t", print_buffer(temp0, 16);//
	for(int r=1;r<8;r+=2)
	{
		//printf("%02x %02x %02x %02x\n", temp0[0][0], temp0[1][1], temp0[2][2], temp0[3][3]);//
		//printf("%08x %08x %08x %08x\n", E0[temp0[0][0]], E1[temp0[1][1]], E2[temp0[2][2]], E3[temp0[3][3]]);//
		*(int*)temp1[0]=*(int*)key[r  ][0]^E0[temp0[0][0]]^E1[temp0[1][1]]^E2[temp0[2][2]]^E3[temp0[3][3]];
		*(int*)temp1[1]=*(int*)key[r  ][1]^E0[temp0[1][0]]^E1[temp0[2][1]]^E2[temp0[3][2]]^E3[temp0[0][3]];
		*(int*)temp1[2]=*(int*)key[r  ][2]^E0[temp0[2][0]]^E1[temp0[3][1]]^E2[temp0[0][2]]^E3[temp0[1][3]];
		*(int*)temp1[3]=*(int*)key[r  ][3]^E0[temp0[3][0]]^E1[temp0[0][1]]^E2[temp0[1][2]]^E3[temp0[2][3]];
		*(int*)temp0[0]=*(int*)key[r+1][0]^E0[temp1[0][0]]^E1[temp1[1][1]]^E2[temp1[2][2]]^E3[temp1[3][3]];
		*(int*)temp0[1]=*(int*)key[r+1][1]^E0[temp1[1][0]]^E1[temp1[2][1]]^E2[temp1[3][2]]^E3[temp1[0][3]];
		*(int*)temp0[2]=*(int*)key[r+1][2]^E0[temp1[2][0]]^E1[temp1[3][1]]^E2[temp1[0][2]]^E3[temp1[1][3]];
		*(int*)temp0[3]=*(int*)key[r+1][3]^E0[temp1[3][0]]^E1[temp1[0][1]]^E2[temp1[1][2]]^E3[temp1[2][3]];
		//std::cout<<"r"<<r<<":\t", print_buffer(temp1, 16);//
		//std::cout<<"r"<<r+1<<":\t", print_buffer(temp0, 16);//
	}
	*(int*)temp1[0]=*(int*)key[9][0]^E0[temp0[0][0]]^E1[temp0[1][1]]^E2[temp0[2][2]]^E3[temp0[3][3]];
	*(int*)temp1[1]=*(int*)key[9][1]^E0[temp0[1][0]]^E1[temp0[2][1]]^E2[temp0[3][2]]^E3[temp0[0][3]];
	*(int*)temp1[2]=*(int*)key[9][2]^E0[temp0[2][0]]^E1[temp0[3][1]]^E2[temp0[0][2]]^E3[temp0[1][3]];
	*(int*)temp1[3]=*(int*)key[9][3]^E0[temp0[3][0]]^E1[temp0[0][1]]^E2[temp0[1][2]]^E3[temp0[2][3]];
	*(int*) text    =*(int*)key[10][0]^s_box[temp1[0][0]]^SBS8[temp1[1][1]]^SBS16[temp1[2][2]]^SBS24[temp1[3][3]];
	*(int*)(text+ 4)=*(int*)key[10][1]^s_box[temp1[1][0]]^SBS8[temp1[2][1]]^SBS16[temp1[3][2]]^SBS24[temp1[0][3]];
	*(int*)(text+ 8)=*(int*)key[10][2]^s_box[temp1[2][0]]^SBS8[temp1[3][1]]^SBS16[temp1[0][2]]^SBS24[temp1[1][3]];
	*(int*)(text+12)=*(int*)key[10][3]^s_box[temp1[3][0]]^SBS8[temp1[0][1]]^SBS16[temp1[1][2]]^SBS24[temp1[2][3]];
	//std::cout<<"r9:\t", print_buffer(temp1, 16);//
	//std::cout<<"ct:\t", print_buffer(text, 16), std::cout<<endl;////*/

/*	__m128i t=_mm_loadu_si128((__m128i*)text);
	__m128i k=_mm_load_si128((__m128i*)key);
	t=_mm_xor_si128(t, k);
	for(int r=1;r<10;++r)
	{
		__m128i t0=_mm_set_epi32(E0[t.m128i_i8[12]], E0[t.m128i_i8[ 8]], E0[t.m128i_i8[4]], E0[t.m128i_i8[0]]);
		__m128i t1=_mm_set_epi32(E1[t.m128i_i8[1]], E1[t.m128i_i8[13]], E1[t.m128i_i8[ 9]], E1[t.m128i_i8[5]]);
		__m128i t2=_mm_set_epi32(E2[t.m128i_i8[6]], E2[t.m128i_i8[2]], E2[t.m128i_i8[14]], E2[t.m128i_i8[10]]);
		t=_mm_set_epi32(E3[t.m128i_i8[11]], E3[t.m128i_i8[7]], E3[t.m128i_i8[3]], E3[t.m128i_i8[15]]);
		t=_mm_xor_si128(t, t2);
		t=_mm_xor_si128(t, t1);
		t=_mm_xor_si128(t, t0);
		k=_mm_load_si128((__m128i*)key+r);
		t=_mm_xor_si128(t, k);
	}
	unsigned char (*temp1)[4]=(unsigned char(*)[4])&t;
	*(int*) text    =*(int*)key[10][0]^s_box[temp1[0][0]]^SBS8[temp1[1][1]]^SBS16[temp1[2][2]]^SBS24[temp1[3][3]];
	*(int*)(text+ 4)=*(int*)key[10][1]^s_box[temp1[1][0]]^SBS8[temp1[2][1]]^SBS16[temp1[3][2]]^SBS24[temp1[0][3]];
	*(int*)(text+ 8)=*(int*)key[10][2]^s_box[temp1[2][0]]^SBS8[temp1[3][1]]^SBS16[temp1[0][2]]^SBS24[temp1[1][3]];
	*(int*)(text+12)=*(int*)key[10][3]^s_box[temp1[3][0]]^SBS8[temp1[0][1]]^SBS16[temp1[1][2]]^SBS24[temp1[2][3]];//*/

	//for(int k=0;k<16;++k)text2[k]=((char*)temp0)[k];

/*	add_round_key(0);//round 0
	for(int k=1;k<10;++k)//rounds 1~9
	{
		substitute_bytes();
		shift_rows();
		mix_columns();
		add_round_key(k);
	}
	substitute_bytes();//round 10
	shift_rows();
	add_round_key(10);*/
}
void		AES::encrypt_sse2(unsigned char text[16])
{
	//unsigned char text2[16];
	//for(int k=0;k<16;++k)
	//	text2[k]=text[k];

	__m128i t=_mm_loadu_si128((__m128i*)text);
	__m128i k=_mm_load_si128((__m128i*)key);
	t=_mm_xor_si128(t, k);
	//std::cout<<"key:\t", print_buffer(key, 16);//
	//std::cout<<"text:\t", print_buffer(text, 16);//
	//std::cout<<"xor k:\t", print_buffer(&t, 16);//
	for(int r=1;r<10;++r)
	{
		//printf("%02x %02x %02x %02x\n", t.m128i_u8[0], t.m128i_u8[5], t.m128i_u8[10], t.m128i_u8[15]);//
		__m128i t0=_mm_set_epi32(E0[t.m128i_u8[12]], E0[t.m128i_u8[ 8]], E0[t.m128i_u8[4]], E0[t.m128i_u8[0]]);
		__m128i t1=_mm_set_epi32(E1[t.m128i_u8[1]], E1[t.m128i_u8[13]], E1[t.m128i_u8[ 9]], E1[t.m128i_u8[5]]);
		__m128i t2=_mm_set_epi32(E2[t.m128i_u8[6]], E2[t.m128i_u8[2]], E2[t.m128i_u8[14]], E2[t.m128i_u8[10]]);
		t=_mm_set_epi32(E3[t.m128i_u8[11]], E3[t.m128i_u8[7]], E3[t.m128i_u8[3]], E3[t.m128i_u8[15]]);
		//printf("%08x %08x %08x %08x\n", t0.m128i_i32[0], t1.m128i_i32[0], t2.m128i_i32[0], t.m128i_i32[0]);//
		t=_mm_xor_si128(t, t2);
		t=_mm_xor_si128(t, t1);
		t=_mm_xor_si128(t, t0);
		k=_mm_load_si128((__m128i*)key+r);
		t=_mm_xor_si128(t, k);
		//std::cout<<"r"<<r<<":\t", print_buffer(&t, 16);//
	}
	unsigned char (*temp1)[4]=(unsigned char(*)[4])&t;
	*(int*) text    =*(int*)key[10][0]^s_box[temp1[0][0]]^SBS8[temp1[1][1]]^SBS16[temp1[2][2]]^SBS24[temp1[3][3]];
	*(int*)(text+ 4)=*(int*)key[10][1]^s_box[temp1[1][0]]^SBS8[temp1[2][1]]^SBS16[temp1[3][2]]^SBS24[temp1[0][3]];
	*(int*)(text+ 8)=*(int*)key[10][2]^s_box[temp1[2][0]]^SBS8[temp1[3][1]]^SBS16[temp1[0][2]]^SBS24[temp1[1][3]];
	*(int*)(text+12)=*(int*)key[10][3]^s_box[temp1[3][0]]^SBS8[temp1[0][1]]^SBS16[temp1[1][2]]^SBS24[temp1[2][3]];
	//std::cout<<"ct:\t", print_buffer(text, 16), std::cout<<endl;//

/*	add_round_key(0);//round 0
	for(int k=1;k<10;++k)//rounds 1~9
	{
		substitute_bytes();
		shift_rows();
		mix_columns();
		add_round_key(k);
	}
	substitute_bytes();//round 10
	shift_rows();
	add_round_key(10);*/
}
void		AES::decrypt(unsigned char *text, unsigned char *key){expand_key(key), dec(text);}
void		AES::decrypt_ni(unsigned char *text)
{
	__m128i m_text=_mm_loadu_si128((__m128i*)text);
	__m128i m_key=_mm_loadu_si128((__m128i*)Dkey);
	m_text=_mm_xor_si128(m_text, m_key);
	
	m_key=_mm_loadu_si128((__m128i*)Dkey+1);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+2);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+3);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+4);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+5);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+6);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+7);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+8);
	m_text=_mm_aesdec_si128(m_text, m_key);
	m_key=_mm_loadu_si128((__m128i*)Dkey+9);
	m_text=_mm_aesdec_si128(m_text, m_key);
	//for(int r=1;r<10;++r)
	//{
	//	m_key=_mm_loadu_si128((__m128i*)Dkey+r);
	//	m_text=_mm_aesdec_si128(m_text, m_key);
	//}
	m_key=_mm_loadu_si128((__m128i*)Dkey+10);
	m_text=_mm_aesdeclast_si128(m_text, m_key);
	_mm_storeu_si128((__m128i*)text, m_text);
}
void		AES::decrypt_ia32(unsigned char *text)
{
	unsigned char temp0[4][4], temp1[4][4];

//	*(long long*)temp[0]=*(long long*)text^*(long long*)key[round][0];
//	*(long long*)temp[2]=*(long long*)(text+8)^*(long long*)key[round][2];

	*(int*)temp0[0]=*(int*) text    ^*(int*)Dkey[0][0];
	*(int*)temp0[1]=*(int*)(text+ 4)^*(int*)Dkey[0][1];
	*(int*)temp0[2]=*(int*)(text+ 8)^*(int*)Dkey[0][2];
	*(int*)temp0[3]=*(int*)(text+12)^*(int*)Dkey[0][3];
	for(int r=1;r<8;r+=2)
	{
		*(int*)temp1[0]=*(int*)Dkey[r  ][0]^D0[temp0[0][0]]^D1[temp0[3][1]]^D2[temp0[2][2]]^D3[temp0[1][3]];
		*(int*)temp1[1]=*(int*)Dkey[r  ][1]^D0[temp0[1][0]]^D1[temp0[0][1]]^D2[temp0[3][2]]^D3[temp0[2][3]];
		*(int*)temp1[2]=*(int*)Dkey[r  ][2]^D0[temp0[2][0]]^D1[temp0[1][1]]^D2[temp0[0][2]]^D3[temp0[3][3]];
		*(int*)temp1[3]=*(int*)Dkey[r  ][3]^D0[temp0[3][0]]^D1[temp0[2][1]]^D2[temp0[1][2]]^D3[temp0[0][3]];
		*(int*)temp0[0]=*(int*)Dkey[r+1][0]^D0[temp1[0][0]]^D1[temp1[3][1]]^D2[temp1[2][2]]^D3[temp1[1][3]];
		*(int*)temp0[1]=*(int*)Dkey[r+1][1]^D0[temp1[1][0]]^D1[temp1[0][1]]^D2[temp1[3][2]]^D3[temp1[2][3]];
		*(int*)temp0[2]=*(int*)Dkey[r+1][2]^D0[temp1[2][0]]^D1[temp1[1][1]]^D2[temp1[0][2]]^D3[temp1[3][3]];
		*(int*)temp0[3]=*(int*)Dkey[r+1][3]^D0[temp1[3][0]]^D1[temp1[2][1]]^D2[temp1[1][2]]^D3[temp1[0][3]];
	}
	*(int*)temp1[0]=*(int*)Dkey[9][0]^D0[temp0[0][0]]^D1[temp0[3][1]]^D2[temp0[2][2]]^D3[temp0[1][3]];
	*(int*)temp1[1]=*(int*)Dkey[9][1]^D0[temp0[1][0]]^D1[temp0[0][1]]^D2[temp0[3][2]]^D3[temp0[2][3]];
	*(int*)temp1[2]=*(int*)Dkey[9][2]^D0[temp0[2][0]]^D1[temp0[1][1]]^D2[temp0[0][2]]^D3[temp0[3][3]];
	*(int*)temp1[3]=*(int*)Dkey[9][3]^D0[temp0[3][0]]^D1[temp0[2][1]]^D2[temp0[1][2]]^D3[temp0[0][3]];
	*(int*) text    =*(int*)Dkey[10][0]^s_box_1[temp1[0][0]]^SB1S8[temp1[3][1]]^SB1S16[temp1[2][2]]^SB1S24[temp1[1][3]];
	*(int*)(text+ 4)=*(int*)Dkey[10][1]^s_box_1[temp1[1][0]]^SB1S8[temp1[0][1]]^SB1S16[temp1[3][2]]^SB1S24[temp1[2][3]];
	*(int*)(text+ 8)=*(int*)Dkey[10][2]^s_box_1[temp1[2][0]]^SB1S8[temp1[1][1]]^SB1S16[temp1[0][2]]^SB1S24[temp1[3][3]];
	*(int*)(text+12)=*(int*)Dkey[10][3]^s_box_1[temp1[3][0]]^SB1S8[temp1[2][1]]^SB1S16[temp1[1][2]]^SB1S24[temp1[0][3]];

/*	*(int*)temp0[0]=*(int*) text    ^*(int*)Dkey[10][0];
	*(int*)temp0[1]=*(int*)(text+ 4)^*(int*)Dkey[10][1];
	*(int*)temp0[2]=*(int*)(text+ 8)^*(int*)Dkey[10][2];
	*(int*)temp0[3]=*(int*)(text+12)^*(int*)Dkey[10][3];
	for(int r=9;r>2;r-=2)
	{
		*(int*)temp1[0]=*(int*)Dkey[r  ][0]^D0[temp0[0][0]]^D1[temp0[3][1]]^D2[temp0[2][2]]^D3[temp0[1][3]];
		*(int*)temp1[1]=*(int*)Dkey[r  ][1]^D0[temp0[1][0]]^D1[temp0[0][1]]^D2[temp0[3][2]]^D3[temp0[2][3]];
		*(int*)temp1[2]=*(int*)Dkey[r  ][2]^D0[temp0[2][0]]^D1[temp0[1][1]]^D2[temp0[0][2]]^D3[temp0[3][3]];
		*(int*)temp1[3]=*(int*)Dkey[r  ][3]^D0[temp0[3][0]]^D1[temp0[2][1]]^D2[temp0[1][2]]^D3[temp0[0][3]];
		*(int*)temp0[0]=*(int*)Dkey[r-1][0]^D0[temp1[0][0]]^D1[temp1[3][1]]^D2[temp1[2][2]]^D3[temp1[1][3]];
		*(int*)temp0[1]=*(int*)Dkey[r-1][1]^D0[temp1[1][0]]^D1[temp1[0][1]]^D2[temp1[3][2]]^D3[temp1[2][3]];
		*(int*)temp0[2]=*(int*)Dkey[r-1][2]^D0[temp1[2][0]]^D1[temp1[1][1]]^D2[temp1[0][2]]^D3[temp1[3][3]];
		*(int*)temp0[3]=*(int*)Dkey[r-1][3]^D0[temp1[3][0]]^D1[temp1[2][1]]^D2[temp1[1][2]]^D3[temp1[0][3]];
	}
	*(int*)temp1[0]=*(int*)Dkey[1][0]^D0[temp0[0][0]]^D1[temp0[3][1]]^D2[temp0[2][2]]^D3[temp0[1][3]];
	*(int*)temp1[1]=*(int*)Dkey[1][1]^D0[temp0[1][0]]^D1[temp0[0][1]]^D2[temp0[3][2]]^D3[temp0[2][3]];
	*(int*)temp1[2]=*(int*)Dkey[1][2]^D0[temp0[2][0]]^D1[temp0[1][1]]^D2[temp0[0][2]]^D3[temp0[3][3]];
	*(int*)temp1[3]=*(int*)Dkey[1][3]^D0[temp0[3][0]]^D1[temp0[2][1]]^D2[temp0[1][2]]^D3[temp0[0][3]];
	*(int*) text    =*(int*)Dkey[0][0]^s_box_1[temp1[0][0]]^SB1S8[temp1[3][1]]^SB1S16[temp1[2][2]]^SB1S24[temp1[1][3]];
	*(int*)(text+ 4)=*(int*)Dkey[0][1]^s_box_1[temp1[1][0]]^SB1S8[temp1[0][1]]^SB1S16[temp1[3][2]]^SB1S24[temp1[2][3]];
	*(int*)(text+ 8)=*(int*)Dkey[0][2]^s_box_1[temp1[2][0]]^SB1S8[temp1[1][1]]^SB1S16[temp1[0][2]]^SB1S24[temp1[3][3]];
	*(int*)(text+12)=*(int*)Dkey[0][3]^s_box_1[temp1[3][0]]^SB1S8[temp1[2][1]]^SB1S16[temp1[1][2]]^SB1S24[temp1[0][3]];*/

/*	add_round_key(10);//round 0
	for(int k=9;k>0;--k)//rounds 1~9
	{
		inverse_shift_rows();
		inverse_sub_bytes();
		add_round_key(k);
		inverse_mix_columns();
	}
	inverse_shift_rows();//round 10
	inverse_sub_bytes();
	add_round_key(0);*/
}