#ifndef		AES_H
#define		AES_H
class		AES
{
	static unsigned char (*key)[4][4], Dkey[11][4][4];
	static unsigned char mult_gf_2_8(unsigned char a, unsigned char b);
	static unsigned char mult_by_4_gf_2_8(unsigned char x){return mult_by_2[mult_by_2[x]];}
	static unsigned char mult_by_8_gf_2_8(unsigned char x){return mult_by_2[mult_by_2[mult_by_2[x]]];}
	static int leftmost_up_bit_pos(int x);
	static int mult_gf_2(int a, int b);
	static int divide_gf_2(int a, int b, int *r);
	static bool mult_inv_gf_2_8(unsigned char x, unsigned char &x_1);
	static void s_box_step_4(unsigned char &x);
	static void s_box_1_step_3(unsigned char &x);
	static void add_round_key(int round);
	static void substitute_bytes();
	static void inverse_sub_bytes();
	static void shift_rows();
	static void inverse_shift_rows();
	static void mix_columns();
	static void inverse_mix_columns();
public:
	static unsigned char s_box[256], s_box_1[256], mult_by_2[256], x_pow_i_4_1[11];
	static unsigned int DK0[256], DK1[256], DK2[256], DK3[256],
		E0[256], E1[256], E2[256], E3[256], SBS8[256], SBS16[256], SBS24[256],
		D0[256], D1[256], D2[256], D3[256], SB1S8[256], SB1S16[256], SB1S24[256];
	static void initiate();
	static void select_AES_NI(bool AES_NI)
	{
		if(AES_NI)
			enc=encrypt_ni, dec=decrypt_ni;
		else
			enc=encrypt_ia32, dec=decrypt_ia32;
	}
	static void expand_key(unsigned char *key);
	static void encrypt_ni(unsigned char *text);
	static void encrypt_ia32(unsigned char *text);
	static void encrypt_sse2(unsigned char *text);//
	static void encrypt(unsigned char *text, unsigned char *key);
	static void (*enc)(unsigned char *text);
	static void decrypt_ni(unsigned char *text);
	static void decrypt_ia32(unsigned char *text);
	static void (*dec)(unsigned char *text);
	static void decrypt(unsigned char *text, unsigned char *key);
};
#endif