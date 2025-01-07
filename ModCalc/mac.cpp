//mac.cpp - main file
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
#include	<stdio.h>
#include	<string>
#include	<iostream>
#ifndef __linux__
#include	<conio.h>
#endif
#include	"mac.h"
const char	file[]=__FILE__;
extern const int	g_buf_size=G_BUF_SIZE;
char				g_buf[G_BUF_SIZE]={0};
char		gfmode=0;//later
static const char helpmsg[]=
//	"01234567890123456789012345678901234567890123456789012345678901234567890123456789"
	"MODCALC: A Modular Arithmetic Calculator\n"
	"Usage:\n"
	"  modcalc [\"expression\"/\"filename\"]\n"
	"\n"
	"Syntax:\n"
	"  Use Matlab-style vectors & matrices, for example:\n"
	"    [a2 a1 a0] or [a2, a1, a0] is a row vector\n"
	"    [a11, a12;  a21, a22]      is a 2x2 matrix\n"
	"  Use same square brackets to access matrix elements, for example:\n"
	"    [1 2; 3 4][0][1] == 2\n"
	"  or use parentheses (index starts at one):\n"
	"    [1 2; 3 4](1,2) == 2\n"
	"  Parentheses of single argument function can be omitted\n"
//	"  For multiline statements:\n"//TODO
//	"    Open a parenthesis at the start of the first command, and close it at the\n"
//	"    end of the last command.\n"
	"  Variable names can be up to 16 characters.\n"//
	"  Multiplication asterix \'*\' is always obligatory.\n"
	"\n"
	"Operators by precedence (first to last):\n"
	"  function call, parentheses, square brackets\n"
	"  \', ^, member access\n"
	"  + - (unary)\n"
	"  * o / \\ % .* ./\n"
	"  + -\n"
	"  :\n"
	"  < <= > >=\n"
	"  == !=\n"
	"  mod\n"
	"  = += -= *= /= \\= %=\n"
	"Notes:\n"
	"  \'     is transpose (eg: M\')\n"
	"  o     is the tensor product\n"
	"  \\     is matrix division from left (first matrix should be square)\n"
	"  .* ./ are element-wise operations\n"
	"  mod   is the same as % but with low precedence\n"
	"\n"
	"Keywords:\n"
	"  help: Prints this info (works in cmd)\n"
	"  clear: Clears all variables from memory\n"
	"  vars: Shows all variables & answer count\n"
	"  open: Opens a text file\n"
	"  benchmark: Toggles CPU cycle counter\n"
	"Crypto:\n"
	"  hist: Histogram of quoted text\n"
	"  sha224/256/384/512: SHA-2 Hash of quoted text/number/variable/matrix\n"
	"  sha3-224/256/384/512: SHA-3\n"
//	"  shake128/256(outbytes, input): SHA-3 extendable output functions\n"
//	"  aes128/256e/d(key, input): Encrypt/decrypt with AES-128/256\n"
	"General functions:\n"
	"  ans(n): The n-th answer\n"
	"  cmd(w, h): Set console buffer size (in characters)\n"
	"  bin/dec/hex: Switch between binary/decimal/hexadecimal\n"
	"Numbers:\n"
	"  min/max/mean/sum/product: Reduce columns then rows (Matlab-style)\n"
	"  rdtsc: Read timestamp counter\n"
	"  rand: 30-bit random number in (seed with RDTSC)\n"
	"  randmod(modulus): Random number in range [0, modulus-1] (seed with RDTSC)\n"
	"  randlog(L): L-bit random number (seed with RDTSC)"
	"  isprime(n): Naive test\n"
	"  mrtest(n, a): Miller-Rabin test of n with witness a\n"
	"  genprime(start): First prime after start\n"
	"  factorize(n): Returns a row vector with prime factors\n"
	"  totient(n): Euler's totient function\n"
	"  carmichael(n): Reduced totient function\n"
	"Modular Arithmetic:\n"
	"  gcd[...]: Greatest common divisor of all elements\n"
	"  lcm[...]: Least common multiple of all elements\n"
	"  invmod(x, n): Inverse of x mod n using extended Euclidean algorithm\n"
	"  ntt/intt(M, [N wx wy anti_cyclic]): Number Theoretic transform\n"
	"    If anti_cyclic is true: provide sqrt(wi) instead of wi\n"
	"  log2/log10(n): floor(log(abs(n))) with the indicated base\n"
	"  proots(n): Primitive roots of n\n"
	"  dlog(x, base, modulus): Discrete log\n"
	"Vectors:\n"
	"  cross(3D vec, 3D vec) -> 3D vec\n"
	"  cross(2D vec, 2D vec) -> scalar\n"
	"  For dot product use transpose (col\'*col or row*row\')\n"
//	"Matrices:\n"
//	"  ref: Row Echelon Form\n"//TODO
//	"  rref: Reduced Row Echelon Form (Gaussian Elimination)\n"
//	"  nullspace: Null space of a matrix\n"
	"Square matrices:\n"
	"  I(n): n-square identity matrix\n"//autosize?
//	"  det: Determinant of square matrix\n"
//	"  inv: Inverse of square matrix\n"
	"  tr: Trace of square matrix\n"
//	"  lu: LU Factorization of square matrix\n"
	"  diag: Make a diagonal matrix\n"
//	"  egval: Eigenvalues of a matrix\n"
//	"  egvec(M, L): Eigenvectors, given matrix M and eigenvalues L\n"
	"Polynomials:\n"
	"  conv: Multiply polynomials\n"
//	"  polpow: Raise a polynomial power an integer\n"
//	"  roots: Find the roots of a polynomial\n"
	"\n"
;
void		print_help()
{
	printf("%s", helpmsg);
}
void		get_str_interactive(std::string &str, const char *cmdstr)
{
	if(gfmode)
		printf("gf");
	if(cmdstr)
		printf("%s", cmdstr);
	printf("> ");
	std::getline(std::cin, str);
}
bool		get_str_from_file(std::string &str)
{
	auto wbuf=open_file_window();
	if(!wbuf)
		return false;
#ifdef __linux__
	FILE *file=fopen(wbuf, "r");
	if(!file)
	{
		printf("%s\n", strerror(errno));
		return false;
	}
#else
	FILE *file=nullptr;
	int ret=_wfopen_s(&file, wbuf, L"r");
	if(ret)
	{
		strerror_s(g_buf, g_buf_size, ret);
		printf("%s\n", g_buf);
		return false;
	}
#endif
	fseek(file, 0, SEEK_END);
	int bytesize=ftell(file);
	fseek(file, 0, SEEK_SET);
	str.resize(bytesize);
	fread(&str[0], 1, bytesize, file);
	fclose(file);
	str.resize(strlen(str.c_str()));//remove extra null terminators at the end
	return true;
}
int			main(int argc, const char **argv)
{
	//set_console_buffer_size(120, 4000);
	printf("MODCALC%s\n\n", argc==2?"":"\t\tCtrl C to exit.");

	std::string str;
	bool quit_prompt=false;
	if(argc>2)
	{
		printf(
			"Usage:\n"
			"  modcalc [\"expression\"/\"filename\"]\n"
			"Please enclose command arguments in doublequotes.\n"
			"Press H for help, X to exit, or any key to continue.\n");
		char c=_getch();
		if((c&0xDF)=='H')
			print_help();
		else if((c&0xDF)=='X')
			return 0;
		get_str_interactive(str, 0);
	}
	else if(argc==2)
	{
		printf("> %s\n\n", argv[1]);
		if(file_is_readablea(argv[1]))
		{
#ifdef __linux__
			FILE *file=fopen(argv[1], "r");
			if(!file)
			{
				printf("%s\n", strerror(errno));
				return EXIT_FAILURE;
			}
#else
			FILE *file=nullptr;
			int ret=fopen_s(&file, argv[1], "r");
			if(ret)
			{
				strerror_s(g_buf, g_buf_size, ret);
				printf("%s\n", g_buf);
				return EXIT_FAILURE;
			}
#endif
			fseek(file, 0, SEEK_END);
			int bytesize=ftell(file);
			fseek(file, 0, SEEK_SET);
			str.resize(bytesize);
			fread(&str[0], 1, bytesize, file);
			fclose(file);
			str.resize(strlen(str.c_str()));
		}
		else
			str=argv[1];
		quit_prompt=true;
	}
	else
	{
		printf(
			"Enter \'help\' for documentation.\n"
			"\n"
			);
		//print_help();//
		get_str_interactive(str, 0);
	}

	for(int result=SOLVE_OK;;)
	{
		while((result=solve(str, result!=SOLVE_OK))==SOLVE_INCOMPLETE)
		{
			std::string str2;
			get_str_interactive(str2, "...");
			str+=str2;
		}
		if(result==SOLVE_PARSE_ERROR)
		{
			printf("\n");
			for(int k=0;k<(int)errors.size();++k)
				printf("%s\n", errors[k].c_str());
			printf("\n");
			errors.clear();
		}
		else if(result==SOLVE_OK_NO_ANS)//success, but don't print answer
		{
			//if(g_answers.size())//clear command clears answers
			//	g_answers.pop_back();
		}
		else//success, print answer
		{
			auto &ans=g_answers.back();
			if(ans.name)
				printf("%s =\n", ans.name);
			else
				printf("ans(%d) =\n", (int)g_answers.size()-1);
			ans.print();
		}
		//if(quit_prompt)
		//	break;
		if(quit_prompt)//
		{
			printf("Quit? [Y/N] ");

			char c=0;
			scanf("%c", &c);
			while(getchar()!='\n');

			//char c=_getche();

			if((c&0xDF)=='Y')
				break;
			quit_prompt=false;
		}

		get_str_interactive(str, 0);
	}
	strings.clear();
	return 0;
}