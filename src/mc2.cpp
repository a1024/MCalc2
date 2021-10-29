//mc2.cpp - main file
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
#include	<stdio.h>
#include	<conio.h>
#include	<string>
#include	<iostream>
extern const int	g_buf_size=G_BUF_SIZE;
char				g_buf[G_BUF_SIZE]={0};
char		gfmode=0;//later
void		print_help()
{
	printf(
	//	"01234567890123456789012345678901234567890123456789012345678901234567890123456789"
		"MCALC2: A Matrix Calculator\n"
		"Usage:\n"
		"  mcalc2 [\"expression\"/filename]\n"
		"\n"
		"Syntax:\n"
		"  Use Matlab-style vectors & matrices.\n"
		"  Elements are separated by commas or spaces, for example:\n"
		"    [a2 a1 a0] or [a2, a1, a0] is a row vector\n"
		"    [a11, a12;  a21, a22]      is a 2x2 matrix\n"
		"  Use same square brackets to access matrix elements, for example:\n"
		"    [1 2; 3 4][0][1] evaluate to 2\n"
	//	"  For multiline statements:\n"
	//	"    Open a parenthesis at the start of the first command, and close it at the\n"
	//	"    end of the last command.\n"
		"  Nested matrices and/or polynomials are not supported yet.\n"
		"  Variable names can be up to 16 characters.\n"//
		"  Multiplication asterix \'*\' is always obligatory.\n"
		"\n"
		"Operators by precedence (first to last):\n"
		"  ^\n"
		"  + - (unary)\n"
		"  \'\n"
		"  / \\ % * o ./ .* - +\n"
		"  = += -= *= /= \\= %%=\n"
		"  ,\n"
		"  ; [ ] ( )\n"
		"Notes:\n"
		"  \'     is transpose (eg: M\')\n"
		"  o     is the tensor product\n"
		"  \\     is for square matrices (eg: A\\B)\n"
	//	"  %%     makes fraction objects (eg: F=num%%den)\n"//X  it's the modulus operator
		"  .* ./ are element-wise operations\n"
		"\n"
		"Keywords:\n"
		"  help: Print this info (works in cmd)\n"
		"  clear: Clears all variables from memory\n"
		"  vars: Shows all variables & answer count\n"
		"  open: Choose a text file to open\n"
		"General functions:\n"
		"  ans(n): The n-th latest answer\n"
		"  cmd(w, h): set console buffer size (in characters)\n"
		"Vectors:\n"
		"  cross(3D vec, 3D vec) -> 3D vec\n"
		"  cross(2D vec, 2D vec) -> scalar\n"
		"  For dot product use transpose (col\'*col or row*row\')\n"
		"Matrices:\n"
		"  ref: Row Echelon Form\n"
		"  rref: Reduced Row Echelon Form (Gaussian Elimination)\n"
		"Square matrices:\n"
		"  I(n): n-square identity matrix\n"//autosize?
		"  det: determinant of square matrix\n"
		"  inv: inverse of square matrix\n"
		"  tr: trace of square matrix\n"
	//	"  diag: diagonal factorization of square matrix\n"//TODO
	//	"  lu: LU factorization of square matrix\n"
	//	"Polynomials:\n"
	//	"  roots: find the roots\n"
	//	"Fraction objects:\n"
	//	"  sample(F): s to z domain\n"
	//	"  ldiv(F, S): long division of fraction F by S steps\n"
	//	"  plot(F, S): long division\n"
	//	"Matrices & polynomials:\n"
	//	"  dft/fft: Discrete Fourier Transform\n"
	//	"  idft/ifft: Inverse Discrete Fourier Transform\n"
		"\n"
		);
}
void		get_str_interactive(std::string &str, const char *cmdstr)
{
	if(gfmode)
		printf("gf");
	if(cmdstr)
		printf(cmdstr);
	printf("> ");
	std::getline(std::cin, str);
}
bool		get_str_from_file(std::string &str)
{
	auto wbuf=open_file_window();
	if(!wbuf)
		return false;
	FILE *file;
	int ret=_wfopen_s(&file, wbuf, L"r");
	if(ret)
	{
		strerror_s(g_buf, g_buf_size, ret);
		printf("%s\n", g_buf);
		return false;
	}
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
	set_console_buffer_size(80, 9001);
	printf("MCALC%s\n\n", argc==2?"":"\t\tCtrl C to exit.");

	print_help();//

	std::string str;
	if(argc>2)
	{
		printf(
			"Usage:\n"
			"  mcalc [\"expression\"/filename]\n"
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
		str=argv[1];
		printf("> %s\n\n", argv[1]);
	}
	else
		get_str_interactive(str, 0);

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
				printf("ans(%d) =\n", g_answers.size()-1);
			ans.print();
		}

		get_str_interactive(str, 0);
	}
}