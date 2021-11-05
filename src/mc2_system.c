//mc2_system.c - OS-dependent code

#include<sys/stat.h>
#ifdef __linux__
#include <termios.h>//https://stackoverflow.com/questions/7469139/what-is-the-equivalent-to-getch-getche-in-linux
static struct termios	old, current;
void	initTermios(int echo)//Initialize new terminal i/o settings
{
	tcgetattr(0, &old);//grab old terminal i/o settings
	current=old;//make new settings same as old settings
	current.c_lflag&=~ICANON;//disable buffered i/o
	if(echo)
		current.c_lflag|=ECHO;//set echo mode
	else
		current.c_lflag&=~ECHO;//set no echo mode
	tcsetattr(0, TCSANOW, &current);//use these new terminal i/o settings now
}
void	resetTermios()//Restore old terminal i/o settings
{
	tcsetattr(0, TCSANOW, &old);
}
char	getch_sel(int echo)//Read 1 character - echo defines echo mode
{
	char ch;
	initTermios(echo);
	ch=getchar();
	resetTermios();
	return ch;
}
char	_getch()//Read 1 character without echo
{
	return getch_sel(0);
}
char	_getche()//Read 1 character with echo
{
	return getch_sel(1);
}
#else
#include	<Windows.h>
#include	<stdio.h>
CONSOLE_SCREEN_BUFFER_INFO csbi;
void		get_console_size(short *w, short *h)
{
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
	*w=csbi.dwSize.X;
	*h=csbi.dwSize.Y;
}
int			set_console_buffer_size(short w, short h)
{
	COORD coords={w, h};
	HANDLE handle=GetStdHandle(STD_OUTPUT_HANDLE);
	int success=SetConsoleScreenBufferSize(handle, coords);
	if(!success)
		printf("Failed to resize console buffer: %d\n\n", GetLastError());
	return success;
}
const int	g_wbuf_size=4096;
wchar_t		g_wbuf[4096]={0};
const wchar_t* open_file_window()
{
	g_wbuf[0]=0;
	HWND hWnd=GetConsoleWindow();
	OPENFILENAMEW ofn=
	{
		sizeof(OPENFILENAMEW),
		hWnd, (HINSTANCE)GetWindowLongW(hWnd, GWL_HINSTANCE),
		L"Text Filex (*.txt)\0*.TXT\0", 0, 0, 1,
		g_wbuf, g_wbuf_size,
		0, 0,//initial filename
		0,
		0,//dialog title
		OFN_CREATEPROMPT|OFN_PATHMUSTEXIST,
		0,//file offset
		0,//extension offset
		L"TXT",//default extension
		0, 0,//data & hook
		0,//template name
		0,//reserved
		0,//reserved
		0,//flags ex
	};
	int success=GetOpenFileNameW(&ofn);
	if(!success)
		return 0;
	memcpy(g_wbuf, ofn.lpstrFile, wcslen(ofn.lpstrFile)*sizeof(wchar_t));
	return g_wbuf;
}
#define	S_ISREG(m)	(((m)&S_IFMT)==S_IFREG)
#endif
int				file_is_readablea(const char *filename)//0: not readable, 1: regular file, 2: folder
{
	struct stat info;
	int error=stat(filename, &info);
	if(!error)
		return 1+!S_ISREG(info.st_mode);
	return 0;
}
int				file_is_readablew(const wchar_t *filename)//0: not readable, 1: regular file, 2: folder
{
	struct _stat32 info;
	int error=_wstat32(filename, &info);
	if(!error)
		return 1+!S_ISREG(info.st_mode);
	return 0;
}