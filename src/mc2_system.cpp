//mc2_system.c - OS-dependent code

#include<sys/stat.h>
#ifdef __linux__
#include	<gtk/gtk.h>
#include	<string>
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
		L"Text Files (*.txt)\0*.TXT\0", 0, 0, 1,
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
#ifndef __linux__
int				file_is_readablew(const wchar_t *filename)//0: not readable, 1: regular file, 2: folder
{
	struct _stat32 info;
	int error=_wstat32(filename, &info);
	if(!error)
		return 1+!S_ISREG(info.st_mode);
	return 0;
}
#endif

#ifdef __linux__
//https://stackoverflow.com/questions/6145910/cross-platform-native-open-save-file-dialogs
//https://github.com/AndrewBelt/osdialog
struct				osdialog_filters
{
	const char *name, **patterns;
	int npatterns;
};
bool				osdialog_file(bool open, const char *default_path, const char *filename, osdialog_filters *filters, int nfilters, std::string &result)
{
	if(!gtk_init_check(NULL, NULL))
		return 0;

	GtkFileChooserAction gtkAction;
	const char* title;
	const char* acceptText;
	//if(action == OSDIALOG_OPEN)
	if(open)
	{
		title = "Open File";
		acceptText = "Open";
		gtkAction = GTK_FILE_CHOOSER_ACTION_OPEN;
	}
	//else if (action == OSDIALOG_OPEN_DIR)
	//{
	//	title = "Open Folder";
	//	acceptText = "Open Folder";
	//	gtkAction = GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER;
	//}
	else
	{
		title = "Save File";
		acceptText = "Save";
		gtkAction = GTK_FILE_CHOOSER_ACTION_SAVE;
	}

	GtkWidget* dialog = gtk_file_chooser_dialog_new(title, NULL, gtkAction,
		"_Cancel", GTK_RESPONSE_CANCEL,
		acceptText, GTK_RESPONSE_ACCEPT, NULL);

	for(int k=0;k<nfilters;++k)
	{
		auto filter=filters+k;
		auto fileFilter=gtk_file_filter_new();
		gtk_file_filter_set_name(fileFilter, filter->name);
		for(int k2=0;k2<filter->npatterns;++k2)
			gtk_file_filter_add_pattern(fileFilter, filter->patterns[k2]);
		gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), fileFilter);
	}

	if(!open)
		gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);

	if(default_path)
		gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), default_path);

	if(!open&&filename)
		gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), filename);

	bool ret=false;
	char *chosen_filename=nullptr;
	if(ret=gtk_dialog_run(GTK_DIALOG(dialog))==GTK_RESPONSE_ACCEPT)
		chosen_filename=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
	gtk_widget_destroy(dialog);

	if (chosen_filename)
		result=chosen_filename;

	while(gtk_events_pending())
		gtk_main_iteration();
	return ret;
}
#define			G_BUF_SIZE		1024
static char		g_buf[G_BUF_SIZE]={0};
const char*		open_file_window()
{
	const char *all_files[]={"*.*"};
	osdialog_filters filters[]=
	{
		{"All Files (*.*)", all_files, 1},
	};
	std::string result;
	if(!osdialog_file(true, nullptr, nullptr, filters, 1, result))
		return nullptr;
	memcpy(g_buf, result.c_str(), result.size()+1);
	return g_buf;
}
#endif