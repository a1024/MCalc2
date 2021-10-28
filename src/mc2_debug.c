//mc2_debug.c - Pointer tracking for debugging
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

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<conio.h>
static void	mem_print(const void *buf, int bytesize)
{
	int k;
	char *pb=(char*)buf;
	if(bytesize)
	{
		printf("%02X", (int)pb[0]&0xFF);
		for(k=1;k<bytesize;++k)
			printf("-%02X", (int)pb[k]&0xFF);
	}
	printf("\n");
}
static void	print_byte_bin(char b)
{
	int k;
	for(k=7;k>=0;--k)
		printf("%d", b>>k&1);
}
static void	mem_print_bin(const void *buf, int bytesize)
{
	int k;
	char *pb=(char*)buf;
	if(bytesize)
	{
		print_byte_bin(pb[0]);
		for(k=1;k<bytesize;++k)
		{
			printf("-");
			print_byte_bin(pb[k]);
			if(!((k+1)%10))
				printf("\n");
		}
	}
}

#define		POINTERS_MAX	0x10000
typedef enum
{
	POINTER_FREED,
	POINTER_RELOCATED,
	POINTER_ACTIVE,
} PointerStatus;
static const char* pointer_status2str(PointerStatus status)
{
	switch(status)
	{
#define	CASE(TYPE)	case TYPE:return #TYPE;
		CASE(POINTER_FREED)
		CASE(POINTER_RELOCATED)
		CASE(POINTER_ACTIVE)
#undef	CASE
	}
	return "<UNDEFINED>";
}
typedef struct
{
	union
	{
		const void *p;
		const char *bytes;
	};
	size_t bytesize;
	PointerStatus status;
} Pointer;
Pointer	pointers[POINTERS_MAX]={0};
int		npointers=0;
static int	pointers_find(const void *p)
{
	int k;

	for(k=0;k<npointers;++k)
		if(pointers[k].p==p)
			return k;
	return npointers;
	//return -1;
}
static void	pointers_action_alloc(const void *p_new, size_t bytesize)
{
	Pointer *pointer;
	int idx=pointers_find(p_new);

	if(idx>=POINTERS_MAX)
	{
		printf("\n\nPOINTERS.ALLOC: Need to track over %d pointers\n\n", POINTERS_MAX);
		_getch();
		abort();
	}

	pointer=pointers+idx;
	pointer->p=p_new;
	pointer->bytesize=bytesize;
	pointer->status=POINTER_ACTIVE;

	npointers+=idx==npointers;
}
static void	pointers_action_realloc(const void *p_old, const void *p_new, size_t bytesize_new)
{
	Pointer *pointer;
	int idx1=pointers_find(p_old), idx2=pointers_find(p_new);
	if(idx1==npointers)
	{
		printf("\n\nPOINTERS.REALLOC: p_old=0x%p not found, p_new=0x%p\n\n", p_old, p_new);
		_getch();
		abort();
	}
	if(p_new==p_old)
		pointers[idx1].bytesize=bytesize_new;
	else
	{
		if(idx2<npointers&&pointers[idx2].status==POINTER_ACTIVE)
		{
			printf("\n\nPOINTERS.REALLOC: p_old=0x%p, p_new=0x%p was active\n\n", p_old, p_new);
			_getch();
			abort();
		}
		if(!p_new)
		{
			printf("\n\nPOINTERS.REALLOC: p_old=0x%p, p_new=0x%p failed\n\n", p_old, p_new);
			_getch();
			abort();
		}
		pointer=pointers+idx1;
		//pointer->bytesize=0;
		pointer->status=POINTER_RELOCATED;

		pointer=pointers+idx2;
		pointer->p=p_new;
		pointer->bytesize=bytesize_new;
		pointer->status=POINTER_ACTIVE;

		npointers+=idx2==npointers;
	}
}
static int	pointers_action_free(const void *p_old)
{
	int idx=pointers_find(p_old);
	if(idx==npointers)
	{
		printf("\n\nPOINTERS.FREE: p_old=0x%p not found\n\n", p_old);
		return -1;
		//_getch();
		//abort();
	}
	if(pointers[idx].status!=POINTER_ACTIVE)
	{
		printf("\n\nPOINTERS.FREE: p_old=0x%p status=%s\n\n", p_old, pointer_status2str(pointers[idx].status));
		return -1;
		//_getch();
		//abort();
	}
	//pointers[idx].bytesize=0;
	pointers[idx].status=POINTER_FREED;
	return 0;
}
static void	pointers_action_access(const void *p)
{
	Pointer *pointer;
	int k, found=0;
	static int access_count=0;
	
	printf("\nPOINTERS.ACCESS: 0x%p\n", p);
	for(k=0;k<npointers;++k)//find in range
	{
		pointer=pointers+k;
		if((char*)p>=pointer->bytes&&(char*)p<pointer->bytes+pointer->bytesize)//found
		{
			found=1;
			printf("\tp=0x%p, %d bytes status=%s\n", pointer->p, pointer->bytesize, pointer_status2str(pointer->status));
			if(pointer->status!=POINTER_ACTIVE)
				_getch();
		}
	}
	if(!found)
		printf("\tnot found\n", p);
	++access_count;
}

const char* get_filename(const char *filename)
{
	int k;
	for(k=strlen(filename)-1;k>=0&&filename[k]!='/'&&filename[k]!='\\';--k);
	return filename+k+1;
}

int			syscall_count=0, emergency_flag=0;
void*		d_alloc		(const char *file, int line, unsigned long bytesize)
{
	void *p;
	++syscall_count;
	printf("%s(%d): #%d malloc %ld", get_filename(file), line, syscall_count, bytesize);
	p=malloc(bytesize);
	printf(" -> 0x%p\n", p);

	pointers_action_alloc(p, bytesize);
	//pointers_add(p);//

	return p;
}
void*		d_realloc	(const char *file, int line, void *p, unsigned long bytesize)
{
	void *p2;
	++syscall_count;
	printf("%s(%d): #%d realloc %ld, 0x%p", get_filename(file), line, syscall_count, bytesize, p);
	p2=realloc(p, bytesize);
	printf(" -> 0x%p\n", p2);

	pointers_action_realloc(p, p2, bytesize);
	//pointers_replace(p, p2);//

	if(p2)
		return p2;
	return p;
}
int			d_free		(const char *file, int line, void *p)
{
	int status=0;

	++syscall_count;
	printf("%s(%d): #%d free 0x%p", get_filename(file), line, syscall_count, p);

	status=pointers_action_free(p);
	if(status<0)
	{
		printf("Expecting a crash at free()\n");
		_getch();
	}

	if(p)
		free(p);

	printf("\n");

	return status;
}
void		d_memset	(const char *file, int line, void *dst, int val, int bytesize)
{
	++syscall_count;
	printf("%s(%d): #%d memset 0x%p := %d, %d before:\n\t", get_filename(file), line, syscall_count, dst, val, bytesize);

	pointers_action_access(dst);//

	mem_print(dst, bytesize);

	memset(dst, val, bytesize);

	printf("\t...after:\n\t");
	mem_print(dst, bytesize);
}
void		d_memcpy	(const char *file, int line, void *dst, const void *src, int bytesize)
{
	++syscall_count;
	printf("%s(%d): #%d memcpy 0x%p := 0x%p, %d before:\n\t", get_filename(file), line, syscall_count, dst, src, bytesize);

	pointers_action_access(src);//
	pointers_action_access(dst);//

	mem_print(dst, bytesize);

	memcpy(dst, src, bytesize);

	printf("\t...after:\n\t");
	mem_print(dst, bytesize);
}
void		d_memmove	(const char *file, int line, void *dst, const void *src, int bytesize)
{
	++syscall_count;
	printf("%s(%d): #%d memmove 0x%p := 0x%p, %d before:\n\t", get_filename(file), line, syscall_count, dst, src, bytesize);

	pointers_action_access(src);//
	pointers_action_access(dst);//

	mem_print(dst, bytesize);

	memmove(dst, src, bytesize);

	printf("\t...after:\n\t");
	mem_print(dst, bytesize);
}