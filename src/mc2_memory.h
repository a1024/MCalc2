//mc2_memory.c - Memory operations
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

#ifndef MC2_CSTUFF
#define MC2_CSTUFF
#include<string.h>


//	#define	DEBUG_MEMORY


#ifdef __linux__
#include<stdlib.h>
#include<assert.h>
//#define		_aligned_malloc(SIZE, ALIGN)	aligned_alloc(ALIGN, SIZE)	//how to realloc?
//#define		_aligned_free(POINTER)			free(POINTER)
#define			GETPTR(PTR)		((void**)(PTR))[-1]
#define			GETSIZE(PTR)	((size_t*)(PTR))[-2]
static void*	_aligned_malloc(size_t size, size_t align)
{
	size_t addr=0, aa=0;
	assert(align>=sizeof(size_t));//
	addr=(size_t)malloc(size+align+2*sizeof(size_t));
	if(!addr)
		return 0;
	aa=addr+2*sizeof(size_t)+align;
	aa-=aa%align;
	GETPTR(aa)=(void*)addr;
	GETSIZE(aa)=size;
	return (void*)aa;
}
static void 	_aligned_free(void *p)
{
	free(GETPTR(p));
}
static void*	_aligned_realloc(void *oldp, size_t size, size_t align)//inherently slow
{
	int copysize=GETSIZE(oldp);
	void *p=_aligned_malloc(size, align);
	if(!p)
		return 0;
	if(copysize>size)
		copysize=size;
	memcpy(p, oldp, copysize);
	_aligned_free(oldp);
	return p;
}
#endif
#ifdef _MSC_VER
#define		scanf	scanf_s
#endif
#ifdef DEBUG_MEMORY
#include	<conio.h>
#define		malloc(SIZE)							d_alloc(file, __LINE__, SIZE)
#define		realloc(POINTER, SIZE)					d_realloc(file, __LINE__, POINTER, SIZE)
#define		free(POINTER)							d_free(file, __LINE__, POINTER)

#define		_aligned_malloc(SIZE, ALIGN)			d_aligned_alloc(file, __LINE__, SIZE, ALIGN)
#define		_aligned_realloc(POINTER, SIZE, ALIGN)	d_aligned_realloc(file, __LINE__, POINTER, SIZE, ALIGN)
#define		_aligned_free(POINTER)					d_aligned_free(file, __LINE__, POINTER)

#define		memcpy(DST, SRC, SIZE)					d_memcpy(file, __LINE__, DST, SRC, SIZE)
#define		memmove(DST, SRC, SIZE)					d_memmove(file, __LINE__, DST, SRC, SIZE)
#define		memset(DST, VAL, SIZE)					d_memset(file, __LINE__, DST, VAL, SIZE)
#ifdef __cplusplus
extern "C"
{
#endif
	extern int	syscall_count, emergency_flag;
	void*	d_alloc				(const char *file, int line, size_t bytesize);
	void*	d_realloc			(const char *file, int line, void *p, size_t bytesize);
	int		d_free				(const char *file, int line, void *p);

	void*	d_aligned_alloc		(const char *file, int line, size_t bytesize, size_t alignment);
	void*	d_aligned_realloc	(const char *file, int line, void *p, size_t bytesize, size_t alignment);
	int		d_aligned_free		(const char *file, int line, void *p);

	void	d_memset			(const char *file, int line, void *dst, int val, size_t bytesize);
	void	d_memcpy			(const char *file, int line, void *dst, const void *src, size_t bytesize);
	void	d_memmove			(const char *file, int line, void *dst, const void *src, size_t bytesize);
#ifdef __cplusplus
}
#endif
#endif
static void	memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
#ifdef DEBUG_MEMORY
	const char file[]=__FILE__;
#endif
	unsigned copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, copied);
	while(copied<<1<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#define		ALLOC(TYPE, NEWP, COUNT)					NEWP=(TYPE*)malloc((COUNT)*sizeof(TYPE))
#define		REALLOC(TYPE, NEWP, OLDP, COUNT)			NEWP=(TYPE*)realloc(OLDP, (COUNT)*sizeof(TYPE))
#define		FREE(OLDP)									free(OLDP)
#define		MEMFILL(TYPE, DST, SRC, DSTCOUNT, SRCCOUNT)	memfill(DST, SRC, (DSTCOUNT)*sizeof(TYPE), (SRCCOUNT)*sizeof(TYPE))
#define		MEMZERO(TYPE, DST, COUNT)					memset(DST, 0, (COUNT)*sizeof(TYPE))
#define		MEMCPY(TYPE, DST, SRC, COUNT)				memcpy(DST, SRC, (COUNT)*sizeof(TYPE))
#define		MEMMOVE(TYPE, DST, SRC, COUNT)				memmove(DST, SRC, (COUNT)*sizeof(TYPE))

#define		DALLOC(NEWP, COUNT)							ALLOC(double, NEWP, COUNT)
#define		DREALLOC(NEWP, OLDP, COUNT)					REALLOC(double, NEWP, OLDP, COUNT)
#define		DFREE(OLDP)									FREE(OLDP)
#define		DMEMFILL(DST, SRC, DSTCOUNT, SRCCOUNT)		MEMFILL(double, DST, SRC, DSTCOUNT, SRCCOUNT)
#define		DMEMZERO(DST, COUNT)						MEMZERO(double, DST, COUNT)
#define		DMEMCPY(DST, SRC, COUNT)					MEMCPY(double, DST, SRC, COUNT)
#define		DMEMMOVE(DST, SRC, COUNT)					MEMMOVE(double, DST, SRC, COUNT)

//for complex data
//typedef		double Comp[2];
typedef struct _Comp
{
	double r, i;
} Comp;
#define		CALLOC(NEWP, COUNT)							NEWP=(Comp*)_aligned_malloc((COUNT)*sizeof(Comp), 16)
#define		CREALLOC(NEWP, OLDP, COUNT)					NEWP=(Comp*)_aligned_realloc(OLDP, (COUNT)*sizeof(Comp), 16)
#define		CFREE(OLDP)									_aligned_free(OLDP)
#define		CMEMFILL(DST, SRC, DSTCOUNT, SRCCOUNT)		memfill(DST, SRD, (DSTCOUNT)*sizeof(Comp), (SRCCOUNT)*sizeof(Comp))
#define		CMEMZERO(DST, COUNT)						memset(DST, 0, (COUNT)*sizeof(Comp))
#define		CMEMCPY(DST, SRC, COUNT)					memcpy(DST, SRC, (COUNT)*sizeof(Comp))
#define		CMEMMOVE(DST, SRC, COUNT)					memmove(DST, SRC, (COUNT)*sizeof(Comp))
#endif