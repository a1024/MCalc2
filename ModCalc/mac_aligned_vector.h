#include<stdlib.h>
#ifdef __linux__
#include<stddef.h>
#include"mc2_memory.h"
#endif
template<typename T, int align>struct AlignedAllocator
{
	typedef AlignedAllocator<T, align> other;

	typedef T value_type;

	typedef value_type *pointer;
	typedef value_type const *const_pointer;
	typedef void *void_pointer;
	typedef const void *const_void_pointer;

	typedef value_type& reference;
	typedef const value_type& const_reference;

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	//typedef false_type propagate_on_container_copy_assignment;
	//typedef false_type propagate_on_container_move_assignment;
	//typedef false_type propagate_on_container_swap;
	AlignedAllocator<T, align> select_on_container_copy_construction()const
	{	// return this allocator
		return *this;
	}

	template<class _Other>struct rebind
	{	// convert this type to allocator<_Other>
		typedef AlignedAllocator<_Other, align> other;
	};

	pointer address(reference _Val) const throw()
	{	// return address of mutable _Val
		return &_Val;
	}

	const_pointer address(const_reference _Val) const throw()
	{	// return address of nonmutable _Val
		return &_Val;
	}

	AlignedAllocator()throw(){}// construct default allocator (do nothing)
	AlignedAllocator(const AlignedAllocator<T, align>&)throw(){}// construct by copying (do nothing)
	template<class _Other>AlignedAllocator(const AlignedAllocator<_Other, align>&)throw(){}// construct from a related allocator (do nothing)

	template<class _Other>AlignedAllocator<T, align>& operator=(const AlignedAllocator<_Other, align>&)
	{	// assign from a related allocator (do nothing)
		return *this;
	}

	void deallocate(pointer _Ptr, size_type)
	{	// deallocate object at _Ptr, ignore size
		_aligned_free(_Ptr);
	}

	pointer allocate(size_type _Count)
	{	// allocate array of _Count elements
		return (T*)_aligned_malloc(_Count*sizeof(T), align);
	}
	pointer allocate(size_type _Count, const void *p)
	{	// allocate array of _Count elements, ignore hint
		return (T*)_aligned_realloc(p, _Count*sizeof(T), align);
	}

	void construct(T *_Ptr)
	{	// default construct object at _Ptr
		::new ((void *)_Ptr) T();
	}
	void construct(T *_Ptr, const T& _Val)
	{	// construct object at _Ptr with value _Val
		::new ((void *)_Ptr) T(_Val);
	}
	//template<class _Objty, class... _Types>void construct(_Objty *_Ptr, _Types&&... _Args)
	//{	// construct _Objty(_Types...) at _Ptr
	//	::new ((void *)_Ptr) _Objty(std::forward<_Types>(_Args)...);
	//}


	template<class _Uty>void destroy(_Uty *_Ptr)
	{	// destroy object at _Ptr
		_Ptr->~_Uty();
	}

	size_t max_size() const throw()
	{	// estimate maximum array size
		return ((size_t)(-1) / sizeof (T));
	}
};
#if 0
//https://stackoverflow.com/questions/12942548/making-stdvector-allocate-aligned-memory
enum class Alignment:size_t
{
    Normal = 8,
    SSE    = 16,
    AVX    = 32,
};
#ifdef _MSC_VER
#define	noexcept
#endif
template<typename T, Alignment Align=Alignment::AVX>class AlignedAllocator;
template<Alignment Align>class AlignedAllocator<void, Align>
{
public:
    typedef void*             pointer;
    typedef const void*       const_pointer;
    typedef void              value_type;

    template <class U> struct rebind { typedef AlignedAllocator<U, Align> other; };
};
template<typename T, Alignment Align>class AlignedAllocator
{
public:
    typedef T         value_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;
    typedef T&        reference;
    typedef const T&  const_reference;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef std::true_type propagate_on_container_move_assignment;

    template<class U>struct rebind { typedef AlignedAllocator<U, Align> other; };

public:
    AlignedAllocator()noexcept{}
    template<class U>AlignedAllocator(const AlignedAllocator<U, Align>&)noexcept{}

    size_type max_size()const noexcept
    {
		return (size_type(~0) - size_type(Align))/sizeof(T);
	}

    pointer address(reference x) const noexcept
    {
		return std::addressof(x);
	}
    const_pointer address(const_reference x) const noexcept
    {
		return std::addressof(x);
	}
    pointer allocate(size_type n, typename AlignedAllocator<void, Align>::const_pointer = 0)
    {
        const size_type alignment = static_cast<size_type>( Align );
        void* ptr=_aligned_malloc(n * sizeof(T), alignment);
        if(!ptr)
            throw std::bad_alloc();
        return reinterpret_cast<pointer>(ptr);
    }
    void deallocate(pointer p, size_type)noexcept
    {
		return _aligned_free(p);
	}
    template <class U, class ...Args>void construct(U* p, Args&&... args)
    {
		::new(p) U(std::forward<Args>(args)...);
	}
    void destroy(pointer p)
    {
		p->~T();
	}
};
template<typename T, Alignment Align>class AlignedAllocator<const T, Align>
{
public:
    typedef T         value_type;
    typedef const T*  pointer;
    typedef const T*  const_pointer;
    typedef const T&  reference;
    typedef const T&  const_reference;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef std::true_type propagate_on_container_move_assignment;

    template <class U>
    struct rebind { typedef AlignedAllocator<U, Align> other; };

public:
    AlignedAllocator()noexcept{}
    template <class U>AlignedAllocator(const AlignedAllocator<U, Align>&)noexcept{}

    size_type max_size()const noexcept
    {
		return (size_type(~0) - size_type(Align)) / sizeof(T);
	}
    const_pointer address(const_reference x) const noexcept
    {
		return std::addressof(x);
	}
    pointer allocate(size_type n, typename AlignedAllocator<void, Align>::const_pointer = 0)
    {
        const size_type alignment = static_cast<size_type>( Align );
        void* ptr = _aligned_malloc(n * sizeof(T), alignment);
        if (ptr == nullptr)
            throw std::bad_alloc();
        return reinterpret_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type)noexcept
    {
		return _aligned_free(p);
	}
    template <class U, class ...Args>void construct(U* p, Args&&... args)
    {
		::new(p) U(std::forward<Args>(args)...);
	}
    void destroy(pointer p)
    {
		p->~T();
	}
};
template <typename T, Alignment TAlign, typename U, Alignment UAlign>inline bool operator==(const AlignedAllocator<T,TAlign>&, const AlignedAllocator<U, UAlign>&)noexcept
{
	return TAlign == UAlign;
}
template <typename T, Alignment TAlign, typename U, Alignment UAlign>inline bool operator!=(const AlignedAllocator<T,TAlign>&, const AlignedAllocator<U, UAlign>&) noexcept
{
	return TAlign != UAlign;
}
#endif