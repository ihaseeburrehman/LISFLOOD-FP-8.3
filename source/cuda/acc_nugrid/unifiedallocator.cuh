#pragma once

template<typename T>
class UnifiedAllocator
{
public:
	T* allocate(std::size_t n);
	void deallocate(T* ptr, std::size_t unused);
};


