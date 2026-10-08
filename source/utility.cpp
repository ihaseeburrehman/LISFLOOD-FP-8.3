#include "utility.h"

static long total_allocated = 0;
static long total_legacy_allocated = 0;

NUMERIC_TYPE* memory_allocate_zero_numeric_legacy(size_t size)
{
	NUMERIC_TYPE* memory = new NUMERIC_TYPE[size]();

	if (memory == NULL)
	{
		printf("memory allocation failed %ldMB, legacy %ld, total %ld\n", total_allocated / 1024 / 1024, total_legacy_allocated / 1024 / 1024, (total_allocated + total_legacy_allocated) / 1024 / 1024);
		exit(-1);
	}
	total_legacy_allocated += size * sizeof(NUMERIC_TYPE);
	return memory;
}

NUMERIC_TYPE* memory_allocate_numeric_legacy(size_t size)
{
	NUMERIC_TYPE* memory = new NUMERIC_TYPE[size];

	if (memory == NULL)
	{
		printf("memory allocation failed %ldMB, legacy %ld, total %ld\n", total_allocated / 1024 / 1024, total_legacy_allocated / 1024 / 1024, (total_allocated + total_legacy_allocated) / 1024 / 1024);
		exit(-1);
	}
	total_legacy_allocated += size * sizeof(NUMERIC_TYPE);
	return memory;
}

void memory_free_legacy(int** memory)
{
	if (*memory != NULL)
	{
		total_legacy_allocated -= sizeof *memory;

		delete[] * memory;
		*memory = NULL;
	}
}

void memory_free_legacy(NUMERIC_TYPE** memory)
{
	if (*memory != NULL)
	{
		total_legacy_allocated -= sizeof *memory;

		delete[] * memory;
		*memory = NULL;
	}
}

void* memory_allocate(size_t size)
{
#if defined(_MSC_VER) || defined (__INTEL_COMPILER)
	void* memory = _mm_malloc(size, 64);
#else
	void* memory = NULL;
	posix_memalign(&memory, 64, size);
#endif
	if (memory == NULL)
	{
		printf("memory allocation failed %ldMB, legacy %ld, total %ld\n", total_allocated / 1024 / 1024, total_legacy_allocated / 1024 / 1024, (total_allocated + total_legacy_allocated) / 1024 / 1024);
		exit(-1);
	}
	total_allocated += size;
	return memory;
}

void memory_free(int** memory)
{
	memory_free((void**)memory);
}

void memory_free(NUMERIC_TYPE** memory)
{
	memory_free((void**)memory);
}

void memory_free(void** memory)
{
#if defined(_MSC_VER) || defined (__INTEL_COMPILER)
	_mm_free(*memory);
#else
	posix_memalign((void**)&memory, sizeof(float) * 100* 128, 64);
#endif
	*memory = NULL;
}

void memory_free(void** memory, size_t size)
{
	memory_free(memory);
	total_allocated += size;
}



// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char *trimwhitespace(char *str)
{
	char *end;

	// Trim leading space
	while (isspace(*str)) str++;

	if (*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while (end > str && isspace(*end)) end--;

	// Write new null terminator
	*(end + 1) = 0;

	return str;
}

void SetArrayValue(int* arr, int value, int length)
{
	for (int j = 0; j < length; j++)
	{
		arr[j] = value;
	}
}
