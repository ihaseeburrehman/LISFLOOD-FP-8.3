#pragma once

#include "lisflood.h"

void *memory_allocate(size_t size);

NUMERIC_TYPE*memory_allocate_zero_numeric_legacy(size_t size);
NUMERIC_TYPE*memory_allocate_numeric_legacy(size_t size);
void memory_free_legacy(void** memory);
void memory_free_legacy(int** memory);
void memory_free_legacy(NUMERIC_TYPE** memory);

void memory_free(void** memory);
void memory_free(int** memory);
void memory_free(NUMERIC_TYPE** memory);


void memory_free(void** memory, size_t size);

char *trimwhitespace(char *str);

void SetArrayValue(int* arr, int value, int length);