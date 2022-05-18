/*
 * File:   ffht.c
 * Author: Canales-Martinez, Isaac Andres
 */

#include <string.h>

#include "ffht.h"

#if defined __AVX__
#include "ffht_avx.h"
#else
#include "ffht_plain.h"
#endif /* __AVX__ */

#define max(a, b)    ((a) > (b) ? (a) : (b))

void ffht(double *v, size_t n)
{
	ffht_part(v, n, 0);
}

void ffht_part(double *v, size_t n, size_t step)
{
	ffht_atom x, y;
	double aux[FFHT_BASE_ELEM];
	size_t i, j;

	if (n < FFHT_BASE_ELEM) {
		if (n == 0)
			return;
		memset(aux, 0, FFHT_BASE_ELEM * sizeof(double));
		memcpy(aux, v, n * sizeof(double));
		ffht_compute_base(aux, FFHT_BASE_ELEM, step);
		memcpy(v, aux, n * sizeof(double));
		return;
	}

	ffht_compute_base(v, n, step);

	step = 1 << step;
	for (step = max(step, FFHT_BASE_ELEM); step < n; step <<= 1) {
		for (i = 0; i < n; i += (step << 1)) {
			for (j = i; j < i + step; j += FFHT_ATOM_ELEM) {
				ffht_atom_load(x, v + j);
				ffht_atom_load(y, v + j + step);
				ffht_atom_store_add(v + j, x, y);
				ffht_atom_store_sub(v + j + step, x, y);
			}
		}
	}

}
