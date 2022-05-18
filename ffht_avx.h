/*
 * File:   ffht_avx.h
 * Author: Canales-Martinez, Isaac Andres
 */

#ifndef FFHT_AVX_H
#define FFHT_AVX_H

#include <stddef.h>
#include <immintrin.h>

/*
 * The fast Fourier-Hadamard transform (FFT) operates on an input vector of real
 * numbers, which are regarded as 64-bit floating-point numbers. ffht_atom is
 * the data type used to store and make computations on these numbers.
 */
#define ffht_atom    __m256d

/*
 * Number of elements (i.e., 64-bit floating-point numbers) stored in a single
 * ffht_atom.
 */
#define FFHT_ATOM_ELEM    4

/*
 * Number of initial iterations of the FFT for which this implementation is
 * optimised. For this purpouse, the function
 *
 *   ffht_compute_base()
 *
 * must be defined. Starting from iteration FFHT_BASE_LOG2 + 1, inclusive, the
 * FFT algorithm must proceed as normal.
 */
#define FFHT_BASE_LOG2    3

/*
 * For convinience, this is set to 2^FFHT_BASE_LOG2.
 */
#define FFHT_BASE_ELEM    8

/*
 * ffht_atom_load()
 *
 *     Set to a the value of (FFHT_ATOM_ELEM elements from) m.
 *
 * @a: ffht_atom to set
 * @m: ffht_atom to load the value from
 */
#define ffht_atom_load(a, m)    \
	((a) = _mm256_loadu_pd((m)))

/*
 * ffht_atom_store_add()
 *
 *     Set to m the value of the element-wise addition a + b.
 *
 * @m: ffht_atom to store the result in
 * @a: addend
 * @b: addend
 */
#define ffht_atom_store_add(m, a, b)    \
	(_mm256_storeu_pd((m), _mm256_add_pd((a), (b))))

/*
 * ffht_atom_store_sub()
 *
 *     Set to m the value of the element-wise subtraction a - b.
 *
 * @m: ffht_atom to store the result in
 * @a: minuend
 * @b: subtrahend
 */
#define ffht_atom_store_sub(m, a, b)    \
	(_mm256_storeu_pd((m), _mm256_sub_pd((a), (b))))

#define ffht_compute_base_0(v, a, b, c, d)    \
	do {\
		(a) = _mm256_set_pd((v)[6], (v)[4], (v)[2], (v)[0]);\
		(b) = _mm256_set_pd((v)[7], (v)[5], (v)[3], (v)[1]);\
		(c) = _mm256_add_pd((a), (b));\
		(d) = _mm256_sub_pd((a), (b));\
		(a) = _mm256_unpacklo_pd((c), (d));\
		(b) = _mm256_unpackhi_pd((c), (d));\
		(c) = _mm256_add_pd((a), (b));\
		(d) = _mm256_sub_pd((a), (b));\
		(a) = _mm256_permute2f128_pd((c), (d), 0x20);\
		(b) = _mm256_permute2f128_pd((c), (d), 0x31);\
		_mm256_storeu_pd((v), _mm256_add_pd((a), (b)));\
		_mm256_storeu_pd(((v) + 4), _mm256_sub_pd((a), (b)));\
	} while (0)

#define ffht_compute_base_1(v, a, b, c, d)    \
	do {\
		(a) = _mm256_set_pd((v)[5], (v)[4], (v)[1], (v)[0]);\
		(b) = _mm256_set_pd((v)[7], (v)[6], (v)[3], (v)[2]);\
		(c) = _mm256_add_pd((a), (b));\
		(d) = _mm256_sub_pd((a), (b));\
		(a) = _mm256_permute2f128_pd((c), (d), 0x20);\
		(b) = _mm256_permute2f128_pd((c), (d), 0x31);\
		_mm256_storeu_pd((v), _mm256_add_pd((a), (b)));\
		_mm256_storeu_pd(((v) + 4), _mm256_sub_pd((a), (b)));\
	} while (0)

#define ffht_compute_base_2(v, a, b, c, d)    \
	do {\
		(a) = _mm256_set_pd((v)[3], (v)[2], (v)[1], (v)[0]);\
		(b) = _mm256_set_pd((v)[7], (v)[6], (v)[5], (v)[4]);\
		_mm256_storeu_pd((v), _mm256_add_pd((a), (b)));\
		_mm256_storeu_pd(((v) + 4), _mm256_sub_pd((a), (b)));\
	} while (0)

/*
 * ffht_compute_base()
 *
 *     Compute the first FFHT_BASE_LOG2 iterations of the fast Fourier-Hadamard
 *     transform on n = 2^m real numbers, given that the first step iterations
 *     have already been computed.
 *
 *     The n input numbers are specified in the array v. The computation is
 *     performed in-place, i.e., the contents of v are overwritten.
 *
 *     It is assumed that
 *       - v is not a NULL pointer;
 *       - n = 2^m, where m >= 0;
 *       - the size of the array pointed to by v is at least n, and
 *       - 0 <= step <= m.
 *
 * @v:       Array of input numbers.
 * @n:       Number of input numbers.
 * @step:    Number of initial iterations already computed.
 */
static void ffht_compute_base(double *v, size_t n, size_t step)
{
	__m256d a, b, c, d;
	size_t i;

	n >>= FFHT_BASE_LOG2;

	switch (step) {
	case 0:
		for (i = 0; i < n; i++, v += FFHT_BASE_ELEM)
			ffht_compute_base_0(v, a, b, c, d);
		break;
	case 1:
		for (i = 0; i < n; i++, v += FFHT_BASE_ELEM)
			ffht_compute_base_1(v, a, b, c, d);
		break;
	case 2:
		for (i = 0; i < n; i++, v += FFHT_BASE_ELEM)
			ffht_compute_base_2(v, a, b, c, d);
		break;
	}
}

#endif /* FFHT_AVX_H */
