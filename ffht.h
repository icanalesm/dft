/*
 * File:   ffht.h
 * Author: Canales-Martinez, Isaac Andres
 */

#ifndef FFHT_H
#define FFHT_H

#include <stddef.h>

/*
 * ffht()
 *
 *     Compute the Fourier-Hadamard spectrum of n = 2^m real numbers using the
 *     fast Fourier-Hadamard transform.
 *
 *     The n input numbers are specified in the array v. The computation is
 *     performed in-place, i.e., the Fourier-Hadamard spectrum is stored in v.
 *
 *     It is assumed that
 *       - v is not a NULL pointer;
 *       - n = 2^m, where m >= 0, and
 *       - the size of the array pointed to by v is at least n.
 *
 * @v:    Array of input numbers.
 * @n:    Number of input numbers.
 */
void ffht(double *v, size_t n);

/*
 * ffht_part()
 *
 *     Compute the Fourier-Hadamard spectrum of n = 2^m real numbers using the
 *     fast Fourier-Hadamard transform, given that the first step iterations
 *     have already been computed.
 *
 *     The n input numbers are specified in the array v. The computation is
 *     performed in-place, i.e., the Fourier-Hadamard spectrum is stored in v.
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
void ffht_part(double *v, size_t n, size_t step);

#endif /* FFHT_H */
