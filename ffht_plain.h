/*
 * File:   ffht_plain.h
 * Author: Canales-Martinez, Isaac Andres
 */

#ifndef FFHT_PLAIN_H
#define FFHT_PLAIN_H

/*
 * The fast Fourier-Hadamard transform (FFT) operates on an input vector of real
 * numbers, which are regarded as 64-bit floating-point numbers. ffht_atom is
 * the data type used to store and make computations on these numbers.
 */
#define ffht_atom    double

/*
 * Number of elements (i.e., 64-bit floating-point numbers) stored in a single
 * ffht_atom.
 */
#define FFHT_ATOM_ELEM    1

/*
 * Number of initial iterations of the FFT for which this implementation is
 * optimised. For this purpouse, the function
 *
 *   ffht_compute_base()
 *
 * must be defined. Starting from iteration FFHT_BASE_LOG2 + 1, inclusive, the
 * FFT algorithm must proceed as normal.
 */
#define FFHT_BASE_LOG2    0

/*
 * For convinience, this is set to 2^FFHT_BASE_LOG2.
 */
#define FFHT_BASE_ELEM    1

/*
 * ffht_atom_load()
 *
 *     Set to a the value of (FFHT_ATOM_ELEM elements from) m.
 *
 * @a: ffht_atom to set
 * @m: ffht_atom to load the value from
 */
#define ffht_atom_load(a, m)    ((a) = *(m))

/*
 * ffht_atom_store_add()
 *
 *     Set to m the value of the element-wise addition a + b.
 *
 * @m: ffht_atom to store the result in
 * @a: addend
 * @b: addend
 */
#define ffht_atom_store_add(m, a, b)    (*(m) = (a) + (b))

/*
 * ffht_atom_store_sub()
 *
 *     Set to m the value of the element-wise subtraction a - b.
 *
 * @m: ffht_atom to store the result in
 * @a: minuend
 * @b: subtrahend
 */
#define ffht_atom_store_sub(m, a, b)    (*(m) = (a) - (b))

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
#define ffht_compute_base(v, n, step)

#endif /* FFHT_PLAIN_H */
