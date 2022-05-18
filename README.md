Implementation of the fast Fourier-Hadamard transform.

## How to use

Copy the three header files and C source file where required. If AVX instructions are supported, specify the required compiler options to use these instructions, otherwise, the non-AVX implementation will be used.

For example, with GCC, to use AVX instructions on a modern CPU
```
gcc -march=native -c -o ffht.o ffht.c
```
To use the non-AVX implementation
```
gcc -c -o ffht.o ffht.c
```
