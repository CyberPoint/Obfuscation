#ifndef __IND_OBFUSCATION__UTILS_H__
#define __IND_OBFUSCATION__UTILS_H__
#include <stddef.h>
#include <gmp.h>
#include<vector>

extern int g_verbose;

double
current_time(void);

int
seed_rng(gmp_randstate_t *rng);

int
load_mpz_scalar(const char *fname, mpz_t x);

int
save_mpz_scalar(const char *fname, const mpz_t x);

int
load_mpz_vector(const char *fname, mpz_t *m, const int len);

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len);

void
mult_mats(mpz_t *result, const mpz_t *left, const mpz_t *right, const mpz_t q,
		  long m, long n, long p);

void
mult_vect_by_mat(mpz_t *v, const mpz_t *m, mpz_t q, int size, mpz_t *tmparray);

void
mult_vect_by_vect(mpz_t out, const mpz_t *m, const mpz_t *v, mpz_t q, int size);

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits);


std::vector<char> mpz_to_vec(mpz_t in);

std::vector<std::vector<char> > mpz_arr_to_vec(mpz_t* in, unsigned long len);

void vec_to_mpz (mpz_t out, std::vector<char> in);

void vec_to_mpz_array (mpz_t* answer, std::vector<std::vector<char> > in);

void broadcast_vector(std::vector<char> * in);

void broadcast_vec_vec(std::vector<std::vector<char> > * in);

#endif
