#ifndef __CLT_MLM_H__
#define __CLT_MLM_H__
#include <stddef.h>
#include <gmpxx.h>
#include <vector>

#include "gperftools-2.4/src/gperftools/profiler.h"

 /* Test for GCC > 4.9.1 */
          #if (__GNUC__ > 4) || (__GNUC__ == 4 && (__GNUC_MINOR__ > 9 ||(__GNUC_MINOR__ == 9 && __GNUC_PATCHLEVEL__ > 1)))
			#ifndef OpenMP_4
				#define OpenMP_4
			#endif
		#endif

struct clt_mlm_state {
	gmp_randstate_t rng;
    unsigned long secparam;
    unsigned long n;
    unsigned long nzs;
    unsigned long rho;
    mpz_t q;
    mpz_t pzt;
    mpz_t *gs;
    mpz_t *crt_coeffs;
    mpz_t *zinvs;
};

int
clt_mlm_setup(struct clt_mlm_state *s, const char *dir, const long *pows,
              long kappa, long size, int verbose);

void
clt_mlm_cleanup(struct clt_mlm_state *s);

void
clt_mlm_encode(struct clt_mlm_state *s, mpz_t out, size_t nins,
			   const mpz_t *ins, unsigned int nzs, const int *indices,
			   const int *pows);

int
clt_mlm_is_zero(const mpz_t c, const mpz_t pzt, const mpz_t q, long nu);

struct clt_mlm_ser_enc{
unsigned long secparam;
unsigned long n;
unsigned long nzs;
unsigned long rho;
std::vector<char> q;
std::vector<char> pzt;
std::vector<std::vector<char> > gs;
std::vector<std::vector<char> > crt_coeffs;
std::vector<std::vector<char> > zinvs;

};

struct clt_mlm_ser_enc * stat_ser(struct clt_mlm_state *s);

struct clt_mlm_state * ser_stat(struct clt_mlm_ser_enc* ser,struct clt_mlm_state * in);

void broadcast_clt_mlm_ser_enc(struct clt_mlm_ser_enc* in);

	#ifdef OpenMP_4
#pragma omp declare reduction (mpz_class_mul_red : mpz_class : omp_out *= omp_in) initializer(omp_priv =1)

#pragma omp declare reduction (mpz_class_add_red : mpz_class : omp_out += omp_in) initializer(omp_priv =0)

#endif

#endif
