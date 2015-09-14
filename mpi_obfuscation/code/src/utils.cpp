#include "utils.h"

#include <fcntl.h>
#include <gmp.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
# include  <mpi.h>


int g_verbose;

// XXX: The use of /dev/urandom is not secure; however, the supercomputer we run
// on doesn't appear to have enough entropy, and blocks for long periods of
// time.  Thus, we use /dev/urandom instead.
#ifndef RANDFILE
#  define RANDFILE "/dev/urandom"
#endif

double
current_time(void)
{
    struct timeval t;
    (void) gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
}

int
seed_rng(gmp_randstate_t *rng)
{
    int file;
    if ((file = open(RANDFILE, O_RDONLY)) == -1) {
        (void) fprintf(stderr, "Error opening %s\n", RANDFILE);
        return 1;
    } else {
        unsigned long seed;
        if (read(file, &seed, sizeof seed) == -1) {
            (void) fprintf(stderr, "Error reading from %s\n", RANDFILE);
            (void) close(file);
            return 1;
        } else {
            if (g_verbose)
                (void) fprintf(stderr, "  Seed: %lu\n", seed);

            gmp_randinit_default(*rng);
            gmp_randseed_ui(*rng, seed);
        }
    }
    if (file != -1)
        (void) close(file);
    return 0;
}

int
load_mpz_scalar(const char *fname, mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    (void) mpz_inp_raw(x, f);
    (void) fclose(f);
    return 0;
}

int
save_mpz_scalar(const char *fname, const mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    if (mpz_out_raw(f, x) == 0) {
        (void) fclose(f);
        return 1;
    }
    (void) fclose(f);
    return 0;
}

int
load_mpz_vector(const char *fname, mpz_t *m, const int len)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        (void) mpz_inp_raw(m[i], f);
    }
    (void) fclose(f);
    return 0;
}

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        if (mpz_out_raw(f, m[i]) == 0) {
            (void) fclose(f);
            return 1;
        }
    }
    (void) fclose(f);
    return 0;
}

void
mult_mats(mpz_t *result, const mpz_t *left, const mpz_t *right, const mpz_t q,
          long m, long n, long p)
{
    mpz_t *tmparray;
    double start, end;

    start = current_time();
    tmparray = (mpz_t *) malloc(sizeof(mpz_t) * m * p);
    for (int i = 0; i < m * p; ++i) {
        mpz_init(tmparray[i]);
    }
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            mpz_t tmp, sum;
            mpz_inits(tmp, sum, NULL);
            for (int k = 0; k < n; ++k) {
                mpz_mul(tmp,
                        left[k * m + (i * m + j) % m],
                        right[k + n * ((i * m + j) / m)]);
                mpz_add(sum, sum, tmp);
                mpz_mod(sum, sum, q);
            }
            mpz_set(tmparray[i * n + j], sum);
            mpz_clears(tmp, sum, NULL);
        }
    }
    for (int i = 0; i < m * p; ++i) {
        mpz_swap(result[i], tmparray[i]);
        mpz_clear(tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, " Multiplying took: %f\n", end - start);
}

void
mult_vect_by_mat(mpz_t *v, const mpz_t *m, mpz_t q, int size, mpz_t *tmparray)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        mpz_t tmp, sum;
        mpz_inits(tmp, sum, NULL);
        for (int j = 0; j < size; ++j) {
            mpz_mul(tmp, v[j], m[i * size + j]);
            mpz_add(sum, sum, tmp);
            mpz_mod(sum, sum, q);
        }
        mpz_set(tmparray[i], sum);
        mpz_clears(tmp, sum, NULL);
    }
    for (int i = 0; i < size; ++i) {
        mpz_swap(v[i], tmparray[i]);
    }
}

void
mult_vect_by_vect(mpz_t out, const mpz_t *v, const mpz_t *u, mpz_t q, int size)
{
    mpz_set_ui(out, 0);
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mul(tmp, v[i], u[i]);
#pragma omp critical
        {
            mpz_add(out, out, tmp);
            mpz_mod(out, out, q);
        }
        mpz_clears(tmp, NULL);
    }
}

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, *rng, nbits);
    mpz_clear(one);
}


std::vector<char> mpz_to_vec(mpz_t in){
	if(mpz_cmp_ui(in,0)==0){
		return *(new std::vector<char>);
	}
	size_t* count = new size_t;
	void* rop = NULL;
	
	char* pack = (char*) mpz_export(rop, count, 1,sizeof(char), 1, 0, in);
	/*size_t realcount = (mpz_sizeinbase(in,2) + 8*sizeof(char) -1)/(8 * sizeof(char));
	fprintf(stderr, "REAL COUNT: %lu\n",realcount);
	fprintf(stderr,"COUNT: %lu\n",*count);
	if (pack==NULL) {fprintf(stderr,"READING FROM NULL?!?!?!?\n");}*/
	std::vector<char> packaged(pack,pack + *count);
	delete count;
	pack = NULL;
	rop = NULL;
	return packaged;
}

std::vector<std::vector<char> > mpz_arr_to_vec(mpz_t* in, unsigned long len){
	std::vector<std::vector<char> > ans;
	for (unsigned long i = 0; i< len; i++){
		ans.push_back(mpz_to_vec(in[i]));
	}
	return ans;
}

void vec_to_mpz (mpz_t out, std::vector<char> in){
	if (in.size()==0){
		mpz_set_ui(out,0);
		return;
	}	
	void* unpack = in.data();
	size_t count = in.size();
	mpz_import(out,count,1,sizeof(char),1,0,unpack);
	
}

void vec_to_mpz_array (mpz_t * answer, std::vector<std::vector<char> > in){
	if (answer!=NULL)	
	delete[] answer;	
	answer=new mpz_t[in.size()];
	for(unsigned long i = 0; i<in.size(); i++){
		mpz_init(answer[i]);
		vec_to_mpz(answer[i],in[i]);
	}		
}


void broadcast_vector(std::vector<char> * in){
	
	int size;
	size = (*in).size();
	MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
	(*in).resize(size);
	MPI_Bcast(&((*in).front()), size, MPI_CHAR, 0, MPI_COMM_WORLD);

}

void broadcast_vec_vec(std::vector<std::vector<char> > * in){
	int sizein = (*in).size();
	MPI_Bcast(&sizein,1,MPI_INT,0,MPI_COMM_WORLD);
	int* locs = new int[sizein+1];
	int running = 0;
	std::vector<char> to_send;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank==0){
	for (int i = 0; i<sizein; i++){
		to_send.insert(to_send.end(),(*in)[i].begin(),(*in)[i].end() );
		locs[i] = running;	//could keep from sending first 0. But probably not significant. 	
		running+=(*in)[i].size();

	}
	locs[sizein] = running;
	}
	MPI_Bcast(locs,sizein+1,MPI_INT,0,MPI_COMM_WORLD);
	broadcast_vector(&to_send);
	(*in).resize(sizein);

	//fprintf(stderr,"MADE IT PAST BCASTS OF VECVEC\n");	

	for(int i = 0; i<sizein; i++)
	{
		std::vector<char> copy;
		int loopsize = locs[i+1]-locs[i];
		copy.resize(loopsize);
		for (int j = 0; j<loopsize; j++){
			copy[j] = to_send[locs[i]+j];
		}
		(*in)[i] = copy;
	}

	

	delete[] locs;

	/*std::vector<char> to_send;	
	int size;
	size = (*in).size();
	MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
	(*in).resize(size);
	for (int i = 0; i<size; i++){
		 to_send = (*in)[i];
		 broadcast_vector(&to_send);
		(*in)[i]=to_send;
	}*/





}

