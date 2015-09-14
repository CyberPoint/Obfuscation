
#include "clt_mlm.h"
#include "utils.h"  
#include <stdarg.h>
#include <mpi.h>
#include <omp.h>

//#include <gperftools/profiler.h>


static void
set_indices_pows(int *indices, int *pows, unsigned int num, ...)
{
    va_list elems;

    va_start(elems, num);
    #pragma omp parallel for
    for (unsigned int i = 0; i < num; ++i) {
        indices[i] = va_arg(elems, int);
        pows[i] = va_arg(elems, int);
    }
}


void try_encoding(struct clt_mlm_state * mlm, int n, int m,int * degs){

mpz_t tmp, c_star;
    mpz_t *alphas, *betas;
 int *indices, *pows;
int idx_set_size;


mpz_inits(c_star, tmp, NULL);

    alphas = (mpz_t *) malloc(sizeof(mpz_t) * n);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        mpz_init(alphas[i]);
        mpz_urandomm(alphas[i], mlm->rng, mlm->gs[1]);
    }
    betas = (mpz_t *) malloc(sizeof(mpz_t) * m);
    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        mpz_init(betas[i]);
        mpz_urandomm(betas[i], mlm->rng, mlm->gs[1]);
    }

/*
std::vector<std::vector<char> > alpha_send = mpz_arr_to_vec(alphas, n);
	std::vector<std::vector<char> > beta_send = mpz_arr_to_vec(betas, m);
//	broadcast_vec_vec(&alpha_send);
//	broadcast_vec_vec(&beta_send);
	vec_to_mpz_array(alphas,alpha_send);
	vec_to_mpz_array(betas,beta_send); */



    // The index set is laid out as follows:
    //   - The first 2 * n entries contain X_i,0 and X_i,1
    //   - The next n entries contain Z_i
    //   - The next n entries contain W_i
    //   - The final entry contains Y
    idx_set_size = 4 * n + 1;

    indices = (int *) malloc(sizeof(int) * idx_set_size);
    pows = (int *) malloc(sizeof(int) * idx_set_size);

    for (int i = 0; i < n; ++i) {
        mpz_t out, elems[2];
        int deg;

        mpz_inits(out, elems[0], elems[1], NULL);

        deg = degs[i];

	fprintf(stderr,"In loop: %d\n",i);

        set_indices_pows(indices, pows, 1, 2 * i, 1);
        mpz_set_ui(elems[0], 0);
        mpz_set(elems[1], alphas[i]);
        clt_mlm_encode(mlm, out, 2, elems, 1, indices, pows);
        //(void) snprintf(fname, fnamelen, "x_%d_0", i);
        //(void) write_element(s->dir, out, fname);

	fprintf(stderr,"After encode 1\n");

        set_indices_pows(indices, pows, 1, 2 * i, 1);
        mpz_set_ui(elems[0], 1);
        mpz_set_ui(elems[1], 1);
        clt_mlm_encode(mlm, out, 2, elems, 1, indices, pows);
        //(void) snprintf(fname, fnamelen, "u_%d_0", i);
        //(void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 2 * i + 1, 1);
        mpz_set_ui(elems[0], 1);
        mpz_set(elems[1], alphas[i]);
        clt_mlm_encode(mlm, out, 2, elems, 1, indices, pows);
        //(void) snprintf(fname, fnamelen, "x_%d_1", i);
        //(void) write_element(s->dir, out, fname);

fprintf(stderr,"After encode 3\n");
        set_indices_pows(indices, pows, 1, 2 * i + 1, 1);
        mpz_set_ui(elems[0], 1);
        mpz_set_ui(elems[1], 1);
        clt_mlm_encode(mlm, out, 2, elems, 1, indices, pows);
        //(void) snprintf(fname, fnamelen, "u_%d_1", i);
        //(void) write_element(s->dir, out, fname);

        mpz_urandomm(elems[0], mlm->rng, mlm->gs[0]);
        mpz_urandomm(elems[1], mlm->rng, mlm->gs[1]);

        set_indices_pows(indices, pows, 3, 2 * i + 1, deg, 2 * n + i, 1,
                         3 * n + i, 1);
        clt_mlm_encode(mlm, out, 2, elems, 3, indices, pows);
        //(void) snprintf(fname, fnamelen, "z_%d_0", i);
        //(void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 3 * n + i, 1);
        mpz_set_ui(elems[0], 0);
        clt_mlm_encode(mlm, out, 2, elems, 1, indices, pows);
        //(void) snprintf(fname, fnamelen, "w_%d_0", i);
        //(void) write_element(s->dir, out, fname);

        mpz_urandomm(elems[0], mlm->rng, mlm->gs[0]);
        mpz_urandomm(elems[1], mlm->rng, mlm->gs[1]);

        set_indices_pows(indices, pows, 3, 2 * i, deg, 2 * n + i, 1, 3 * n + i, 1);
        clt_mlm_encode(mlm, out, 2, elems, 3, indices, pows);
        //(void) snprintf(fname, fnamelen, "z_%d_1", i);
        //(void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 3 * n + i, 1);
        mpz_set_ui(elems[0], 0);
        clt_mlm_encode(mlm, out, 2, elems, 1, indices, pows);
        //(void) snprintf(fname, fnamelen, "w_%d_1", i);
        //(void) write_element(s->dir, out, fname);

        mpz_clears(out, elems[0], elems[1], NULL);
    }

}


bool clt_struct_equal(struct clt_mlm_state * in1, struct clt_mlm_state * in2){
bool answer = true;

	answer &= (mpz_cmp(in1->q,in2->q)==0);
	answer &= (mpz_cmp(in1->pzt,in2->pzt)==0);
	for (int i = 0; i< in1->n; i++){
		fprintf(stderr,"Comparing gs %d \n",i);
		answer &= (mpz_cmp(in1->gs[i],in2->gs[i])==0);
	}


return answer;
}






int main(int argc, char** argv){

MPI_Init(&argc,&argv);



struct clt_mlm_state *  mlm = new struct clt_mlm_state;
/*mlm->secparam = 12;
mlm->nzs = 33;
long int pows [33] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,8 };
*/

mlm->secparam = 8;
mlm->nzs = 17;
long int pows [17] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,4 };


clt_mlm_setup (mlm, "driven", pows, 17,17, 1);
 int degs[4] = {1,1,1,1};


fprintf(stderr,"Before, has g[0] = %lu\n",mpz_get_ui(mlm->gs[0]));
fprintf(stderr,"Before, has q = %lu\n",mpz_get_ui(mlm->q));
struct clt_mlm_ser_enc* sent_mlm = stat_ser(mlm);
	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!1Serialized okay!!!!!!!!!!!!!!!!!!!!!\n");
	broadcast_clt_mlm_ser_enc(sent_mlm);
	fprintf(stderr,"*********************Broadcast okay***********************\n");
	/*mpz_t checker;
mpz_init(checker);
vec_to_mpz(checker,sent_mlm->gs[0]);
fprintf(stderr,"After send BUT before unser, gs[0]= %lu\n",mpz_get_ui(checker));*/
(ser_stat(sent_mlm,mlm));
	fprintf(stderr,"*********************DECODED okay***********************\n");
fprintf(stderr,"After, has g[0] = %lu\n",mpz_get_ui(mlm->gs[0]));
fprintf(stderr,"After, has q = %lu\n",mpz_get_ui(mlm->q));
//fprintf(stderr,"Process has pzt = %lu\n",mpz_get_ui(mlm->pzt));
//if (clt_struct_equal(mlm,a_mlm))fprintf(stderr,"\nEQUAL!\n");
//ProfilerStart("encoding_4_8.prof");
try_encoding(mlm,4,4,degs);
//ProfilerStop();

delete mlm;

MPI_Finalize();

return 0;
}
