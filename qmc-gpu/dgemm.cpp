/*
In case you're wondering, dgemm stands for Double-precision, GEneral
Matrix-Matrix multiplication.
*/
#include "dgemm.h"
#include <stdlib.h>
const char* dgemm_desc = "tuned dgemm";

DGEMM_ROUTINE_FP dgemm_routines[8][8][1] = {
	{{basic_dgemm_1x1x1},{basic_dgemm_1x2x1},{basic_dgemm_1x3x1},{basic_dgemm_1x4x1},{basic_dgemm_1x5x1},{basic_dgemm_1x6x1},{basic_dgemm_1x7x1},{basic_dgemm_1x8x1}},
	{{basic_dgemm_2x1x1},{basic_dgemm_2x2x1},{basic_dgemm_2x3x1},{basic_dgemm_2x4x1},{basic_dgemm_2x5x1},{basic_dgemm_2x6x1},{basic_dgemm_2x7x1},{basic_dgemm_2x8x1}},
	{{basic_dgemm_3x1x1},{basic_dgemm_3x2x1},{basic_dgemm_3x3x1},{basic_dgemm_3x4x1},{basic_dgemm_3x5x1},{basic_dgemm_3x6x1},{basic_dgemm_3x7x1},{basic_dgemm_3x8x1}},
	{{basic_dgemm_4x1x1},{basic_dgemm_4x2x1},{basic_dgemm_4x3x1},{basic_dgemm_4x4x1},{basic_dgemm_4x5x1},{basic_dgemm_4x6x1},{basic_dgemm_4x7x1},{basic_dgemm_4x8x1}},
	{{basic_dgemm_5x1x1},{basic_dgemm_5x2x1},{basic_dgemm_5x3x1},{basic_dgemm_5x4x1},{basic_dgemm_5x5x1},{basic_dgemm_5x6x1},{basic_dgemm_5x7x1},{basic_dgemm_5x8x1}},
	{{basic_dgemm_6x1x1},{basic_dgemm_6x2x1},{basic_dgemm_6x3x1},{basic_dgemm_6x4x1},{basic_dgemm_6x5x1},{basic_dgemm_6x6x1},{basic_dgemm_6x7x1},{basic_dgemm_6x8x1}},
	{{basic_dgemm_7x1x1},{basic_dgemm_7x2x1},{basic_dgemm_7x3x1},{basic_dgemm_7x4x1},{basic_dgemm_7x5x1},{basic_dgemm_7x6x1},{basic_dgemm_7x7x1},{basic_dgemm_7x8x1}},
	{{basic_dgemm_8x1x1},{basic_dgemm_8x2x1},{basic_dgemm_8x3x1},{basic_dgemm_8x4x1},{basic_dgemm_8x5x1},{basic_dgemm_8x6x1},{basic_dgemm_8x7x1},{basic_dgemm_8x8x1}}
};

/* You'll definitely change this... */
#if !defined(BLOCK_SIZE)
int BLOCK_SIZE = 512;
int BLOCK_SIZE_2=128;
#endif

extern double *AA;
int reg_r=4;
int reg_c=8;

static double A_temp[1000*1000];

void do_block_lv1 (const int lda,
			  const double *A, const double *B, double *C,
			  const int i, const int j, const int k);
void do_block_lv2 (const int lda, const int M, const int N, const int K, 
			  const double *A, const double *B, double *C,
			  const int i, const int j, const int k);
void fill_list_rec(int start, int M, int curr_m, int *lst, int lst_idx);

void fill_list(int M, int *lst);

void blocked_dgemm(const int lda, const int M, const int N, const int K,
			  const double *A, const double *B, double *C) ;

void basic_dgemm (const int lda, const int M, const int N, const int K,
			 const double *A, const double *B, double *C)
{
	int i, j, k;

	/*
	To optimize this, think about loop unrolling and software
	pipelining.  Hint:  For the majority of the matmuls, you
	know exactly how many iterations there are (the block size)...
	*/

	for (i = 0; i < M; ++i) {
		const register double *Ai_ = A + i*lda;
		for (j = 0; j < N; ++j) {
			const register double *B_j = B + j*lda;

			register double cij = C [j*lda + i];

			for (k = 0; k < K; ++k) {
				cij += Ai_[k] * B_j[k];
			}

			C[j*lda + i] = cij;
		}
	}
}


/*
A is M-by-K
B is K-by-N
C is M-by-N

lda is the leading dimension of the matrix (the M of square_dgemm).
*/


void
do_block_lv1 (const int lda,
			  const double *A, const double *B, double *C,
			  const int i, const int j, const int k)
{

	const int M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
	const int N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
	const int K = (k+BLOCK_SIZE > lda? lda-k : BLOCK_SIZE);

	//basic_dgemm (lda, M, N, K,
	//       A + i + k*lda, B + k + j*lda, C + i + j*lda);
	blocked_dgemm (lda, M, N, K,
		A + k + i*lda, B + k + j*lda, C + i + j*lda);
}
void
do_block_lv2 (const int lda, const int M, const int N, const int K, 
			  const double *A, const double *B, double *C,
			  const int i, const int j, const int k)
{
	DGEMM_ROUTINE_FP fp;
	const int Mp = (i+BLOCK_SIZE_2 > M? M-i : BLOCK_SIZE_2);
	const int Np = (j+BLOCK_SIZE_2 > N? N-j : BLOCK_SIZE_2);
	const int Kp = (k+BLOCK_SIZE_2 > K? K-k : BLOCK_SIZE_2);

	fp = dgemm_routines[reg_r-1][reg_c-1][0];

	(*fp)(lda, Mp, Np, Kp,
		A + k + i*lda, B + k + j*lda, C + i + j*lda);

}
/*Method for filling in order of list with spiral pattern.. uses recursion*/

void fill_list_rec(int start, int M, int curr_m, int *lst, int lst_idx) {
	int i;
	int idx=lst_idx;



	if(curr_m == 1) {
		lst[lst_idx] = start;
		return;
	}
	if(curr_m == 2) {
		lst[lst_idx] = start;
		lst[lst_idx+1] = start+1;
		lst[lst_idx+2] = start+1+M;
		lst[lst_idx+3] = start+M;
		return;
	}


	for (i=0; i<curr_m; i++) {
		lst[idx] = i+start;
		idx++;
	}

	for(i=1; i<curr_m; i++) {
		lst[idx] = lst[idx-1]+M;
		idx++;
	}

	for(i=curr_m-2; i>=0; i--) {
		lst[idx] = lst[idx-1]-1;
		idx++;
	}

	for(i=curr_m-2; i>=1; i--) {
		lst[idx] = lst[idx-1]-M;
		idx++;
	}
	fill_list_rec(start+M+1, M, curr_m-2, lst, idx);

}
void fill_list(int M, int *lst) {
	fill_list_rec(0, M, M, lst, 0);
}

void transpose_major_block(const int lda, int i, int j,
						   double *A, double *AA){

	const int M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
	const int N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
	int ii,jj;

	for(ii=0; ii<M; ii++){
		for (jj=0; jj<N; jj++){
			AA[ii*lda+jj] = A[i+jj*lda];
		}
	}
}

void col_maj_2_row_maj(const int M,
	const double *A, double *AA){

	/*   int lda=M; */
	/*   int bi, bj; */
	/*   const int n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0); */

	/*   for(bi=0; bi<n_blocks; bi++) { */
	/*     const int i = bi*BLOCK_SIZE; */
	/*     for(bj=0; bj<n_blocks; bj++) { */
	/*       const int j = bj*BLOCK_SIZE; */
	/*       transpose_major_block(M, i, j, A +  i + j*lda, AA + i*lda + j); */
	/*     } */
	/*   } */

	double *col, *col_s;
	double* AA_s;
	int ii,jj;
	col_s = (double*) A;
	AA_s = AA;

	for(jj=0; jj<M; jj++){
		col = col_s;
		AA = AA_s;
		for(ii=0; ii<M; ii++){
			*AA = *col;
			col++;
			AA += M;
		}
		col_s += M;
		AA_s ++;
	}
	AA = AA_s - M;

	}

void square_dgemm (const int M, const double *A, const double *B, double *C)
{
	const int n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0);
	int bi, bj, bk, bx;
	int * block_order = new int[n_blocks*n_blocks];
	//int block_order[n_blocks*n_blocks];
	double *AA;

	fill_list(n_blocks, block_order);
	if(M<=1000)
		AA = A_temp;
	else {
		AA = (double*)malloc(sizeof(double)*M*M);
	}
	col_maj_2_row_maj(M,A,AA);



	/*    block_order = (int*) malloc(sizeof(int)*n_blocks*n_blocks); */
	/*   fill_list(n_blocks, block_order); */



	for(bi =0; bi < n_blocks; ++bi) {
		const int i = bi * BLOCK_SIZE;

		for(bj = 0; bj < n_blocks; ++bj) {
			const int j = (bi&0x1) ?  ((n_blocks-1-bj)*BLOCK_SIZE) : (bj * BLOCK_SIZE);

			for (bk = 0; bk < n_blocks; ++bk) {
				const int k = ((bi+bj)&0x1) ? ((n_blocks-1-bk)*BLOCK_SIZE) : (bk * BLOCK_SIZE);
				do_block_lv1 (M, AA, B, C, i, j, k);
			}
		}
	}
	/*   for (bx = 0; bx < n_blocks*n_blocks; bx++) { */
	/*     const int i = block_order[bx]%n_blocks * BLOCK_SIZE; */
	/*     const int j = block_order[bx]/n_blocks * BLOCK_SIZE; */

	/*     for (bk = 0; bk < n_blocks; ++bk) { */
	/*       const int k = bk * BLOCK_SIZE; */
	/*       do_block_lv1 (M, A, B, C, i, j, k); */
	/*     } */
	/*   } */
	if(M>1000)
		free(AA);

}

void 
blocked_dgemm(const int lda, const int M, const int N, const int K,
			  const double *A, const double *B, double *C) 
{

	int bi, bj, bk;
	int *block_order;
	const int m_blocks = M / BLOCK_SIZE_2 + (M%BLOCK_SIZE_2? 1 : 0);
	const int n_blocks = N / BLOCK_SIZE_2 + (N%BLOCK_SIZE_2? 1 : 0);
	const int k_blocks = K / BLOCK_SIZE_2 + (K%BLOCK_SIZE_2? 1 : 0);

	for (bi=0; bi<m_blocks; bi++) {
		const int i = bi*BLOCK_SIZE_2;

		for (bj=0; bj<n_blocks; bj++) {
			//const int j = bj*BLOCK_SIZE_2;
			const int j = (bi&0x1) ?  ((n_blocks-1-bj)*BLOCK_SIZE_2) : (bj * BLOCK_SIZE_2);

			for (bk=0; bk<k_blocks; bk++) {
				//const int k = bk*BLOCK_SIZE_2;
				const int k = ((bi+bj)&0x1) ? ((k_blocks-1-bk)*BLOCK_SIZE_2) : (bk * BLOCK_SIZE_2);	 
				do_block_lv2(lda, M,N,K,A,B,C, i,j,k);
			}
		}
	}
	/*   for (bi=0; bi<m_blocks; bi++) { */
	/*     const int i = bi*BLOCK_SIZE_2; */

	/*     for (bj=0; bj<n_blocks; bj++) { */
	/*       const int j = bj*BLOCK_SIZE_2; */

	/*       for (bk=0; bk<k_blocks; bk++) { */
	/* 	const int k = bk*BLOCK_SIZE_2; */
	/* 	do_block_lv2(lda, M,N,K,A,B,C, i,j,k); */
	/*       } */
	/*     } */
	/*   } */
}

void set_l_block_size(int bls) {
	BLOCK_SIZE = bls;
}

void set_s_block_size(int bls) {
	BLOCK_SIZE_2 = bls;
}

void set_r(int r) {
	reg_r = r;
}

void set_c(int c){
	reg_c = c;
}



void basic_dgemm_1x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double * Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			register double cij_3 = C[((j+2)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
			C[((j+2)*lda)+i+0] = cij_3 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			register double cij_3 = C[((j+2)*lda)+i+0];
			register double cij_4 = C[((j+3)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
			C[((j+2)*lda)+i+0] = cij_3 ;
			C[((j+3)*lda)+i+0] = cij_4 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			register double cij_3 = C[((j+2)*lda)+i+0];
			register double cij_4 = C[((j+3)*lda)+i+0];
			register double cij_5 = C[((j+4)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
			C[((j+2)*lda)+i+0] = cij_3 ;
			C[((j+3)*lda)+i+0] = cij_4 ;
			C[((j+4)*lda)+i+0] = cij_5 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			register double cij_3 = C[((j+2)*lda)+i+0];
			register double cij_4 = C[((j+3)*lda)+i+0];
			register double cij_5 = C[((j+4)*lda)+i+0];
			register double cij_6 = C[((j+5)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
			C[((j+2)*lda)+i+0] = cij_3 ;
			C[((j+3)*lda)+i+0] = cij_4 ;
			C[((j+4)*lda)+i+0] = cij_5 ;
			C[((j+5)*lda)+i+0] = cij_6 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			register double cij_3 = C[((j+2)*lda)+i+0];
			register double cij_4 = C[((j+3)*lda)+i+0];
			register double cij_5 = C[((j+4)*lda)+i+0];
			register double cij_6 = C[((j+5)*lda)+i+0];
			register double cij_7 = C[((j+6)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
			C[((j+2)*lda)+i+0] = cij_3 ;
			C[((j+3)*lda)+i+0] = cij_4 ;
			C[((j+4)*lda)+i+0] = cij_5 ;
			C[((j+5)*lda)+i+0] = cij_6 ;
			C[((j+6)*lda)+i+0] = cij_7 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_1x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 1;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+1)*lda)+i+0];
			register double cij_3 = C[((j+2)*lda)+i+0];
			register double cij_4 = C[((j+3)*lda)+i+0];
			register double cij_5 = C[((j+4)*lda)+i+0];
			register double cij_6 = C[((j+5)*lda)+i+0];
			register double cij_7 = C[((j+6)*lda)+i+0];
			register double cij_8 = C[((j+7)*lda)+i+0];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+1)*lda)+i+0] = cij_2 ;
			C[((j+2)*lda)+i+0] = cij_3 ;
			C[((j+3)*lda)+i+0] = cij_4 ;
			C[((j+4)*lda)+i+0] = cij_5 ;
			C[((j+5)*lda)+i+0] = cij_6 ;
			C[((j+6)*lda)+i+0] = cij_7 ;
			C[((j+7)*lda)+i+0] = cij_8 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			register double cij_5 = C[((j+2)*lda)+i+0];
			register double cij_6 = C[((j+2)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
			C[((j+2)*lda)+i+0] = cij_5 ;
			C[((j+2)*lda)+i+1] = cij_6 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			register double cij_5 = C[((j+2)*lda)+i+0];
			register double cij_6 = C[((j+2)*lda)+i+1];
			register double cij_7 = C[((j+3)*lda)+i+0];
			register double cij_8 = C[((j+3)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
			C[((j+2)*lda)+i+0] = cij_5 ;
			C[((j+2)*lda)+i+1] = cij_6 ;
			C[((j+3)*lda)+i+0] = cij_7 ;
			C[((j+3)*lda)+i+1] = cij_8 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			register double cij_5 = C[((j+2)*lda)+i+0];
			register double cij_6 = C[((j+2)*lda)+i+1];
			register double cij_7 = C[((j+3)*lda)+i+0];
			register double cij_8 = C[((j+3)*lda)+i+1];
			register double cij_9 = C[((j+4)*lda)+i+0];
			register double cij_10 = C[((j+4)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
			C[((j+2)*lda)+i+0] = cij_5 ;
			C[((j+2)*lda)+i+1] = cij_6 ;
			C[((j+3)*lda)+i+0] = cij_7 ;
			C[((j+3)*lda)+i+1] = cij_8 ;
			C[((j+4)*lda)+i+0] = cij_9 ;
			C[((j+4)*lda)+i+1] = cij_10 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			register double cij_5 = C[((j+2)*lda)+i+0];
			register double cij_6 = C[((j+2)*lda)+i+1];
			register double cij_7 = C[((j+3)*lda)+i+0];
			register double cij_8 = C[((j+3)*lda)+i+1];
			register double cij_9 = C[((j+4)*lda)+i+0];
			register double cij_10 = C[((j+4)*lda)+i+1];
			register double cij_11 = C[((j+5)*lda)+i+0];
			register double cij_12 = C[((j+5)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
			C[((j+2)*lda)+i+0] = cij_5 ;
			C[((j+2)*lda)+i+1] = cij_6 ;
			C[((j+3)*lda)+i+0] = cij_7 ;
			C[((j+3)*lda)+i+1] = cij_8 ;
			C[((j+4)*lda)+i+0] = cij_9 ;
			C[((j+4)*lda)+i+1] = cij_10 ;
			C[((j+5)*lda)+i+0] = cij_11 ;
			C[((j+5)*lda)+i+1] = cij_12 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			register double cij_5 = C[((j+2)*lda)+i+0];
			register double cij_6 = C[((j+2)*lda)+i+1];
			register double cij_7 = C[((j+3)*lda)+i+0];
			register double cij_8 = C[((j+3)*lda)+i+1];
			register double cij_9 = C[((j+4)*lda)+i+0];
			register double cij_10 = C[((j+4)*lda)+i+1];
			register double cij_11 = C[((j+5)*lda)+i+0];
			register double cij_12 = C[((j+5)*lda)+i+1];
			register double cij_13 = C[((j+6)*lda)+i+0];
			register double cij_14 = C[((j+6)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
			C[((j+2)*lda)+i+0] = cij_5 ;
			C[((j+2)*lda)+i+1] = cij_6 ;
			C[((j+3)*lda)+i+0] = cij_7 ;
			C[((j+3)*lda)+i+1] = cij_8 ;
			C[((j+4)*lda)+i+0] = cij_9 ;
			C[((j+4)*lda)+i+1] = cij_10 ;
			C[((j+5)*lda)+i+0] = cij_11 ;
			C[((j+5)*lda)+i+1] = cij_12 ;
			C[((j+6)*lda)+i+0] = cij_13 ;
			C[((j+6)*lda)+i+1] = cij_14 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_2x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 2;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+1)*lda)+i+0];
			register double cij_4 = C[((j+1)*lda)+i+1];
			register double cij_5 = C[((j+2)*lda)+i+0];
			register double cij_6 = C[((j+2)*lda)+i+1];
			register double cij_7 = C[((j+3)*lda)+i+0];
			register double cij_8 = C[((j+3)*lda)+i+1];
			register double cij_9 = C[((j+4)*lda)+i+0];
			register double cij_10 = C[((j+4)*lda)+i+1];
			register double cij_11 = C[((j+5)*lda)+i+0];
			register double cij_12 = C[((j+5)*lda)+i+1];
			register double cij_13 = C[((j+6)*lda)+i+0];
			register double cij_14 = C[((j+6)*lda)+i+1];
			register double cij_15 = C[((j+7)*lda)+i+0];
			register double cij_16 = C[((j+7)*lda)+i+1];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_4 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+1)*lda)+i+0] = cij_3 ;
			C[((j+1)*lda)+i+1] = cij_4 ;
			C[((j+2)*lda)+i+0] = cij_5 ;
			C[((j+2)*lda)+i+1] = cij_6 ;
			C[((j+3)*lda)+i+0] = cij_7 ;
			C[((j+3)*lda)+i+1] = cij_8 ;
			C[((j+4)*lda)+i+0] = cij_9 ;
			C[((j+4)*lda)+i+1] = cij_10 ;
			C[((j+5)*lda)+i+0] = cij_11 ;
			C[((j+5)*lda)+i+1] = cij_12 ;
			C[((j+6)*lda)+i+0] = cij_13 ;
			C[((j+6)*lda)+i+1] = cij_14 ;
			C[((j+7)*lda)+i+0] = cij_15 ;
			C[((j+7)*lda)+i+1] = cij_16 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			register double cij_7 = C[((j+2)*lda)+i+0];
			register double cij_8 = C[((j+2)*lda)+i+1];
			register double cij_9 = C[((j+2)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
			C[((j+2)*lda)+i+0] = cij_7 ;
			C[((j+2)*lda)+i+1] = cij_8 ;
			C[((j+2)*lda)+i+2] = cij_9 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			register double cij_7 = C[((j+2)*lda)+i+0];
			register double cij_8 = C[((j+2)*lda)+i+1];
			register double cij_9 = C[((j+2)*lda)+i+2];
			register double cij_10 = C[((j+3)*lda)+i+0];
			register double cij_11 = C[((j+3)*lda)+i+1];
			register double cij_12 = C[((j+3)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_11 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_12 += Ai_[k+2*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
			C[((j+2)*lda)+i+0] = cij_7 ;
			C[((j+2)*lda)+i+1] = cij_8 ;
			C[((j+2)*lda)+i+2] = cij_9 ;
			C[((j+3)*lda)+i+0] = cij_10 ;
			C[((j+3)*lda)+i+1] = cij_11 ;
			C[((j+3)*lda)+i+2] = cij_12 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			register double cij_7 = C[((j+2)*lda)+i+0];
			register double cij_8 = C[((j+2)*lda)+i+1];
			register double cij_9 = C[((j+2)*lda)+i+2];
			register double cij_10 = C[((j+3)*lda)+i+0];
			register double cij_11 = C[((j+3)*lda)+i+1];
			register double cij_12 = C[((j+3)*lda)+i+2];
			register double cij_13 = C[((j+4)*lda)+i+0];
			register double cij_14 = C[((j+4)*lda)+i+1];
			register double cij_15 = C[((j+4)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_11 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_12 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
			C[((j+2)*lda)+i+0] = cij_7 ;
			C[((j+2)*lda)+i+1] = cij_8 ;
			C[((j+2)*lda)+i+2] = cij_9 ;
			C[((j+3)*lda)+i+0] = cij_10 ;
			C[((j+3)*lda)+i+1] = cij_11 ;
			C[((j+3)*lda)+i+2] = cij_12 ;
			C[((j+4)*lda)+i+0] = cij_13 ;
			C[((j+4)*lda)+i+1] = cij_14 ;
			C[((j+4)*lda)+i+2] = cij_15 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			register double cij_7 = C[((j+2)*lda)+i+0];
			register double cij_8 = C[((j+2)*lda)+i+1];
			register double cij_9 = C[((j+2)*lda)+i+2];
			register double cij_10 = C[((j+3)*lda)+i+0];
			register double cij_11 = C[((j+3)*lda)+i+1];
			register double cij_12 = C[((j+3)*lda)+i+2];
			register double cij_13 = C[((j+4)*lda)+i+0];
			register double cij_14 = C[((j+4)*lda)+i+1];
			register double cij_15 = C[((j+4)*lda)+i+2];
			register double cij_16 = C[((j+5)*lda)+i+0];
			register double cij_17 = C[((j+5)*lda)+i+1];
			register double cij_18 = C[((j+5)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_11 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_12 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
			C[((j+2)*lda)+i+0] = cij_7 ;
			C[((j+2)*lda)+i+1] = cij_8 ;
			C[((j+2)*lda)+i+2] = cij_9 ;
			C[((j+3)*lda)+i+0] = cij_10 ;
			C[((j+3)*lda)+i+1] = cij_11 ;
			C[((j+3)*lda)+i+2] = cij_12 ;
			C[((j+4)*lda)+i+0] = cij_13 ;
			C[((j+4)*lda)+i+1] = cij_14 ;
			C[((j+4)*lda)+i+2] = cij_15 ;
			C[((j+5)*lda)+i+0] = cij_16 ;
			C[((j+5)*lda)+i+1] = cij_17 ;
			C[((j+5)*lda)+i+2] = cij_18 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			register double cij_7 = C[((j+2)*lda)+i+0];
			register double cij_8 = C[((j+2)*lda)+i+1];
			register double cij_9 = C[((j+2)*lda)+i+2];
			register double cij_10 = C[((j+3)*lda)+i+0];
			register double cij_11 = C[((j+3)*lda)+i+1];
			register double cij_12 = C[((j+3)*lda)+i+2];
			register double cij_13 = C[((j+4)*lda)+i+0];
			register double cij_14 = C[((j+4)*lda)+i+1];
			register double cij_15 = C[((j+4)*lda)+i+2];
			register double cij_16 = C[((j+5)*lda)+i+0];
			register double cij_17 = C[((j+5)*lda)+i+1];
			register double cij_18 = C[((j+5)*lda)+i+2];
			register double cij_19 = C[((j+6)*lda)+i+0];
			register double cij_20 = C[((j+6)*lda)+i+1];
			register double cij_21 = C[((j+6)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_11 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_12 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
			C[((j+2)*lda)+i+0] = cij_7 ;
			C[((j+2)*lda)+i+1] = cij_8 ;
			C[((j+2)*lda)+i+2] = cij_9 ;
			C[((j+3)*lda)+i+0] = cij_10 ;
			C[((j+3)*lda)+i+1] = cij_11 ;
			C[((j+3)*lda)+i+2] = cij_12 ;
			C[((j+4)*lda)+i+0] = cij_13 ;
			C[((j+4)*lda)+i+1] = cij_14 ;
			C[((j+4)*lda)+i+2] = cij_15 ;
			C[((j+5)*lda)+i+0] = cij_16 ;
			C[((j+5)*lda)+i+1] = cij_17 ;
			C[((j+5)*lda)+i+2] = cij_18 ;
			C[((j+6)*lda)+i+0] = cij_19 ;
			C[((j+6)*lda)+i+1] = cij_20 ;
			C[((j+6)*lda)+i+2] = cij_21 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_3x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 3;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+1)*lda)+i+0];
			register double cij_5 = C[((j+1)*lda)+i+1];
			register double cij_6 = C[((j+1)*lda)+i+2];
			register double cij_7 = C[((j+2)*lda)+i+0];
			register double cij_8 = C[((j+2)*lda)+i+1];
			register double cij_9 = C[((j+2)*lda)+i+2];
			register double cij_10 = C[((j+3)*lda)+i+0];
			register double cij_11 = C[((j+3)*lda)+i+1];
			register double cij_12 = C[((j+3)*lda)+i+2];
			register double cij_13 = C[((j+4)*lda)+i+0];
			register double cij_14 = C[((j+4)*lda)+i+1];
			register double cij_15 = C[((j+4)*lda)+i+2];
			register double cij_16 = C[((j+5)*lda)+i+0];
			register double cij_17 = C[((j+5)*lda)+i+1];
			register double cij_18 = C[((j+5)*lda)+i+2];
			register double cij_19 = C[((j+6)*lda)+i+0];
			register double cij_20 = C[((j+6)*lda)+i+1];
			register double cij_21 = C[((j+6)*lda)+i+2];
			register double cij_22 = C[((j+7)*lda)+i+0];
			register double cij_23 = C[((j+7)*lda)+i+1];
			register double cij_24 = C[((j+7)*lda)+i+2];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_5 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_11 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_12 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_22 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_23 += Ai_[k+1*lda]*B_j[k+7*lda];
				cij_24 += Ai_[k+2*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+1)*lda)+i+0] = cij_4 ;
			C[((j+1)*lda)+i+1] = cij_5 ;
			C[((j+1)*lda)+i+2] = cij_6 ;
			C[((j+2)*lda)+i+0] = cij_7 ;
			C[((j+2)*lda)+i+1] = cij_8 ;
			C[((j+2)*lda)+i+2] = cij_9 ;
			C[((j+3)*lda)+i+0] = cij_10 ;
			C[((j+3)*lda)+i+1] = cij_11 ;
			C[((j+3)*lda)+i+2] = cij_12 ;
			C[((j+4)*lda)+i+0] = cij_13 ;
			C[((j+4)*lda)+i+1] = cij_14 ;
			C[((j+4)*lda)+i+2] = cij_15 ;
			C[((j+5)*lda)+i+0] = cij_16 ;
			C[((j+5)*lda)+i+1] = cij_17 ;
			C[((j+5)*lda)+i+2] = cij_18 ;
			C[((j+6)*lda)+i+0] = cij_19 ;
			C[((j+6)*lda)+i+1] = cij_20 ;
			C[((j+6)*lda)+i+2] = cij_21 ;
			C[((j+7)*lda)+i+0] = cij_22 ;
			C[((j+7)*lda)+i+1] = cij_23 ;
			C[((j+7)*lda)+i+2] = cij_24 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			register double cij_9 = C[((j+2)*lda)+i+0];
			register double cij_10 = C[((j+2)*lda)+i+1];
			register double cij_11 = C[((j+2)*lda)+i+2];
			register double cij_12 = C[((j+2)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
			C[((j+2)*lda)+i+0] = cij_9 ;
			C[((j+2)*lda)+i+1] = cij_10 ;
			C[((j+2)*lda)+i+2] = cij_11 ;
			C[((j+2)*lda)+i+3] = cij_12 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			register double cij_9 = C[((j+2)*lda)+i+0];
			register double cij_10 = C[((j+2)*lda)+i+1];
			register double cij_11 = C[((j+2)*lda)+i+2];
			register double cij_12 = C[((j+2)*lda)+i+3];
			register double cij_13 = C[((j+3)*lda)+i+0];
			register double cij_14 = C[((j+3)*lda)+i+1];
			register double cij_15 = C[((j+3)*lda)+i+2];
			register double cij_16 = C[((j+3)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
			C[((j+2)*lda)+i+0] = cij_9 ;
			C[((j+2)*lda)+i+1] = cij_10 ;
			C[((j+2)*lda)+i+2] = cij_11 ;
			C[((j+2)*lda)+i+3] = cij_12 ;
			C[((j+3)*lda)+i+0] = cij_13 ;
			C[((j+3)*lda)+i+1] = cij_14 ;
			C[((j+3)*lda)+i+2] = cij_15 ;
			C[((j+3)*lda)+i+3] = cij_16 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			register double cij_9 = C[((j+2)*lda)+i+0];
			register double cij_10 = C[((j+2)*lda)+i+1];
			register double cij_11 = C[((j+2)*lda)+i+2];
			register double cij_12 = C[((j+2)*lda)+i+3];
			register double cij_13 = C[((j+3)*lda)+i+0];
			register double cij_14 = C[((j+3)*lda)+i+1];
			register double cij_15 = C[((j+3)*lda)+i+2];
			register double cij_16 = C[((j+3)*lda)+i+3];
			register double cij_17 = C[((j+4)*lda)+i+0];
			register double cij_18 = C[((j+4)*lda)+i+1];
			register double cij_19 = C[((j+4)*lda)+i+2];
			register double cij_20 = C[((j+4)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
			C[((j+2)*lda)+i+0] = cij_9 ;
			C[((j+2)*lda)+i+1] = cij_10 ;
			C[((j+2)*lda)+i+2] = cij_11 ;
			C[((j+2)*lda)+i+3] = cij_12 ;
			C[((j+3)*lda)+i+0] = cij_13 ;
			C[((j+3)*lda)+i+1] = cij_14 ;
			C[((j+3)*lda)+i+2] = cij_15 ;
			C[((j+3)*lda)+i+3] = cij_16 ;
			C[((j+4)*lda)+i+0] = cij_17 ;
			C[((j+4)*lda)+i+1] = cij_18 ;
			C[((j+4)*lda)+i+2] = cij_19 ;
			C[((j+4)*lda)+i+3] = cij_20 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			register double cij_9 = C[((j+2)*lda)+i+0];
			register double cij_10 = C[((j+2)*lda)+i+1];
			register double cij_11 = C[((j+2)*lda)+i+2];
			register double cij_12 = C[((j+2)*lda)+i+3];
			register double cij_13 = C[((j+3)*lda)+i+0];
			register double cij_14 = C[((j+3)*lda)+i+1];
			register double cij_15 = C[((j+3)*lda)+i+2];
			register double cij_16 = C[((j+3)*lda)+i+3];
			register double cij_17 = C[((j+4)*lda)+i+0];
			register double cij_18 = C[((j+4)*lda)+i+1];
			register double cij_19 = C[((j+4)*lda)+i+2];
			register double cij_20 = C[((j+4)*lda)+i+3];
			register double cij_21 = C[((j+5)*lda)+i+0];
			register double cij_22 = C[((j+5)*lda)+i+1];
			register double cij_23 = C[((j+5)*lda)+i+2];
			register double cij_24 = C[((j+5)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
			C[((j+2)*lda)+i+0] = cij_9 ;
			C[((j+2)*lda)+i+1] = cij_10 ;
			C[((j+2)*lda)+i+2] = cij_11 ;
			C[((j+2)*lda)+i+3] = cij_12 ;
			C[((j+3)*lda)+i+0] = cij_13 ;
			C[((j+3)*lda)+i+1] = cij_14 ;
			C[((j+3)*lda)+i+2] = cij_15 ;
			C[((j+3)*lda)+i+3] = cij_16 ;
			C[((j+4)*lda)+i+0] = cij_17 ;
			C[((j+4)*lda)+i+1] = cij_18 ;
			C[((j+4)*lda)+i+2] = cij_19 ;
			C[((j+4)*lda)+i+3] = cij_20 ;
			C[((j+5)*lda)+i+0] = cij_21 ;
			C[((j+5)*lda)+i+1] = cij_22 ;
			C[((j+5)*lda)+i+2] = cij_23 ;
			C[((j+5)*lda)+i+3] = cij_24 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			register double cij_9 = C[((j+2)*lda)+i+0];
			register double cij_10 = C[((j+2)*lda)+i+1];
			register double cij_11 = C[((j+2)*lda)+i+2];
			register double cij_12 = C[((j+2)*lda)+i+3];
			register double cij_13 = C[((j+3)*lda)+i+0];
			register double cij_14 = C[((j+3)*lda)+i+1];
			register double cij_15 = C[((j+3)*lda)+i+2];
			register double cij_16 = C[((j+3)*lda)+i+3];
			register double cij_17 = C[((j+4)*lda)+i+0];
			register double cij_18 = C[((j+4)*lda)+i+1];
			register double cij_19 = C[((j+4)*lda)+i+2];
			register double cij_20 = C[((j+4)*lda)+i+3];
			register double cij_21 = C[((j+5)*lda)+i+0];
			register double cij_22 = C[((j+5)*lda)+i+1];
			register double cij_23 = C[((j+5)*lda)+i+2];
			register double cij_24 = C[((j+5)*lda)+i+3];
			register double cij_25 = C[((j+6)*lda)+i+0];
			register double cij_26 = C[((j+6)*lda)+i+1];
			register double cij_27 = C[((j+6)*lda)+i+2];
			register double cij_28 = C[((j+6)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
			C[((j+2)*lda)+i+0] = cij_9 ;
			C[((j+2)*lda)+i+1] = cij_10 ;
			C[((j+2)*lda)+i+2] = cij_11 ;
			C[((j+2)*lda)+i+3] = cij_12 ;
			C[((j+3)*lda)+i+0] = cij_13 ;
			C[((j+3)*lda)+i+1] = cij_14 ;
			C[((j+3)*lda)+i+2] = cij_15 ;
			C[((j+3)*lda)+i+3] = cij_16 ;
			C[((j+4)*lda)+i+0] = cij_17 ;
			C[((j+4)*lda)+i+1] = cij_18 ;
			C[((j+4)*lda)+i+2] = cij_19 ;
			C[((j+4)*lda)+i+3] = cij_20 ;
			C[((j+5)*lda)+i+0] = cij_21 ;
			C[((j+5)*lda)+i+1] = cij_22 ;
			C[((j+5)*lda)+i+2] = cij_23 ;
			C[((j+5)*lda)+i+3] = cij_24 ;
			C[((j+6)*lda)+i+0] = cij_25 ;
			C[((j+6)*lda)+i+1] = cij_26 ;
			C[((j+6)*lda)+i+2] = cij_27 ;
			C[((j+6)*lda)+i+3] = cij_28 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_4x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 4;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+1)*lda)+i+0];
			register double cij_6 = C[((j+1)*lda)+i+1];
			register double cij_7 = C[((j+1)*lda)+i+2];
			register double cij_8 = C[((j+1)*lda)+i+3];
			register double cij_9 = C[((j+2)*lda)+i+0];
			register double cij_10 = C[((j+2)*lda)+i+1];
			register double cij_11 = C[((j+2)*lda)+i+2];
			register double cij_12 = C[((j+2)*lda)+i+3];
			register double cij_13 = C[((j+3)*lda)+i+0];
			register double cij_14 = C[((j+3)*lda)+i+1];
			register double cij_15 = C[((j+3)*lda)+i+2];
			register double cij_16 = C[((j+3)*lda)+i+3];
			register double cij_17 = C[((j+4)*lda)+i+0];
			register double cij_18 = C[((j+4)*lda)+i+1];
			register double cij_19 = C[((j+4)*lda)+i+2];
			register double cij_20 = C[((j+4)*lda)+i+3];
			register double cij_21 = C[((j+5)*lda)+i+0];
			register double cij_22 = C[((j+5)*lda)+i+1];
			register double cij_23 = C[((j+5)*lda)+i+2];
			register double cij_24 = C[((j+5)*lda)+i+3];
			register double cij_25 = C[((j+6)*lda)+i+0];
			register double cij_26 = C[((j+6)*lda)+i+1];
			register double cij_27 = C[((j+6)*lda)+i+2];
			register double cij_28 = C[((j+6)*lda)+i+3];
			register double cij_29 = C[((j+7)*lda)+i+0];
			register double cij_30 = C[((j+7)*lda)+i+1];
			register double cij_31 = C[((j+7)*lda)+i+2];
			register double cij_32 = C[((j+7)*lda)+i+3];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_6 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_29 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_30 += Ai_[k+1*lda]*B_j[k+7*lda];
				cij_31 += Ai_[k+2*lda]*B_j[k+7*lda];
				cij_32 += Ai_[k+3*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+1)*lda)+i+0] = cij_5 ;
			C[((j+1)*lda)+i+1] = cij_6 ;
			C[((j+1)*lda)+i+2] = cij_7 ;
			C[((j+1)*lda)+i+3] = cij_8 ;
			C[((j+2)*lda)+i+0] = cij_9 ;
			C[((j+2)*lda)+i+1] = cij_10 ;
			C[((j+2)*lda)+i+2] = cij_11 ;
			C[((j+2)*lda)+i+3] = cij_12 ;
			C[((j+3)*lda)+i+0] = cij_13 ;
			C[((j+3)*lda)+i+1] = cij_14 ;
			C[((j+3)*lda)+i+2] = cij_15 ;
			C[((j+3)*lda)+i+3] = cij_16 ;
			C[((j+4)*lda)+i+0] = cij_17 ;
			C[((j+4)*lda)+i+1] = cij_18 ;
			C[((j+4)*lda)+i+2] = cij_19 ;
			C[((j+4)*lda)+i+3] = cij_20 ;
			C[((j+5)*lda)+i+0] = cij_21 ;
			C[((j+5)*lda)+i+1] = cij_22 ;
			C[((j+5)*lda)+i+2] = cij_23 ;
			C[((j+5)*lda)+i+3] = cij_24 ;
			C[((j+6)*lda)+i+0] = cij_25 ;
			C[((j+6)*lda)+i+1] = cij_26 ;
			C[((j+6)*lda)+i+2] = cij_27 ;
			C[((j+6)*lda)+i+3] = cij_28 ;
			C[((j+7)*lda)+i+0] = cij_29 ;
			C[((j+7)*lda)+i+1] = cij_30 ;
			C[((j+7)*lda)+i+2] = cij_31 ;
			C[((j+7)*lda)+i+3] = cij_32 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			register double cij_11 = C[((j+2)*lda)+i+0];
			register double cij_12 = C[((j+2)*lda)+i+1];
			register double cij_13 = C[((j+2)*lda)+i+2];
			register double cij_14 = C[((j+2)*lda)+i+3];
			register double cij_15 = C[((j+2)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+4*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
			C[((j+2)*lda)+i+0] = cij_11 ;
			C[((j+2)*lda)+i+1] = cij_12 ;
			C[((j+2)*lda)+i+2] = cij_13 ;
			C[((j+2)*lda)+i+3] = cij_14 ;
			C[((j+2)*lda)+i+4] = cij_15 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			register double cij_11 = C[((j+2)*lda)+i+0];
			register double cij_12 = C[((j+2)*lda)+i+1];
			register double cij_13 = C[((j+2)*lda)+i+2];
			register double cij_14 = C[((j+2)*lda)+i+3];
			register double cij_15 = C[((j+2)*lda)+i+4];
			register double cij_16 = C[((j+3)*lda)+i+0];
			register double cij_17 = C[((j+3)*lda)+i+1];
			register double cij_18 = C[((j+3)*lda)+i+2];
			register double cij_19 = C[((j+3)*lda)+i+3];
			register double cij_20 = C[((j+3)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_19 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+4*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
			C[((j+2)*lda)+i+0] = cij_11 ;
			C[((j+2)*lda)+i+1] = cij_12 ;
			C[((j+2)*lda)+i+2] = cij_13 ;
			C[((j+2)*lda)+i+3] = cij_14 ;
			C[((j+2)*lda)+i+4] = cij_15 ;
			C[((j+3)*lda)+i+0] = cij_16 ;
			C[((j+3)*lda)+i+1] = cij_17 ;
			C[((j+3)*lda)+i+2] = cij_18 ;
			C[((j+3)*lda)+i+3] = cij_19 ;
			C[((j+3)*lda)+i+4] = cij_20 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			register double cij_11 = C[((j+2)*lda)+i+0];
			register double cij_12 = C[((j+2)*lda)+i+1];
			register double cij_13 = C[((j+2)*lda)+i+2];
			register double cij_14 = C[((j+2)*lda)+i+3];
			register double cij_15 = C[((j+2)*lda)+i+4];
			register double cij_16 = C[((j+3)*lda)+i+0];
			register double cij_17 = C[((j+3)*lda)+i+1];
			register double cij_18 = C[((j+3)*lda)+i+2];
			register double cij_19 = C[((j+3)*lda)+i+3];
			register double cij_20 = C[((j+3)*lda)+i+4];
			register double cij_21 = C[((j+4)*lda)+i+0];
			register double cij_22 = C[((j+4)*lda)+i+1];
			register double cij_23 = C[((j+4)*lda)+i+2];
			register double cij_24 = C[((j+4)*lda)+i+3];
			register double cij_25 = C[((j+4)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_19 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_25 += Ai_[k+4*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
			C[((j+2)*lda)+i+0] = cij_11 ;
			C[((j+2)*lda)+i+1] = cij_12 ;
			C[((j+2)*lda)+i+2] = cij_13 ;
			C[((j+2)*lda)+i+3] = cij_14 ;
			C[((j+2)*lda)+i+4] = cij_15 ;
			C[((j+3)*lda)+i+0] = cij_16 ;
			C[((j+3)*lda)+i+1] = cij_17 ;
			C[((j+3)*lda)+i+2] = cij_18 ;
			C[((j+3)*lda)+i+3] = cij_19 ;
			C[((j+3)*lda)+i+4] = cij_20 ;
			C[((j+4)*lda)+i+0] = cij_21 ;
			C[((j+4)*lda)+i+1] = cij_22 ;
			C[((j+4)*lda)+i+2] = cij_23 ;
			C[((j+4)*lda)+i+3] = cij_24 ;
			C[((j+4)*lda)+i+4] = cij_25 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			register double cij_11 = C[((j+2)*lda)+i+0];
			register double cij_12 = C[((j+2)*lda)+i+1];
			register double cij_13 = C[((j+2)*lda)+i+2];
			register double cij_14 = C[((j+2)*lda)+i+3];
			register double cij_15 = C[((j+2)*lda)+i+4];
			register double cij_16 = C[((j+3)*lda)+i+0];
			register double cij_17 = C[((j+3)*lda)+i+1];
			register double cij_18 = C[((j+3)*lda)+i+2];
			register double cij_19 = C[((j+3)*lda)+i+3];
			register double cij_20 = C[((j+3)*lda)+i+4];
			register double cij_21 = C[((j+4)*lda)+i+0];
			register double cij_22 = C[((j+4)*lda)+i+1];
			register double cij_23 = C[((j+4)*lda)+i+2];
			register double cij_24 = C[((j+4)*lda)+i+3];
			register double cij_25 = C[((j+4)*lda)+i+4];
			register double cij_26 = C[((j+5)*lda)+i+0];
			register double cij_27 = C[((j+5)*lda)+i+1];
			register double cij_28 = C[((j+5)*lda)+i+2];
			register double cij_29 = C[((j+5)*lda)+i+3];
			register double cij_30 = C[((j+5)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_19 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_25 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_27 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_28 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_29 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_30 += Ai_[k+4*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
			C[((j+2)*lda)+i+0] = cij_11 ;
			C[((j+2)*lda)+i+1] = cij_12 ;
			C[((j+2)*lda)+i+2] = cij_13 ;
			C[((j+2)*lda)+i+3] = cij_14 ;
			C[((j+2)*lda)+i+4] = cij_15 ;
			C[((j+3)*lda)+i+0] = cij_16 ;
			C[((j+3)*lda)+i+1] = cij_17 ;
			C[((j+3)*lda)+i+2] = cij_18 ;
			C[((j+3)*lda)+i+3] = cij_19 ;
			C[((j+3)*lda)+i+4] = cij_20 ;
			C[((j+4)*lda)+i+0] = cij_21 ;
			C[((j+4)*lda)+i+1] = cij_22 ;
			C[((j+4)*lda)+i+2] = cij_23 ;
			C[((j+4)*lda)+i+3] = cij_24 ;
			C[((j+4)*lda)+i+4] = cij_25 ;
			C[((j+5)*lda)+i+0] = cij_26 ;
			C[((j+5)*lda)+i+1] = cij_27 ;
			C[((j+5)*lda)+i+2] = cij_28 ;
			C[((j+5)*lda)+i+3] = cij_29 ;
			C[((j+5)*lda)+i+4] = cij_30 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			register double cij_11 = C[((j+2)*lda)+i+0];
			register double cij_12 = C[((j+2)*lda)+i+1];
			register double cij_13 = C[((j+2)*lda)+i+2];
			register double cij_14 = C[((j+2)*lda)+i+3];
			register double cij_15 = C[((j+2)*lda)+i+4];
			register double cij_16 = C[((j+3)*lda)+i+0];
			register double cij_17 = C[((j+3)*lda)+i+1];
			register double cij_18 = C[((j+3)*lda)+i+2];
			register double cij_19 = C[((j+3)*lda)+i+3];
			register double cij_20 = C[((j+3)*lda)+i+4];
			register double cij_21 = C[((j+4)*lda)+i+0];
			register double cij_22 = C[((j+4)*lda)+i+1];
			register double cij_23 = C[((j+4)*lda)+i+2];
			register double cij_24 = C[((j+4)*lda)+i+3];
			register double cij_25 = C[((j+4)*lda)+i+4];
			register double cij_26 = C[((j+5)*lda)+i+0];
			register double cij_27 = C[((j+5)*lda)+i+1];
			register double cij_28 = C[((j+5)*lda)+i+2];
			register double cij_29 = C[((j+5)*lda)+i+3];
			register double cij_30 = C[((j+5)*lda)+i+4];
			register double cij_31 = C[((j+6)*lda)+i+0];
			register double cij_32 = C[((j+6)*lda)+i+1];
			register double cij_33 = C[((j+6)*lda)+i+2];
			register double cij_34 = C[((j+6)*lda)+i+3];
			register double cij_35 = C[((j+6)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_19 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_25 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_27 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_28 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_29 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_30 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_31 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_32 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_33 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_34 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_35 += Ai_[k+4*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
			C[((j+2)*lda)+i+0] = cij_11 ;
			C[((j+2)*lda)+i+1] = cij_12 ;
			C[((j+2)*lda)+i+2] = cij_13 ;
			C[((j+2)*lda)+i+3] = cij_14 ;
			C[((j+2)*lda)+i+4] = cij_15 ;
			C[((j+3)*lda)+i+0] = cij_16 ;
			C[((j+3)*lda)+i+1] = cij_17 ;
			C[((j+3)*lda)+i+2] = cij_18 ;
			C[((j+3)*lda)+i+3] = cij_19 ;
			C[((j+3)*lda)+i+4] = cij_20 ;
			C[((j+4)*lda)+i+0] = cij_21 ;
			C[((j+4)*lda)+i+1] = cij_22 ;
			C[((j+4)*lda)+i+2] = cij_23 ;
			C[((j+4)*lda)+i+3] = cij_24 ;
			C[((j+4)*lda)+i+4] = cij_25 ;
			C[((j+5)*lda)+i+0] = cij_26 ;
			C[((j+5)*lda)+i+1] = cij_27 ;
			C[((j+5)*lda)+i+2] = cij_28 ;
			C[((j+5)*lda)+i+3] = cij_29 ;
			C[((j+5)*lda)+i+4] = cij_30 ;
			C[((j+6)*lda)+i+0] = cij_31 ;
			C[((j+6)*lda)+i+1] = cij_32 ;
			C[((j+6)*lda)+i+2] = cij_33 ;
			C[((j+6)*lda)+i+3] = cij_34 ;
			C[((j+6)*lda)+i+4] = cij_35 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_5x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 5;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+1)*lda)+i+0];
			register double cij_7 = C[((j+1)*lda)+i+1];
			register double cij_8 = C[((j+1)*lda)+i+2];
			register double cij_9 = C[((j+1)*lda)+i+3];
			register double cij_10 = C[((j+1)*lda)+i+4];
			register double cij_11 = C[((j+2)*lda)+i+0];
			register double cij_12 = C[((j+2)*lda)+i+1];
			register double cij_13 = C[((j+2)*lda)+i+2];
			register double cij_14 = C[((j+2)*lda)+i+3];
			register double cij_15 = C[((j+2)*lda)+i+4];
			register double cij_16 = C[((j+3)*lda)+i+0];
			register double cij_17 = C[((j+3)*lda)+i+1];
			register double cij_18 = C[((j+3)*lda)+i+2];
			register double cij_19 = C[((j+3)*lda)+i+3];
			register double cij_20 = C[((j+3)*lda)+i+4];
			register double cij_21 = C[((j+4)*lda)+i+0];
			register double cij_22 = C[((j+4)*lda)+i+1];
			register double cij_23 = C[((j+4)*lda)+i+2];
			register double cij_24 = C[((j+4)*lda)+i+3];
			register double cij_25 = C[((j+4)*lda)+i+4];
			register double cij_26 = C[((j+5)*lda)+i+0];
			register double cij_27 = C[((j+5)*lda)+i+1];
			register double cij_28 = C[((j+5)*lda)+i+2];
			register double cij_29 = C[((j+5)*lda)+i+3];
			register double cij_30 = C[((j+5)*lda)+i+4];
			register double cij_31 = C[((j+6)*lda)+i+0];
			register double cij_32 = C[((j+6)*lda)+i+1];
			register double cij_33 = C[((j+6)*lda)+i+2];
			register double cij_34 = C[((j+6)*lda)+i+3];
			register double cij_35 = C[((j+6)*lda)+i+4];
			register double cij_36 = C[((j+7)*lda)+i+0];
			register double cij_37 = C[((j+7)*lda)+i+1];
			register double cij_38 = C[((j+7)*lda)+i+2];
			register double cij_39 = C[((j+7)*lda)+i+3];
			register double cij_40 = C[((j+7)*lda)+i+4];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_7 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_12 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_13 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_17 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_18 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_19 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_22 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_23 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_24 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_25 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_27 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_28 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_29 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_30 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_31 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_32 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_33 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_34 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_35 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_36 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_37 += Ai_[k+1*lda]*B_j[k+7*lda];
				cij_38 += Ai_[k+2*lda]*B_j[k+7*lda];
				cij_39 += Ai_[k+3*lda]*B_j[k+7*lda];
				cij_40 += Ai_[k+4*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+1)*lda)+i+0] = cij_6 ;
			C[((j+1)*lda)+i+1] = cij_7 ;
			C[((j+1)*lda)+i+2] = cij_8 ;
			C[((j+1)*lda)+i+3] = cij_9 ;
			C[((j+1)*lda)+i+4] = cij_10 ;
			C[((j+2)*lda)+i+0] = cij_11 ;
			C[((j+2)*lda)+i+1] = cij_12 ;
			C[((j+2)*lda)+i+2] = cij_13 ;
			C[((j+2)*lda)+i+3] = cij_14 ;
			C[((j+2)*lda)+i+4] = cij_15 ;
			C[((j+3)*lda)+i+0] = cij_16 ;
			C[((j+3)*lda)+i+1] = cij_17 ;
			C[((j+3)*lda)+i+2] = cij_18 ;
			C[((j+3)*lda)+i+3] = cij_19 ;
			C[((j+3)*lda)+i+4] = cij_20 ;
			C[((j+4)*lda)+i+0] = cij_21 ;
			C[((j+4)*lda)+i+1] = cij_22 ;
			C[((j+4)*lda)+i+2] = cij_23 ;
			C[((j+4)*lda)+i+3] = cij_24 ;
			C[((j+4)*lda)+i+4] = cij_25 ;
			C[((j+5)*lda)+i+0] = cij_26 ;
			C[((j+5)*lda)+i+1] = cij_27 ;
			C[((j+5)*lda)+i+2] = cij_28 ;
			C[((j+5)*lda)+i+3] = cij_29 ;
			C[((j+5)*lda)+i+4] = cij_30 ;
			C[((j+6)*lda)+i+0] = cij_31 ;
			C[((j+6)*lda)+i+1] = cij_32 ;
			C[((j+6)*lda)+i+2] = cij_33 ;
			C[((j+6)*lda)+i+3] = cij_34 ;
			C[((j+6)*lda)+i+4] = cij_35 ;
			C[((j+7)*lda)+i+0] = cij_36 ;
			C[((j+7)*lda)+i+1] = cij_37 ;
			C[((j+7)*lda)+i+2] = cij_38 ;
			C[((j+7)*lda)+i+3] = cij_39 ;
			C[((j+7)*lda)+i+4] = cij_40 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			register double cij_13 = C[((j+2)*lda)+i+0];
			register double cij_14 = C[((j+2)*lda)+i+1];
			register double cij_15 = C[((j+2)*lda)+i+2];
			register double cij_16 = C[((j+2)*lda)+i+3];
			register double cij_17 = C[((j+2)*lda)+i+4];
			register double cij_18 = C[((j+2)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+5*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
			C[((j+2)*lda)+i+0] = cij_13 ;
			C[((j+2)*lda)+i+1] = cij_14 ;
			C[((j+2)*lda)+i+2] = cij_15 ;
			C[((j+2)*lda)+i+3] = cij_16 ;
			C[((j+2)*lda)+i+4] = cij_17 ;
			C[((j+2)*lda)+i+5] = cij_18 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			register double cij_13 = C[((j+2)*lda)+i+0];
			register double cij_14 = C[((j+2)*lda)+i+1];
			register double cij_15 = C[((j+2)*lda)+i+2];
			register double cij_16 = C[((j+2)*lda)+i+3];
			register double cij_17 = C[((j+2)*lda)+i+4];
			register double cij_18 = C[((j+2)*lda)+i+5];
			register double cij_19 = C[((j+3)*lda)+i+0];
			register double cij_20 = C[((j+3)*lda)+i+1];
			register double cij_21 = C[((j+3)*lda)+i+2];
			register double cij_22 = C[((j+3)*lda)+i+3];
			register double cij_23 = C[((j+3)*lda)+i+4];
			register double cij_24 = C[((j+3)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_22 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+5*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
			C[((j+2)*lda)+i+0] = cij_13 ;
			C[((j+2)*lda)+i+1] = cij_14 ;
			C[((j+2)*lda)+i+2] = cij_15 ;
			C[((j+2)*lda)+i+3] = cij_16 ;
			C[((j+2)*lda)+i+4] = cij_17 ;
			C[((j+2)*lda)+i+5] = cij_18 ;
			C[((j+3)*lda)+i+0] = cij_19 ;
			C[((j+3)*lda)+i+1] = cij_20 ;
			C[((j+3)*lda)+i+2] = cij_21 ;
			C[((j+3)*lda)+i+3] = cij_22 ;
			C[((j+3)*lda)+i+4] = cij_23 ;
			C[((j+3)*lda)+i+5] = cij_24 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			register double cij_13 = C[((j+2)*lda)+i+0];
			register double cij_14 = C[((j+2)*lda)+i+1];
			register double cij_15 = C[((j+2)*lda)+i+2];
			register double cij_16 = C[((j+2)*lda)+i+3];
			register double cij_17 = C[((j+2)*lda)+i+4];
			register double cij_18 = C[((j+2)*lda)+i+5];
			register double cij_19 = C[((j+3)*lda)+i+0];
			register double cij_20 = C[((j+3)*lda)+i+1];
			register double cij_21 = C[((j+3)*lda)+i+2];
			register double cij_22 = C[((j+3)*lda)+i+3];
			register double cij_23 = C[((j+3)*lda)+i+4];
			register double cij_24 = C[((j+3)*lda)+i+5];
			register double cij_25 = C[((j+4)*lda)+i+0];
			register double cij_26 = C[((j+4)*lda)+i+1];
			register double cij_27 = C[((j+4)*lda)+i+2];
			register double cij_28 = C[((j+4)*lda)+i+3];
			register double cij_29 = C[((j+4)*lda)+i+4];
			register double cij_30 = C[((j+4)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_22 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
			C[((j+2)*lda)+i+0] = cij_13 ;
			C[((j+2)*lda)+i+1] = cij_14 ;
			C[((j+2)*lda)+i+2] = cij_15 ;
			C[((j+2)*lda)+i+3] = cij_16 ;
			C[((j+2)*lda)+i+4] = cij_17 ;
			C[((j+2)*lda)+i+5] = cij_18 ;
			C[((j+3)*lda)+i+0] = cij_19 ;
			C[((j+3)*lda)+i+1] = cij_20 ;
			C[((j+3)*lda)+i+2] = cij_21 ;
			C[((j+3)*lda)+i+3] = cij_22 ;
			C[((j+3)*lda)+i+4] = cij_23 ;
			C[((j+3)*lda)+i+5] = cij_24 ;
			C[((j+4)*lda)+i+0] = cij_25 ;
			C[((j+4)*lda)+i+1] = cij_26 ;
			C[((j+4)*lda)+i+2] = cij_27 ;
			C[((j+4)*lda)+i+3] = cij_28 ;
			C[((j+4)*lda)+i+4] = cij_29 ;
			C[((j+4)*lda)+i+5] = cij_30 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			register double cij_13 = C[((j+2)*lda)+i+0];
			register double cij_14 = C[((j+2)*lda)+i+1];
			register double cij_15 = C[((j+2)*lda)+i+2];
			register double cij_16 = C[((j+2)*lda)+i+3];
			register double cij_17 = C[((j+2)*lda)+i+4];
			register double cij_18 = C[((j+2)*lda)+i+5];
			register double cij_19 = C[((j+3)*lda)+i+0];
			register double cij_20 = C[((j+3)*lda)+i+1];
			register double cij_21 = C[((j+3)*lda)+i+2];
			register double cij_22 = C[((j+3)*lda)+i+3];
			register double cij_23 = C[((j+3)*lda)+i+4];
			register double cij_24 = C[((j+3)*lda)+i+5];
			register double cij_25 = C[((j+4)*lda)+i+0];
			register double cij_26 = C[((j+4)*lda)+i+1];
			register double cij_27 = C[((j+4)*lda)+i+2];
			register double cij_28 = C[((j+4)*lda)+i+3];
			register double cij_29 = C[((j+4)*lda)+i+4];
			register double cij_30 = C[((j+4)*lda)+i+5];
			register double cij_31 = C[((j+5)*lda)+i+0];
			register double cij_32 = C[((j+5)*lda)+i+1];
			register double cij_33 = C[((j+5)*lda)+i+2];
			register double cij_34 = C[((j+5)*lda)+i+3];
			register double cij_35 = C[((j+5)*lda)+i+4];
			register double cij_36 = C[((j+5)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_22 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_32 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_33 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_34 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_35 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_36 += Ai_[k+5*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
			C[((j+2)*lda)+i+0] = cij_13 ;
			C[((j+2)*lda)+i+1] = cij_14 ;
			C[((j+2)*lda)+i+2] = cij_15 ;
			C[((j+2)*lda)+i+3] = cij_16 ;
			C[((j+2)*lda)+i+4] = cij_17 ;
			C[((j+2)*lda)+i+5] = cij_18 ;
			C[((j+3)*lda)+i+0] = cij_19 ;
			C[((j+3)*lda)+i+1] = cij_20 ;
			C[((j+3)*lda)+i+2] = cij_21 ;
			C[((j+3)*lda)+i+3] = cij_22 ;
			C[((j+3)*lda)+i+4] = cij_23 ;
			C[((j+3)*lda)+i+5] = cij_24 ;
			C[((j+4)*lda)+i+0] = cij_25 ;
			C[((j+4)*lda)+i+1] = cij_26 ;
			C[((j+4)*lda)+i+2] = cij_27 ;
			C[((j+4)*lda)+i+3] = cij_28 ;
			C[((j+4)*lda)+i+4] = cij_29 ;
			C[((j+4)*lda)+i+5] = cij_30 ;
			C[((j+5)*lda)+i+0] = cij_31 ;
			C[((j+5)*lda)+i+1] = cij_32 ;
			C[((j+5)*lda)+i+2] = cij_33 ;
			C[((j+5)*lda)+i+3] = cij_34 ;
			C[((j+5)*lda)+i+4] = cij_35 ;
			C[((j+5)*lda)+i+5] = cij_36 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			register double cij_13 = C[((j+2)*lda)+i+0];
			register double cij_14 = C[((j+2)*lda)+i+1];
			register double cij_15 = C[((j+2)*lda)+i+2];
			register double cij_16 = C[((j+2)*lda)+i+3];
			register double cij_17 = C[((j+2)*lda)+i+4];
			register double cij_18 = C[((j+2)*lda)+i+5];
			register double cij_19 = C[((j+3)*lda)+i+0];
			register double cij_20 = C[((j+3)*lda)+i+1];
			register double cij_21 = C[((j+3)*lda)+i+2];
			register double cij_22 = C[((j+3)*lda)+i+3];
			register double cij_23 = C[((j+3)*lda)+i+4];
			register double cij_24 = C[((j+3)*lda)+i+5];
			register double cij_25 = C[((j+4)*lda)+i+0];
			register double cij_26 = C[((j+4)*lda)+i+1];
			register double cij_27 = C[((j+4)*lda)+i+2];
			register double cij_28 = C[((j+4)*lda)+i+3];
			register double cij_29 = C[((j+4)*lda)+i+4];
			register double cij_30 = C[((j+4)*lda)+i+5];
			register double cij_31 = C[((j+5)*lda)+i+0];
			register double cij_32 = C[((j+5)*lda)+i+1];
			register double cij_33 = C[((j+5)*lda)+i+2];
			register double cij_34 = C[((j+5)*lda)+i+3];
			register double cij_35 = C[((j+5)*lda)+i+4];
			register double cij_36 = C[((j+5)*lda)+i+5];
			register double cij_37 = C[((j+6)*lda)+i+0];
			register double cij_38 = C[((j+6)*lda)+i+1];
			register double cij_39 = C[((j+6)*lda)+i+2];
			register double cij_40 = C[((j+6)*lda)+i+3];
			register double cij_41 = C[((j+6)*lda)+i+4];
			register double cij_42 = C[((j+6)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_22 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_32 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_33 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_34 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_35 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_36 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_37 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_38 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_39 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_40 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_41 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_42 += Ai_[k+5*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
			C[((j+2)*lda)+i+0] = cij_13 ;
			C[((j+2)*lda)+i+1] = cij_14 ;
			C[((j+2)*lda)+i+2] = cij_15 ;
			C[((j+2)*lda)+i+3] = cij_16 ;
			C[((j+2)*lda)+i+4] = cij_17 ;
			C[((j+2)*lda)+i+5] = cij_18 ;
			C[((j+3)*lda)+i+0] = cij_19 ;
			C[((j+3)*lda)+i+1] = cij_20 ;
			C[((j+3)*lda)+i+2] = cij_21 ;
			C[((j+3)*lda)+i+3] = cij_22 ;
			C[((j+3)*lda)+i+4] = cij_23 ;
			C[((j+3)*lda)+i+5] = cij_24 ;
			C[((j+4)*lda)+i+0] = cij_25 ;
			C[((j+4)*lda)+i+1] = cij_26 ;
			C[((j+4)*lda)+i+2] = cij_27 ;
			C[((j+4)*lda)+i+3] = cij_28 ;
			C[((j+4)*lda)+i+4] = cij_29 ;
			C[((j+4)*lda)+i+5] = cij_30 ;
			C[((j+5)*lda)+i+0] = cij_31 ;
			C[((j+5)*lda)+i+1] = cij_32 ;
			C[((j+5)*lda)+i+2] = cij_33 ;
			C[((j+5)*lda)+i+3] = cij_34 ;
			C[((j+5)*lda)+i+4] = cij_35 ;
			C[((j+5)*lda)+i+5] = cij_36 ;
			C[((j+6)*lda)+i+0] = cij_37 ;
			C[((j+6)*lda)+i+1] = cij_38 ;
			C[((j+6)*lda)+i+2] = cij_39 ;
			C[((j+6)*lda)+i+3] = cij_40 ;
			C[((j+6)*lda)+i+4] = cij_41 ;
			C[((j+6)*lda)+i+5] = cij_42 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_6x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 6;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+1)*lda)+i+0];
			register double cij_8 = C[((j+1)*lda)+i+1];
			register double cij_9 = C[((j+1)*lda)+i+2];
			register double cij_10 = C[((j+1)*lda)+i+3];
			register double cij_11 = C[((j+1)*lda)+i+4];
			register double cij_12 = C[((j+1)*lda)+i+5];
			register double cij_13 = C[((j+2)*lda)+i+0];
			register double cij_14 = C[((j+2)*lda)+i+1];
			register double cij_15 = C[((j+2)*lda)+i+2];
			register double cij_16 = C[((j+2)*lda)+i+3];
			register double cij_17 = C[((j+2)*lda)+i+4];
			register double cij_18 = C[((j+2)*lda)+i+5];
			register double cij_19 = C[((j+3)*lda)+i+0];
			register double cij_20 = C[((j+3)*lda)+i+1];
			register double cij_21 = C[((j+3)*lda)+i+2];
			register double cij_22 = C[((j+3)*lda)+i+3];
			register double cij_23 = C[((j+3)*lda)+i+4];
			register double cij_24 = C[((j+3)*lda)+i+5];
			register double cij_25 = C[((j+4)*lda)+i+0];
			register double cij_26 = C[((j+4)*lda)+i+1];
			register double cij_27 = C[((j+4)*lda)+i+2];
			register double cij_28 = C[((j+4)*lda)+i+3];
			register double cij_29 = C[((j+4)*lda)+i+4];
			register double cij_30 = C[((j+4)*lda)+i+5];
			register double cij_31 = C[((j+5)*lda)+i+0];
			register double cij_32 = C[((j+5)*lda)+i+1];
			register double cij_33 = C[((j+5)*lda)+i+2];
			register double cij_34 = C[((j+5)*lda)+i+3];
			register double cij_35 = C[((j+5)*lda)+i+4];
			register double cij_36 = C[((j+5)*lda)+i+5];
			register double cij_37 = C[((j+6)*lda)+i+0];
			register double cij_38 = C[((j+6)*lda)+i+1];
			register double cij_39 = C[((j+6)*lda)+i+2];
			register double cij_40 = C[((j+6)*lda)+i+3];
			register double cij_41 = C[((j+6)*lda)+i+4];
			register double cij_42 = C[((j+6)*lda)+i+5];
			register double cij_43 = C[((j+7)*lda)+i+0];
			register double cij_44 = C[((j+7)*lda)+i+1];
			register double cij_45 = C[((j+7)*lda)+i+2];
			register double cij_46 = C[((j+7)*lda)+i+3];
			register double cij_47 = C[((j+7)*lda)+i+4];
			register double cij_48 = C[((j+7)*lda)+i+5];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_8 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_14 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_15 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_20 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_21 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_22 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_32 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_33 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_34 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_35 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_36 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_37 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_38 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_39 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_40 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_41 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_42 += Ai_[k+5*lda]*B_j[k+6*lda];
				cij_43 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_44 += Ai_[k+1*lda]*B_j[k+7*lda];
				cij_45 += Ai_[k+2*lda]*B_j[k+7*lda];
				cij_46 += Ai_[k+3*lda]*B_j[k+7*lda];
				cij_47 += Ai_[k+4*lda]*B_j[k+7*lda];
				cij_48 += Ai_[k+5*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+1)*lda)+i+0] = cij_7 ;
			C[((j+1)*lda)+i+1] = cij_8 ;
			C[((j+1)*lda)+i+2] = cij_9 ;
			C[((j+1)*lda)+i+3] = cij_10 ;
			C[((j+1)*lda)+i+4] = cij_11 ;
			C[((j+1)*lda)+i+5] = cij_12 ;
			C[((j+2)*lda)+i+0] = cij_13 ;
			C[((j+2)*lda)+i+1] = cij_14 ;
			C[((j+2)*lda)+i+2] = cij_15 ;
			C[((j+2)*lda)+i+3] = cij_16 ;
			C[((j+2)*lda)+i+4] = cij_17 ;
			C[((j+2)*lda)+i+5] = cij_18 ;
			C[((j+3)*lda)+i+0] = cij_19 ;
			C[((j+3)*lda)+i+1] = cij_20 ;
			C[((j+3)*lda)+i+2] = cij_21 ;
			C[((j+3)*lda)+i+3] = cij_22 ;
			C[((j+3)*lda)+i+4] = cij_23 ;
			C[((j+3)*lda)+i+5] = cij_24 ;
			C[((j+4)*lda)+i+0] = cij_25 ;
			C[((j+4)*lda)+i+1] = cij_26 ;
			C[((j+4)*lda)+i+2] = cij_27 ;
			C[((j+4)*lda)+i+3] = cij_28 ;
			C[((j+4)*lda)+i+4] = cij_29 ;
			C[((j+4)*lda)+i+5] = cij_30 ;
			C[((j+5)*lda)+i+0] = cij_31 ;
			C[((j+5)*lda)+i+1] = cij_32 ;
			C[((j+5)*lda)+i+2] = cij_33 ;
			C[((j+5)*lda)+i+3] = cij_34 ;
			C[((j+5)*lda)+i+4] = cij_35 ;
			C[((j+5)*lda)+i+5] = cij_36 ;
			C[((j+6)*lda)+i+0] = cij_37 ;
			C[((j+6)*lda)+i+1] = cij_38 ;
			C[((j+6)*lda)+i+2] = cij_39 ;
			C[((j+6)*lda)+i+3] = cij_40 ;
			C[((j+6)*lda)+i+4] = cij_41 ;
			C[((j+6)*lda)+i+5] = cij_42 ;
			C[((j+7)*lda)+i+0] = cij_43 ;
			C[((j+7)*lda)+i+1] = cij_44 ;
			C[((j+7)*lda)+i+2] = cij_45 ;
			C[((j+7)*lda)+i+3] = cij_46 ;
			C[((j+7)*lda)+i+4] = cij_47 ;
			C[((j+7)*lda)+i+5] = cij_48 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			register double cij_15 = C[((j+2)*lda)+i+0];
			register double cij_16 = C[((j+2)*lda)+i+1];
			register double cij_17 = C[((j+2)*lda)+i+2];
			register double cij_18 = C[((j+2)*lda)+i+3];
			register double cij_19 = C[((j+2)*lda)+i+4];
			register double cij_20 = C[((j+2)*lda)+i+5];
			register double cij_21 = C[((j+2)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+6*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
			C[((j+2)*lda)+i+0] = cij_15 ;
			C[((j+2)*lda)+i+1] = cij_16 ;
			C[((j+2)*lda)+i+2] = cij_17 ;
			C[((j+2)*lda)+i+3] = cij_18 ;
			C[((j+2)*lda)+i+4] = cij_19 ;
			C[((j+2)*lda)+i+5] = cij_20 ;
			C[((j+2)*lda)+i+6] = cij_21 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			register double cij_15 = C[((j+2)*lda)+i+0];
			register double cij_16 = C[((j+2)*lda)+i+1];
			register double cij_17 = C[((j+2)*lda)+i+2];
			register double cij_18 = C[((j+2)*lda)+i+3];
			register double cij_19 = C[((j+2)*lda)+i+4];
			register double cij_20 = C[((j+2)*lda)+i+5];
			register double cij_21 = C[((j+2)*lda)+i+6];
			register double cij_22 = C[((j+3)*lda)+i+0];
			register double cij_23 = C[((j+3)*lda)+i+1];
			register double cij_24 = C[((j+3)*lda)+i+2];
			register double cij_25 = C[((j+3)*lda)+i+3];
			register double cij_26 = C[((j+3)*lda)+i+4];
			register double cij_27 = C[((j+3)*lda)+i+5];
			register double cij_28 = C[((j+3)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+6*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
			C[((j+2)*lda)+i+0] = cij_15 ;
			C[((j+2)*lda)+i+1] = cij_16 ;
			C[((j+2)*lda)+i+2] = cij_17 ;
			C[((j+2)*lda)+i+3] = cij_18 ;
			C[((j+2)*lda)+i+4] = cij_19 ;
			C[((j+2)*lda)+i+5] = cij_20 ;
			C[((j+2)*lda)+i+6] = cij_21 ;
			C[((j+3)*lda)+i+0] = cij_22 ;
			C[((j+3)*lda)+i+1] = cij_23 ;
			C[((j+3)*lda)+i+2] = cij_24 ;
			C[((j+3)*lda)+i+3] = cij_25 ;
			C[((j+3)*lda)+i+4] = cij_26 ;
			C[((j+3)*lda)+i+5] = cij_27 ;
			C[((j+3)*lda)+i+6] = cij_28 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			register double cij_15 = C[((j+2)*lda)+i+0];
			register double cij_16 = C[((j+2)*lda)+i+1];
			register double cij_17 = C[((j+2)*lda)+i+2];
			register double cij_18 = C[((j+2)*lda)+i+3];
			register double cij_19 = C[((j+2)*lda)+i+4];
			register double cij_20 = C[((j+2)*lda)+i+5];
			register double cij_21 = C[((j+2)*lda)+i+6];
			register double cij_22 = C[((j+3)*lda)+i+0];
			register double cij_23 = C[((j+3)*lda)+i+1];
			register double cij_24 = C[((j+3)*lda)+i+2];
			register double cij_25 = C[((j+3)*lda)+i+3];
			register double cij_26 = C[((j+3)*lda)+i+4];
			register double cij_27 = C[((j+3)*lda)+i+5];
			register double cij_28 = C[((j+3)*lda)+i+6];
			register double cij_29 = C[((j+4)*lda)+i+0];
			register double cij_30 = C[((j+4)*lda)+i+1];
			register double cij_31 = C[((j+4)*lda)+i+2];
			register double cij_32 = C[((j+4)*lda)+i+3];
			register double cij_33 = C[((j+4)*lda)+i+4];
			register double cij_34 = C[((j+4)*lda)+i+5];
			register double cij_35 = C[((j+4)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_32 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_33 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+6*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
			C[((j+2)*lda)+i+0] = cij_15 ;
			C[((j+2)*lda)+i+1] = cij_16 ;
			C[((j+2)*lda)+i+2] = cij_17 ;
			C[((j+2)*lda)+i+3] = cij_18 ;
			C[((j+2)*lda)+i+4] = cij_19 ;
			C[((j+2)*lda)+i+5] = cij_20 ;
			C[((j+2)*lda)+i+6] = cij_21 ;
			C[((j+3)*lda)+i+0] = cij_22 ;
			C[((j+3)*lda)+i+1] = cij_23 ;
			C[((j+3)*lda)+i+2] = cij_24 ;
			C[((j+3)*lda)+i+3] = cij_25 ;
			C[((j+3)*lda)+i+4] = cij_26 ;
			C[((j+3)*lda)+i+5] = cij_27 ;
			C[((j+3)*lda)+i+6] = cij_28 ;
			C[((j+4)*lda)+i+0] = cij_29 ;
			C[((j+4)*lda)+i+1] = cij_30 ;
			C[((j+4)*lda)+i+2] = cij_31 ;
			C[((j+4)*lda)+i+3] = cij_32 ;
			C[((j+4)*lda)+i+4] = cij_33 ;
			C[((j+4)*lda)+i+5] = cij_34 ;
			C[((j+4)*lda)+i+6] = cij_35 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			register double cij_15 = C[((j+2)*lda)+i+0];
			register double cij_16 = C[((j+2)*lda)+i+1];
			register double cij_17 = C[((j+2)*lda)+i+2];
			register double cij_18 = C[((j+2)*lda)+i+3];
			register double cij_19 = C[((j+2)*lda)+i+4];
			register double cij_20 = C[((j+2)*lda)+i+5];
			register double cij_21 = C[((j+2)*lda)+i+6];
			register double cij_22 = C[((j+3)*lda)+i+0];
			register double cij_23 = C[((j+3)*lda)+i+1];
			register double cij_24 = C[((j+3)*lda)+i+2];
			register double cij_25 = C[((j+3)*lda)+i+3];
			register double cij_26 = C[((j+3)*lda)+i+4];
			register double cij_27 = C[((j+3)*lda)+i+5];
			register double cij_28 = C[((j+3)*lda)+i+6];
			register double cij_29 = C[((j+4)*lda)+i+0];
			register double cij_30 = C[((j+4)*lda)+i+1];
			register double cij_31 = C[((j+4)*lda)+i+2];
			register double cij_32 = C[((j+4)*lda)+i+3];
			register double cij_33 = C[((j+4)*lda)+i+4];
			register double cij_34 = C[((j+4)*lda)+i+5];
			register double cij_35 = C[((j+4)*lda)+i+6];
			register double cij_36 = C[((j+5)*lda)+i+0];
			register double cij_37 = C[((j+5)*lda)+i+1];
			register double cij_38 = C[((j+5)*lda)+i+2];
			register double cij_39 = C[((j+5)*lda)+i+3];
			register double cij_40 = C[((j+5)*lda)+i+4];
			register double cij_41 = C[((j+5)*lda)+i+5];
			register double cij_42 = C[((j+5)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_32 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_33 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_37 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_38 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_39 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_40 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_41 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_42 += Ai_[k+6*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
			C[((j+2)*lda)+i+0] = cij_15 ;
			C[((j+2)*lda)+i+1] = cij_16 ;
			C[((j+2)*lda)+i+2] = cij_17 ;
			C[((j+2)*lda)+i+3] = cij_18 ;
			C[((j+2)*lda)+i+4] = cij_19 ;
			C[((j+2)*lda)+i+5] = cij_20 ;
			C[((j+2)*lda)+i+6] = cij_21 ;
			C[((j+3)*lda)+i+0] = cij_22 ;
			C[((j+3)*lda)+i+1] = cij_23 ;
			C[((j+3)*lda)+i+2] = cij_24 ;
			C[((j+3)*lda)+i+3] = cij_25 ;
			C[((j+3)*lda)+i+4] = cij_26 ;
			C[((j+3)*lda)+i+5] = cij_27 ;
			C[((j+3)*lda)+i+6] = cij_28 ;
			C[((j+4)*lda)+i+0] = cij_29 ;
			C[((j+4)*lda)+i+1] = cij_30 ;
			C[((j+4)*lda)+i+2] = cij_31 ;
			C[((j+4)*lda)+i+3] = cij_32 ;
			C[((j+4)*lda)+i+4] = cij_33 ;
			C[((j+4)*lda)+i+5] = cij_34 ;
			C[((j+4)*lda)+i+6] = cij_35 ;
			C[((j+5)*lda)+i+0] = cij_36 ;
			C[((j+5)*lda)+i+1] = cij_37 ;
			C[((j+5)*lda)+i+2] = cij_38 ;
			C[((j+5)*lda)+i+3] = cij_39 ;
			C[((j+5)*lda)+i+4] = cij_40 ;
			C[((j+5)*lda)+i+5] = cij_41 ;
			C[((j+5)*lda)+i+6] = cij_42 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			register double cij_15 = C[((j+2)*lda)+i+0];
			register double cij_16 = C[((j+2)*lda)+i+1];
			register double cij_17 = C[((j+2)*lda)+i+2];
			register double cij_18 = C[((j+2)*lda)+i+3];
			register double cij_19 = C[((j+2)*lda)+i+4];
			register double cij_20 = C[((j+2)*lda)+i+5];
			register double cij_21 = C[((j+2)*lda)+i+6];
			register double cij_22 = C[((j+3)*lda)+i+0];
			register double cij_23 = C[((j+3)*lda)+i+1];
			register double cij_24 = C[((j+3)*lda)+i+2];
			register double cij_25 = C[((j+3)*lda)+i+3];
			register double cij_26 = C[((j+3)*lda)+i+4];
			register double cij_27 = C[((j+3)*lda)+i+5];
			register double cij_28 = C[((j+3)*lda)+i+6];
			register double cij_29 = C[((j+4)*lda)+i+0];
			register double cij_30 = C[((j+4)*lda)+i+1];
			register double cij_31 = C[((j+4)*lda)+i+2];
			register double cij_32 = C[((j+4)*lda)+i+3];
			register double cij_33 = C[((j+4)*lda)+i+4];
			register double cij_34 = C[((j+4)*lda)+i+5];
			register double cij_35 = C[((j+4)*lda)+i+6];
			register double cij_36 = C[((j+5)*lda)+i+0];
			register double cij_37 = C[((j+5)*lda)+i+1];
			register double cij_38 = C[((j+5)*lda)+i+2];
			register double cij_39 = C[((j+5)*lda)+i+3];
			register double cij_40 = C[((j+5)*lda)+i+4];
			register double cij_41 = C[((j+5)*lda)+i+5];
			register double cij_42 = C[((j+5)*lda)+i+6];
			register double cij_43 = C[((j+6)*lda)+i+0];
			register double cij_44 = C[((j+6)*lda)+i+1];
			register double cij_45 = C[((j+6)*lda)+i+2];
			register double cij_46 = C[((j+6)*lda)+i+3];
			register double cij_47 = C[((j+6)*lda)+i+4];
			register double cij_48 = C[((j+6)*lda)+i+5];
			register double cij_49 = C[((j+6)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_32 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_33 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_37 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_38 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_39 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_40 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_41 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_42 += Ai_[k+6*lda]*B_j[k+5*lda];
				cij_43 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_44 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_45 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_46 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_47 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_48 += Ai_[k+5*lda]*B_j[k+6*lda];
				cij_49 += Ai_[k+6*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
			C[((j+2)*lda)+i+0] = cij_15 ;
			C[((j+2)*lda)+i+1] = cij_16 ;
			C[((j+2)*lda)+i+2] = cij_17 ;
			C[((j+2)*lda)+i+3] = cij_18 ;
			C[((j+2)*lda)+i+4] = cij_19 ;
			C[((j+2)*lda)+i+5] = cij_20 ;
			C[((j+2)*lda)+i+6] = cij_21 ;
			C[((j+3)*lda)+i+0] = cij_22 ;
			C[((j+3)*lda)+i+1] = cij_23 ;
			C[((j+3)*lda)+i+2] = cij_24 ;
			C[((j+3)*lda)+i+3] = cij_25 ;
			C[((j+3)*lda)+i+4] = cij_26 ;
			C[((j+3)*lda)+i+5] = cij_27 ;
			C[((j+3)*lda)+i+6] = cij_28 ;
			C[((j+4)*lda)+i+0] = cij_29 ;
			C[((j+4)*lda)+i+1] = cij_30 ;
			C[((j+4)*lda)+i+2] = cij_31 ;
			C[((j+4)*lda)+i+3] = cij_32 ;
			C[((j+4)*lda)+i+4] = cij_33 ;
			C[((j+4)*lda)+i+5] = cij_34 ;
			C[((j+4)*lda)+i+6] = cij_35 ;
			C[((j+5)*lda)+i+0] = cij_36 ;
			C[((j+5)*lda)+i+1] = cij_37 ;
			C[((j+5)*lda)+i+2] = cij_38 ;
			C[((j+5)*lda)+i+3] = cij_39 ;
			C[((j+5)*lda)+i+4] = cij_40 ;
			C[((j+5)*lda)+i+5] = cij_41 ;
			C[((j+5)*lda)+i+6] = cij_42 ;
			C[((j+6)*lda)+i+0] = cij_43 ;
			C[((j+6)*lda)+i+1] = cij_44 ;
			C[((j+6)*lda)+i+2] = cij_45 ;
			C[((j+6)*lda)+i+3] = cij_46 ;
			C[((j+6)*lda)+i+4] = cij_47 ;
			C[((j+6)*lda)+i+5] = cij_48 ;
			C[((j+6)*lda)+i+6] = cij_49 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_7x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 7;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+1)*lda)+i+0];
			register double cij_9 = C[((j+1)*lda)+i+1];
			register double cij_10 = C[((j+1)*lda)+i+2];
			register double cij_11 = C[((j+1)*lda)+i+3];
			register double cij_12 = C[((j+1)*lda)+i+4];
			register double cij_13 = C[((j+1)*lda)+i+5];
			register double cij_14 = C[((j+1)*lda)+i+6];
			register double cij_15 = C[((j+2)*lda)+i+0];
			register double cij_16 = C[((j+2)*lda)+i+1];
			register double cij_17 = C[((j+2)*lda)+i+2];
			register double cij_18 = C[((j+2)*lda)+i+3];
			register double cij_19 = C[((j+2)*lda)+i+4];
			register double cij_20 = C[((j+2)*lda)+i+5];
			register double cij_21 = C[((j+2)*lda)+i+6];
			register double cij_22 = C[((j+3)*lda)+i+0];
			register double cij_23 = C[((j+3)*lda)+i+1];
			register double cij_24 = C[((j+3)*lda)+i+2];
			register double cij_25 = C[((j+3)*lda)+i+3];
			register double cij_26 = C[((j+3)*lda)+i+4];
			register double cij_27 = C[((j+3)*lda)+i+5];
			register double cij_28 = C[((j+3)*lda)+i+6];
			register double cij_29 = C[((j+4)*lda)+i+0];
			register double cij_30 = C[((j+4)*lda)+i+1];
			register double cij_31 = C[((j+4)*lda)+i+2];
			register double cij_32 = C[((j+4)*lda)+i+3];
			register double cij_33 = C[((j+4)*lda)+i+4];
			register double cij_34 = C[((j+4)*lda)+i+5];
			register double cij_35 = C[((j+4)*lda)+i+6];
			register double cij_36 = C[((j+5)*lda)+i+0];
			register double cij_37 = C[((j+5)*lda)+i+1];
			register double cij_38 = C[((j+5)*lda)+i+2];
			register double cij_39 = C[((j+5)*lda)+i+3];
			register double cij_40 = C[((j+5)*lda)+i+4];
			register double cij_41 = C[((j+5)*lda)+i+5];
			register double cij_42 = C[((j+5)*lda)+i+6];
			register double cij_43 = C[((j+6)*lda)+i+0];
			register double cij_44 = C[((j+6)*lda)+i+1];
			register double cij_45 = C[((j+6)*lda)+i+2];
			register double cij_46 = C[((j+6)*lda)+i+3];
			register double cij_47 = C[((j+6)*lda)+i+4];
			register double cij_48 = C[((j+6)*lda)+i+5];
			register double cij_49 = C[((j+6)*lda)+i+6];
			register double cij_50 = C[((j+7)*lda)+i+0];
			register double cij_51 = C[((j+7)*lda)+i+1];
			register double cij_52 = C[((j+7)*lda)+i+2];
			register double cij_53 = C[((j+7)*lda)+i+3];
			register double cij_54 = C[((j+7)*lda)+i+4];
			register double cij_55 = C[((j+7)*lda)+i+5];
			register double cij_56 = C[((j+7)*lda)+i+6];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_9 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_16 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_17 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_23 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_24 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_25 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_30 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_31 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_32 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_33 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_37 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_38 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_39 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_40 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_41 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_42 += Ai_[k+6*lda]*B_j[k+5*lda];
				cij_43 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_44 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_45 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_46 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_47 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_48 += Ai_[k+5*lda]*B_j[k+6*lda];
				cij_49 += Ai_[k+6*lda]*B_j[k+6*lda];
				cij_50 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_51 += Ai_[k+1*lda]*B_j[k+7*lda];
				cij_52 += Ai_[k+2*lda]*B_j[k+7*lda];
				cij_53 += Ai_[k+3*lda]*B_j[k+7*lda];
				cij_54 += Ai_[k+4*lda]*B_j[k+7*lda];
				cij_55 += Ai_[k+5*lda]*B_j[k+7*lda];
				cij_56 += Ai_[k+6*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+1)*lda)+i+0] = cij_8 ;
			C[((j+1)*lda)+i+1] = cij_9 ;
			C[((j+1)*lda)+i+2] = cij_10 ;
			C[((j+1)*lda)+i+3] = cij_11 ;
			C[((j+1)*lda)+i+4] = cij_12 ;
			C[((j+1)*lda)+i+5] = cij_13 ;
			C[((j+1)*lda)+i+6] = cij_14 ;
			C[((j+2)*lda)+i+0] = cij_15 ;
			C[((j+2)*lda)+i+1] = cij_16 ;
			C[((j+2)*lda)+i+2] = cij_17 ;
			C[((j+2)*lda)+i+3] = cij_18 ;
			C[((j+2)*lda)+i+4] = cij_19 ;
			C[((j+2)*lda)+i+5] = cij_20 ;
			C[((j+2)*lda)+i+6] = cij_21 ;
			C[((j+3)*lda)+i+0] = cij_22 ;
			C[((j+3)*lda)+i+1] = cij_23 ;
			C[((j+3)*lda)+i+2] = cij_24 ;
			C[((j+3)*lda)+i+3] = cij_25 ;
			C[((j+3)*lda)+i+4] = cij_26 ;
			C[((j+3)*lda)+i+5] = cij_27 ;
			C[((j+3)*lda)+i+6] = cij_28 ;
			C[((j+4)*lda)+i+0] = cij_29 ;
			C[((j+4)*lda)+i+1] = cij_30 ;
			C[((j+4)*lda)+i+2] = cij_31 ;
			C[((j+4)*lda)+i+3] = cij_32 ;
			C[((j+4)*lda)+i+4] = cij_33 ;
			C[((j+4)*lda)+i+5] = cij_34 ;
			C[((j+4)*lda)+i+6] = cij_35 ;
			C[((j+5)*lda)+i+0] = cij_36 ;
			C[((j+5)*lda)+i+1] = cij_37 ;
			C[((j+5)*lda)+i+2] = cij_38 ;
			C[((j+5)*lda)+i+3] = cij_39 ;
			C[((j+5)*lda)+i+4] = cij_40 ;
			C[((j+5)*lda)+i+5] = cij_41 ;
			C[((j+5)*lda)+i+6] = cij_42 ;
			C[((j+6)*lda)+i+0] = cij_43 ;
			C[((j+6)*lda)+i+1] = cij_44 ;
			C[((j+6)*lda)+i+2] = cij_45 ;
			C[((j+6)*lda)+i+3] = cij_46 ;
			C[((j+6)*lda)+i+4] = cij_47 ;
			C[((j+6)*lda)+i+5] = cij_48 ;
			C[((j+6)*lda)+i+6] = cij_49 ;
			C[((j+7)*lda)+i+0] = cij_50 ;
			C[((j+7)*lda)+i+1] = cij_51 ;
			C[((j+7)*lda)+i+2] = cij_52 ;
			C[((j+7)*lda)+i+3] = cij_53 ;
			C[((j+7)*lda)+i+4] = cij_54 ;
			C[((j+7)*lda)+i+5] = cij_55 ;
			C[((j+7)*lda)+i+6] = cij_56 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x1x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 1;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x2x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 2;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x3x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 3;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			register double cij_17 = C[((j+2)*lda)+i+0];
			register double cij_18 = C[((j+2)*lda)+i+1];
			register double cij_19 = C[((j+2)*lda)+i+2];
			register double cij_20 = C[((j+2)*lda)+i+3];
			register double cij_21 = C[((j+2)*lda)+i+4];
			register double cij_22 = C[((j+2)*lda)+i+5];
			register double cij_23 = C[((j+2)*lda)+i+6];
			register double cij_24 = C[((j+2)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_23 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_24 += Ai_[k+7*lda]*B_j[k+2*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
			C[((j+2)*lda)+i+0] = cij_17 ;
			C[((j+2)*lda)+i+1] = cij_18 ;
			C[((j+2)*lda)+i+2] = cij_19 ;
			C[((j+2)*lda)+i+3] = cij_20 ;
			C[((j+2)*lda)+i+4] = cij_21 ;
			C[((j+2)*lda)+i+5] = cij_22 ;
			C[((j+2)*lda)+i+6] = cij_23 ;
			C[((j+2)*lda)+i+7] = cij_24 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x4x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 4;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			register double cij_17 = C[((j+2)*lda)+i+0];
			register double cij_18 = C[((j+2)*lda)+i+1];
			register double cij_19 = C[((j+2)*lda)+i+2];
			register double cij_20 = C[((j+2)*lda)+i+3];
			register double cij_21 = C[((j+2)*lda)+i+4];
			register double cij_22 = C[((j+2)*lda)+i+5];
			register double cij_23 = C[((j+2)*lda)+i+6];
			register double cij_24 = C[((j+2)*lda)+i+7];
			register double cij_25 = C[((j+3)*lda)+i+0];
			register double cij_26 = C[((j+3)*lda)+i+1];
			register double cij_27 = C[((j+3)*lda)+i+2];
			register double cij_28 = C[((j+3)*lda)+i+3];
			register double cij_29 = C[((j+3)*lda)+i+4];
			register double cij_30 = C[((j+3)*lda)+i+5];
			register double cij_31 = C[((j+3)*lda)+i+6];
			register double cij_32 = C[((j+3)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_23 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_24 += Ai_[k+7*lda]*B_j[k+2*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_31 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_32 += Ai_[k+7*lda]*B_j[k+3*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
			C[((j+2)*lda)+i+0] = cij_17 ;
			C[((j+2)*lda)+i+1] = cij_18 ;
			C[((j+2)*lda)+i+2] = cij_19 ;
			C[((j+2)*lda)+i+3] = cij_20 ;
			C[((j+2)*lda)+i+4] = cij_21 ;
			C[((j+2)*lda)+i+5] = cij_22 ;
			C[((j+2)*lda)+i+6] = cij_23 ;
			C[((j+2)*lda)+i+7] = cij_24 ;
			C[((j+3)*lda)+i+0] = cij_25 ;
			C[((j+3)*lda)+i+1] = cij_26 ;
			C[((j+3)*lda)+i+2] = cij_27 ;
			C[((j+3)*lda)+i+3] = cij_28 ;
			C[((j+3)*lda)+i+4] = cij_29 ;
			C[((j+3)*lda)+i+5] = cij_30 ;
			C[((j+3)*lda)+i+6] = cij_31 ;
			C[((j+3)*lda)+i+7] = cij_32 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x5x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 5;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			register double cij_17 = C[((j+2)*lda)+i+0];
			register double cij_18 = C[((j+2)*lda)+i+1];
			register double cij_19 = C[((j+2)*lda)+i+2];
			register double cij_20 = C[((j+2)*lda)+i+3];
			register double cij_21 = C[((j+2)*lda)+i+4];
			register double cij_22 = C[((j+2)*lda)+i+5];
			register double cij_23 = C[((j+2)*lda)+i+6];
			register double cij_24 = C[((j+2)*lda)+i+7];
			register double cij_25 = C[((j+3)*lda)+i+0];
			register double cij_26 = C[((j+3)*lda)+i+1];
			register double cij_27 = C[((j+3)*lda)+i+2];
			register double cij_28 = C[((j+3)*lda)+i+3];
			register double cij_29 = C[((j+3)*lda)+i+4];
			register double cij_30 = C[((j+3)*lda)+i+5];
			register double cij_31 = C[((j+3)*lda)+i+6];
			register double cij_32 = C[((j+3)*lda)+i+7];
			register double cij_33 = C[((j+4)*lda)+i+0];
			register double cij_34 = C[((j+4)*lda)+i+1];
			register double cij_35 = C[((j+4)*lda)+i+2];
			register double cij_36 = C[((j+4)*lda)+i+3];
			register double cij_37 = C[((j+4)*lda)+i+4];
			register double cij_38 = C[((j+4)*lda)+i+5];
			register double cij_39 = C[((j+4)*lda)+i+6];
			register double cij_40 = C[((j+4)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_23 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_24 += Ai_[k+7*lda]*B_j[k+2*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_31 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_32 += Ai_[k+7*lda]*B_j[k+3*lda];
				cij_33 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_37 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_38 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_39 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_40 += Ai_[k+7*lda]*B_j[k+4*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
			C[((j+2)*lda)+i+0] = cij_17 ;
			C[((j+2)*lda)+i+1] = cij_18 ;
			C[((j+2)*lda)+i+2] = cij_19 ;
			C[((j+2)*lda)+i+3] = cij_20 ;
			C[((j+2)*lda)+i+4] = cij_21 ;
			C[((j+2)*lda)+i+5] = cij_22 ;
			C[((j+2)*lda)+i+6] = cij_23 ;
			C[((j+2)*lda)+i+7] = cij_24 ;
			C[((j+3)*lda)+i+0] = cij_25 ;
			C[((j+3)*lda)+i+1] = cij_26 ;
			C[((j+3)*lda)+i+2] = cij_27 ;
			C[((j+3)*lda)+i+3] = cij_28 ;
			C[((j+3)*lda)+i+4] = cij_29 ;
			C[((j+3)*lda)+i+5] = cij_30 ;
			C[((j+3)*lda)+i+6] = cij_31 ;
			C[((j+3)*lda)+i+7] = cij_32 ;
			C[((j+4)*lda)+i+0] = cij_33 ;
			C[((j+4)*lda)+i+1] = cij_34 ;
			C[((j+4)*lda)+i+2] = cij_35 ;
			C[((j+4)*lda)+i+3] = cij_36 ;
			C[((j+4)*lda)+i+4] = cij_37 ;
			C[((j+4)*lda)+i+5] = cij_38 ;
			C[((j+4)*lda)+i+6] = cij_39 ;
			C[((j+4)*lda)+i+7] = cij_40 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x6x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 6;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			register double cij_17 = C[((j+2)*lda)+i+0];
			register double cij_18 = C[((j+2)*lda)+i+1];
			register double cij_19 = C[((j+2)*lda)+i+2];
			register double cij_20 = C[((j+2)*lda)+i+3];
			register double cij_21 = C[((j+2)*lda)+i+4];
			register double cij_22 = C[((j+2)*lda)+i+5];
			register double cij_23 = C[((j+2)*lda)+i+6];
			register double cij_24 = C[((j+2)*lda)+i+7];
			register double cij_25 = C[((j+3)*lda)+i+0];
			register double cij_26 = C[((j+3)*lda)+i+1];
			register double cij_27 = C[((j+3)*lda)+i+2];
			register double cij_28 = C[((j+3)*lda)+i+3];
			register double cij_29 = C[((j+3)*lda)+i+4];
			register double cij_30 = C[((j+3)*lda)+i+5];
			register double cij_31 = C[((j+3)*lda)+i+6];
			register double cij_32 = C[((j+3)*lda)+i+7];
			register double cij_33 = C[((j+4)*lda)+i+0];
			register double cij_34 = C[((j+4)*lda)+i+1];
			register double cij_35 = C[((j+4)*lda)+i+2];
			register double cij_36 = C[((j+4)*lda)+i+3];
			register double cij_37 = C[((j+4)*lda)+i+4];
			register double cij_38 = C[((j+4)*lda)+i+5];
			register double cij_39 = C[((j+4)*lda)+i+6];
			register double cij_40 = C[((j+4)*lda)+i+7];
			register double cij_41 = C[((j+5)*lda)+i+0];
			register double cij_42 = C[((j+5)*lda)+i+1];
			register double cij_43 = C[((j+5)*lda)+i+2];
			register double cij_44 = C[((j+5)*lda)+i+3];
			register double cij_45 = C[((j+5)*lda)+i+4];
			register double cij_46 = C[((j+5)*lda)+i+5];
			register double cij_47 = C[((j+5)*lda)+i+6];
			register double cij_48 = C[((j+5)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_23 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_24 += Ai_[k+7*lda]*B_j[k+2*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_31 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_32 += Ai_[k+7*lda]*B_j[k+3*lda];
				cij_33 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_37 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_38 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_39 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_40 += Ai_[k+7*lda]*B_j[k+4*lda];
				cij_41 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_42 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_43 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_44 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_45 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_46 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_47 += Ai_[k+6*lda]*B_j[k+5*lda];
				cij_48 += Ai_[k+7*lda]*B_j[k+5*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
			C[((j+2)*lda)+i+0] = cij_17 ;
			C[((j+2)*lda)+i+1] = cij_18 ;
			C[((j+2)*lda)+i+2] = cij_19 ;
			C[((j+2)*lda)+i+3] = cij_20 ;
			C[((j+2)*lda)+i+4] = cij_21 ;
			C[((j+2)*lda)+i+5] = cij_22 ;
			C[((j+2)*lda)+i+6] = cij_23 ;
			C[((j+2)*lda)+i+7] = cij_24 ;
			C[((j+3)*lda)+i+0] = cij_25 ;
			C[((j+3)*lda)+i+1] = cij_26 ;
			C[((j+3)*lda)+i+2] = cij_27 ;
			C[((j+3)*lda)+i+3] = cij_28 ;
			C[((j+3)*lda)+i+4] = cij_29 ;
			C[((j+3)*lda)+i+5] = cij_30 ;
			C[((j+3)*lda)+i+6] = cij_31 ;
			C[((j+3)*lda)+i+7] = cij_32 ;
			C[((j+4)*lda)+i+0] = cij_33 ;
			C[((j+4)*lda)+i+1] = cij_34 ;
			C[((j+4)*lda)+i+2] = cij_35 ;
			C[((j+4)*lda)+i+3] = cij_36 ;
			C[((j+4)*lda)+i+4] = cij_37 ;
			C[((j+4)*lda)+i+5] = cij_38 ;
			C[((j+4)*lda)+i+6] = cij_39 ;
			C[((j+4)*lda)+i+7] = cij_40 ;
			C[((j+5)*lda)+i+0] = cij_41 ;
			C[((j+5)*lda)+i+1] = cij_42 ;
			C[((j+5)*lda)+i+2] = cij_43 ;
			C[((j+5)*lda)+i+3] = cij_44 ;
			C[((j+5)*lda)+i+4] = cij_45 ;
			C[((j+5)*lda)+i+5] = cij_46 ;
			C[((j+5)*lda)+i+6] = cij_47 ;
			C[((j+5)*lda)+i+7] = cij_48 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x7x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 7;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			register double cij_17 = C[((j+2)*lda)+i+0];
			register double cij_18 = C[((j+2)*lda)+i+1];
			register double cij_19 = C[((j+2)*lda)+i+2];
			register double cij_20 = C[((j+2)*lda)+i+3];
			register double cij_21 = C[((j+2)*lda)+i+4];
			register double cij_22 = C[((j+2)*lda)+i+5];
			register double cij_23 = C[((j+2)*lda)+i+6];
			register double cij_24 = C[((j+2)*lda)+i+7];
			register double cij_25 = C[((j+3)*lda)+i+0];
			register double cij_26 = C[((j+3)*lda)+i+1];
			register double cij_27 = C[((j+3)*lda)+i+2];
			register double cij_28 = C[((j+3)*lda)+i+3];
			register double cij_29 = C[((j+3)*lda)+i+4];
			register double cij_30 = C[((j+3)*lda)+i+5];
			register double cij_31 = C[((j+3)*lda)+i+6];
			register double cij_32 = C[((j+3)*lda)+i+7];
			register double cij_33 = C[((j+4)*lda)+i+0];
			register double cij_34 = C[((j+4)*lda)+i+1];
			register double cij_35 = C[((j+4)*lda)+i+2];
			register double cij_36 = C[((j+4)*lda)+i+3];
			register double cij_37 = C[((j+4)*lda)+i+4];
			register double cij_38 = C[((j+4)*lda)+i+5];
			register double cij_39 = C[((j+4)*lda)+i+6];
			register double cij_40 = C[((j+4)*lda)+i+7];
			register double cij_41 = C[((j+5)*lda)+i+0];
			register double cij_42 = C[((j+5)*lda)+i+1];
			register double cij_43 = C[((j+5)*lda)+i+2];
			register double cij_44 = C[((j+5)*lda)+i+3];
			register double cij_45 = C[((j+5)*lda)+i+4];
			register double cij_46 = C[((j+5)*lda)+i+5];
			register double cij_47 = C[((j+5)*lda)+i+6];
			register double cij_48 = C[((j+5)*lda)+i+7];
			register double cij_49 = C[((j+6)*lda)+i+0];
			register double cij_50 = C[((j+6)*lda)+i+1];
			register double cij_51 = C[((j+6)*lda)+i+2];
			register double cij_52 = C[((j+6)*lda)+i+3];
			register double cij_53 = C[((j+6)*lda)+i+4];
			register double cij_54 = C[((j+6)*lda)+i+5];
			register double cij_55 = C[((j+6)*lda)+i+6];
			register double cij_56 = C[((j+6)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_23 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_24 += Ai_[k+7*lda]*B_j[k+2*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_31 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_32 += Ai_[k+7*lda]*B_j[k+3*lda];
				cij_33 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_37 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_38 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_39 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_40 += Ai_[k+7*lda]*B_j[k+4*lda];
				cij_41 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_42 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_43 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_44 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_45 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_46 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_47 += Ai_[k+6*lda]*B_j[k+5*lda];
				cij_48 += Ai_[k+7*lda]*B_j[k+5*lda];
				cij_49 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_50 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_51 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_52 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_53 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_54 += Ai_[k+5*lda]*B_j[k+6*lda];
				cij_55 += Ai_[k+6*lda]*B_j[k+6*lda];
				cij_56 += Ai_[k+7*lda]*B_j[k+6*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
			C[((j+2)*lda)+i+0] = cij_17 ;
			C[((j+2)*lda)+i+1] = cij_18 ;
			C[((j+2)*lda)+i+2] = cij_19 ;
			C[((j+2)*lda)+i+3] = cij_20 ;
			C[((j+2)*lda)+i+4] = cij_21 ;
			C[((j+2)*lda)+i+5] = cij_22 ;
			C[((j+2)*lda)+i+6] = cij_23 ;
			C[((j+2)*lda)+i+7] = cij_24 ;
			C[((j+3)*lda)+i+0] = cij_25 ;
			C[((j+3)*lda)+i+1] = cij_26 ;
			C[((j+3)*lda)+i+2] = cij_27 ;
			C[((j+3)*lda)+i+3] = cij_28 ;
			C[((j+3)*lda)+i+4] = cij_29 ;
			C[((j+3)*lda)+i+5] = cij_30 ;
			C[((j+3)*lda)+i+6] = cij_31 ;
			C[((j+3)*lda)+i+7] = cij_32 ;
			C[((j+4)*lda)+i+0] = cij_33 ;
			C[((j+4)*lda)+i+1] = cij_34 ;
			C[((j+4)*lda)+i+2] = cij_35 ;
			C[((j+4)*lda)+i+3] = cij_36 ;
			C[((j+4)*lda)+i+4] = cij_37 ;
			C[((j+4)*lda)+i+5] = cij_38 ;
			C[((j+4)*lda)+i+6] = cij_39 ;
			C[((j+4)*lda)+i+7] = cij_40 ;
			C[((j+5)*lda)+i+0] = cij_41 ;
			C[((j+5)*lda)+i+1] = cij_42 ;
			C[((j+5)*lda)+i+2] = cij_43 ;
			C[((j+5)*lda)+i+3] = cij_44 ;
			C[((j+5)*lda)+i+4] = cij_45 ;
			C[((j+5)*lda)+i+5] = cij_46 ;
			C[((j+5)*lda)+i+6] = cij_47 ;
			C[((j+5)*lda)+i+7] = cij_48 ;
			C[((j+6)*lda)+i+0] = cij_49 ;
			C[((j+6)*lda)+i+1] = cij_50 ;
			C[((j+6)*lda)+i+2] = cij_51 ;
			C[((j+6)*lda)+i+3] = cij_52 ;
			C[((j+6)*lda)+i+4] = cij_53 ;
			C[((j+6)*lda)+i+5] = cij_54 ;
			C[((j+6)*lda)+i+6] = cij_55 ;
			C[((j+6)*lda)+i+7] = cij_56 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}
void basic_dgemm_8x8x1 (const int lda, const int M, const int N, const int K, const double *A, const double *B, double *C) {
	int i, j, k;
	int bi, bj, bk;
	int num_m, num_n, num_k;
	int r = 8;
	int c = 8;
	num_m = M/r;
	num_n = N/c;
	num_k = K;
	for (bi=0; bi<num_m; bi++) {
		i = bi*r;
		const register double *Ai_ = A + i*lda;
		for (bj=0; bj<num_n; bj++) {
			j = bj*c;
			const register double *B_j = B + j*lda;
			register double cij_1 = C[((j+0)*lda)+i+0];
			register double cij_2 = C[((j+0)*lda)+i+1];
			register double cij_3 = C[((j+0)*lda)+i+2];
			register double cij_4 = C[((j+0)*lda)+i+3];
			register double cij_5 = C[((j+0)*lda)+i+4];
			register double cij_6 = C[((j+0)*lda)+i+5];
			register double cij_7 = C[((j+0)*lda)+i+6];
			register double cij_8 = C[((j+0)*lda)+i+7];
			register double cij_9 = C[((j+1)*lda)+i+0];
			register double cij_10 = C[((j+1)*lda)+i+1];
			register double cij_11 = C[((j+1)*lda)+i+2];
			register double cij_12 = C[((j+1)*lda)+i+3];
			register double cij_13 = C[((j+1)*lda)+i+4];
			register double cij_14 = C[((j+1)*lda)+i+5];
			register double cij_15 = C[((j+1)*lda)+i+6];
			register double cij_16 = C[((j+1)*lda)+i+7];
			register double cij_17 = C[((j+2)*lda)+i+0];
			register double cij_18 = C[((j+2)*lda)+i+1];
			register double cij_19 = C[((j+2)*lda)+i+2];
			register double cij_20 = C[((j+2)*lda)+i+3];
			register double cij_21 = C[((j+2)*lda)+i+4];
			register double cij_22 = C[((j+2)*lda)+i+5];
			register double cij_23 = C[((j+2)*lda)+i+6];
			register double cij_24 = C[((j+2)*lda)+i+7];
			register double cij_25 = C[((j+3)*lda)+i+0];
			register double cij_26 = C[((j+3)*lda)+i+1];
			register double cij_27 = C[((j+3)*lda)+i+2];
			register double cij_28 = C[((j+3)*lda)+i+3];
			register double cij_29 = C[((j+3)*lda)+i+4];
			register double cij_30 = C[((j+3)*lda)+i+5];
			register double cij_31 = C[((j+3)*lda)+i+6];
			register double cij_32 = C[((j+3)*lda)+i+7];
			register double cij_33 = C[((j+4)*lda)+i+0];
			register double cij_34 = C[((j+4)*lda)+i+1];
			register double cij_35 = C[((j+4)*lda)+i+2];
			register double cij_36 = C[((j+4)*lda)+i+3];
			register double cij_37 = C[((j+4)*lda)+i+4];
			register double cij_38 = C[((j+4)*lda)+i+5];
			register double cij_39 = C[((j+4)*lda)+i+6];
			register double cij_40 = C[((j+4)*lda)+i+7];
			register double cij_41 = C[((j+5)*lda)+i+0];
			register double cij_42 = C[((j+5)*lda)+i+1];
			register double cij_43 = C[((j+5)*lda)+i+2];
			register double cij_44 = C[((j+5)*lda)+i+3];
			register double cij_45 = C[((j+5)*lda)+i+4];
			register double cij_46 = C[((j+5)*lda)+i+5];
			register double cij_47 = C[((j+5)*lda)+i+6];
			register double cij_48 = C[((j+5)*lda)+i+7];
			register double cij_49 = C[((j+6)*lda)+i+0];
			register double cij_50 = C[((j+6)*lda)+i+1];
			register double cij_51 = C[((j+6)*lda)+i+2];
			register double cij_52 = C[((j+6)*lda)+i+3];
			register double cij_53 = C[((j+6)*lda)+i+4];
			register double cij_54 = C[((j+6)*lda)+i+5];
			register double cij_55 = C[((j+6)*lda)+i+6];
			register double cij_56 = C[((j+6)*lda)+i+7];
			register double cij_57 = C[((j+7)*lda)+i+0];
			register double cij_58 = C[((j+7)*lda)+i+1];
			register double cij_59 = C[((j+7)*lda)+i+2];
			register double cij_60 = C[((j+7)*lda)+i+3];
			register double cij_61 = C[((j+7)*lda)+i+4];
			register double cij_62 = C[((j+7)*lda)+i+5];
			register double cij_63 = C[((j+7)*lda)+i+6];
			register double cij_64 = C[((j+7)*lda)+i+7];
			for (k=0; k<K; ++k) {
				cij_1 += Ai_[k+0*lda]*B_j[k+0*lda];
				cij_2 += Ai_[k+1*lda]*B_j[k+0*lda];
				cij_3 += Ai_[k+2*lda]*B_j[k+0*lda];
				cij_4 += Ai_[k+3*lda]*B_j[k+0*lda];
				cij_5 += Ai_[k+4*lda]*B_j[k+0*lda];
				cij_6 += Ai_[k+5*lda]*B_j[k+0*lda];
				cij_7 += Ai_[k+6*lda]*B_j[k+0*lda];
				cij_8 += Ai_[k+7*lda]*B_j[k+0*lda];
				cij_9 += Ai_[k+0*lda]*B_j[k+1*lda];
				cij_10 += Ai_[k+1*lda]*B_j[k+1*lda];
				cij_11 += Ai_[k+2*lda]*B_j[k+1*lda];
				cij_12 += Ai_[k+3*lda]*B_j[k+1*lda];
				cij_13 += Ai_[k+4*lda]*B_j[k+1*lda];
				cij_14 += Ai_[k+5*lda]*B_j[k+1*lda];
				cij_15 += Ai_[k+6*lda]*B_j[k+1*lda];
				cij_16 += Ai_[k+7*lda]*B_j[k+1*lda];
				cij_17 += Ai_[k+0*lda]*B_j[k+2*lda];
				cij_18 += Ai_[k+1*lda]*B_j[k+2*lda];
				cij_19 += Ai_[k+2*lda]*B_j[k+2*lda];
				cij_20 += Ai_[k+3*lda]*B_j[k+2*lda];
				cij_21 += Ai_[k+4*lda]*B_j[k+2*lda];
				cij_22 += Ai_[k+5*lda]*B_j[k+2*lda];
				cij_23 += Ai_[k+6*lda]*B_j[k+2*lda];
				cij_24 += Ai_[k+7*lda]*B_j[k+2*lda];
				cij_25 += Ai_[k+0*lda]*B_j[k+3*lda];
				cij_26 += Ai_[k+1*lda]*B_j[k+3*lda];
				cij_27 += Ai_[k+2*lda]*B_j[k+3*lda];
				cij_28 += Ai_[k+3*lda]*B_j[k+3*lda];
				cij_29 += Ai_[k+4*lda]*B_j[k+3*lda];
				cij_30 += Ai_[k+5*lda]*B_j[k+3*lda];
				cij_31 += Ai_[k+6*lda]*B_j[k+3*lda];
				cij_32 += Ai_[k+7*lda]*B_j[k+3*lda];
				cij_33 += Ai_[k+0*lda]*B_j[k+4*lda];
				cij_34 += Ai_[k+1*lda]*B_j[k+4*lda];
				cij_35 += Ai_[k+2*lda]*B_j[k+4*lda];
				cij_36 += Ai_[k+3*lda]*B_j[k+4*lda];
				cij_37 += Ai_[k+4*lda]*B_j[k+4*lda];
				cij_38 += Ai_[k+5*lda]*B_j[k+4*lda];
				cij_39 += Ai_[k+6*lda]*B_j[k+4*lda];
				cij_40 += Ai_[k+7*lda]*B_j[k+4*lda];
				cij_41 += Ai_[k+0*lda]*B_j[k+5*lda];
				cij_42 += Ai_[k+1*lda]*B_j[k+5*lda];
				cij_43 += Ai_[k+2*lda]*B_j[k+5*lda];
				cij_44 += Ai_[k+3*lda]*B_j[k+5*lda];
				cij_45 += Ai_[k+4*lda]*B_j[k+5*lda];
				cij_46 += Ai_[k+5*lda]*B_j[k+5*lda];
				cij_47 += Ai_[k+6*lda]*B_j[k+5*lda];
				cij_48 += Ai_[k+7*lda]*B_j[k+5*lda];
				cij_49 += Ai_[k+0*lda]*B_j[k+6*lda];
				cij_50 += Ai_[k+1*lda]*B_j[k+6*lda];
				cij_51 += Ai_[k+2*lda]*B_j[k+6*lda];
				cij_52 += Ai_[k+3*lda]*B_j[k+6*lda];
				cij_53 += Ai_[k+4*lda]*B_j[k+6*lda];
				cij_54 += Ai_[k+5*lda]*B_j[k+6*lda];
				cij_55 += Ai_[k+6*lda]*B_j[k+6*lda];
				cij_56 += Ai_[k+7*lda]*B_j[k+6*lda];
				cij_57 += Ai_[k+0*lda]*B_j[k+7*lda];
				cij_58 += Ai_[k+1*lda]*B_j[k+7*lda];
				cij_59 += Ai_[k+2*lda]*B_j[k+7*lda];
				cij_60 += Ai_[k+3*lda]*B_j[k+7*lda];
				cij_61 += Ai_[k+4*lda]*B_j[k+7*lda];
				cij_62 += Ai_[k+5*lda]*B_j[k+7*lda];
				cij_63 += Ai_[k+6*lda]*B_j[k+7*lda];
				cij_64 += Ai_[k+7*lda]*B_j[k+7*lda];
			}
			C[((j+0)*lda)+i+0] = cij_1 ;
			C[((j+0)*lda)+i+1] = cij_2 ;
			C[((j+0)*lda)+i+2] = cij_3 ;
			C[((j+0)*lda)+i+3] = cij_4 ;
			C[((j+0)*lda)+i+4] = cij_5 ;
			C[((j+0)*lda)+i+5] = cij_6 ;
			C[((j+0)*lda)+i+6] = cij_7 ;
			C[((j+0)*lda)+i+7] = cij_8 ;
			C[((j+1)*lda)+i+0] = cij_9 ;
			C[((j+1)*lda)+i+1] = cij_10 ;
			C[((j+1)*lda)+i+2] = cij_11 ;
			C[((j+1)*lda)+i+3] = cij_12 ;
			C[((j+1)*lda)+i+4] = cij_13 ;
			C[((j+1)*lda)+i+5] = cij_14 ;
			C[((j+1)*lda)+i+6] = cij_15 ;
			C[((j+1)*lda)+i+7] = cij_16 ;
			C[((j+2)*lda)+i+0] = cij_17 ;
			C[((j+2)*lda)+i+1] = cij_18 ;
			C[((j+2)*lda)+i+2] = cij_19 ;
			C[((j+2)*lda)+i+3] = cij_20 ;
			C[((j+2)*lda)+i+4] = cij_21 ;
			C[((j+2)*lda)+i+5] = cij_22 ;
			C[((j+2)*lda)+i+6] = cij_23 ;
			C[((j+2)*lda)+i+7] = cij_24 ;
			C[((j+3)*lda)+i+0] = cij_25 ;
			C[((j+3)*lda)+i+1] = cij_26 ;
			C[((j+3)*lda)+i+2] = cij_27 ;
			C[((j+3)*lda)+i+3] = cij_28 ;
			C[((j+3)*lda)+i+4] = cij_29 ;
			C[((j+3)*lda)+i+5] = cij_30 ;
			C[((j+3)*lda)+i+6] = cij_31 ;
			C[((j+3)*lda)+i+7] = cij_32 ;
			C[((j+4)*lda)+i+0] = cij_33 ;
			C[((j+4)*lda)+i+1] = cij_34 ;
			C[((j+4)*lda)+i+2] = cij_35 ;
			C[((j+4)*lda)+i+3] = cij_36 ;
			C[((j+4)*lda)+i+4] = cij_37 ;
			C[((j+4)*lda)+i+5] = cij_38 ;
			C[((j+4)*lda)+i+6] = cij_39 ;
			C[((j+4)*lda)+i+7] = cij_40 ;
			C[((j+5)*lda)+i+0] = cij_41 ;
			C[((j+5)*lda)+i+1] = cij_42 ;
			C[((j+5)*lda)+i+2] = cij_43 ;
			C[((j+5)*lda)+i+3] = cij_44 ;
			C[((j+5)*lda)+i+4] = cij_45 ;
			C[((j+5)*lda)+i+5] = cij_46 ;
			C[((j+5)*lda)+i+6] = cij_47 ;
			C[((j+5)*lda)+i+7] = cij_48 ;
			C[((j+6)*lda)+i+0] = cij_49 ;
			C[((j+6)*lda)+i+1] = cij_50 ;
			C[((j+6)*lda)+i+2] = cij_51 ;
			C[((j+6)*lda)+i+3] = cij_52 ;
			C[((j+6)*lda)+i+4] = cij_53 ;
			C[((j+6)*lda)+i+5] = cij_54 ;
			C[((j+6)*lda)+i+6] = cij_55 ;
			C[((j+6)*lda)+i+7] = cij_56 ;
			C[((j+7)*lda)+i+0] = cij_57 ;
			C[((j+7)*lda)+i+1] = cij_58 ;
			C[((j+7)*lda)+i+2] = cij_59 ;
			C[((j+7)*lda)+i+3] = cij_60 ;
			C[((j+7)*lda)+i+4] = cij_61 ;
			C[((j+7)*lda)+i+5] = cij_62 ;
			C[((j+7)*lda)+i+6] = cij_63 ;
			C[((j+7)*lda)+i+7] = cij_64 ;
		}
	}
	if(M%r !=0 || N%c !=0) {
		basic_dgemm_1x1x1(lda, M, N-num_n*c, K,A, B+num_n*c*lda, C+num_n*c*lda);
		basic_dgemm_1x1x1(lda, M-num_m*r, num_n*c, K, A+num_m*r*lda, B, C+num_m*r);
	}
}


