/*
 * SSYMM01.h
 *
 *  Created on: 2021Äê10ÔÂ4ÈÕ
 *      Author: windows
 */

#ifndef SSYMM01_H_
#define SSYMM01_H_

#include<cblas.h>
#include<cassert>

#include<iostream>
using namespace std;

void SSYMM01(OPENBLAS_CONST enum CBLAS_ORDER Order,
				OPENBLAS_CONST enum CBLAS_SIDE Side,
				OPENBLAS_CONST enum CBLAS_UPLO Uplo,
				OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
				OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,
				OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B,
				OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta,
				float *C, OPENBLAS_CONST blasint ldc){

	int Limit=30;

	if(Side==CblasLeft){
		if(M==1){
		    for(int i=0;i<N;i++)Order==CblasRowMajor?C[i]+=A[0]*B[i]:C[i*ldc]+=A[0]*B[i*ldb];
		}
		else if(M<=Limit){
			float* cpyA=new float[M*M];
			int row;
			for(int i=0;i<M;i++)
				for(int j=0;j<M;j++){
					Uplo==CblasUpper?(row=i<j?i:j):(row=i>j?i:j);
					if(Order==CblasRowMajor)cpyA[i*M+j]=A[row*lda+i+j-row];
					else cpyA[j*M+i]=A[(i+j-row)*lda+row];
				}
			cblas_sgemm(Order, CblasNoTrans, CblasNoTrans, M, N, M, alpha, cpyA, M, B, ldb, beta, C, ldc);
			delete[] cpyA;
		}
		else if(N==1){
			int k1=M%2?M/2+1:M/2,k2=M/2;
			const float* A11=A,*A22=A+k1*(lda+1);
			const float* B1=B,*B2=B+k1*ldb;
			if(Order==CblasColMajor){
				B2=B+k1;
			}
			float* C1=C,*C2=C+k1*ldc;
			if(Order==CblasColMajor){
				C2=C+k1;
			}

			const float* edgeA=A+k1;
			if(Order==CblasColMajor)edgeA=A+k1*lda;
			enum CBLAS_TRANSPOSE TranC1=CblasNoTrans;
			enum CBLAS_TRANSPOSE TranC2=CblasTrans;
			if(Uplo==CblasLower){
				edgeA=A+k1*lda;
				if(Order==CblasColMajor)edgeA=A+k1;
				TranC2=CblasNoTrans;
				TranC1=CblasTrans;
			}

			//C1
			SSYMM01(Order, Side, Uplo, k1, 1, alpha, A11, lda, B1, ldb, beta, C1, ldc);
			cblas_sgemm(Order, TranC1, CblasNoTrans, k1, 1, k2, alpha, edgeA, lda, B2, ldb, beta, C1, ldc);
			//C2
			cblas_sgemm(Order, TranC2, CblasNoTrans, k2, 1, k1, alpha, edgeA, lda, B1, ldb, beta, C2, ldc);
			SSYMM01(Order, Side, Uplo, k2, 1, alpha, A22, lda, B2, ldb, beta, C2, ldc);

		}
		else{
			assert(M>Limit);
			//odd and even
			int k1=M%2?M/2+1:M/2, k2=M/2;
			int j1=N%2?N/2+1:N/2, j2=N/2;
			const float* A11=A,*A22=A+k1*(lda+1);
			const float* B11=B,*B12=B+j1,*B21=B+k1*ldb,*B22=B+k1*ldb+j1;
			if(Order==CblasColMajor){
				B12=B+j1*ldb;
				B21=B+k1;
				B22=B+j1*ldb+k1;
			}
			float* C11=C,*C12=C+j1,*C21=C+k1*ldc,*C22=C+k1*ldc+j1;
			if(Order==CblasColMajor){
				C12=C+j1*ldc;
				C21=C+k1;
				C22=C+j1*ldc+k1;
			}
			const float* edgeA=A+k1;
			if(Order==CblasColMajor)edgeA=A+k1*lda;
			enum CBLAS_TRANSPOSE TranC1=CblasNoTrans;
			enum CBLAS_TRANSPOSE TranC2=CblasTrans;
			if(Uplo==CblasLower){
				edgeA=A+k1*lda;
				if(Order==CblasColMajor)edgeA=A+k1;
				TranC2=CblasNoTrans;
				TranC1=CblasTrans;
			}

			//C11
			SSYMM01(Order, Side, Uplo, k1, j1, alpha, A11, lda, B11, ldb, beta, C11, ldc);
			cblas_sgemm(Order, TranC1, CblasNoTrans, k1, j1, k2, alpha, edgeA, lda, B21, ldb, beta, C11, ldc);
			//C12
			SSYMM01(Order, Side, Uplo, k1, j2, alpha, A11, lda, B12, ldb, beta, C12, ldc);
			cblas_sgemm(Order, TranC1, CblasNoTrans, k1, j2, k2, alpha, edgeA, lda, B22, ldb, beta, C12, ldc);
			//C21,A21 need to transfer
			cblas_sgemm(Order, TranC2, CblasNoTrans, k2, j1, k1, alpha, edgeA, lda, B11, ldb, beta, C21, ldc);
			SSYMM01(Order, Side, Uplo, k2, j1, alpha, A22, lda, B21, ldb, beta, C21, ldc);
			//C22,A21 need to transfer
			cblas_sgemm(Order, TranC2, CblasNoTrans, k2, j2, k1, alpha, edgeA, lda, B12, ldb, beta, C22, ldc);
			SSYMM01(Order, Side, Uplo, k2, j2, alpha, A22, lda, B22, ldb, beta, C22, ldc);
		}
	}
	else{
		assert(Side==CblasRight);

		if(N==1){
		    for(int i=0;i<M;i++)Order==CblasRowMajor?C[i*ldc]+=B[i*ldb]*A[0]:C[i]+=B[i]*A[0];
		}
		else if(N<=Limit){
			float* cpyA=new float[N*N];
			int row;
			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++){
					Uplo==CblasUpper?(row=i<j?i:j):(row=i>j?i:j);
					if(Order==CblasRowMajor)cpyA[i*N+j]=A[row*lda+i+j-row];
					else cpyA[j*N+i]=A[(i+j-row)*lda+row];
				}
			cblas_sgemm(Order, CblasNoTrans, CblasNoTrans, M, N, N, alpha, B, ldb, cpyA, N, beta, C, ldc);
			delete[] cpyA;
		}
		else if(M==1){
			int k1=N%2?N/2+1:N/2,k2=N/2;
			const float* A11=A,*A22=A+k1*(lda+1);
			const float* B1=B,*B2=B+k1;
			if(Order==CblasColMajor){
				B2=B+k1*ldb;
			}
			float* C1=C,*C2=C+k1;
			if(Order==CblasColMajor){
				C2=C+k1*ldc;
			}

			const float* edgeA=A+k1;
			if(Order==CblasColMajor)edgeA=A+k1*lda;
			enum CBLAS_TRANSPOSE TranC1=CblasTrans;
			enum CBLAS_TRANSPOSE TranC2=CblasNoTrans;
			if(Uplo==CblasLower){
				edgeA=A+k1*lda;
				if(Order==CblasColMajor)edgeA=A+k1;
				TranC2=CblasTrans;
				TranC1=CblasNoTrans;
			}

			//C1
			SSYMM01(Order, Side, Uplo, 1, k1, alpha, A11, lda, B1, ldb, beta, C1, ldc);
			cblas_sgemm(Order, CblasNoTrans, TranC1, 1, k1, k2, alpha, B2, ldb, edgeA, lda, beta, C1, ldc);
			//C2
			cblas_sgemm(Order, CblasNoTrans, TranC2, 1, k2, k1, alpha, B1, ldb, edgeA, lda, beta, C2, ldc);
			SSYMM01(Order, Side, Uplo, 1, k2, alpha, A22, lda, B2, ldb, beta, C2, ldc);

		}
		else{
			assert(N>Limit);
			//odd and even
			int k1=N%2?N/2+1:N/2, k2=N/2;
			int j1=M%2?M/2+1:M/2, j2=M/2;
			const float* A11=A,*A22=A+k1*(lda+1);
			const float* B11=B,*B12=B+k1,*B21=B+j1*ldb,*B22=B+j1*ldb+k1;
			if(Order==CblasColMajor){
				B12=B+k1*ldb;
				B21=B+j1;
				B22=B+k1*ldb+j1;
			}
			float* C11=C,*C12=C+k1,*C21=C+j1*ldc,*C22=C+j1*ldc+k1;
			if(Order==CblasColMajor){
				C12=C+k1*ldc;
				C21=C+j1;
				C22=C+k1*ldc+j1;
			}
			const float* edgeA=A+k1;
			if(Order==CblasColMajor)edgeA=A+k1*lda;
			enum CBLAS_TRANSPOSE TranC1=CblasTrans;
			enum CBLAS_TRANSPOSE TranC2=CblasNoTrans;
			if(Uplo==CblasLower){
				edgeA=A+k1*lda;
				if(Order==CblasColMajor)edgeA=A+k1;
				TranC2=CblasTrans;
				TranC1=CblasNoTrans;
			}

			//C11
			SSYMM01(Order, Side, Uplo, j1, k1, alpha, A11, lda, B11, ldb, beta, C11, ldc);
			cblas_sgemm(Order, CblasNoTrans, TranC1, j1, k1, k2, alpha, B12, ldb, edgeA, lda, beta, C11, ldc);
			//C21
			SSYMM01(Order, Side, Uplo, j2, k1, alpha, A11, lda, B21, ldb, beta, C21, ldc);
			cblas_sgemm(Order, CblasNoTrans, TranC1, j2, k1, k2, alpha, B22, ldb, edgeA, lda, beta, C21, ldc);
			//C12
			cblas_sgemm(Order, CblasNoTrans, TranC2, j1, k2, k1, alpha, B11, ldb, edgeA, lda, beta, C12, ldc);
			SSYMM01(Order, Side, Uplo, j1, k2, alpha, A22, lda, B12, ldb, beta, C12, ldc);
			//C22,A21 need to transfer
			cblas_sgemm(Order, CblasNoTrans, TranC2, j2, k2, k1, alpha, B21, ldb, edgeA, lda, beta, C22, ldc);
			SSYMM01(Order, Side, Uplo, j2, k2, alpha, A22, lda, B22, ldb, beta, C22, ldc);
		}

	}
}

#endif /* SSYMM01_H_ */

