/*
 * main.cpp
 *
 *  Created on: 2021年10月4日
 *      Author: windows
 */

#include"SSYMM01.h"
#include<cblas.h>
#include<string.h>
#include<iostream>
#include<ctime>

using namespace std;

void show(const float* A,int M,int N);
int main(){
	int flag=1;
	char ifshow1='y',ifshow2='y';
	enum CBLAS_ORDER Order=CblasRowMajor;
	enum CBLAS_UPLO Uplo=CblasUpper;
	enum CBLAS_SIDE Side=CblasLeft;
	cout<<"Press number to choose your mode"<<endl;
	cout<<"Row\\Up\\Left(default:1 or other)\t""Row\\Up\\Right(2)\t"
			"Row\\Low\\Left(3)\t""Row\\Low\\Right(4)\n"
			"Col\\Up\\Left(5)\t""Col\\Up\\Right(6)\t"
			"Col\\Low\\Left(7)\t""Col\\Low\\Right(8)"
			<<endl;
	cin>>flag;
	switch(flag){
	case 2:
		Side=CblasRight;
		break;
	case 3:
		Uplo=CblasLower;
		break;
	case 4:
		Uplo=CblasLower;
		Side=CblasRight;
		break;
	case 5:
		Order=CblasColMajor;
		break;
	case 6:
		Order=CblasColMajor;
		Side=CblasRight;
		break;
	case 7:
		Order=CblasColMajor;
		Uplo=CblasLower;
		break;
	case 8:
		Order=CblasColMajor;
		Uplo=CblasLower;
		Side=CblasRight;
		break;
	default:
		flag=1;
	}

	cout<<"Then you will start a MxM*MxN(left)OR MxN*NxN(right) matrices operation\n";
	cout<<"Choose the dimension\nM:";
	int M, N, K, lda, ldb ,ldc;
	cin>>M;
	cout<<"N:";
	cin>>N;

	cout<<"To show matrices Atest,A,B enter 'y'(otherwise 'n')\nifshow1:";
	cin>>ifshow1;
	cout<<"To show the result matrices C(=Atest*B),D(=A*B) enter 'y'(otherwise 'n')\nifshow2:";
	cin>>ifshow2;

	Side==CblasLeft?K=M:K=N;
	lda=K;
	Order==CblasRowMajor?(ldb=N,ldc=N):(ldb=M,ldc=M);
	//而且若Order==CblasColMajor，则AL代表上三角，AU。。。
	float alpha=1.0;
	float beta=1.0;
	float* A=new float[K*K];
	float* B=new float[M*N];
	float* C=new float[M*N];
	float* D=new float[M*N];
	float* Atest=new float[K*K];
	memset(C,0,sizeof(float)*M*N);
	memset(D,0,sizeof(float)*M*N);
	memset(A,0,sizeof(float)*K*K);
	memset(B,0,sizeof(float)*M*N);
	memset(Atest,0,sizeof(float)*K*K);

	for(int i=0;i<K;i++)//Symmerisz
		for(int j=0;j<K;j++){
			if(j>=i)A[i*K+j]=rand()%20+1;
			else A[i*K+j]=A[j*K+i];
			if(3<=flag&&flag<=6&&i>=j)Atest[i*K+j]=A[i*K+j];
			else if((flag<3||flag>6)&&i<=j)Atest[i*K+j]=A[i*K+j];
		}
	for(int i=0;i<M*N;i++){
		B[i]=rand()%20+1;
	}
	if(ifshow1!='n'){
		cout<<"Atest:"<<endl;
		show(Atest,K,K);
		cout<<"A:"<<endl;
		show(A,K,K);
		cout<<"B:"<<endl;
		show(B,M,N);
	}
	clock_t start,end1,end2;

	float* AorB1,*AorB2,ld1,ld2;
	Side==CblasLeft?(AorB1=A,ld1=lda,AorB2=B,ld2=ldb):(AorB1=B,ld1=ldb,AorB2=A,ld2=lda);
	start=clock();
	cblas_sgemm(Order, CblasNoTrans, CblasNoTrans, M, N, K, alpha, AorB1, ld1, AorB2, ld2, beta, D, ldc);
	end1=clock();
	SSYMM01(Order, Side, Uplo, M, N, alpha, Atest, lda, B, ldb, beta, C, ldc);
	end2=clock();
	if(ifshow2!='n'){
		cout<<"C:"<<endl;
		show(C,M,N);
		cout<<"D:"<<endl;
		show(D,M,N);
	}

	cout<<"cblas_sgemm Spend:"<<end1-start<<"ms"<<endl;
	cout<<"SSYMM01 Spend:"<<end2-end1<<"ms"<<endl;
	//比较C和D，若都一样，则Pass！
	if(memcmp(C,D,sizeof(float)*M*N)==0)cout<<"Pass!"<<endl;
	else cout<<"Fail!"<<endl;

	cout<<"Enter any number to end program"<<endl;
	cin>>flag;

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] D;
	delete[] Atest;
}

void show(const float* A,int M,int N){
	for(int i=0;i<M;i++){
		cout<<"[ "<<i+1<<" ]=";
		for(int j=0;j<N;j++)
			cout<<A[i*N+j]<<',';
		cout<<endl;
	}
}
