/*
 * median.h
 *
 *  Created on: 2021年10月23日
 *      Author: windows
 */

#ifndef MEDIAN_H_
#define MEDIAN_H_



/*
 * QSort.h
 *
 *  Created on: 2021年10月23日
 *      Author: windows
 */

#ifndef QSORT_H_
#define QSORT_H_

#include<iostream>
#include<ctime>
#include<cassert>


template<class T>
class Array{
	int len;
	T* store;

	//The range is[p,r) !!
	T med(int p,int r,int K){
		assert(r-p>=0);
		int q=partition(p,r);
		if(q-p+1==K)return store[q];
		else if(q-p+1>K)return med(p,q,K);
		else return med(q+1,r,K-(q-p+1));
	}
	T rdmed(int p,int r,int K){
		assert(r-p>=0);
		int q=randpartition(p,r);
		if(q-p+1==K)return store[q];
		else if(q-p+1>K)return rdmed(p,q,K);
		else return rdmed(q+1,r,K-(q-p+1));

	}
	void exchange(int i,int j){
		T temp=store[i];
		store[i]=store[j];
		store[j]=temp;
	}
	int partition(int p,int r){
		assert(r-p>0);
		int tail=p;
		T pivot=store[r-1];
		for(int i=p;i<r-1;i++){
			if(store[i]<pivot){
				exchange(tail,i);
				tail++;
			}
		}
		exchange(tail,r-1);
		return tail;
	}
	int randpartition(int p,int r){
		int q=(rand()*1000+rand())%(r-p)+p;
		exchange(q,r-1);
		return partition(p,r);
	}
	void rdqsort(int p,int r){
		assert(r-p>=0);
		if(r-p==1||r-p==0)return;
		else{
			int q=randpartition(p,r);
			//int q=partition(p,r);
			rdqsort(p,q);
			rdqsort(q+1,r);
		}
	}

public:
	Array(int n):len(n){
		assert(n>0);
		store=new T[len];
		srand(time(0));
		for(int i=0;i<len;i++){
			store[i]=T((rand()*10000+rand())%10000000);
			//store[i]=T(rand()%100);
		}
		//std::cout<<"store:"<<std::endl;
		//Show();
	}

	Array(const Array& A){
		len=A.len;
		store=new T[len];
		srand(time(0));
		for(int i=0;i<len;i++){
			store[i]=A.store[i];
		}
		//std::cout<<"storeCpy:"<<std::endl;
		//Show();
	}
	void Show(){
		for(int i=0;i<len;i++){
			if(!(i%100)&&i)std::cout<<std::endl;
			std::cout<<store[i]<<" ";
		}
		std::cout<<std::endl;
	}

	T Median(){
		return med(0,len,(len+1)/2);
	}
	T RDMedian(){
		return rdmed(0,len,(len+1)/2);
	}
	void RDQSort(){
		rdqsort(0,len);
		std::cout<<"RDQSort:"<<std::endl;
		//Show();
	}
};

















#endif /* QSORT_H_ */













#endif /* MEDIAN_H_ */
