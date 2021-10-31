/*
 * main.cpp
 *
 *  Created on: 2021年10月23日
 *      Author: windows
 */




#include<iostream>
#include<windows.h>
#include<psapi.h>
#include"median.h"
using namespace std;

template <class T>
void CountT1(string name,T* ptr,void (T::*func)()){
	clock_t start,end;
	start=clock();
	(ptr->*func)();
	end=clock();
	cout<<name<<":"<<end-start<<" ms"<<endl;
}

template < template <class>class T,class G>
void CountT2(string name,T<G>* ptr,G (T<G>::*func)()){
	clock_t start,end;
	G temp;
	start=clock();
	temp=(ptr->*func)();
	end=clock();
	cout<<name<<":"<<temp<<endl;
	cout<<name<<":"<<end-start<<" ms"<<endl;
}

int main(){
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	size_t s0=pmc.WorkingSetSize;
	size_t s1,s2,sk;
	//size_t s3;

	int n=100000000;
	Array<int> A(n);//,B(A);

	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	sk=pmc.WorkingSetSize;
	//CountT1("RDQSort",&A,&Array<int>::RDQSort);


	//A.RDQSort();
	//B.RDQSort();


	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	s1=pmc.WorkingSetSize;

	CountT2("RDMedian",&A,&Array<int>::RDMedian);
	//CountT2("Median",&B,&Array<int>::Median);
	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	s2=pmc.WorkingSetSize;

//	CountT2("Median",&B,&Array<int>::Median);
//	//CountT2("RDMedian",&A,&Array<int>::RDMedian);
//	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
//	s3=pmc.WorkingSetSize;

	cout<<"Space Consumed(RDS):"<<s1-sk<<endl;
	cout<<"Space Consumed(RDM):"<<s2-s1<<endl;
	cout<<"Allocation:"<<sk-s0<<endl;

//	cout<<"Space Consumed:"<<s3-s2<<endl;
//	cout<<"Space Consumed(Stack):"<<s3-sk<<endl;
//	cout<<"Space Consumed(Total):"<<s3-s0<<endl;


}

