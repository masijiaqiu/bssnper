//
//  main.cpp
//  bssnper
//
//  Created by masijiaqiu on 13/1/2017.
//  Copyright © 2017 masijiaqiu. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <math.h>
using namespace std;

float minhetfreq =0.1;
float minhomfreq =0.85;
int minquali =15;
int mincover =10;
int minread2 =2;
int maxcover =1000;
float errorate =0.02;
int mapvalue =20;

double c_logged_factorial[151]; 

void InitializeLoggedFactorialArray();
double LoggedFactorial(int n);
char* GetGenotype(int n);
void GetBayes();
int GetFactorial(int n);
void PrintGenotype();


int main(int argc, const char * argv[]) {
	// if (argc<2){
	// 	cerr << "Example: bssnper.o --fa hg19.fa --input BSMAP.sort.bam --output snp.candidate.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>SNP.log" << std::endl;
	// 	return 1;
	// }

	// cout << "Hello, world!\n";
	// if (system("sh ./echo.sh")){
	// 	return 1;
	// }
	
	cout<<c_logged_factorial[2]<<endl;
	ifstream snp_file("../sed.sed.candidate");
	if(!snp_file) {
		cout << "Cannot open input file.\n";
		return 1;
	}

	InitializeLoggedFactorialArray();
	cout<< c_logged_factorial[2]<<endl;
	char snpline[255];
	char* snp;
	while(snp_file) {
	    snp_file.getline(snpline, 255);  // delim defaults to '\n'
	    snp = &snpline[0];
	    cout << snpline[0] << endl;
	}
	
	cerr << "SNPs & genotype finish.\n";
	return 0;
}

void InitializeLoggedFactorialArray(){
	for (int n=0; n<101; n++){
		c_logged_factorial[n] = LoggedFactorial(n);
	}
}

double LoggedFactorial(int n){
	if (n>1){
		return log(n)+LoggedFactorial(n-1);
	}else
		return 0;
}

double GetFactorialPoly(int a, int t, int c, int g){
	return c_logged_factorial[a+t+c+g]-c_logged_factorial[a]-c_logged_factorial[t]-c_logged_factorial[c]-c_logged_factorial[g];
}

// void GetBayes(){

// };

// char* GetGenotype(int n){
// 	char* pStr=new char[n+1];//last one for '\0'
// 	pStr[n]='\0';
// 	int i;
// 	for(i=0;i<n;i++)
// 		pStr[i]=i+97;
// 	return pStr;
// }



// float GetFactorial(int n){
// 	if (n<=1) {
// 		return 0;
// 	}else{
// 		return 
// 	}

// }

