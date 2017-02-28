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
#include <vector>
#include <cstdlib>
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
const float prior_a_aa = log(0.985);
const float prior_a_tt = log(0.000083);
const float prior_a_cc = log(0.000083);
const float prior_a_gg = log(0.00033);
const float prior_a_at = log(0.00017);
const float prior_a_ac = log(0.00017);
const float prior_a_ag = log(0.000667);
const float prior_a_ct = (log(2.78) - 8*log(10));
const float prior_a_gt = (log(1.1) - 7*log(10));
const float prior_a_cg = (log(1.1) - 7*log(10));

const float prior_t_aa = log(0.000083);
const float prior_t_tt = log(0.985);
const float prior_t_cc = log(0.00033);
const float prior_t_gg = log(0.000083);
const float prior_t_at = log(0.00017);
const float prior_t_ac = (log(1.1) - 7*log(10));
const float prior_t_ag = (log(2.78) - 8*log(10));
const float prior_t_ct = log(0.000667);
const float prior_t_gt = log(0.00017);
const float prior_t_cg = (log(1.1) - 7*log(10));

const float prior_c_aa = log(0.000083);
const float prior_c_tt = log(0.00033);
const float prior_c_cc = log(0.985);
const float prior_c_gg = log(0.000083);
const float prior_c_at = (log(1.1) - 7*log(10));
const float prior_c_ac = log(0.00017);
const float prior_c_ag = (log(2.78) - 8*log(10));
const float prior_c_ct = log(0.000667);
const float prior_c_gt = (log(1.1) - 7*log(10));
const float prior_c_cg = log(0.00017);

const float prior_g_aa = log(0.00033);
const float prior_g_tt = log(0.000083);
const float prior_g_cc = log(0.000083);
const float prior_g_gg = log(0.9985);
const float prior_g_at = (log(1.1) - 7*log(10));
const float prior_g_ac = (log(1.1) - 7*log(10));
const float prior_g_ag = (log(6.67) - 4*log(10));
const float prior_g_ct = (log(2.78) - 8*log(10));
const float prior_g_gt = (log(1.67) - 4*log(10));
const float prior_g_cg = (log(1.67) - 4*log(10));

void InitializeLoggedFactorialArray();
double LoggedFactorial(int n);
double GetFactorialPoly(int a, int t, int c, int g);
string CalculateBayes(string str,char c);
vector<string> split(const string &s, const string &seperator);
vector<int> split2int(const string &s, const string &seperator);
string PrintGenotype(string str);
string GetGenotype(vector<int> v, string str);
// void PrintSNPOutput(string chr, string genotype, string atcg, string qatcg);
// void PrintGenotype();

int main(int argc, const char * argv[]) {
	// if (argc<2){
	// 	cerr << "Example: bssnper.o --fa hg19.fa --input BSMAP.sort.bam --output snp.candidate.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>SNP.log" << std::endl;
	// 	return 1;
	// }

	// cout << "Hello, world!\n";
	// if (system("sh ./echo.sh")){
	// 	return 1;
	// }
	InitializeLoggedFactorialArray();              //initialize the logged gactorial array

	ifstream snp_file("/Users/choumasijia/Projects/bssnper/sed10.sed.candidate");
	if(!snp_file) {
		cout << "Cannot open input file.\n";
		return 1;
	}
    ofstream out_file("/Users/choumasijia/Projects/bssnper/bssnper/foo.out");
    char snpline[255];
    char* snp;
    string strline;
    snp_file.getline(snpline, 255);
    while(getline(snp_file, strline)){
        cout<<strline<<endl;
        vector<string> original_line = split(strline, "\t");


        string refbase = original_line[2];
        vector<int> num_watson = split2int(original_line[3], ",");
        vector<int> num_crick = split2int(original_line[4], ",");
        vector<int> qual_watson = split2int(original_line[5], ",");
        vector<int> qual_crick = split2int(original_line[6], ",");
        
     //    string bayes;
    	// bayes = CalculateBayes(strline, 'G');
    	// out_file<< strline<<endl;
    	// cout<<bayes;
//    	cout<<PrintGenotype(bayes);
    }
//	while(snp_file) {
//	    snp_file.getline(snpline, 255);  // delim defaults to '\n'
//	    snp = &snpline[0];
        
//	    cout << snpline << endl;
//	}
	
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

vector<string> split(const string &s, const string &seperator){
    vector<string> result;
    string::size_type pos1, pos2;
    pos2 = s.find(seperator);
    pos1 = 0;
    while (string::npos != pos2)
    {
        result.push_back(s.substr(pos1, pos2 - pos1));
        pos1 = pos2 + 1;
        pos2 = s.find(seperator, pos1);
    }
    result.push_back(s.substr(pos1));
    return result;
}

vector<int> split2int(const string &s, const string &seperator){
    vector<int> result;
    string::size_type pos1, pos2;
    pos2 = s.find(seperator);
    pos1 = 0;
    while (string::npos != pos2)
    {
        result.push_back(stoi(s.substr(pos1, pos2 - pos1)));
        pos1 = pos2 + 1;
        pos2 = s.find(seperator, pos1);
    }
    result.push_back(stoi(s.substr(pos1)));
    return result;
}

string GetGenotype(vector<int> v, string str){

   return "nn";
}

string PrintGenotype(string str){
    int watson [4]={0,4,0,10};
    int crick[4]={0,0,0,0};
    int watsonq[4]={37};
    int crickq[4]={37};
    int totaldepth=1;
    float baseq_watsonA, baseq_watsonT, baseq_watsonC, baseq_watsonG = 0.00019952623149688793;
    float baseq_crickA, baseq_crickT,baseq_crickC, baseq_crickG = 0.00019952623149688793;
    string genotype="AA";
    char reference = 'G';
    //what if "totaldepth == 0 ???"
    float qvalue;
    int depth;
    int var;
    if(genotype=="AA"){
    	switch(reference){
    		case 'A':
    			return "REF";
    		case 'T':
            {
                qvalue= (watsonq[0]>crickq[0])?watsonq[0]:crickq[0];
				depth=watson[0]+crick[0]+watson[1]+crick[1];
				var=watson[0]+crick[0];
				float T2A_weight = var/totaldepth;
				if(depth >= mincover && qvalue >= minquali && var >= minread2 && T2A_weight >= minhomfreq){
					cout<< "$lines[0]\t$lines[1]\t.\t$lines[2]\tA\t$genoqual\tPASS\tAA\t$T2A\t\n"<<endl;
				}else{
					cout<< "$lines[0]\t$lines[1]\t.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$T2A\t"<<endl;
                }
            }
				return "T>AA";
            case 'C':
            {
				qvalue= (watsonq[0]>crickq[0])?watsonq[0]:crickq[0];
				depth=watson[0]+crick[0]+watson[1]+crick[1];
				var=watson[0]+crick[0];
				float C2A_weight = var/totaldepth;
				if(depth >= mincover && qvalue >= minquali && var >= minread2 && C2A_weight >= minhomfreq){
					cout<< "$lines[0]\t$lines[1]\t.\t$lines[2]\tA\t$genoqual\tPASS\tAA\t$C2A\t\n"<<endl;
				}else{
					cout<< "$lines[0]\t$lines[1]\t.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$C2A\t"<<endl;
                }
            }
				return "C>AA";
            case 'G':
            {
				qvalue= watsonq[0];
				depth=watson[0]+watson[3];
				var=watson[0];
				float G2A_weight; 
				if(depth){
					G2A_weight = var/totaldepth;
				}else{
					G2A_weight = 0;
				}
				if(depth >= mincover && qvalue >= minquali && var >= minread2 && G2A_weight >= minhomfreq){
					cout<< "$lines[0]\t$lines[1]\t.\t$lines[2]\tA\t$genoqual\tPASS\tAA\t$G2A\t\n"<<endl;
				}else{
					cout<< "$lines[0]\t$lines[1]\t.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$G2A\t"<<endl;
                }
            }
				return "G>AA";
    	}
    	cout<<"compare yes!"<<endl;
    }

    return "NN\n";
}


string CalculateBayes(string candidate_line, char refbase){
    float aa=0., at=0., ac=0., ag=0., cc=0., cg=0., ct=0., gg=0., gt=0., tt=0.;
    int watson [4]={0,4,0,10};
    int crick[4]={0,0,0,0};
    
    int watsonq[4]={37};
    int crickq[4]={37};
    
    float baseq_watsonA, baseq_watsonT, baseq_watsonC, baseq_watsonG = 0.00019952623149688793;
    float baseq_crickA, baseq_crickT,baseq_crickC, baseq_crickG = 0.00019952623149688793;
    
    //if all of watsonq and crickq is 0, genotypemaybe = NN, qualerr = 0
    
    double nn = GetFactorialPoly(watson[0], crick[1], crick[2], watson[3]);
    if(watsonq[0]>0){      //A>0
        float a_aa = 1- baseq_watsonA;
        float other = baseq_watsonA/3.;
        aa = nn + watson[0]*log(a_aa) + (crick[1]+watson[1]+crick[2]+watson[2]+crick[3]+watson[3])*log(other);
        if(crick[1]>0){  //AT
            float a_at=(1-(baseq_watsonA+baseq_crickT)/2)/2;      //provid genotype is AT, the probility to find A.
            other=(baseq_watsonA+baseq_crickT)/4;   //the probability of other 2 types opear if genotype is AT.
            at = nn+ (watson[0]+crick[1])*log(a_at) + (crick[2]+watson[3]+watson[2]+crick[3])*log(other);
        }
        if(crick[2]>0 || watson[2]>0){  //AC
            float a_ac= (1-baseq_watsonA)/2;
            other=baseq_watsonA/3;
            ac = nn+ (watson[0]+watson[2]+crick[2])*log(a_ac) + (crick[1]+watson[3]+crick[3])*log(other);
        }
        if(crick[3]>0 || watson[3]>0){ //AG
            float a_ag=(1-baseq_watsonA)/2;
            other= baseq_watsonA/3;
            ag = nn+ (watson[0]+watson[3]+crick[3])*log(a_ag) + (crick[2]+crick[1]+watson[2])*log(other);
        }
    }
    if(crickq[1]>0 ){  //filter wsq[1]>0 but crq[1]==0
        float t_tt=1-baseq_crickT;
        float other = baseq_crickT/3;
        //type TT
        tt= nn + crick[1]*log(t_tt)+ (watson[0]+crick[2]+watson[3]+crick[3]+watson[2])*log(other);
        if(watson[2]>0 || crick[2]>0){
            float t_ct=(1-baseq_crickT)/2;
            other=baseq_crickT/3;
            //type CT
            ct= nn + (watson[2]+crick[1]+crick[2])*log(t_ct) + (watson[0]+watson[3]+crick[3])*log(other);
        }
        if(watson[3]>0 || crick[3]>0){
            float t_gt=(1-baseq_crickT)/2;
            other=baseq_crickT/3;
            gt= nn + (crick[1]+watson[3]+crick[3]) * log(t_gt) + (watson[0]+crick[2]+watson[2])*log(other);
        }
    }
    if(crickq[2]>0 || watsonq[2]>0){  //CC CG
        float baseq_C=(baseq_crickC>=baseq_watsonC)?baseq_watsonC:baseq_crickC;
        float c_cc=1-baseq_C;
        float other = baseq_C/3;
        //type CC
        cc = nn + (crick[2]+watson[2])*log(c_cc) + (watson[0]+crick[0]+crick[1]+watson[3]+crick[3])*log(other);
        if(watson[3]>0 || crick[3]>0){ //CG
            float baseq_G= (baseq_watsonG>=baseq_crickG)?baseq_crickG:baseq_watsonG;
            /////////////   baseq_G ? baseq_C ?   ////////////
            float c_cg=(1-baseq_C)/2;
            other=(baseq_C)/2;
            cg = nn + (crick[2]+watson[3]+watson[2]+crick[3])*log(c_cg) + (watson[0]+crick[1])*log(other);
        }
    }
    if(watsonq[3]>0 || crickq[3]>0 ){
        float baseq_G = (baseq_watsonG>=baseq_crickG)?baseq_crickG:baseq_watsonG;
        float g_gg=1-baseq_G;
        float other = baseq_G/3;
        gg = nn + (watson[3]+crick[3])*log(g_gg) + (watson[0]+crick[1]+watson[1]+crick[2]+watson[2])*log(other);
    }
    cout<< candidate_line<<endl;
    
    if(refbase == 'A'){
        aa+=prior_a_aa;
        tt+=prior_a_tt;
        cc+=prior_a_cc;
        gg+=prior_a_gg;
        at+=prior_a_at;
        ac+=prior_a_ac;
        ag+=prior_a_ag;
        ct+=prior_a_ct;
        gt+=prior_a_gt;
        cg+=prior_a_cg;
    }else if(refbase == 'T'){
        aa+=prior_t_aa;
        tt+=prior_t_tt;
        cc+=prior_t_cc;
        gg+=prior_t_gg;
        at+=prior_t_at;
        ac+=prior_t_ac;
        ag+=prior_t_ag;
        ct+=prior_t_ct;
        gt+=prior_t_gt;
        cg+=prior_t_cg;
    }else if(refbase == 'C'){
        aa+=prior_c_aa;
        tt+=prior_c_tt;
        cc+=prior_c_cc;
        gg+=prior_c_gg;
        at+=prior_c_at;
        ac+=prior_c_ac;
        ag+=prior_c_ag;
        ct+=prior_c_ct;
        gt+=prior_c_gt;
        cg+=prior_c_cg;
    }else if(refbase == 'G'){
        aa+=prior_g_aa;
        tt+=prior_g_tt;
        cc+=prior_g_cc;
        gg+=prior_g_gg;
        at+=prior_g_at;
        ac+=prior_g_ac;
        ag+=prior_g_ag;
        ct+=prior_g_ct;
        gt+=prior_g_gt;
        cg+=prior_g_cg;
    }
    
    float prob=0., qual=0.;
    return "NN";
    
};

// char* foo(int n){
// 	char* pStr=new char[n+1];//last one for '\0'
// 	pStr[n]='\0';
// 	int i;
// 	for(i=0;i<n;i++)
// 		pStr[i]=i+97;
// 	return pStr;
// }

