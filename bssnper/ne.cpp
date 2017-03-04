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
#include <map>
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
//const double logindex = 0.1*log(0.1);
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
float LoggedFactorial(int n);
float GetFactorialPoly(int a, int t, int c, int g);
string CalculateBayes(string str,char c);
vector<string> split(const string &s, const string &seperator);
vector<int> split2int(const string &s, const string &seperator);
std::vector<float> qual2baseq(std::vector<int> intvec, int &f);
string PrintGenotype(string str);
string GetGenotype(string ref, std::vector<int> w, std::vector<int> c, std::vector<int> wq, std::vector<int> cq);
// void PrintSNPOutput(string chr, string genotype, string atcg, string qatcg);
// void PrintGenotype();

int main(int argc, const char * argv[]) {
    if(argc<2){
        cerr<<"wrong argc."<<endl;
        return 1;
    }
    

    // if (argc<2){
    //     cerr << "Example: bssnper.o --fa hg19.fa --input BSMAP.sort.bam --output snp.candidate.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>SNP.log" << std::endl;
    //     return 1;
    // }

    // cout << "Hello, world!\n";
    // if (system("sh ./echo.sh")){
    //     return 1;
    // }
    InitializeLoggedFactorialArray();              //initialize the logged gactorial array

//    ifstream snp_file("/Users/choumasijia/Projects/bssnper/sed10.sed.candidate");
    ifstream snp_file(argv[1]);
    if(!snp_file) {
        cout << "Cannot open input file.\n";
        return 1;
    }
    ofstream out_file("/Users/choumasijia/Projects/bssnper/bssnper/foo.out");
    char snpline[255];
    string strline;
    snp_file.getline(snpline, 255);
    while(getline(snp_file, strline)){
        
        vector<string> original_line = split(strline, "\t");
        string refbase = original_line[2];
        vector<int> num_watson = split2int(original_line[3], ",");
        vector<int> num_crick = split2int(original_line[4], ",");
        vector<int> qual_watson = split2int(original_line[5], ",");
        vector<int> qual_crick = split2int(original_line[6], ",");
        string genotype_info = GetGenotype(refbase, num_watson, num_crick, qual_watson, qual_crick);
        if(genotype_info == "R" || genotype_info.empty()){
            continue;
        }else{
            cout<< strline << "\t"<<genotype_info<<endl;
        }
    }

//    while(snp_file) {
//        snp_file.getline(snpline, 255);  // delim defaults to '\n'
//        snp = &snpline[0];
        
//        cout << snpline << endl;
//    }
    
    ///////////////////////////////////////////////////
//    map<float, string>genotypes;
//    genotypes[-0.444]="AA";
//    genotypes[-7.3]="AT";
//    genotypes[-34.444]="AG";
//    genotypes[0.3]="AC";
//    genotypes[5.6]="CC";
//    map<float,string>::iterator genotypes_it;
//    for(genotypes_it=genotypes.begin(); genotypes_it!=genotypes.end();++genotypes_it){
//        cout<<"key: "<<genotypes_it->first <<" value: "<<genotypes_it->second<<endl;
//    }
    ///////////////////////////////////////////////
    cerr << "SNPs & genotype finish.\n";
    return 0;
}

void InitializeLoggedFactorialArray(){
    for (int n=0; n<101; n++){
        c_logged_factorial[n] = LoggedFactorial(n);
    }
}

float LoggedFactorial(int n){
    if (n>1){
        return log(n)+LoggedFactorial(n-1);
    }else
        return 0;
}

float GetFactorialPoly(int a, int t, int c, int g){
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

 vector<float> qual2baseq(vector<int> intvec, int &f){
     std::vector<float> result;
     for(std::vector<int>::iterator it = intvec.begin(); it != intvec.end(); ++it) {
         result.push_back( pow(0.1, 0.1* (*it)));
         if (*it == 0){
             f++;
         }
     }
    return result;
 }

string GetGenotype(string ref, vector<int> watson, vector<int> crick, vector<int> qual_watson, vector<int> qual_crick){
    int flag_zero = 0;
    vector<float> baseq_watson = qual2baseq(qual_watson, flag_zero);
    vector<float> baseq_crick = qual2baseq(qual_crick, flag_zero);
    if (flag_zero == 8){
        return "R";
    }

    float nn = GetFactorialPoly(watson[0], crick[1], crick[2], watson[3]);
    float aa=0., at=0., ac=0., ag=0., cc=0., cg=0., ct=0., gg=0., gt=0., tt=0.;
    if(qual_watson[0]>0){      //A>0
        float a_aa = 1- baseq_watson[0];
        float other = baseq_watson[0]/3.;
        aa = nn + watson[0]*log(a_aa) + (crick[1]+watson[1]+crick[2]+watson[2]+crick[3]+watson[3])*log(other);
        if(crick[1]>0){  //AT
            float a_at=(1-(baseq_watson[0]+baseq_crick[1])/2)/2;      //provid genotype is AT, the probility to find A.
            other=(baseq_watson[0]+baseq_crick[1])/4;   //the probability of other 2 types opear if genotype is AT.
            at = nn+ (watson[0]+crick[1])*log(a_at) + (crick[2]+watson[3]+watson[2]+crick[3])*log(other);
        }
        if(crick[2]>0 || watson[2]>0){  //AC
            float a_ac= (1-baseq_watson[0])/2;
            other=baseq_watson[0]/3;
            ac = nn+ (watson[0]+watson[2]+crick[2])*log(a_ac) + (crick[1]+watson[3]+crick[3])*log(other);
        }
        if(crick[3]>0 || watson[3]>0){ //AG
            float a_ag=(1-baseq_watson[0])/2;
            other= baseq_watson[0]/3;
            ag = nn+ (watson[0]+watson[3]+crick[3])*log(a_ag) + (crick[2]+crick[1]+watson[2])*log(other);
        }
    }
    if(qual_crick[1]>0 ){  //filter wsq[1]>0 but crq[1]==0
        float t_tt=1-baseq_crick[1];
        float other = baseq_crick[1]/3;
        //type TT
        tt= nn + crick[1]*log(t_tt)+ (watson[0]+crick[2]+watson[3]+crick[3]+watson[2])*log(other);
        if(watson[2]>0 || crick[2]>0){
            float t_ct=(1-baseq_crick[1])/2;
            other=baseq_crick[1]/3;
            //type CT
            ct= nn + (watson[2]+crick[1]+crick[2])*log(t_ct) + (watson[0]+watson[3]+crick[3])*log(other);
        }
        if(watson[3]>0 || crick[3]>0){
            float t_gt=(1-baseq_crick[1])/2;
            other=baseq_crick[1]/3;
            gt= nn + (crick[1]+watson[3]+crick[3]) * log(t_gt) + (watson[0]+crick[2]+watson[2])*log(other);
        }
    }
    if(qual_crick[2]>0 || qual_watson[2]>0){  //CC CG
        float baseq_C=(baseq_crick[2]>=baseq_watson[2])?baseq_watson[2]:baseq_crick[2];
        float c_cc=1-baseq_C;
        float other = baseq_C/3;
        //type CC
        cc = nn + (crick[2]+watson[2])*log(c_cc) + (watson[0]+crick[0]+crick[1]+watson[3]+crick[3])*log(other);
        if(watson[3]>0 || crick[3]>0){ //CG
            float baseq_G= (baseq_watson[3]>=baseq_crick[3])?baseq_crick[3]:baseq_watson[3];
            /////////////   baseq_G ? baseq_C ?   ////////////
            float c_cg=(1-baseq_C)/2;
            other=(baseq_C)/2;
            cg = nn + (crick[2]+watson[3]+watson[2]+crick[3])*log(c_cg) + (watson[0]+crick[1])*log(other);
        }
    }
    if(qual_watson[3]>0 || qual_crick[3]>0 ){
        float baseq_G = (baseq_watson[3]>=baseq_crick[3])?baseq_crick[3]:baseq_watson[3];
        float g_gg=1-baseq_G;
        float other = baseq_G/3;
        gg = nn + (watson[3]+crick[3])*log(g_gg) + (watson[0]+crick[1]+watson[1]+crick[2]+watson[2])*log(other);
    }
    
    map<float, string> genotypes;
    
     if(ref == "A"){
         if(aa) genotypes[(aa+=prior_a_aa)]= "AA";
         if(tt) genotypes[(tt+=prior_a_tt)]= "TT";
         if(cc) genotypes[(cc+=prior_a_cc)]= "CC";
         if(gg) genotypes[(gg+=prior_a_gg)]= "GG";
         if(at) genotypes[(at+=prior_a_at)]= "AT";
         if(ac) genotypes[(ac+=prior_a_ac)]= "AC";
         if(ag) genotypes[(ag+=prior_a_ag)]= "AG";
         if(ct) genotypes[(ct+=prior_a_ct)]= "CT";
         if(gt) genotypes[(gt+=prior_a_gt)]= "GT";
         if(cg) genotypes[(cg+=prior_a_cg)]= "CG";
     }else if(ref == "T"){
         if(aa) genotypes[(aa+=prior_t_aa)]= "AA";
         if(tt) genotypes[(tt+=prior_t_tt)]= "TT";
         if(cc) genotypes[(cc+=prior_t_cc)]= "CC";
         if(gg) genotypes[(gg+=prior_t_gg)]= "GG";
         if(at) genotypes[(at+=prior_t_at)]= "AT";
         if(ac) genotypes[(ac+=prior_t_ac)]= "AC";
         if(ag) genotypes[(ag+=prior_t_ag)]= "AG";
         if(ct) genotypes[(ct+=prior_t_ct)]= "CT";
         if(gt) genotypes[(gt+=prior_t_gt)]= "GT";
         if(cg) genotypes[(cg+=prior_t_cg)]= "CG";
     }else if(ref == "C"){
         if(aa) genotypes[(aa+=prior_c_aa)]= "AA";
         if(tt) genotypes[(tt+=prior_c_tt)]= "TT";
         if(cc) genotypes[(cc+=prior_c_cc)]= "CC";
         if(gg) genotypes[(gg+=prior_c_gg)]= "GG";
         if(at) genotypes[(at+=prior_c_at)]= "AT";
         if(ac) genotypes[(ac+=prior_c_ac)]= "AC";
         if(ag) genotypes[(ag+=prior_c_ag)]= "AG";
         if(ct) genotypes[(ct+=prior_c_ct)]= "CT";
         if(gt) genotypes[(gt+=prior_c_gt)]= "GT";
         if(cg) genotypes[(cg+=prior_c_cg)]= "CG";
     }else if(ref == "G"){
         if(aa) genotypes[(aa+=prior_g_aa)]= "AA";
         if(tt) genotypes[(tt+=prior_g_tt)]= "TT";
         if(cc) genotypes[(cc+=prior_g_cc)]= "CC";
         if(gg) genotypes[(gg+=prior_g_gg)]= "GG";
         if(at) genotypes[(at+=prior_g_at)]= "AT";
         if(ac) genotypes[(ac+=prior_g_ac)]= "AC";
         if(ag) genotypes[(ag+=prior_g_ag)]= "AG";
         if(ct) genotypes[(ct+=prior_g_ct)]= "CT";
         if(gt) genotypes[(gt+=prior_g_gt)]= "GT";
         if(cg) genotypes[(cg+=prior_g_cg)]= "CG";
     }

    float sum_qual=0.;
    float max_qual=0.;
    float genotype_qual;
    float genotype_prob;
    string genotype="R";
    map<float,string>::iterator genotypes_it;
    for(genotypes_it=genotypes.begin(); genotypes_it!=genotypes.end();++genotypes_it){
        genotype = genotypes_it->second;
        max_qual = pow(2.7, genotypes_it->first);
        sum_qual += max_qual;
    }

    int genotypes_size = genotypes.size();
      if (genotypes_size==0){
        return "R";
    }else if(genotypes_size==1){
        if(genotype == "AA"){
            genotype_prob = 1-1/(1+pow(0.5, watson[0]));
        }else if(genotype == "TT"){
            genotype_prob = 1-1/(1+pow(0.5, crick[1]));
        }else if(genotype == "CC"){
            genotype_prob = 1-1/(1+pow(0.5, (watson[2]+crick[2])));
        }else if(genotype == "GG"){
            genotype_prob = 1-1/(1+pow(0.5, (watson[3]+crick[3])));
        }else{
            genotype_prob = 1;
        }
        if(genotype_prob){
            genotype_qual= -10*log(genotype_prob)/log(10);
        }else{
            genotype_qual = 1000;
        }
    }else if(sum_qual && (sum_qual != max_qual)){
        genotype_prob = 1- max_qual/sum_qual;
        genotype_qual = -10*log(genotype_prob)/log(10);
    }else{
        genotype_qual = 1000;
    }
    char genotype_qual_appr[250];
    sprintf(genotype_qual_appr,"%.0f\t",genotype_qual);



    float totaldepth = static_cast<float>(watson[0]+watson[1]+watson[2]+watson[3]+crick[0]+crick[1]+crick[2]+crick[3]);
    float depth;
    int qvalue,qvalueA,qvalueT,qvalueC,qvalueG;
    string alt,filter;
    char genotype_freq[250];
    // AA
    if(genotype == "AA"){
        if(ref == "A"){
            return "R";
        }else if(ref == "T"){    //T>AA
            qvalue=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            depth=static_cast<float>(watson[0]+crick[0]+watson[1]+crick[1]);
            int var=watson[0]+crick[0];
            float T2A = var/totaldepth;
            sprintf(genotype_freq,"%.3f\t",T2A);
            alt = "A\t";
            if(depth >= mincover  && qvalue >= minquali && var >=minread2 && T2A>=minhomfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "C"){    //C>AA
            qvalue=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            depth=static_cast<float>(watson[2]+crick[2]+watson[0]+crick[0]);
            int var=watson[0]+crick[0];
            float C2A = var/totaldepth;
            sprintf(genotype_freq,"%.3f\t",C2A);
            alt = "A\t";
            if(depth >= mincover  && qvalue >= minquali && var >=minread2 && C2A>=minhomfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "G"){   // G>AA
            qvalue=qual_watson[0];
            depth=static_cast<float>(watson[3]+watson[0]);
            int var=watson[0];
            alt = "A\t";
            float G2A =0;
            if(depth>0){
                G2A = var/totaldepth;
            }
            sprintf(genotype_freq,"%.3f\t",G2A);
            if(depth >= mincover  && qvalue >= minquali && var >=minread2 && G2A>=minhomfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }   
    }
    //AT
    else if(genotype =="AT"){
        if(ref == "A"){     //A>AT    
            qvalue=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];;
            depth=static_cast<float>(watson[0]+watson[1]+crick[0]+crick[1]);
            int var=watson[1]+crick[1];
            float A2T = var/totaldepth;
            sprintf(genotype_freq,"%.3f\t",A2T);
            alt = "T\t";
            if(depth >= mincover  && qvalue >= minquali && var >=minread2 && A2T>=minhomfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "T"){    //T>AT
            qvalue=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            depth=static_cast<float>(watson[0]+crick[0]+watson[1]+crick[1]);
            int var=watson[0]+crick[0];
            float T2A = var/totaldepth;
            sprintf(genotype_freq,"%.3f\t",T2A);
            alt = "A\t";
            if(depth >= mincover  && qvalue >= minquali && var >=minread2 && T2A>=minhomfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "C"){    //C>AT
            qvalueA=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            qvalueT=qual_crick[1];
            int varA=watson[0]+crick[0];
            int varT=crick[1];
            depth=static_cast<float>(watson[0]+crick[0]+crick[1]);
            float C2A = varA/totaldepth;
            float C2T = varT/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f\t",C2A,C2T);
            alt = "AT\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueT>=minquali && varA >=minread2 && varT >=minread2 && C2A>=minhetfreq && C2T>= minhetfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "G"){
            qvalueA=qual_watson[0];
            qvalueT=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
            int varA=watson[0];
            int varT=crick[1]+watson[1];
            depth=static_cast<float>(watson[0]+crick[1]+watson[1]);
            float G2A = varA/totaldepth;
            float G2T = varT/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f\t",G2A,G2T);
            alt = "AT\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueT>=minquali && varA >=minread2 && varT >=minread2 && G2A>=minhetfreq && G2T>= minhetfreq){
                filter = "Pass\t";
            }else{
                filter = "Low\t";
            }
        }
    }    
    //AC
    else if(genotype =="AC"){
        if(ref =="A"){
            qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            int varC=watson[2]+crick[2];
            depth=static_cast<float>(watson[0]+crick[2]+watson[2]);
            float A2C=varC/totaldepth;
            sprintf(genotype_freq,"%.3f\t",A2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 && A2C>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref == "T"){
            qvalueA=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            int varA=watson[0]+crick[0];
            int varC=crick[2]+watson[2];
            depth=static_cast<float>(watson[0]+crick[2]+watson[2]+crick[0]+crick[1]);
            float T2A=varA/totaldepth;
            float T2C=varC/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f\t",T2A,T2C);
            alt = "AC\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueC>=minquali && varA >=minread2 && varC>=minread2 && T2A>=minhetfreq && T2C>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "C"){
            qvalueA=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            int varA=watson[0]+crick[0];
            depth=static_cast<float>(watson[0]+crick[2]+watson[2]);
            float C2A=varA/totaldepth;
            sprintf(genotype_freq, "%.3f\t", C2A);
            alt = "A\t";
            if(depth >= mincover  && qvalueA >= minquali  && varA >=minread2 && C2A>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "G"){
            qvalueA=qual_watson[0];
            qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            int varA=watson[0];
            int varC=crick[2]+watson[2];
            depth=static_cast<float>(watson[0]+crick[2]+watson[2]+crick[0]);
            float G2A=varA/totaldepth;
            float G2C=varC/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f\t",G2A,G2C);
            alt = "AC\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueC>=minquali && varA>=minread2 && varC>=minread2 && G2A>=minhetfreq && G2C>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }
        }
    }
    else if(genotype =="AG"){
        if(ref == "A"){
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            int varG=watson[3]+crick[3];
            depth=static_cast<float>(watson[0]+crick[3]+watson[3]+watson[1]+crick[1]+watson[2]+crick[2]);
            float A2G=varG/depth;
            sprintf(genotype_freq,"%.3f\t",A2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 && A2G>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }
        }else if(ref == "T"){         //T>AG
            qvalueA=qual_watson[0];
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            int varA=watson[0];
            int varG=crick[3]+watson[3];
            depth=static_cast<float>(watson[0]+crick[3]+watson[3]+crick[0]+crick[1]+watson[1]);
            float T2A=varA/totaldepth;
            float T2G=varG/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f\t",T2A,T2G);
            alt = "AG\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueG>=minquali && varA >=minread2 && varG>=minread2 && T2A>=minhetfreq && T2G>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }
        }else if (ref == "C")      //C>AG
        {
            qvalueA=qual_watson[0];
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            int varA=watson[0];
            int varG=crick[3]+watson[3];
            depth=static_cast<float>(watson[0]+crick[2]+watson[2]+crick[0]+crick[1]+watson[1]);
            float C2A=varA/totaldepth;
            float C2G=varG/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f\t",C2A,C2G);
            alt = "AG\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueG>=minquali && varA >=minread2 && varG>=minread2 && C2A>=minhetfreq && C2G>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }            
        }else if (ref == "G"){
            /*G>AG*/
            qvalueA=qual_watson[0];
               int varA=watson[0];
            depth=static_cast<float>(crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2]);
            float G2A=varA/depth;
            sprintf(genotype_freq,"%.3f\t",G2A);
            alt = "A\t";
            if(depth >= mincover  && qvalueA >= minquali  && varA >=minread2 && G2A>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter = "Low\t";
            }
        }
    }
    else if(genotype =="TT"){
        if(ref =="A"){
            qvalueT=qual_crick[1];
            depth=static_cast<float>(watson[0]+crick[0]+crick[1]);
            int varT=crick[1];
            float A2T = varT/totaldepth;
            sprintf(genotype_freq,"%.3f\t",A2T);
            alt = "T\t";
            if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 && A2T>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref == "T"){
            return "R";
        }else if(ref == "C"){
            qvalueT=qual_crick[1];
            depth=static_cast<float>(crick[1]+crick[3]);
            int varT=crick[1];
            float C2T=0.;
            if(depth>0){
                C2T=varT/depth;
            }
            sprintf(genotype_freq,"%.3f\t",C2T);
            alt = "T\t";
            if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 && C2T>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref == "G"){
            qvalueT=qual_crick[1];
            depth=static_cast<float>(watson[3]+crick[3]+crick[1]);
            int varT=crick[1]+watson[1];
            float G2T=varT/totaldepth;
            sprintf(genotype_freq,"%.3f\t",G2T);
            alt = "T\t";
            if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 && G2T>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }
    }
    else if(genotype =="CC"){
        if(ref =="A"){
            qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
            depth=static_cast<float>(watson[0]+crick[0]+crick[2]+watson[2]+watson[3]+watson[3]+crick[1]);
            int varC=watson[2]+crick[2];
            float A2C = varC/depth;
            sprintf(genotype_freq,"%.3f\t",A2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 && A2C>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref =="T"){
            qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
            depth=static_cast<float>(crick[1]+crick[2]+watson[2]+watson[0]+crick[0]+watson[3]+crick[3]);
            int varC=watson[2]+crick[2];
            float T2C = varC/depth;
            sprintf(genotype_freq,"%.3f\t",T2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 && T2C>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref =="C"){
            return "R";
        }else if(ref =="G"){
            qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
            depth=static_cast<float>(watson[3]+crick[3]+crick[2]+watson[2]+watson[0]+crick[0]+crick[1]);
            int varC=watson[2]+crick[2];
            float G2C = varC/depth;
            sprintf(genotype_freq,"%.3f\t",G2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 && G2C>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }            
       } 
       else if(genotype =="GG"){
        if(ref =="A"){
            qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
            depth=static_cast<float>(watson[0]+crick[3]+watson[3]+watson[1]+crick[1]+watson[2]+crick[2]);
            int varG=watson[3]+crick[3];
            float A2G = varG/depth;
            sprintf(genotype_freq,"%.3f\t",A2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 && A2G>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref =="T"){
            qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
            depth=static_cast<float>(crick[1]+watson[1]+crick[3]+watson[3]+watson[0]+watson[1]+crick[1]);
            int varG=watson[3]+crick[3];
            float T2G = varG/depth;
            sprintf(genotype_freq,"%.3f\t",T2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 && T2G>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref =="C"){
            qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
            depth=static_cast<float>(crick[1]+watson[1]+crick[3]+watson[3]+watson[0]+watson[1]+crick[1]);
            int varG=watson[3]+crick[3];
            float C2G = varG/depth;
            sprintf(genotype_freq,"%.3f\t",C2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 && C2G>=minhomfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref== "G"){
            return "R";
        }                
       } 
       else if(genotype =="CT"){
         if(ref =="A"){
            qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            qvalueT=qual_crick[1];
            int varC=watson[2]+crick[2];
            int varT=crick[1];
            depth=static_cast<float>(totaldepth-watson[1]);
            float A2C = varC/depth;
            float A2T = varC/depth;
            sprintf(genotype_freq,"%.3f,%.3f\t",A2C,A2T);
            alt = "CT\t";
            if(depth >= mincover  && qvalueC >= minquali && qvalueT >= minquali && varC >=minread2 && varT >=minread2 && A2C>=minhetfreq && A2T>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref == "T"){
             qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            int varC=watson[2]+crick[2];
            depth=static_cast<float>(totaldepth-watson[1]);
            float T2C = varC/depth;
            sprintf(genotype_freq,"%.3f\t",T2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 && T2C>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref == "C"){
             qvalueT=qual_crick[1];
            int varT=crick[1];
            depth=static_cast<float>(totaldepth-watson[1]);
            float C2T=varT/depth;
            sprintf(genotype_freq,"%.3f\t",C2T);
            alt = "T\t";
            if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 && C2T>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref == "G"){
             qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            qvalueT=qual_crick[1];
            int varC=watson[2]+crick[2];
            int varT=crick[1];
            depth=static_cast<float>(totaldepth-watson[1]);
            float G2C=varC/depth;
            float G2T=varT/depth;
            sprintf(genotype_freq,"%.3f,%.3f\t",G2C, G2T);
            alt = "CT\t";
            if(depth >= mincover  && qvalueC >= minquali && qvalueT>=minquali && varC >=minread2 && varT>=minread2 && G2C>=minhetfreq && G2T>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }
    }
    else if(genotype =="GT"){
        if(ref =="A"){
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueT=qual_crick[1];
            int varG=watson[3]+crick[3];
            int varT=crick[1];
            depth=static_cast<float>(totaldepth-watson[1]-crick[0]);
            float A2G = varG/depth;
            float A2T = varT/depth;
            sprintf(genotype_freq,"%.3f,%.3f\t",A2G,A2T);
            alt = "GT\t";
            if(depth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 && A2G>=minhetfreq && A2T>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
        }else if(ref =="T"){
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            int varG=watson[3]+crick[3];
            depth=static_cast<float>(watson[0]+watson[1]+watson[2]+watson[3]+crick[1]+crick[2]+crick[3]);
            float T2G = varG/depth;
            sprintf(genotype_freq,"%.3f\t",T2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 && T2G>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref =="C"){
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueT=qual_crick[1];
            int varG=watson[3]+crick[3];
            int varT=crick[1];
            depth=static_cast<float>(totaldepth-watson[1]-crick[0]);
            float C2G = varG/depth;
            float C2T = varT/depth;
            sprintf(genotype_freq,"%.3f,%.3f\t",C2G,C2T);
            alt = "GT\t";
            if(depth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 && C2G>=minhetfreq && C2T>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref =="G"){
            qvalueT=qual_crick[1];
            int varT=crick[1];
            float G2T = varT/totaldepth;
            sprintf(genotype_freq,"%.3f\t",G2T);
            alt = "T\t";
            if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 && G2T>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }
    }
    else if(genotype =="CG"){
        if(ref =="A"){
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueC=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
            int varG=watson[3]+crick[3];
            int varC=crick[1]+watson[1];
            depth=totaldepth-static_cast<float>(watson[1]+crick[0]);
            float A2G = varG/depth;
            float A2C = varC/depth;
            sprintf(genotype_freq,"%.3f,%.3f\t",A2C,A2G);
            alt = "CG\t";
            if(depth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 && A2G>=minhetfreq && A2C>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref =="T"){
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueC=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
            int varG=watson[3]+crick[3];
            int varC=crick[1]+watson[1];
            depth=static_cast<float>(totaldepth-watson[1]-crick[0]);
            float T2G = varG/depth;
            float T2C = varC/depth;
            sprintf(genotype_freq,"%.3f,%.3f\t",T2C,T2G);
            alt = "CG\t";
            if(depth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 && T2G>=minhetfreq && T2C>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref =="C"){
            qvalueG=qual_watson[3];
            int varG=watson[3];
            depth=totaldepth-static_cast<float>(watson[1]+crick[0]);
            float C2G = varG/depth;
            sprintf(genotype_freq,"%.3f\t",C2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 && C2G>=minhetfreq){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }else if(ref =="G"){
            qvalueC=qual_crick[2];
            int varC=crick[2];
            depth=totaldepth-static_cast<float>(watson[1]+crick[0]);
            float G2C = varC/depth;
            sprintf(genotype_freq,"%.3f\t",G2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2){
                filter = "PASS\t";
            }else{
                filter ="Low\t";
            }
         }
    }
    else{
        return "[ACGT]\t";
    }

    string geno_info = alt.append(filter).append(genotype_qual_appr).append(genotype).append("\t").append(genotype_freq);

    return geno_info;
}


