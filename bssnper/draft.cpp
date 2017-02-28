#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <math.h>
#include <vector>
#include <cstdlib>
using namespace std;

void sw(int a, int b){
        int c;
        c =a;
        a =b;
        b =c;
        cout<<b<<endl;
        return;
}

int main(int argc, const char * argv[]){
    int s = 1;
    int t = 3;
    sw(s,t);
    cout<<t<<endl;
    return 0;
}


{

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
}