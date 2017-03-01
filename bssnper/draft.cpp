//include <iostream>
//include <fstream>
//include <string>
//include <typeinfo>
//include <math.h>
//include <vector>
//include <cstdlib>
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

// char* foo(int n){
// 	char* pStr=new char[n+1];//last one for '\0'
// 	pStr[n]='\0';
// 	int i;
// 	for(i=0;i<n;i++)
// 		pStr[i]=i+97;
// 	return pStr;
// }


{

#TT
	if(genotypemaybe eq "TT"){
                if(lines[2] =~/A/i){#A>TT
			 qvalueT=qual_crick[1];
                         depth=watson[0]+crick[0]+crick[1];
                         varT=crick[1];

                        if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 ){
                                 A2T=sprintf("%.3f",varT/totaldepth);
                                if(A2T>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tPASS\tTT\tA2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tTT\tA2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 A2T=sprintf("%.3f",varT/totaldepth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tTT\tA2T\t".join("\t",@lines[3..6])."\n";
                        }
                }
                if(lines[2] =~/T/i){
			return "REF"; 
                }
                if(lines[2] =~/C/i){#C>TT
			 qvalueT=qual_crick[1];
                         depth=crick[1]+crick[3];
                         varT=crick[1];
			 C2T;
			if(depth>0){
				C2T=sprintf("%.3f",varT/depth);

			}else{
				C2T=0;
				
			}
			
                        if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 ){
                                if(C2T>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tPASS\tTT\tC2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tTT\tC2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tTT\tC2T\t".join("\t",@lines[3..6])."\n";
                        }
	
                }
                if(lines[2] =~/G/i){#G>TT
			 qvalueT=qual_crick[1];
                         depth=watson[3]+crick[3]+crick[1];
                         varT=crick[1]+watson[1];

                        if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 ){
                                 G2T=sprintf("%.3f",varT/totaldepth);
                                if(G2T>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tPASS\tTT\tG2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tTT\tG2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 G2T=sprintf("%.3f",varT/totaldepth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tTT\tG2T\t".join("\t",@lines[3..6])."\n";
                        }
	
                }
        }  

#CC     
        if(genotypemaybe eq "CC"){
                if(lines[2] =~/A/i){#A>CC
			 qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
                         depth=watson[0]+crick[0]+crick[2]+watson[2]+watson[3]+watson[3]+crick[1];
                         varC=watson[2]+crick[2];

                        if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 ){
                                 A2C=sprintf("%.3f",varC/depth);
                                if(A2C>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tPASS\tCC\tA2C\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCC\tA2C\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 A2C=sprintf("%.3f",varC/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCC\tA2C\t".join("\t",@lines[3..6])."\n";
                        }	
                }   
                if(lines[2] =~/T/i){#T>CC
			 qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
                         depth=crick[1]+crick[2]+watson[2]+watson[0]+crick[0]+watson[3]+crick[3];
                         varC=watson[2]+crick[2];

                        if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 ){
                                 T2C=sprintf("%.3f",varC/depth);
                                if(T2C>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tPASS\tCC\tT2C\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCC\tT2C\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 T2C=sprintf("%.3f",varC/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCC\tT2C\t".join("\t",@lines[3..6])."\n";
                        }

                }   
                if(lines[2] =~/C/i){#C>CC  
			return "REF";
                }   
                if(lines[2] =~/G/i){#G>CC
			 qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
                         depth=watson[3]+crick[3]+crick[2]+watson[2]+watson[0]+crick[0]+crick[1];
                         varC=watson[2]+crick[2];

                        if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 ){
                                 G2C=sprintf("%.3f",varC/depth);
                                if(G2C>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tPASS\tCC\tG2C\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCC\tG2C\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 G2C=sprintf("%.3f",varC/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCC\tG2C\t".join("\t",@lines[3..6])."\n";
                        }
                }   
        }   
#GG    
        if(genotypemaybe eq "GG"){
                if(lines[2] =~/A/i){#A>GG
			 qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
                         depth=watson[0]+crick[3]+watson[3]+watson[1]+crick[1]+watson[2]+crick[2];
                         varG=watson[3]+crick[3];

                        if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 ){
                                 A2G=sprintf("%.3f",varG/depth);
                                if(A2G>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tPASS\tGG\tA2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGG\tA2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 A2G=sprintf("%.3f",varG/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGG\tA2G\t".join("\t",@lines[3..6])."\n";
                        }
			
                }   
                if(lines[2] =~/T/i){#T>GG
			 qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
                         depth=crick[1]+watson[1]+crick[3]+watson[3]+watson[0]+watson[1]+crick[1];
                         varG=watson[3]+crick[3];

                        if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 ){
                                 T2G=sprintf("%.3f",varG/depth);
                                if(T2G>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tPASS\tGG\tT2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGG\tT2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 T2G=sprintf("%.3f",varG/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGG\tT2G\t".join("\t",@lines[3..6])."\n";
                        }
	
                } 

  
                if(lines[2] =~/C/i){#C>GG
			 qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
                         depth=crick[1]+watson[1]+crick[3]+watson[3]+watson[0]+watson[1]+crick[1];
                         varG=watson[3]+crick[3];

                        if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 ){
                                 C2G=sprintf("%.3f",varG/depth);
                                if(C2G>=minhomfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tPASS\tGG\tC2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGG\tC2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 C2G=sprintf("%.3f",varG/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGG\tC2G\t".join("\t",@lines[3..6])."\n";
                        }
	
                }   
                if(lines[2] =~/G/i){
			return "REF";
                }   
        }    
#CT	
	if(genotypemaybe eq "CT"){
                if(lines[2] =~/A/i){#A>CT
			 qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
                         qvalueT=qual_crick[1];
                         varC=watson[2]+crick[2];
                         varT=crick[1];
                        # depth=watson[0]+crick[0]+crick[1];
			 depth=totaldepth-watson[1];
                        if(depth >= mincover  && qvalueC >= minquali && qvalueT>=minquali && varC >=minread2 && varT>=minread2 ){
                                 A2C=sprintf("%.3f",varC/depth);
                                 A2T=sprintf("%.3f",varT/depth);
                                if(A2C>=minhetfreq && A2T>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCT\tgenoqual\tPASS\tCT\tA2C\,A2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCT\tgenoqual\tLow\tCT\tA2C\,A2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 A2C= sprintf("%.3f",varC/depth);
                                 A2T= sprintf("%.3f",varT/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tCT\tgenoqual\tLow\tCT\tA2C\,A2T\t".join("\t",@lines[3..6])."\n";
                        }

                }
                if(lines[2] =~/T/i){ #T>CT
			 qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
                         varC=watson[2]+crick[2];
                        # depth=crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2];
			 depth=totaldepth-watson[1];
                         T2C= sprintf("%.3f",varC/depth);
                        if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 ){
                                if(T2C>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tPASS\tCT\tT2C\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCT\tT2C\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tCT\tT2C\t".join("\t",@lines[3..6])."\n";
                        }
                }
                if(lines[2] =~/C/i){ #C>CT
			 qvalueT=qual_crick[1];
                         varT=crick[1];
                        # depth=crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2];
                         depth=totaldepth-watson[1];
                         C2T=sprintf("%.3f",varT/depth);
                        if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 ){
                                if(C2T>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tPASS\tCT\tC2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tCT\tC2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tCT\tC2T\t".join("\t",@lines[3..6])."\n";
                        }

                }
                if(lines[2] =~/G/i){#G>CT
			 qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
                         qvalueT=qual_crick[1];
                         varC=watson[2]+crick[2];
                         varT=crick[1];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=totaldepth-watson[1];
                        if(depth >= mincover  && qvalueC >= minquali && qvalueT>=minquali && varC >=minread2 && varT>=minread2 ){
                                 G2C=sprintf("%.3f",varC/depth);
                                 G2T=sprintf("%.3f",varT/depth);
                                if(G2C>=minhetfreq && G2T>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCT\tgenoqual\tPASS\tCT\tG2C\,G2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCT\tgenoqual\tLow\tCT\tG2C\,G2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 G2C= sprintf("%.3f",varC/depth);
                                 G2T= sprintf("%.3f",varT/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tCT\tgenoqual\tLow\tCT\tG2C\,G2T\t".join("\t",@lines[3..6])."\n";
                        }			
                }
        }   
#GT
	if(genotypemaybe eq "GT"){
                if(lines[2] =~/A/i){ #A>GT
			 qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
                         qvalueT=qual_crick[1];
                         varG=watson[3]+crick[3];
                         varT=crick[1];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=totaldepth-watson[1]-crick[0];
                        if(depth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 ){
                                 A2G= sprintf("%.3f",varG/depth);
                                 A2T= sprintf("%.3f",varT/depth);
                                if(A2G>=minhetfreq && A2T>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tGT\tgenoqual\tPASS\tGT\tA2G\,A2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tGT\tgenoqual\tLow\tGT\tA2G\,A2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 A2G= sprintf("%.3f",varG/depth);
                                 A2T= sprintf("%.3f",varT/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tGT\tgenoqual\tLow\tGT\tA2G\,A2T\t".join("\t",@lines[3..6])."\n";
                        }

                }
                if(lines[2] =~/T/i){#T>G
			 qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
                         varG=watson[3]+crick[3];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=watson[0]+watson[1]+watson[2]+watson[3]+crick[1]+crick[2]+crick[3];
                        if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 ){
                                 T2G= sprintf("%.3f",varG/depth);
                                if(T2G>=minhetfreq ){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tPASS\tGT\tT2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGT\tT2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 T2G= sprintf("%.3f",varG/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tGT\tT2G\t".join("\t",@lines[3..6])."\n";
                        }	
                }
                if(lines[2] =~/C/i){#C>GT
			 qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
                         qvalueT=qual_crick[1];
                         varG=watson[3]+crick[3];
                         varT=crick[1];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=totaldepth-watson[1]-crick[0];
                        if(depth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 ){
                                 C2G= sprintf("%.3f",varG/depth);
                                 C2T= sprintf("%.3f",varT/depth);
                                if(C2G>=minhetfreq && C2T>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tGT\tgenoqual\tPASS\tGT\tC2G\,C2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tGT\tgenoqual\tLow\tGT\tC2G\,C2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 C2G= sprintf("%.3f",varG/depth);
                                 C2T= sprintf("%.3f",varT/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tGT\tgenoqual\tLow\tGT\tC2G\,C2T\t".join("\t",@lines[3..6])."\n";
                        }

                }
                if(lines[2] =~/G/i){#G>GT
			 qvalueT=qual_crick[1];
                         varT=crick[1];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=totaldepth;
                        if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 ){
                                 G2T= sprintf("%.3f",varT/depth);
                                if(G2T>=minhetfreq ){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tPASS\tGT\tG2T\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tGT\tG2T\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 G2T= sprintf("%.3f",varT/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tT\tgenoqual\tLow\tGT\tG2T\t".join("\t",@lines[3..6])."\n";
                        }
                }
        }   
#CG
	if(genotypemaybe eq "CG"){
                if(lines[2] =~/A/i){#A>CG
			 qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
                         qvalueC=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
                         varG=watson[3]+crick[3];
                         varC=crick[1]+watson[1];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=totaldepth-watson[1]-crick[0];
                        if(depth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 ){
                                 A2G= sprintf("%.3f",varG/depth);
                                 A2C= sprintf("%.3f",varC/depth);
                                if(A2G>=minhetfreq && A2C>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCG\tgenoqual\tPASS\tCG\tA2C\,A2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCG\tgenoqual\tLow\tCG\tA2C\,A2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 A2G= sprintf("%.3f",varG/depth);
                                 A2C= sprintf("%.3f",varC/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tCG\tgenoqual\tLow\tCG\tA2G\,A2C\t".join("\t",@lines[3..6])."\n";
                        }
	
                }
                if(lines[2] =~/T/i){#T>CG
			 qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
                         qvalueC=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
                         varG=watson[3]+crick[3];
                         varC=crick[1]+watson[1];
                        # depth=watson[0]+crick[0]+crick[1];
                         depth=totaldepth-watson[1]-crick[0];
                        if(depth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 ){
                                 T2G= sprintf("%.3f",varG/depth);
                                 T2C= sprintf("%.3f",varC/depth);
                                if(T2G>=minhetfreq && T2C>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCG\tgenoqual\tPASS\tCG\tT2C\,T2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tCG\tgenoqual\tLow\tCG\tT2C\,T2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                 T2G= sprintf("%.3f",varG/depth);
                                 T2C= sprintf("%.3f",varC/depth);
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tCG\tgenoqual\tLow\tCG\tT2G\,T2C\t".join("\t",@lines[3..6])."\n";
                        }
                }
                if(lines[2] =~/C/i){#C>CG
			 qvalueG=qual_watson[3];
                         varG=watson[3];
                        # depth=crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2];
                         depth=totaldepth-watson[1]-crick[0];
                         C2G= sprintf("%.3f",varG/depth);
                        if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 ){
                                if(C2G>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tPASS\tCG\tC2G\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tCG\tC2G\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tG\tgenoqual\tLow\tCG\tC2G\t".join("\t",@lines[3..6])."\n";
                        }			
                }
                if(lines[2] =~/G/i){#G>CG
			 qvalueC=qual_crick[2];
                         varC=crick[2];
                        # depth=crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2];
                         depth=totaldepth-watson[1]-crick[0];
                         G2C= sprintf("%.3f",varC/depth);
                        if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 ){
                                if(G2C>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tPASS\tCG\tG2C\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCG\tG2C\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCG\tG2C\t".join("\t",@lines[3..6])."\n";
                        }
	
                }
        }   
	unless(lines[2]=~/[ACGT]/i){
        print "[ACGT]\n";
		return 0;
		#print "lines[0]\tlines[1]\t\.\tlines[2]\t0\tSuper\tNN\t\.\t".join("\t",@lines[3..6])."\n";
	}

	
	
}
string snp_info;
{	// AA
	//AC
	if(genotype =="AC"){
		if(ref =="A"){
			qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            varC=watson[2]+crick[2];
            depth=watson[0]+crick[2]+watson[2];
			A2C=varC/totaldepth;
			sprintf(genotype_freq,"%.3f",A2C);
			alt = "C\t";
			if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 && A2C>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}else if(ref == "T"){
			qvalueA=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            varA=watson[0]+crick[0];
            varC=crick[2]+watson[2];
            depth=watson[0]+crick[2]+watson[2]+crick[0]+crick[1];
			T2A=varA/totaldepth;
            T2C=varC/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f",T2A,T2C);
            alt = "AC\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueC>=minquali && varA >=minread2 && varC>=minread2 && T2A>=minhetfreq && T2C>=minhetfreq){
            	filter = "PASS\t";
            }else{
            	filter = "Low\t";
            }
		}else if(ref == "C"){
			qvalueA=(qual_watson[0]>qual_crick[0])?qual_watson[0]:qual_crick[0];
            varA=watson[0]+crick[0];
            depth=watson[0]+crick[2]+watson[2];
            C2A=varA/totaldepth;
            sprintf(genotype_freq, "%.3f", C2A);
            alt = "A\t";
            if(depth >= mincover  && qvalueA >= minquali  && varA >=minread2 && C2A>=minhetfreq){
            	filter = "PASS\t";
            }else{
            	filter = "Low\t";
            }
		}else if(ref == "G"){
			qvalueA=qual_watson[0];
            qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            varA=watson[0];
            varC=crick[2]+watson[2];
            depth=watson[0]+crick[2]+watson[2]+crick[0];
            G2A=varA/totaldepth;
            G2C=varC/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f",G2A,G2C);
            alt = "AC\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueC>=minquali && varA>=minread2 && varC>=minread2 && G2A>=minhetfreq && G2C>=minhetfreq){
            	filter = "PASS\t";
            }else{
            	filter = "Low\t";
            }
		}
	}

	if(genotype =="AG"){
		if(ref == "A"){
			qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            varG=watson[3]+crick[3];
            depth=watson[0]+crick[3]+watson[3]+watson[1]+crick[1]+watson[2]+crick[2];
            A2G=varG/depth;
            sprintf(genotype_freq,"%.3f",A2G);
            alt = "G\t";
            if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 && A2G>=minhetfreq){
            	filter = "PASS\t";
            }else{
            	filter = "Low\t";
            }
		}else if(ref == "T"){         //T>AG
			qvalueA=qual_watson[0];
            qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            varA=watson[0];
            varG=crick[3]+watson[3];
            depth=watson[0]+crick[3]+watson[3]+crick[0]+crick[1]+watson[1];
            T2A=varA/totaldepth;
            T2G=varG/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f",T2A,T2G);
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
            varA=watson[0];
            varG=crick[3]+watson[3];
            depth=watson[0]+crick[2]+watson[2]+crick[0]+crick[1]+watson[1];
            C2A=varA/totaldepth;
            C2G=varG/totaldepth;
            sprintf(genotype_freq,"%.3f,%.3f",C2A,C2G);
            alt = "AG\t";
            if(depth >= mincover  && qvalueA >= minquali && qvalueG>=minquali && varA >=minread2 && varG>=minread2 && C2A>=minhetfreq && C2G>=minhetfreq){
            	filter = "PASS\t";
            }else{
            	filter = "Low\t";
            }            
		}else if (ref == "G"){
			/*G>AG*/
			qvalueA=qual_watson[0];
           	varA=watson[0];
            depth=crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2];
            G2A=varA/depth;
            sprintf(genotype_freq,"%.3f",G2A);
            alt = "A\t";
            if(depth >= mincover  && qvalueA >= minquali  && varA >=minread2 && G2A>=minhetfreq){
            	filter = "PASS\t";
            }else{
            	filter = "Low\t";
            }
		}
	}

	if(genotype =="TT"){
		if(ref =="A"){
			qvalueT=qual_crick[1];
            depth=watson[0]+crick[0]+crick[1];
            varT=crick[1];
            A2T = varT/totaldepth;
            sprintf(genotype_freq,"%.3f",A2T);
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
            depth=crick[1]+crick[3];
            varT=crick[1];
            C2T=0.;
            if(depth>0){
				C2T=varT/depth;
			}
			sprintf(genotype_freq,"%.3f",C2T);
			alt = "T\t";
			if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 && C2T>=minhomfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}else if(ref == "G"){
			qvalueT=qual_crick[1];
            depth=watson[3]+crick[3]+crick[1];
            varT=crick[1]+watson[1];
            G2T=varT/totaldepth;
			sprintf(genotype_freq,"%.3f",G2T);
			alt = "T\t";
			if(depth >= mincover  && qvalueT >= minquali && varT >=minread2 && G2T>=minhomfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}
	}

	if(genotype =="CC"){
		if(ref =="A"){
			qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
            depth=watson[0]+crick[0]+crick[2]+watson[2]+watson[3]+watson[3]+crick[1];
            varC=watson[2]+crick[2];
            A2C = varC/depth;
            sprintf(genotype_freq,"%.3f",A2C);
			alt = "C\t";
			if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 && A2C>=minhomfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}else if(ref =="T"){
			qvalueC=(qual_crick[2]>qual_watson[2])?qual_crick[2]:qual_watson[2];
            depth=crick[1]+crick[2]+watson[2]+watson[0]+crick[0]+watson[3]+crick[3];
            varC=watson[2]+crick[2];
            T2C = varC/depth;
            sprintf(genotype_freq,"%.3f",T2C);
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
            depth=watson[3]+crick[3]+crick[2]+watson[2]+watson[0]+crick[0]+crick[1];
            varC=watson[2]+crick[2];
            G2C = varC/depth;
            sprintf(genotype_freq,"%.3f",G2C);
			alt = "C\t";
			if(depth >= mincover  && qvalueC >= minquali && varC >=minread2 && G2C>=minhomfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}            
   	}   

   	if(genotype =="GG"){
		if(ref =="A"){
			qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
            depth=watson[0]+crick[3]+watson[3]+watson[1]+crick[1]+watson[2]+crick[2];
            varG=watson[3]+crick[3];
            A2G = varG/depth;
            sprintf(genotype_freq,"%.3f",A2G);
			alt = "G\t";
			if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 && A2G>=minhomfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}else if(ref =="T"){
			qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
            depth=crick[1]+watson[1]+crick[3]+watson[3]+watson[0]+watson[1]+crick[1];
            varG=watson[3]+crick[3];
            T2G = varG/depth;
            sprintf(genotype_freq,"%.3f",T2G);
			alt = "G\t";
			if(depth >= mincover  && qvalueG >= minquali && varG >=minread2 && T2G>=minhomfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}else if(ref =="C"){
			qvalueG=(qual_crick[3]>qual_watson[3])?qual_crick[3]:qual_watson[3];
            depth=crick[1]+watson[1]+crick[3]+watson[3]+watson[0]+watson[1]+crick[1];
            varG=watson[3]+crick[3];
            C2G = varG/depth;
            sprintf(genotype_freq,"%.3f",C2G);
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
 	
 	if(genotype =="CT"){
 		if(ref =="A"){
			qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            qvalueT=qual_crick[1];
            varC=watson[2]+crick[2];
            varT=crick[1];
			depth=totaldepth-watson[1];
            A2C = varC/depth;
            A2T = varC/depth;
            sprintf(genotype_freq,"%.3f,%.3f",A2C,A2T);
			alt = "CT\t";
			if(depth >= mincover  && qvalueC >= minquali && qvalueT >= minquali && varC >=minread2 && varT >=minread2 && A2C>=minhetfreq && A2T>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}else if(ref == "T"){
 			qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            varC=watson[2]+crick[2];
            depth=totaldepth-watson[1];
            T2C = varC/depth;
            sprintf(genotype_freq,"%.3f",T2C);
            alt = "C\t";
            if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 && T2C>=minhetfreq){
            	filter = "PASS\t";
			}else{
				filter ="Low\t";
            }
 		}else if(ref == "C"){
 			qvalueT=qual_crick[1];
            varT=crick[1];
            depth=totaldepth-watson[1];
            C2T=varT/depth;
            sprintf(genotype_freq,"%.3f",C2T);
            alt = "T\t";
            if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 && C2T>=minhetfreq){
            	filter = "PASS\t";
			}else{
				filter ="Low\t";
            }
 		}else if(ref == "G"){
 			qvalueC=(qual_watson[2]>qual_crick[2])?qual_watson[2]:qual_crick[2];
            qvalueT=qual_crick[1];
            varC=watson[2]+crick[2];
            varT=crick[1];
            depth=totaldepth-watson[1];
            G2C=varC/depth;
            G2T=varT/depth;
            sprintf(genotype_freq,"%.3f,%.3f",G2C, G2T);
            alt = "CT\t";
            if(depth >= mincover  && qvalueC >= minquali && qvalueT>=minquali && varC >=minread2 && varT>=minread2 && G2C>=minhetfreq && G2T>=minhetfreq){
            	filter = "PASS\t";
			}else{
				filter ="Low\t";
            }
 		}
	}

	if(genotype =="GT"){
		if(ref =="A"){
			qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueT=qual_crick[1];
            varG=watson[3]+crick[3];
            varT=crick[1];
            depth=totaldepth-watson[1]-crick[0];
            A2G = varG/depth;
            A2T = varT/depth;
            sprintf(genotype_freq,"%.3f,%.3f",A2G,A2T);
			alt = "GT\t";
			if(depth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 && A2G>=minhetfreq && A2T>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
		}else if(ref =="T"){
			qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            varG=watson[3]+crick[3];
            depth=watson[0]+watson[1]+watson[2]+watson[3]+crick[1]+crick[2]+crick[3];
            T2G = varG/depth;
            sprintf(genotype_freq,"%.3f",T2G);
			alt = "G\t";
			if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 && T2G>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}else if(ref =="C"){
			qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
			qvalueT=qual_crick[1];
            varG=watson[3]+crick[3];
            varT=crick[1];
            depth=totaldepth-watson[1]-crick[0];
            C2G = varG/depth;
            C2T = varT/depth;
            sprintf(genotype_freq,"%.3f,%.3f",C2G,C2T);
			alt = "GT\t";
			if(depth >= mincover  && qvalueG >= minquali && qvalueT>=minquali && varG >=minread2 && varT>=minread2 && C2G>=minhetfreq && C2T>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}else if(ref =="G"){
			qvalueT=qual_crick[1];
			varT=crick[1];
			depth=totaldepth;
            G2T = varT/depth;
            sprintf(genotype_freq,"%.3f",T2G);
			alt = "T\t";
			if(depth >= mincover  && qvalueT >= minquali  && varT >=minread2 && G2T>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}
	}

	if(genotype =="CG"){
		if(ref =="A"){
			qqvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueC=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
            varG=watson[3]+crick[3];
            varC=crick[1]+watson[1];
            depth=totaldepth-watson[1]-crick[0];
            A2G = varG/depth;
            A2C = varC/depth;
            sprintf(genotype_freq,"%.3f,%.3f",A2C,A2G);
			alt = "CG\t";
			if(depth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 && A2G>=minhetfreq && A2C>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}else if(ref =="T"){
			qvalueG=(qual_watson[3]>qual_crick[3])?qual_watson[3]:qual_crick[3];
            qvalueC=(qual_watson[1]>qual_crick[1])?qual_watson[1]:qual_crick[1];
            varG=watson[3]+crick[3];
            varC=crick[1]+watson[1];
            depth=totaldepth-watson[1]-crick[0];
            T2G = varG/depth;
            T2C = varC/depth;
            sprintf(genotype_freq,"%.3f,%.3f",T2C,T2G);
			alt = "CG\t";
			if(depth >= mincover  && qvalueG >= minquali && qvalueC>=minquali && varG >=minread2 && varC>=minread2 && T2G>=minhetfreq && T2C>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}else if(ref =="C"){
			qvalueG=qual_watson[3];
            varG=watson[3];
            depth=totaldepth-watson[1]-crick[0];
            C2G = varG/depth;
            sprintf(genotype_freq,"%.3f",C2G);
			alt = "G\t";
			if(depth >= mincover  && qvalueG >= minquali  && varG >=minread2 && C2G>=minhetfreq){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}else if(ref =="G"){
			qvalueC=qual_crick[2];
            varC=crick[2];
            depth=totaldepth-watson[1]-crick[0];
            G2C = varC/depth;
            sprintf(genotype_freq,"%.3f",G2C);
			alt = "C\t";
			if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2){
				filter = "PASS\t";
			}else{
				filter ="Low\t";
			}
 		}
	}


#CG
	if(genotypemaybe eq "CG"){

                if(lines[2] =~/G/i){#G>CG
			 qvalueC=qual_crick[2];
                         varC=crick[2];
                        # depth=crick[3]+watson[3]+watson[0]+watson[1]+crick[1]+watson[2]+crick[2];
                         depth=totaldepth-watson[1]-crick[0];
                         G2C= sprintf("%.3f",varC/depth);
                        if(depth >= mincover  && qvalueC >= minquali  && varC >=minread2 ){
                                if(G2C>=minhetfreq){
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tPASS\tCG\tG2C\t".join("\t",@lines[3..6])."\n";
                                }else{
                                        print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCG\tG2C\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
                                print "lines[0]\tlines[1]\t\.\tlines[2]\tC\tgenoqual\tLow\tCG\tG2C\t".join("\t",@lines[3..6])."\n";
                        }
	
                }
        }   
	unless(lines[2]=~/[ACGT]/i){
        print "[ACGT]\n";
		return 0;
		#print "lines[0]\tlines[1]\t\.\tlines[2]\t0\tSuper\tNN\t\.\t".join("\t",@lines[3..6])."\n";
	}

	
	