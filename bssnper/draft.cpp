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