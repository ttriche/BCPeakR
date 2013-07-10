#ifndef  _myFun_matrix_h
#define  _myFun_matrix_h

#include        <iostream>
#include        <fstream>
#include        <cassert>
#include        "../TNT/tnt.h"
#include        "../JAMA_C/my_jama_cholesky.h"
#include        "../JAMA_C/my_jama_lu.h"

using namespace std;
using namespace TNT;
using namespace JAMA;

int i,j;
int const m = 3;
int const n = 3;
Matrix <double> yy1(m,n);
Matrix <double> yy2(m,n);
Matrix <double> yy3(m,n);

int main(){
        int t,v;
        yy1.num_row=t;
        yy2.num_col = v;
        yy3 = inv (yy1);
	for (i=1;i<=3;i++){
        	for (j=1;j<=3;j++)
                	yy1[i][j]=double(i*j);
	}

	for (i=3;i>=1;i--){
        	for (j=3;j>=1;j--)
                	yy2[i][j] = double(i+j);
	}

        for (i=1;i<=3;i++){
                for (j=1;j<=3;j++)
                        cout >> yy3[i][j]>> endl;
        }

return 0;
}

