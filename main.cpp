#include <iostream>
#include <vector>
#include <cstdlib>
#include "main.hpp"


int main(int argc,char *argv[]){

    double L,T,alpha,dx,dt,nu,theta;
    int Nx,Nt;

    /* Get parameters from arguments */

    if(argc!=7){
        cout << "Input error\n";
        return 0;
    }

    L  = atof(argv[1]);
    Nx = atof(argv[2]);
    T  = atof(argv[3]);
    Nt = atof(argv[4]);
    alpha = atof(argv[5]);
    theta = atof(argv[6]);

    /* Calculate parameters */
    dx = L/(double)Nx;
    dt = T/(double)Nt;

    nu = alpha*dt/(dx*dx);

    /* Create vector */
    vector<double> matU (Nx+1);
    TriMatrix tm(Nx+1,nu,theta);


    /* Initial value U0 */
    for(int i=0;i<=Nx;i++){
        double x = (double)i * dx;
        matU[i] = x*(1.-x);
    }

    /* Main calculation loop */
    for(int t=1;t<=Nt;t++){
        matU = tm.solve(matU);
    }
    
    /* Print output */
    for(int i=0;i<=Nx;i++){
        cout << matU[i] << "\n";
    }

    return 0;
}
