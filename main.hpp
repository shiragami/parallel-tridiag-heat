/* Header for main3.cpp */
#include <iostream>
#include <vector>
#include <cstdlib>
#include <pthread.h>
#include <thread>

using namespace std;


class TriMatrix{
        int N;
        int Nmid;

        /* TDM1 */
        vector<double> Dbot;
        vector<double> Dmid;
        vector<double> Dtop;

        /* TDM1 */
        vector<double> D2bot;
        vector<double> D2mid;
        vector<double> D2top;

    public:
        TriMatrix(int,double,double);
        vector<double> cross(vector<double>);
        vector<double> solve(vector<double>);

        void parallelcross(int start,int end,double *pX,double *pY);
        void parallelsolve_first(int process,double *m1,double *dm);
        void parallelsolve_second(int process,double *m1,double *m2,double *dm);

        void print();
};


/* Multi-threading matrix solver : First part */
void TriMatrix::parallelsolve_first(int process,double *m1,double *dm){
    if(process==1){
        /* Process 1 : top-bot */
        for(int i=1;i<=Nmid;i++){
            double m  = Dbot[i]/ *(dm+i-1);
            *(dm+i)  -= Dtop[i-1]*m;
            *(m1+i)  -= *(m1+i-1)*m;
        }
    }else{
        /* Process 2 : bot-top */
        for(int i=N-2;i>Nmid;i--){
            double m  = Dtop[i]/ *(dm+i+1);
            *(dm+i)  -= Dbot[i+1]*m;
            *(m1+i)  -= *(m1+i+1)*m;
        }
    }
}

/* Multi-threading matrix solver : Second part */
void TriMatrix::parallelsolve_second(int process,double *m1,double *m2,double *dm){
    if(process==1){
        /* Process 1 : bot-top */
        for(int i=Nmid;i>=0;i--){
            *(m2+i) = ( *(m1+i) - *(m2+i+1)*Dtop[i] ) / *(dm+i);
        }
    }else{
        /* Process 2 : top-bot */
        for(int i=Nmid+2;i<N;i++){
            *(m2+i) = ( *(m1+i) -  *(m2+i-1)*Dbot[i] ) / *(dm+i);
        }
    }
}

/* Multi-threading matrix multiplication*/
void TriMatrix::parallelcross(int start,int end,double *pX,double *pY){
    for(int i=start;i<end;i++){
        *(pY+i-start) = *(pX+i-1)*D2bot[i] + *(pX+i)*D2mid[i] + *(pX+i+1)*D2top[i];
    }
}

/* Solve the inverse matrix */
vector<double> TriMatrix::solve(vector<double> matX){

    /* Calculate the cross of TDM2 and matX */
    vector<double> matU2 (N);
    vector<double> matU1 = cross(matX);
    vector<double> mid = Dmid;

    double *m1 = matU1.data();
    double *m2 = matU2.data();
    double *dm = mid.data();

    /* Parallel TDMA algorithm */

    /* Run the first half */
    thread thread1a (&TriMatrix::parallelsolve_first,this,1,m1,dm);
    thread thread2a (&TriMatrix::parallelsolve_first,this,2,m1,dm);

    /* Wait both threads to finish */
    thread1a.join();
    thread2a.join();

    /* Solve middle point */

    double k1 = mid[Nmid];
    double k2 = Dtop[Nmid];
    double k3 = Dbot[Nmid+1];
    double k4 = mid[Nmid+1];
    double det = k1*k4 - k2*k3;

    matU2[Nmid]   = ( k4*matU1[Nmid] - k2*matU1[Nmid+1])/det;
    matU2[Nmid+1] = (-k3*matU1[Nmid] + k1*matU1[Nmid+1])/det;

    thread thread1b (&TriMatrix::parallelsolve_second,this,1,m1,m2,dm);
    thread thread2b (&TriMatrix::parallelsolve_second,this,2,m1,m2,dm);

    thread1b.join();
    thread2b.join();

    return matU2;
}


/* Print content of matT for debugging */
void TriMatrix::print(){
    cout << Dmid[0] << endl;
    for(int i=1;i<N;i++){
        cout << Dbot[i-1] << " " << Dmid[i] << " " << Dtop[i-1] << endl;
    }
    cout << Dmid[N-1] << "\n";
}


/* Hard coded tri-diag matrix multiplication */
vector<double> TriMatrix::cross(vector<double> matX){

    /* Y = D*X */
    vector<double> matY (N);
    vector<double> matY1 (Nmid);
    vector<double> matY2 (Nmid+1);

    /* Assign pointers to vectors */
    double *pY1 = matY1.data();
    double *pY2 = matY2.data();
    double *pX  = matX.data();

    /* First element */
    matY[0] = matX[0]*D2mid[0] + matX[1]*D2top[0];

    /* Run calculation in 2 parallel threads */
    thread thread1 (&TriMatrix::parallelcross,this,1,Nmid,pX,pY1);
    thread thread2 (&TriMatrix::parallelcross,this,Nmid,N-1,pX,pY2);

    /* Wait both threads to finish */
    thread1.join();
    thread2.join();
 
    /* Merge element */
    for(int i=1;i<Nmid;i++) matY[i] = matY1[i-1];
    for(int i=Nmid;i<N-1;i++) matY[i] = matY2[i-Nmid];

    /* Last element */
    matY[N-1] = matX[N-1]*D2mid[N-1] + matX[N-2]*D2bot[N-1];

    return matY;
}

/* Init */
TriMatrix::TriMatrix(int Nx,double nu,double theta){
    N = Nx;
    Nmid = (N-1)/2;    

    /* Resize matrix */
    Dbot.resize(N);
    Dmid.resize(N);
    Dtop.resize(N);

    D2bot.resize(N);
    D2mid.resize(N);
    D2top.resize(N);

    /* Initialize matrix 1 */
    Dmid[0] = 1.; Dtop[0] = 0.;

    for(int i=1;i<N-1;i++){
        Dbot[i] = -nu*theta;
        Dmid[i] = 1. + 2.*nu*theta;
        Dtop[i] = -nu*theta;
    }

    Dmid[N-1] = 1.; Dbot[N-1] = 0.;


    /* Initialize matrix 2*/
    D2mid[0] = 1.; D2top[0] = 0.;

    for(int i=1;i<N-1;i++){
        D2bot[i] = nu*(1.-theta);
        D2mid[i] = 1. - 2.*nu*(1.-theta);
        D2top[i] = nu*(1.-theta);
    }
    D2mid[N-1] = 1.; D2bot[N-1] = 0.;

}
