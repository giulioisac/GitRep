#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>
using namespace arma;
using namespace std;
/*
Solving an eigenvalue problem with jacobi's method
Particles in central potential field
2particles: omega= 0.01, 0.5, 1, 5;
*/
void maximum(mat A, int N, int &col, int &row);
void jacobi(mat &A, mat &R, int N, int &col, int &row);
void Matrix(mat &A, double h, int N);
void jacobi_matrix(mat &A, mat &R, int N, int &col, int &row);
void minimumEigenvalues(mat A,int N);
void TEST_EIGENVALUES(mat &A, double h, int N);
void exit(mat &A, mat &R, double h, int N, double tolerance,int max);
void for_rotation(mat &A, mat &R, double h, int N, double tolerance,int max);
void for_Algo(mat &A, mat &R, double h, int N, double tolerance,int max);
void TEST_ORTHOGONALITY(mat &R,int N);
void Matrix2(mat &A, double h, int N);
void exit_test(mat &A, mat &R, double h, int N);
void TEST_MAX(mat A, int N);
int main()
{
    int n = 2*pow(10,2);
    int max = 100000;
    double tolerance = pow(10,-4);
    double rmax = 6.5;
    double h = rmax/n;
    int N = n-1;
    mat A(N,N);
    mat R(N,N);

    //DIFFERENT FUNCTIONS TO CALL IN ORDER TO GET WHAT YOU WANT

    //for_rotation(A,R,h,N,tolerance,max);
    for_Algo(A,R,h,N,tolerance,max);
    TEST_ORTHOGONALITY(R,N);
    TEST_EIGENVALUES(A,h,N);
    //exit(A,R,h,N,tolerance,max);
    //exit_test(A,R,h,N);
    return 0;
}
void exit_test(mat &A, mat &R, double h, int N){
    Matrix2(A,h,N);
    R.eye();
    vec eigval(N);
    Matrix2(A,h,N);
    eig_sym(eigval, R, A);
    A=diagmat(eigval);
    //cout<<A<<endl;
    //position of first eigenvalue
    //out of ro
    double ro=h;
    for (int i=0;i<N;i++){
        cout<<ro<<endl;
        ro+=h;
    }
    cout<<"change"<<endl;
    //out of first eigenvector
    ro=h;
    for(int i=0; i<N;i++){
        cout<<R(i,0)*R(i,0)<<endl;
        ro+=h;
    }
    cout<<"change"<<endl;
    //out second eigenvector;
    ro=h;
    for(int i=0; i<N;i++){
        cout<<R(i,1)*R(i,1)<<endl;
        ro+=h;
    }

    cout<<"change"<<endl;
    //out second eigenvector;
    ro=h;
    for(int i=0; i<N;i++){
        cout<<R(i,2)*R(i,2)<<endl;
        ro+=h;
    }

    cout<<"change"<<endl;
    cout<<"break"<<endl;
}
void Matrix2(mat &A, double h, int N){
    vec v(N);
    A.zeros();
    //potential
    double omega=1.0;
    for(int i=0;i<N;i++){
        v(i)=omega*omega*(i+1)*h*(i+1)*h+1/((i+1)*h);
    }
    for(int i=1; i<N-1;i++){
        A(i,i)= 2/(h*h)+ v(i);
        A(i,i+1)= -1/(h*h);
        A(i,i-1)= -1/(h*h);
    }
    A(0,0)= 2/(h*h)+ v(0);
    A(0,1)= -1/(h*h);
    A(N-1,N-1)= 2/(h*h)+ v(N-1);
    A(N-1,N-2)= -1/(h*h);
}
void maximum(mat A, int N, int &col, int &row){
    for(int i=0; i<N ; i++){
        for(int j=i+1; j<N ; j++){
            if (fabs(A(i,j))>fabs(A(row,col))){
                row=i;
                col=j;
            }
        }
    }
}
void jacobi(mat &A, mat &R, int N, int &col, int &row){

    double s,c;
//    if(A(row,col)!=0){
//    double tau,theta;
//    tau=(A(col,col)-A(row,row))/(2*A(row,col));
//    theta=atan(1.0/tau)/2.0;
//    c=cos(theta);
//    s=sin(theta);
//    }
//    else {
//        c = 1.0;
//        s = 0.0;
//      }
    if ( A(row,col) != 0.0 ) {
    double t, tau;
    tau = (A(col,col) - A(row,row))/(2*A(row,col));
    if ( tau > 0 ) {
    t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
    t = -1.0/( -tau + sqrt(1.0 + tau*tau));
    }
    c = 1/sqrt(1+t*t);
    s = c*t;
    } else {
    c = 1.0;
    s = 0.0;
    }


    //ROTATIONS!

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;

    //diagonal elements
    a_kk= A(row,row);
    a_ll= A(col,col);
    A(row,row) = c*c*a_kk - 2.0*c*s*A(row,col) + s*s*a_ll;
    A(col,col) = s*s*a_kk + 2.0*c*s*A(row,col) + c*c*a_ll;

    // hard-coding non-diagonal elements by hand
    A(row,col) = 0.0;
    A(col,row) = 0.0;
    // other elements

    for ( int i = 0; i < N; i++ ) {
    if ( i != row && i != col ) {
      a_ik = A(i,row);
      a_il = A(i,col);

      A(i,row) = c*a_ik - s*a_il;
      A(row,i) = A(i,row);
      A(i,col) = c*a_il + s*a_ik;
      A(col,i) = A(i,col);
      }

    //  And finally the new eigenvectors
        r_ik = R(i,row);
        r_il = R(i,col);

        R(i,row) = c*r_ik - s*r_il;
        R(i,col) = c*r_il + s*r_ik;
     }
}
void Matrix(mat &A, double h, int N){
    vec v(N);
    A.zeros();
    //potential
    for(int i=0;i<N;i++){
        v(i)=(i+1)*h*(i+1)*h;
    }
    for(int i=1; i<N-1;i++){
        A(i,i)= 2/(h*h)+ v(i);
        A(i,i+1)= -1/(h*h);
        A(i,i-1)= -1/(h*h);
    }
    A(0,0)= 2/(h*h)+ v(0);
    A(0,1)= -1/(h*h);
    A(N-1,N-1)= 2/(h*h)+ v(N-1);
    A(N-1,N-2)= -1/(h*h);
}
void jacobi_matrix(mat &A, mat &R, int N, int &col, int &row){

    //Find Angle of rotation
    double tau,theta,cose, sen;

    tau=(A(col,col)-A(row,row))/(2*A(row,col));
    theta=atan(1.0/tau)/2.0;
    cose=cos(theta);
    sen=sin(theta);

    //build  matrix of rotation
    mat J(N,N);
    J.eye();
    J(row,col)=sen;
    J(col,row)=-sen;
    J(col,col)=cose;
    J(row,row)=cose;

    //ROTATE!
    A=J.t()*A*J;
    R=J.t()*R;
}
void minimumEigenvalues(mat A,int N){
    int j=N-1;
    cout<<"EIGENVALUES AFTER JACOBI METHOD: "<<endl;
    for(int i=0;i<N;i++){
        if(fabs(A(i,i))<fabs(A(j,j))){
            j=i;
        }
    }
    cout<<A(j,j)<<endl;
    int k=N-1;
    for(int i=0;i<N;i++){
        if((fabs(A(i,i))<fabs(A(k,k))) && (i!=j)){
            k=i;
        }
    }
    cout<<A(k,k)<<endl;
    int l=N-1;
    for(int i=0;i<N;i++){
        if(fabs(A(i,i))<fabs(A(l,l)) && i!=j && i!=k){
            l=i;
        }
    }
    cout<<A(l,l)<<endl;
}
void TEST_EIGENVALUES(mat &A, double h, int N){
    cout<<"-------------------TEST EIGENVALUES------------------"<<endl;
    clock_t start, finish;
    vec eigval(N);
    Matrix2(A,h,N);
    start = clock();
    eig_sym(eigval,A);
    finish = clock();
    float time=((float)finish - (float)start) / CLOCKS_PER_SEC;
    cout<< time <<" seconds."<<endl;
    cout<<"test eigenvalues  matrix: "<<endl;
    cout<<eigval(0)<<endl;
    cout<<eigval(1)<<endl;
    cout<<eigval(2)<<endl;
}
void exit(mat &A, mat &R, double h, int N, double tolerance,int max){
    int counter=0;
    int col = 0;
    int row = 1;
    Matrix2(A,h,N);
    R.eye();
    maximum(A,N,col,row);
    while (fabs(A(row,col))>tolerance && counter<max){
    jacobi(A,R,N,col,row);
    counter++;
    maximum(A,N,col,row);
    }
//    vec eigval(N);
//    Matrix2(A,h,N);
//    eig_sym(eigval, R, A);
    //position of first eigenvalue
    int m=0;
    for (int i=0;i<N;i++){
        if (fabs(A(i,i))<fabs(A(m,m))){m=i;}
    }
    //out of ro
    double ro=h;
    for (int i=0;i<N;i++){
        cout<<ro<<endl;
        ro+=h;
    }
    cout<<"change"<<endl;
    //out of first eigenvector
    ro=h;
    for(int i=0; i<N;i++){
        cout<<R(i,m)*R(i,m)<<endl;
        ro+=h;
    }
    cout<<"change"<<endl;
    //position second eigenvalue
    int k=0;
    for (int i=0;i<N;i++){
        if (fabs(A(i,i))<fabs(A(k,k))&& A(i,i)!=A(m,m)){k=i;}
    }
    //out second eigenvector;
    ro=h;
    for(int i=0; i<N;i++){
        cout<<R(i,k)*R(i,k)<<endl;
        ro+=h;
    }

    cout<<"change"<<endl;
    //position third eigenvalue
    int l=0;
    for (int i=0;i<N;i++){
        if (fabs(A(i,i))<fabs(A(l,l))&& A(i,i)!=A(m,m)&&A(i,i)!=A(k,k)){l=i;}
    }
    //out second eigenvector;
    ro=h;
    for(int i=0; i<N;i++){
        cout<<R(i,l)*R(i,l)<<endl;
        ro+=h;
    }

    cout<<"change"<<endl;
    cout<<"break"<<endl;
}
void for_rotation(mat &A, mat &R, double h, int N, double tolerance,int max){
    int counter=0;
    int col = N-1;
    int row = N-2;
    Matrix(A,h,N);
    R.eye();
    cout<<"--------------JACOBI ROTATION MATRIX--------------"<<endl;
    clock_t start, finish;
    start = clock();
    maximum(A,N,col,row);
    while (fabs(A(row,col))>tolerance && counter<max){
    jacobi_matrix(A,R,N,col,row);
    counter++;
    maximum(A,N,col,row);
    }
    finish = clock();
    float time=((float)finish - (float)start) / CLOCKS_PER_SEC;
    cout<<counter<<" ripetitions and  "<< time <<" seconds."<<endl;
    minimumEigenvalues(A,N);
    return;
}
void for_Algo(mat &A, mat &R, double h, int N, double tolerance,int max){
    int counter=0;
    int col = 0;
    int row = 1;

    Matrix2(A,h,N);
    R.eye();
    cout<<"--------------JACOBI ROTATION ALGORITHM--------------"<<endl;
    clock_t start, finish;
    start = clock();
    maximum(A,N,col,row);
    while (fabs(A(row,col))>tolerance && counter<max){
    jacobi(A,R,N,col,row);
    counter++;
    maximum(A,N,col,row);
    //TEST_MAX(A,N);
    }
    finish = clock();
    float time=((float)finish - (float)start) / CLOCKS_PER_SEC;
    cout<<counter<<" ripetitions and  "<< time <<" seconds."<<endl;
    minimumEigenvalues(A,N);
    //    cout<<"post-algorithm"<<endl;
    //    cout<<A<<endl;
    return;
}
void TEST_ORTHOGONALITY(mat &R,int N){
    bool err=true;
    cout<<endl<<"-------------------TEST ORTHOGONALITY----------------"<<endl;
    double tol=pow(10,-3);
    for(int i=0;i<N;i++){
        if(norm(R.col(i))>1+tol || norm(R.col(i))<1-tol){
            cout<<"errore: colonna "<<i<<" non ha norma 1"<<endl;
            err=false;
        }
    }
    for(int i=0;i<N;i++){
        for(int j=1;j<N;j++){
            if(i!=j && (dot( R.col(i) , R.col(j) )> 0.0 +tol || dot( R.col(i) , R.col(j) )<0.0-tol )){
                cout<<"column "<<i<<" and "<<j<<" are not orthogonal"<<endl;
                err=false;
            }
        }
    }

    if(err){cout<<"The matrix of eigenvectors is orthogonal."<<endl;}

    cout<<endl;
}
void TEST_MAX(mat A, int N){
    int col=N-2 ,row=N-3;
    maximum(A,N,col,row);
    A=A-diagmat(A);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            A(i,j)=abs(A(i,j));
        }
    }
    if(abs(A(row,col))!=max(max(A))){
        cout<<"wrong maximum"<<endl;
    }
}
