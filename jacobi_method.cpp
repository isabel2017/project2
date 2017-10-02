#include "jacobi_method.h"



// ***********************find the max nondiag element(value and index)
//A is a symmetrical matrix, only upper triangular elements are compared here
double max_nondiag(double** A,int* k, int* l,int n)
{
    double max=0;
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            if(fabs(A[i][j])>max)
           {
                max=fabs(A[i][j]);
                *k=i;
                *l=j;
            }
        }
    }
    return max;
}
//********************jacobi rotation
void jacobi_rotation(double** A,int k, int l,int n)
{
    double a_kk=A[k][k];
    double a_ll=A[l][l];
    double a_kl=A[k][l];
    double a_ik,a_il;
    double t,s,c;

//calculate rotation angel
    if(a_kl!=0){

        double tau=(a_ll-a_kk)/(2*a_kl);
        //to make fabs(angel)<=pi/4,should choose smaller fabs(t)
        if(tau>0){
            t=-tau+sqrt(1+pow(tau,2));
        }else{
            t=-tau-sqrt(1+pow(tau,2));
        }
            c=1/sqrt(1+pow(t,2));
            s=t*c;
    }
// make rotation to A
    A[k][k]=a_kk*c*c-2*a_kl*c*s+a_ll*s*s;
    A[l][l]=a_ll*c*c+2*a_kl*c*s+a_kk*s*s;
    A[k][l]=0;
    A[l][k]=0;
    for(int i=0;i<n;i++){
        if (i!=k&&i!=l){
            a_ik=A[i][k];
            a_il=A[i][l];
            A[i][k]=a_ik*c-a_il*s;
            A[i][l]=a_il*c+a_ik*s;
            A[k][i]=A[i][k];
            A[l][i]=A[i][l];
        }
    }
return;
}



void jacobi_method(double** A, int n)
{
    int k,l;
    double epsilon=1.0e-8;
    int iteration_count=0;
    int max_iteration=(double)n*(double)n*(double)n;
    double max_nondiag_element=max_nondiag(A,&k, & l,n);
    while(max_nondiag_element>epsilon && (double)iteration_count<max_iteration){
        max_nondiag_element=max_nondiag(A,&k, &l,n);
        jacobi_rotation(A,k,l,n);
        iteration_count++;

    }
    cout<<"iteration_count="<<iteration_count<<endl;


}


