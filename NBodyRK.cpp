#include<iostream>
#include<cmath>
#include <string>
#include<fstream>
#include<sstream>
#include <boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include<boost/shared_ptr.hpp>
#include"Array.hpp"
#include <ctime>

using namespace boost::numeric::ublas;
using namespace boost;

#define G 4*(M_PI)*(M_PI)
#define R_shock 0.000036
#define Eps_H 0.00000000000005

typedef Array<vector<double> > arvec;

shared_ptr<vector<double> > F(double t, arvec R, Array<double> m, Array<double> norm_r, unsigned j){
    shared_ptr<vector<double> > Derivative;

    Derivative.reset(new vector<double>(R[j].size()));
    for(int i=0;i<3;i++) (*Derivative)(i)=R[j](i+3);
    for(int i=0;i<3;i++){
        (*Derivative)(i+3)=0;
        for(unsigned k=0;k<j;++k) 
            (*Derivative)(i+3)+=G*m[k]*(R[k](i)-R[j](i))/(norm_r[k]*norm_r[k]*norm_r[k]);
        for(unsigned k=j+1;k<R.size();++k) 
            (*Derivative)(i+3)+=G*m[k]*(R[k](i)-R[j](i))/(norm_r[k]*norm_r[k]*norm_r[k]);
    }
    return Derivative;
}

Array<Array<double> >norm_r(arvec R){
    Array<Array<double> > norm_R(R.size(),Array<double>(R.size()));

   for(unsigned i=0;i<R.size();++i){
        for(unsigned j=0;j<i;++j)
            norm_R[i][j]=norm_R[j][i];
        for(unsigned j=i+1;j<R.size();++j)
            norm_R[i][j]=sqrt((R[i](0)-R[j](0))*(R[i](0)-R[j](0))+(R[i](1)-R[j](1))*(R[i](1)-R[j](1))+
                                                +(R[i](2)-R[j](2))*(R[i](2)-R[j](2)));
    }
    return norm_R;
}


double EnergyIntegral(arvec R, Array<Array<double> > norm_r, Array<double> m){
    double E=0;
    for(unsigned i=0;i<m.size();++i){
        for(unsigned j=i+1;j<m.size();++j)
            E+=-G*m[j]*m[i]/(norm_r[i][j]);
        E+=m[i]*(R[i](3)*R[i](3)+R[i](4)*R[i](4)+R[i](5)*R[i](5))/2;
    }    
    return E;
}

Array<arvec > RungeKutta(double& a, arvec X0, Array<double> m, unsigned &n, 
    shared_ptr<vector<double> > F(double, arvec, Array<double>, Array<double>, unsigned)){
    Array<arvec> A(n+1,arvec(X0.size(),vector<double>(X0[0].size())));
    Array<Array<double> > norm_A;
    double H0,H1;
    unsigned i=0;
    arvec k1(X0.size(),vector<double>(X0[0].size())),k2(X0.size(),vector<double>(X0[0].size())),
    k3(X0.size(),vector<double>(X0[0].size())),k4(X0.size(),vector<double>(X0[0].size()));
    double t=0.,h=a/n;
    A[0]=X0;
    H0=EnergyIntegral(A[0], norm_r(A[0]),m);
    unsigned t1=clock();
    while(i<n){
        norm_A=norm_r(A[i]);
        for(unsigned j=0; j<X0.size();++j)
            k1[j]=h*(*F(t,A[i],m,norm_A[j],j));

        norm_A=norm_r(A[i]+k1/2.);
        for(unsigned j=0; j<X0.size();++j)    
            k2[j]=h*(*F(t+h/2,A[i]+k1/2.,m,norm_A[j],j));
        
        norm_A=norm_r(A[i]+k2/2.);
        for(unsigned j=0; j<X0.size();++j)    
            k3[j]=h*(*F(t+h/2,A[i]+k2/2.,m,norm_A[j],j));
        
        norm_A=norm_r(A[i]+k3);
        for(unsigned j=0; j<X0.size();++j)    
            k4[j]=h*(*F(t+h,A[i]+k3,m,norm_A[j],j));
        
        A[i+1]=A[i] + (k1+2.*k2+2.*k3+k4)/6.;
        H1=EnergyIntegral(A[i+1], norm_r(A[i+1]),m);
        if(fabs(H1-H0)>Eps_H){
            A.subjoin(n-i,arvec(X0.size()));
            n=2*n-i; h=h/2;
        }else{
            if(sqrt((A[i+1][0](0)-A[i+1][1](0))*(A[i+1][0](0)-A[i+1][1](0))+(A[i+1][0](1)-A[i+1][1](1))*
            (A[i+1][0](1)-A[i+1][1](1))+(A[i+1][0](2)-A[i+1][1](2))*(A[i+1][0](2)-A[i+1][1](2))) < R_shock){
                std::cout << sqrt((A[i+1][0](0)-A[i+1][1](0))*(A[i+1][0](0)-A[i+1][1](0))+(A[i+1][0](1)-A[i+1][1](1))*
                (A[i+1][0](1)-A[i+1][1](1))+(A[i+1][0](2)-A[i+1][1](2))*(A[i+1][0](2)-A[i+1][1](2))) << std::endl;
                std::cout << t << std::endl;
                a=t; n=i+1;
                unsigned t2=clock();
                std::cout << "t2-t1= " << (t2-t1)/1000. << std::endl; 
                return A;
            }
            t+=h;
            i++;
        }
        H0=H1;    
    }
    unsigned t2=clock();
    std::cout << "t2-t1= " << (t2-t1)/1000. << std::endl;
    return A;
}

const char* filenamestr(const char* s1, unsigned j, const char* s2){
    std::stringstream ss;
    ss << s1 << j << s2;
    return (ss.str()).c_str();
}



int main(){
    unsigned t0=clock();
    Array<arvec> A;
    std::ofstream out,energy;
    std::ifstream file_size,file_mass,file_initial;
    arvec X0;
    Array<double> m;
    unsigned N,M;
    double B;


    file_size.open("file_size.dat",std::ios_base::in);
    file_size >> B >> N >> M;
    file_size.close();
    
    X0=arvec(M,vector<double>(6));
    m=Array<double>(M);

    file_mass.open("file_mass.dat",std::ios_base::in);
    for(unsigned i=0;i<m.size();++i){
        file_mass >> m[i];
    }    
    file_mass.close();
 
    file_initial.open("file_initial.dat",std::ios_base::in);
    for(unsigned i=0;i<M;++i)
        for(int j=0;j<6;++j)
            file_initial >> X0[i](j);
    file_initial.close();

    A=RungeKutta(B,X0,m,N,F);

    unsigned t3=clock();
    std::cout << "t3-t0= " << (t3-t0)/1000. << std::endl;
    
    for(unsigned j=0;j<M;++j){
        out.open(filenamestr("NBodyRK_",j,".dat"),std::ios_base::trunc);
        for(unsigned i=0;i<=N;++i){
            for(unsigned k=0;k<A[i][j].size();++k)
                out << A[i][j](k) << " ";
            out << std::endl;
        }        
        out.close();
    }
    
    
    energy.open("NBodyRK_energy.dat",std::ios_base::trunc);
    for(unsigned i=0;i<=N;++i)
        energy <<  "h(" << i << ")= " << EnergyIntegral(A[i], norm_r(A[i]),m) << std::endl;
    std::cout << fabs(EnergyIntegral(A[0], norm_r(A[0]),m)-EnergyIntegral(A[N], norm_r(A[N]),m)) << std::endl;
    energy.close();
    std::cout << "N=" << N << std::endl;
    return 0;
}