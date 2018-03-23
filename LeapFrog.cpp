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
    for(int i=0;i<3;i++){
        (*Derivative)(i)=0;
        for(unsigned k=0;k<j;++k) 
            (*Derivative)(i)+=G*m[k]*(R[k](i)-R[j](i))/(norm_r[k]*norm_r[k]*norm_r[k]);
        for(unsigned k=j+1;k<R.size();++k) 
            (*Derivative)(i)+=G*m[k]*(R[k](i)-R[j](i))/(norm_r[k]*norm_r[k]*norm_r[k]);
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


double EnergyIntegral(arvec V, Array<Array<double> > norm_r, Array<double> m){
    double E=0;
    for(unsigned i=0;i<m.size();++i){
        for(unsigned j=i+1;j<m.size();++j)
            E+=-G*m[j]*m[i]/(norm_r[i][j]);
        E+=m[i]*(V[i](0)*V[i](0)+V[i](1)*V[i](1)+V[i](2)*V[i](2))/2;
    }    
    return E;
}

template <typename T>
void flatten(const T& value, std::ostream& out=std::cout){
    out << value << " ";
}

template <typename T>
void flatten(const Array<T>& array, std::ostream& out=std::cout)
{ 
    for(unsigned i=0; i<array.size(); ++i)
        flatten(array[i],out);
}

void LeapFrog(double& a, Array<double> m, unsigned &n, 
    shared_ptr<vector<double> > F(double, arvec, Array<double>, Array<double>, unsigned),
    Array<arvec> X, Array<arvec> V){
    Array<Array<double> > norm_X;
    double t=0.,h=a/n;
    unsigned t1=clock();
    for(unsigned i=0;i<n;++i){
        X[i+1]=X[i]+h*V[i];
        norm_X=norm_r(X[i+1]);
        for(unsigned k=0;k<V[0].size();++k)
            V[i+1][k]=V[i][k]+h*(*F(t,X[i+1],m,norm_X[k],k));
        t+=h;
    }
    unsigned t2=clock();
    std::cout << "t2-t1= " << (t2-t1)/1000. << std::endl;
}

const char* filenamestr(const char* s1, unsigned j, const char* s2){
    std::stringstream ss;
    ss << s1 << j << s2;
    return (ss.str()).c_str();
}



int main(){
    unsigned t0=clock();
    Array<arvec> X,V;
    std::ofstream out,energy;
    std::ifstream file_size,file_mass,file_initial;
    Array<double> m;
    arvec V0;
    Array<Array<double> > norm_X0;
    unsigned N,M;
    double B;


    file_size.open("file_size.dat",std::ios_base::in);
    file_size >> B >> N >> M;
    file_size.close();
    
    m=Array<double>(M);
    V0=arvec(M,vector<double>(3));

    file_mass.open("file_mass.dat",std::ios_base::in);
    for(unsigned i=0;i<m.size();++i){
        file_mass >> m[i];
    }    
    file_mass.close();
 
    X=Array<arvec>(N+1,arvec(M,vector<double>(3)));
    V=Array<arvec>(N+1,arvec(M,vector<double>(3)));
    
    
    file_initial.open("file_initial.dat",std::ios_base::in);
    for(unsigned i=0;i<M;++i){
        for(int j=0;j<3;++j)
            file_initial >> X[0][i](j);
        for(int j=0;j<3;++j)
            file_initial >> V0[i](j);
    }        
    file_initial.close();
    norm_X0=norm_r(X[0]);
    for(unsigned k=0;k<M;++k)
        V[0][k]=V0[k]+B/N*(*F(0,X[0],m,norm_X0[k],k));
    
    LeapFrog(B,m,N,F,X,V);

    unsigned t3=clock();
    std::cout << "t3-t0= " << (t3-t0)/1000. << std::endl;
    
    for(unsigned j=0;j<M;++j){
        out.open(filenamestr("NBodyLP_",j,".dat"),std::ios_base::trunc);
        for(unsigned i=0;i<=N;++i){
            for(unsigned k=0;k<X[i][j].size();++k)
                out << X[i][j](k) << " ";
            out << std::endl;
        }        
        out.close();
    }
    
    
    energy.open("NBodyLP_energy.dat",std::ios_base::trunc);
    energy <<  "h(" << 0 << ")= " << EnergyIntegral((V0+V[0])/2, norm_r(X[0]),m) << std::endl;
    for(unsigned i=1;i<=N;++i)
        energy <<  "h(" << i << ")= " << EnergyIntegral((V[i-1]+V[i])/2, norm_r(X[i]),m) << std::endl;
    std::cout << fabs(EnergyIntegral((V0+V[0])/2, norm_r(X[0]),m)-
    EnergyIntegral((V[N-1]+V[N])/2, norm_r(X[N]),m)) << std::endl;
    energy.close();
    std::cout << "N=" << N << std::endl;
    return 0;
}