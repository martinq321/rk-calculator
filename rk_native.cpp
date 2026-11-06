// Native C++ evaluation of the same PRECOMPUTED quadrature as the browser.
// This is not the previous FLINT/Arb certified evaluator.
// Compile: g++ -O3 -std=c++17 rk_native.cpp -o rk_native
// Run from this directory: ./rk_native 2 1000000
#include "kernels_config.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

static std::vector<double> load_kernel(const std::string& file, std::size_t n) {
    std::vector<double> a((KMAX+1)*n);
    std::ifstream f(file,std::ios::binary);
    if(!f) throw std::runtime_error("Cannot open "+file+"; run from the package directory.");
    f.read(reinterpret_cast<char*>(a.data()),a.size()*sizeof(double));
    if(!f || f.peek()!=EOF) throw std::runtime_error("Unexpected kernel size: "+file);
    // Kernel data are little-endian IEEE float64 (usual x86/ARM/Wasm hosts).
    return a;
}
static double dot(const std::vector<double>& a,std::size_t n,int k,double x,double lo=0,bool below=false) {
    double sum=0, correction=0, logx=std::log(x);
    for(std::size_t i=0;i<n;++i) {
        if(below ? a[i]>lo : a[i]<=lo) continue;
        double term=std::exp(logx*(a[i]-1.0))*a[std::size_t(k)*n+i];
        double y=term-correction, t=sum+y;
        correction=(t-sum)-y; sum=t;
    }
    return sum*x;
}
static double riemann_r(double x) {
    double y=std::log(x),term=1.0/x,sum=term;
    int max_n=static_cast<int>(y+20*std::sqrt(y+1)+251);
    for(int n=1;n<=max_n;++n) {
        term*=y/n;
        double add=term/(n<=200 ? GRAM_DEN[n-1] : static_cast<double>(n)); sum+=add;
        if(n>y+10 && std::fabs(add)<1e-16*std::fabs(sum)) return sum*x;
    }
    throw std::runtime_error("Gram series stopping rule not reached.");
}
int main(int argc,char**argv) {
    try {
        if(argc!=3) {std::cerr<<"Usage: rk_native k x (k=1..6, 1<x<=1e300)\n";return 1;}
        std::size_t pos=0;
        int k=std::stoi(argv[1],&pos); if(pos!=std::string(argv[1]).size()) throw std::runtime_error("Invalid k");
        double x=std::stod(argv[2],&pos); if(pos!=std::string(argv[2]).size()) throw std::runtime_error("Invalid x");
        if(k<1||k>int(KMAX)||!std::isfinite(x)||x<=1||x>1e300) throw std::runtime_error("Inputs outside demo range");
        std::cout<<std::setprecision(16);
        if(k==1) {std::cout<<"R1="<<riemann_r(x)<<"\nMETHOD=GRAM_FLOAT64\n";return 0;}
        auto a=load_kernel("kernels_fine.bin",NFINE), b=load_kernel("kernels_mesh.bin",NMESH);
        double v=dot(a,NFINE,k,x), small=dot(a,NFINE,k,x,1.0/65.0,true), coarse=dot(b,NMESH,k,x);
        std::cout<<"R"<<k<<"_APPROX="<<v<<"\nQ=128\nLOWER_CUT=0..1/129 (omitted)\n"
                 <<"CUTOFF_CHANGE="<<std::fabs(small)<<"\nMESH_CHANGE="<<std::fabs(v-coarse)
                 <<"\nSTATUS=UNCERTIFIED_FINITE_CUT\n";
    } catch(const std::exception&e) {std::cerr<<e.what()<<"\n";return 2;}
}
