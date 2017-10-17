#ifndef  _VEC_H_    /* only process this file once */
#define  _VEC_H_
#include <cmath>
struct vec{
	double a[2];
	vec(double* arr){
		for(int i=0; i<2; i++) a[i] = arr[i];
	}
	vec(){
		for(int i=0; i<2; i++) a[i] = 0.0;
	}
	vec(double m, double n){
		a[0] = m;
		a[1] = n;
	}
	double norm(){
		return sqrt(a[0]*a[0] + a[1]*a[1]);
	}
	double dot(const vec& o){
		return a[0]*o.a[0] + a[1]*o.a[1];
	}
	vec add(const vec& o){
		double b[2];
		for(int i=0; i<2; i++) b[i] = a[i] + o.a[i];
		return vec(b);
	}
	vec minus(const vec& o){
		double b[2];
		for(int i=0; i<2; i++) b[i] = a[i] - o.a[i];
		return vec(b);
	}
	vec times(double c){
		double b[2];
		for(int i=0; i<2; i++) b[i] = c*a[i];
		return vec(b);
	}
	void print(){
		std::cout<<a[0]<<" "<<a[1]<<std::endl;
	}
	
};
#endif
