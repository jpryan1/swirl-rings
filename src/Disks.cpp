#include "Disks.h"
#include <assert.h>
#include <omp.h>
Animation* Disks::animation;
//TODO currentCollisions probably doesn't need to be vector!


void Disks::initialize(int N){
	this->num_of_disks = N;
	disks = new Disk[N];

	
	cells = new std::vector<int>[30*30];
	nbinx = 30;
	//Initialize boundary variables
	this->boundpos[0] = 0;
	this->boundpos[1] = 0;
	this->boundvel[0] = 2;//0.5;
	this->boundvel[1] = 0;
	this->swirl_time=0;
	
	//Get disks coordinates from input file
	std::ifstream input("input.txt", std::ifstream::in);
	for(int i=0; i<N; i++){
		input >> disks[i].pos[0] >> disks[i].pos[1] >> disks[i].vel[0] >> disks[i].vel[1];
		if(GAS) disks[i].vel[1]+=1;
		disks[i].ID = i;
		disks[i].ang = 0;
		disks[i].ang_vel = 0;
		
		
	}
	input.close();
	
	
	if(animation){
		animation->setDisks(disks, boundpos, boundvel);
		animation->notReady = false;
	}
	
	inner_ang=0;
	inner_vel=0;
	//stats.initialize(N);
	//stats.update(disks, boundpos,boundvel, -1);
}

double Disks::getAngVel(){
	
	
	double len;
	
	double angvel = 0;
	vec CoM = centerOfMass();
	vec CoM_vel = centerOfMassVel();
	for(int i=0; i<num_of_disks; i++) {
		double r1 = disks[i].pos[0] - CoM.a[0];
		double r2 = disks[i].pos[1] - CoM.a[1];
		double v1 = disks[i].vel[0] - CoM_vel.a[0];
		double v2 = disks[i].vel[1] - CoM_vel.a[1];
		double norm = r1*r1 + r2*r2;
		angvel +=  (r1*v2 - r2*v1)/norm;
	}
	return angvel / num_of_disks;
	
}

vec Disks::centerOfMass(){
	double a = 0;
	double b = 0;
	for(int i=0; i<num_of_disks; i++){
		a += disks[i].pos[0];
		b += disks[i].pos[1];
	}
	a = a / num_of_disks;
	b = b / num_of_disks;
	return vec(a,b);
}


vec Disks::centerOfMassVel(){
	double a = 0;
	double b = 0;
	for(int i=0; i<num_of_disks; i++){
		a += disks[i].vel[0];
		b += disks[i].vel[1];
	}
	a = a / num_of_disks;
	b = b / num_of_disks;
	return vec(a,b);
}void Disks::printStats(){
	
	//stats.printHeatMap();
	//Insert your favorite print function here.
}





