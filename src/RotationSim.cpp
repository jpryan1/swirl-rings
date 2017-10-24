
#include "animation.h"
#include <thread>
#include <iostream>
#include <vector>
#include "Disks.h"
#include "Collision.h"
#include <ctime>

#define NUM_OF_ITERATIONS 20//1000000.0

void simulation(int N);


int main(int argc, char** argv){
	
	//By default, we use 55 disks. The max that we can use with our input file is 56
	int NUM_OF_DISKS = 55;
	std::cout.precision(16);
	
	
	if(argc<2){
		std::cout<<"Include arg - a for animate or r for run (without animation)"<<std::endl;
		return 0;
	}
	
	
	/*
	 
	 We accept one of three options. 'r' means just run the simulation on NUM_OF_DISKS disks,
	 and output the average angular velocity.
	 'a' is the same, but with an animation included
	 't' means 'test', and runs the simulation for N from 5 to 56, outputting stats results
	
	 */
	
	//We allow command line arg input of NUM_OF_DISKS
	if(argc > 2) NUM_OF_DISKS = std::stoi(argv[2], NULL);
	
	
	
	if(*(argv[1])=='r'){
		
		simulation(NUM_OF_DISKS);
		return 0;
		
	}
	
	else if(*(argv[1])=='a'){
		
		Animation animation(NUM_OF_DISKS);
		animation.setup();
		Disks::animation = &animation;

		//Set the simulation running
		std::thread drawer(simulation, NUM_OF_DISKS);
		//Begin the animation
		animation.draw();
		//Wait for the simulation to finish
		drawer.join();
		return 0;
		
	}
	
	else if(*(argv[1])=='t'){
		
		for(NUM_OF_DISKS=5; NUM_OF_DISKS<= 56; NUM_OF_DISKS+=1){
			simulation(NUM_OF_DISKS);
		}
		
	}
	
	else{
		
		std::cout<<"Option must be a for animate or r for run (without animation)"<<std::endl;
		return 0;
		
	}

	
}




void simulation(int NUM_OF_DISKS){
	std::vector<Collision> currentCollisions;
	Disks disks;
	disks.initialize(NUM_OF_DISKS);

	double total_ang_vel = 0;
	double total_time=0;
	for(int iterations = 0; iterations<NUM_OF_ITERATIONS; iterations++){
		disks.nextCollisions(currentCollisions);

		//fill the currentCollisions vector with the next collisions.
		disks.updatePositions(currentCollisions[0].getTime());
		total_time +=currentCollisions[0].getTime();
		//move everyone to their next position, record time
	
		
		for(int i=0; i<currentCollisions.size(); i++){
			disks.processCollision(currentCollisions[i]);
			
		}



		//update relevant trajectories - bonk!
		
		
		double meanAV = disks.getAngVel();
		total_ang_vel += meanAV;
		//record angvel, add to average
	}
	
//	disks.printStats();
	std::cout<<(total_ang_vel/NUM_OF_ITERATIONS)<<std::endl;

}










