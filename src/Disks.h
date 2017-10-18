#ifndef  _DISKS_H_    /* only process this file once */
#define  _DISKS_H_


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "vec.h"
#include "type.h"
#include "disk.h"
#include "stats.h"
#include "animation.h"
#include "Collision.h"

#define DO_ANNULUS 0
#define GAS 0
class Collision;


class Disks{
	
	public:
	
	//Constants
		//Time interval between "swirl collisions"
		const double swirl_interval = 1;
		//Radius of boundary
		const double boundrad = 20;
		//Angle of swirl collision. For example, pi/6 means
		//the boundary traces a dodecagonal path
		const double swirl_angle = 3.14159265359 / 6;
		//Marble/Marble friction parameter
	const double mu = 10;
		//Marble/Wall friction parameter
	const double wmu = 10;
	
	
	//Set by the main function in RotationSim.cpp depending on whether we're animating
		static Animation* animation;
	
	
	
		Disks(){}
		~Disks(){
			if(disks){
				delete [] disks;
			}
		}
	
		//Sets up disks, reads positions/velocities from file
		void initialize(int N);
		//Update positions based on trajectories and time passed
		void updatePositions(double collisionTime);
		//Used in animation
		void draw();
		//Populate currentCollisions vector with next collisions to take place
		void nextCollisions(std::vector<Collision>& currentCollisions);
	
		//Returns next collision between a,b if there is one
		Collision nextDiskCollision(Disk& a, Disk& b);
		//Returns next collision of disk with wall
		Collision nextWallCollision(Disk& disk);
	
		//Iterates through all disks, adds next disk/disk collision
		void checkDiskCollisions( std::vector<Collision>& currentCollisions );
		void fillCells();
		//Iterates through all disks, adds next disk/wall collision
		void checkWallCollisions( std::vector<Collision>& currentCollisions );
		//Adds a swirl collision
		void checkSwirlCollision( std::vector<Collision>& currentCollisions );

		//This function is used by the above three, it adds "collision" to the vector only if
		//it is the next to occur, otherwise it leaves the vector unchanged.
		void addCollision(std::vector<Collision>& currentCollisions, Collision& collision);

	//Update trajectories for the various collisions. This is where the magic happens
		void processCollision(Collision& collision);
		void processNormalCollision(Collision& collision);
		void processWallCollision(Collision& collision);
		void processInnerWallCollision(Collision& collision);
		void swirl();
	
	
	//Various stats functions.
	void printStats();
	double getAngVel();
	vec centerOfMass();
	vec centerOfMassVel();

	double squareSum();
	double getAngVelVariance(double mean);
	
	//The disks!
	Disk* disks;
	
	private:
		double inner_ang;
		double inner_vel;

		std::vector<int>* cells;
		int nbinx;

		int num_of_disks;
	
		Stats stats;
		double boundpos[2];
		double boundvel[2];
		double swirl_time;
};


#endif
