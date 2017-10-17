#ifndef  _COLLISION_H_    /* only process this file once */
#define  _COLLISION_H_


#include <stdlib.h>
#include <iostream>     
#include <fstream>
#include <cmath>
#include "type.h"
#include "disk.h"

class Collision{
	
	public:
		Collision(){}
	//For a "swirl collision", when the boundary's trajectory changes
		Collision(double t) : time(t), type(SWIRL) {}
	
	//For a NORMAL collisions, between two marbles
		Collision(double t, Disk& a, Disk& b){
			time = t;
			type = NORMAL;
			disks[0] = &a;
			disks[1] = &b;
		}
	
	//For a collision of one marble with the wall
		Collision(double t, Disk& a){
			time = t;
			type = WALL;
			disks[0] = &a;
			disks[1] = NULL;
		}
	
	//Copy constructor
		Collision(const Collision &other);
	
	
		double getTime();
	
		Type getType();
	
	//returns value between 0 and 2pi, corresponding to loc
	double getWallHitPosition();
	
	
		Disk* disks[2];
	
	
	private:
		Type type;
		double time;
		
	
};


#endif
