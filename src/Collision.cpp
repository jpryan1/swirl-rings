#include "Collision.h"


Collision::Collision(const Collision &other){
	type = other.type;
	time = other.time;
	disks[0] = other.disks[0];
	disks[1] = other.disks[1];
}


double Collision::getTime(){
	return this->time;
}


Type Collision::getType(){
	return this->type;
}
