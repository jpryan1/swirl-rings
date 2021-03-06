#include "Disks.h"
#include <assert.h>
#include <omp.h>
//TODO currentCollisions probably doesn't need to be vector!



void Disks::updatePositions(double time){


	if(animation){
		animation->setDisks(disks, boundpos, boundvel, inner_ang, inner_vel);
		animation->moveDisks(time);
	}


	boundpos[0] += time*boundvel[0];
	boundpos[1] += time*boundvel[1];
	for(int i=0; i<num_of_disks; i++){
		for(int j=0; j<2; j++){
			disks[i].pos[j] += time*disks[i].vel[j];
			
		}
		disks[i].ang += time*disks[i].ang_vel;
		
	}
	inner_ang += time*inner_vel;
	//stats.update(disks, boundpos,boundvel, time);
}

/*
 The following 4 functions are the processing of events in the simulation, 
 one function per type of event, and one function to rule them all.
 */

void Disks::processCollision(Collision& collision){
	//Obligation here is just to change the trajectories of the affected disks
	double angvel_before, angvel_after;
	int which;
	
	angvel_before = getAngVel();
	
	switch(collision.getType()){
		case NORMAL:
			
			which = 0;
			processNormalCollision(collision);
	
			break;
			
		case WALL:
			which = 1;
			processWallCollision(collision);
			
			break;
			
		case SWIRL:
			which = 2;
			swirl();
			
			break;
	}
	
	angvel_after = getAngVel();
	
	//stats.updateContributions(angvel_after-angvel_before, which);
	
	
}


void Disks::processNormalCollision(Collision& collision){
	Disk s1 = *(collision.disks[0]);
	Disk s2 = *(collision.disks[1]);
	vec v1 = vec(collision.disks[0]->vel);
	vec v2 = vec(collision.disks[1]->vel);
	vec x1 = vec(collision.disks[0]->pos);
	vec x2 = vec(collision.disks[1]->pos);

	vec unit = x2.minus(x1).times(0.5);
	double angv1 = s1.ang_vel;
	double angv2 = s2.ang_vel;
	
	//Process as if elastic
	vec old_normal = unit.times(v1.dot(unit));
	double dot1 = v1.minus(v2).dot(x1.minus(x2));
	dot1 /= pow(x1.minus(x2).norm(), 2);
	vec s1vel = v1.minus(x1.minus(x2).times(dot1));
	double dot2 = v2.minus(v1).dot(x2.minus(x1));
	dot2 /= pow( x2.minus(x1).norm(), 2);
	vec s2vel = v2.minus(x2.minus(x1).times(dot2));
	
	
	//Part 3a - Break vectors into components
	vec newPerp1 = unit.times(s1vel.dot(unit));
	vec par1 = s1vel.minus(newPerp1);
	if(fabs(par1.norm())<1e-13){
		par1 = vec();
	}
	unit = unit.times(-1); //points towards x1
	vec newPerp2 = unit.times(s2vel.dot(unit));
	vec par2 = s2vel.minus(newPerp2);
	if(fabs(par2.norm())<1e-13){
		par2 = vec();
	}
	
	vec delv1perp = newPerp1.minus(old_normal);
	double DELTA = delv1perp.norm();
	//So named because it is the same for both disks
	
	//Part 3b - r0 is x2vel minus x1vel
	vec r0unit(unit.a[1], -unit.a[0]);
	
	vec r0 = par1.minus(par2).add(r0unit.times(angv1+angv2));
	int s0 = (r0.dot(r0unit) > 0 ? 1 : -1);
	
	//Part 3c - apply sticking rules
	double t = (r0.norm())/(mu*DELTA*2*((1.0/DISK_MASS)+(1.0/DISK_MOMENT)));
	if(mu==0){
		t=0;
	}
	t = fmin(t, 1);
	double impulse = mu*s0*DELTA*t;
	
	vec newPar1 = par1.minus(r0unit.times(impulse/DISK_MASS));
	vec newPar2 = par2.add(r0unit.times(impulse/DISK_MASS));
	
	memcpy(collision.disks[0]->vel, newPar1.add(newPerp1).a, 2*sizeof(double));
	memcpy(collision.disks[1]->vel, newPar2.add(newPerp2).a, 2*sizeof(double));
	angv1 -= impulse/DISK_MOMENT;
	angv2 -= impulse/DISK_MOMENT;
	collision.disks[0]->ang_vel = angv1;
	collision.disks[1]->ang_vel = angv2;
	
}



void Disks::processWallCollision(Collision& collision){



	Disk s = *(collision.disks[0]);


	vec rad(s.pos[0]-boundpos[0],s.pos[1]-boundpos[1]);
	double radNorm = rad.norm();
	double wall_moment;
	double wall_vel;
	bool inner;
	if(radNorm>boundrad-1.1){
		inner = false;
		wall_vel = OUTER_WALL_VEL;
		wall_moment = OUTER_WALL_MOMENT;
	}else{
		inner=true;
		wall_vel = inner_vel;
		wall_moment = INNER_WALL_MOMENT;
	}

	rad = rad.times(1/radNorm); //rad is normal.
	
	double angv = s.ang_vel;
	vec bv(boundvel);
	vec newVel(s.vel);
	newVel = newVel.minus(bv);  //now we're in stationary boundary frame of reference
	
	vec perp = rad.times(newVel.dot(rad));
	vec par = newVel.minus(perp);
	perp = perp.times(-1); //bonk

	double DELTA = 2*perp.norm();


	if(inner){
		
		vec r0unit=vec(rad.a[1], -rad.a[0]);
		
		vec r0 = par.add(r0unit.times(angv+wall_vel*20));
		int s0 = (r0.dot(r0unit) > 0 ? 1 : -1);
		
		double t = (r0.norm())/(wmu*DELTA*((1.0/DISK_MASS)+(1.0/DISK_MOMENT)+(20.0/wall_moment)));
		if(wmu==0){
			t=0;
		}
		t=fmin(t,1);
		double impulse = wmu*s0*DELTA*t;
		
		
		vec newPar = par.minus(r0unit.times(impulse/DISK_MASS));
		angv -= impulse/DISK_MOMENT;
		inner_vel -= impulse/INNER_WALL_MOMENT;
		
		vec result = perp.add(newPar);
		result = result.add(bv);
		
		for(int i=0; i<2; i++) collision.disks[0]->vel[i] = result.a[i];
		collision.disks[0]->ang_vel = angv;

		

	}
	else{
		vec r0unit=vec(-rad.a[1], rad.a[0]);
	
		
		vec r0 = par.add(r0unit.times(angv-wall_vel*30));
		int s0 = (r0.dot(r0unit) > 0 ? 1 : -1);
		
		double t = (r0.norm())/(wmu*DELTA*((1.0/DISK_MASS)+(1.0/DISK_MOMENT)+(30.0/wall_moment)));
		if(wmu==0){
			t=0;
		}
		t=fmin(t,1);
		double impulse = wmu*s0*DELTA*t;
		
		
		vec newPar = par.minus(r0unit.times(impulse/DISK_MASS));
		angv -= impulse/DISK_MOMENT;
		vec result = perp.add(newPar);
		result = result.add(bv);
		
		for(int i=0; i<2; i++) collision.disks[0]->vel[i] = result.a[i];
		collision.disks[0]->ang_vel = angv;

		
		
	}
	
}


void Disks::swirl(){
	double v0 = boundvel[0];
	boundvel[0] = cos(swirl_angle) * boundvel[0] - sin(swirl_angle) * boundvel[1];
	boundvel[1] = sin(swirl_angle) * v0 + cos(swirl_angle) * boundvel[1];
}


