#include "Disks.h"


Animation* Disks::animation;
//TODO currentCollisions probably doesn't need to be vector!


void Disks::initialize(int N){
	this->num_of_disks = N;
	disks = new Disk[N];
	//Initialize boundary variables
	this->boundpos[0] = 0;
	this->boundpos[1] = 0;
	this->boundvel[0] = 0.5;
	this->boundvel[1] = 0;
	this->swirl_time=0;
	
	//Get disks coordinates from input file
	std::ifstream input("input.txt", std::ifstream::in);
	for(int i=0; i<N; i++){
		input >> disks[i].pos[0] >> disks[i].pos[1] >> disks[i].vel[0] >> disks[i].vel[1];
		disks[i].ID = i;
		disks[i].ang = 0;
		disks[i].ang_vel = 0;
		
		
	}
	input.close();
	
	
	if(animation){
		animation->setDisks(disks, boundpos, boundvel);
		animation->notReady = false;
	}
	
	
	//stats.initialize(N);
	//stats.update(disks, boundpos,boundvel, -1);
}


void Disks::updatePositions(double time){

	boundpos[0] += time*boundvel[0];
	boundpos[1] += time*boundvel[1];
	for(int i=0; i<num_of_disks; i++){
		for(int j=0; j<2; j++){
			disks[i].pos[j] += time*disks[i].vel[j];
			
		}
		disks[i].ang += time*disks[i].ang_vel;
		
	}
	
	if(animation){
		animation->moveDisks(time);
		animation->setDisks(disks, boundpos, boundvel);
		
	}
	//stats.update(disks, boundpos,boundvel, time);
}

/*
 The following two functions calculate the next time of collision between either a pair of marbles or a marble and the wall
 */

Collision Disks::nextDiskCollision(Disk& a, Disk& b){
	double time;
	
	double dv0 = a.vel[0]-b.vel[0];
	double dv1 = a.vel[1]-b.vel[1];
	double dp0 = a.pos[0]-b.pos[0];
	double dp1 = a.pos[1]-b.pos[1];
	
	
	double A = pow(dv0,2) + pow(dv1,2);
	double B = 2 * ( dp0* dv0 + dp1 *dv1);
	double C = pow(dp0,2) + pow(dp1,2) -4;
	//ball radius is 1, dist between two centers is 2, 2^2=4
	
	double det = B*B-4*A*C;
	if(det<0){
		time = -1;
	}
	else{
		double t = (0.0 - B - sqrt(det))/(2*A);
		if(t>1e-13){
			time = t;
			
		}else{
			time = -1;
		}
	}
	return Collision( time, a, b);
	
	
}

Collision Disks::nextWallCollision(Disk& a){
	
	double time;
	double dv0 = a.vel[0]-boundvel[0];
	double dv1 = a.vel[1]-boundvel[1];
	double dp0 = a.pos[0]-boundpos[0];
	double dp1 = a.pos[1]-boundpos[1];

	
	double A = pow(dv0,2) + pow(dv1,2);
	double B = 2 * ( dp0* dv0 + dp1 *dv1);
	double C = pow(dp0,2) + pow(dp1,2) - pow(boundrad-1, 2);
	double det = B*B-4*A*C;
	if(det<0){
		time = -1;
	}else{
		double t = (0.0 - B + sqrt(det))/(2*A);
		if(t>1e-13){
			time = t;
		}
		else{
			time = -1;
		}
	}
	
//	double firsttime = time;
//	
//	
//	C = pow(dp0,2) + pow(dp1,2) - pow(20+1, 2);
//	det = B*B-4*A*C;
//	
//	if(det<0){
//		return Collision(firsttime, a);
//	}else{
//		double t = (0.0 - B + sqrt(det))/(2*A);
//		if(t>1e-13){
//			time = t;
//		}
//		else{
//			return Collision(firsttime, a);
//		}
//	}
//
//	if(firsttime==-1) return Collision(time, a);
//	if(firsttime>time) std::cout<<"inner in "<<time<<std::endl;
	return Collision( time, a);
	
}


/*The following four functions look for the next collision to occur by comparing the times
 of all collisions that will happen in the future
 */
void Disks::nextCollisions(std::vector<Collision>& currentCollisions){
	currentCollisions.clear();
	checkDiskCollisions( currentCollisions );
	checkWallCollisions( currentCollisions );
	if(boundvel[0] || boundvel[1] )  checkSwirlCollision( currentCollisions );
	
}

void Disks::checkDiskCollisions( std::vector<Collision>& currentCollisions ){
	Collision collision;
	
	for(int i=0; i<num_of_disks-1; i++){
		for(int j=i+1; j<num_of_disks; j++){
			collision = nextDiskCollision(disks[i], disks[j]);
			if(collision.getTime() == -1) continue;
			addCollision(currentCollisions, collision);
		}
	}
}

void Disks::checkWallCollisions( std::vector<Collision>& currentCollisions ){
	
	Collision collision;
	for(int i=0; i<num_of_disks; i++){
		collision = nextWallCollision(disks[i]);
		if(collision.getTime() ==-1) continue;
		addCollision(currentCollisions, collision);
	}
}

void Disks::checkSwirlCollision( std::vector<Collision>& currentCollisions ){
	
	Collision collision(swirl_interval - swirl_time);
	if(swirl_interval - swirl_time - currentCollisions[0].getTime() < 1e-13){
		//this will evaluate to true only if the swirl will be added
		swirl_time = 0;
	}else{
		swirl_time += currentCollisions[0].getTime();
	}
	addCollision( currentCollisions, collision);
	
	
}



/*
 This function deals with checking whether an event is the next to occur.
 It checks the event's time against the times of other events that are to take place shortly, 
 and adds the event if it is at least as soon as the soonest next event.
 Recall that we are using a vector of events because some events may happen at exactly the same time.
 */
void Disks::addCollision(std::vector<Collision>& currentCollisions, Collision& collision){
	
	
	//if currentcollisions is empty, just add
	//Time of vector is found by first collision. If we're within 1e13, we add.
	//If we're 1e13 faster than first one, we delete all that aren't as fast
	
	
	if(currentCollisions.size()==0){
		currentCollisions.push_back(collision);
		return;
	}
	
	double diff = currentCollisions[0].getTime() - collision.getTime();
	
	if( diff > 1e-13){
		currentCollisions.clear();
		currentCollisions.push_back(collision);
	}
	else if(fabs(diff) < 1e-13){
		currentCollisions.push_back(collision);
	}
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
	double r0 = par1.norm() - par2.norm() + angv1 + angv2;
	int s0 = (r0>0) ? 1 : -1;
	double t = (r0)/(4*mu*s0*DELTA);
	if(mu==0){
		t=0;
	}
	
	
	//Part 3c - apply sticking rules
	
	t = fmin(t, 1);
	double impulse = mu*s0*DELTA*t;
	vec unitPar1, unitPar2, newPar1, newPar2;
	
	if(par1.norm()!=0||par2.norm()!=0){
		if(par1.norm()==0){
			unitPar2 = par2.times(1.0/par2.norm());
			unitPar1 = vec(unitPar2.a);
		}
		else if(par2.norm()==0){
			unitPar1 = par1.times(1.0/par1.norm());
			unitPar2 = vec(unitPar1.a);
		}
		else{
			unitPar1 = par1.times(1.0/par1.norm());
			unitPar2 = par2.times(1.0/par2.norm());
		}
		newPar1 = unitPar1.times(par1.norm()-impulse);
		newPar2 = unitPar2.times(par2.norm()+impulse);
		memcpy(collision.disks[0]->vel, newPar1.add(newPerp1).a, 2*sizeof(double));
		memcpy(collision.disks[1]->vel, newPar2.add(newPerp2).a, 2*sizeof(double));
	}else{
		memcpy(collision.disks[0]->vel, newPerp1.a, 2*sizeof(double));
		memcpy(collision.disks[1]->vel, newPerp2.a, 2*sizeof(double));
	}
	angv1 -= impulse;
	angv2 -= impulse;
	collision.disks[0]->ang_vel = angv1;
	collision.disks[1]->ang_vel = angv2;
	
}
void Disks::processWallCollision(Collision& collision){

	Disk s = *(collision.disks[0]);
	double angv = s.ang_vel;
	vec bv(boundvel);
	vec rad(s.pos[0]-boundpos[0],s.pos[1]-boundpos[1]);
	double radNorm = rad.norm();
	rad = rad.times(1/radNorm); //rad is normal.
	
	
	vec newVel(s.vel);
	newVel = newVel.minus(bv);  //now we're in stationary boundary frame of reference
	
	vec perp = rad.times(newVel.dot(rad));
	vec par = newVel.minus(perp);
	perp = perp.times(-1); //bonk
	
	
	double damp = 2*wmu*perp.norm();
	double parnorm = par.norm();
	vec normPar;
	if(fabs(parnorm)<1e-13){
		normPar=vec();
	}else{
		normPar = par.times(1.0/parnorm);
	}
	double t = (parnorm+angv)/(2*damp);
	if(wmu==0){
		t=0;
	}
	t = fmin(t, 1);
	vec newPar = par.minus(normPar.times(t*damp));
	angv -= t*damp;
	perp = perp.add(newPar);
	perp = perp.add(bv);
	for(int i=0; i<2; i++) collision.disks[0]->vel[i] = perp.a[i];
	collision.disks[0]->ang_vel = angv;
	//std::cout<<collision.disks[0]->ang_vel<<std::endl;
}


void Disks::swirl(){
	double v0 = boundvel[0];
	boundvel[0] = cos(swirl_angle) * boundvel[0] - sin(swirl_angle) * boundvel[1];
	boundvel[1] = sin(swirl_angle) * v0 + cos(swirl_angle) * boundvel[1];
}





/*
 The below functions are solely used for collecting data and generating plots
 
 */


double Disks::squareSum(){
	
	
	//This method isn't used anymore, angvel is calculated differently now
	double sum = 0;
	
	for(int i=0; i<num_of_disks; i++){
		double a = disks[i].pos[0] - boundpos[0];
		double b = disks[i].pos[1] - boundpos[1];
		
		sum += a*a + b*b;
	}
	return sum;
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
}


double Disks::getAngVelVariance(double mean){
	return 0;
	// This needs rewriting based on new angvel formula
//	vec c = centerOfMass();
//	double sum = squareSum();
//	double angvelvar = 0;
//	for(int i=0; i<num_of_disks; i++) {
//		double a = disks[i].pos[0] - c.a[0];
//		double b = disks[i].pos[1] - c.a[1];
//		double dist = pow(a,2) + pow(b,2);
//		double cross = a * disks[i].vel[1] - b * disks[i].vel[0];
//		angvelvar += pow(cross,2)/dist;
//	}
//	return sqrt(  (angvelvar / sum) - pow(mean,2) );
}

void Disks::printStats(){
	
	//stats.printHeatMap();
	//Insert your favorite print function here.
}





