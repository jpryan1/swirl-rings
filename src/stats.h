
#ifndef  _STATS_H_    /* only process this file once */
#define  _STATS_H_
#include "disk.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string.h>

#define HEATMAP_BINS 100
#define QUIVER_BINS 20
#define DENSITY_BINS 100
#define COLORWHEEL_BINS 100

class Stats{
	public:
		//Initialized the class, fills the arrays with zeroes
		void initialize(int N);
	
		//Update ALL stats based on disk pos/vels and boundary pos/vels and time passed
		void update(Disk* d, double* b,double* bv, double time);
		void updateColorwheel(double dif, int disk_ID);
	void updateContributions(double dif, int which);
	//Print stats
		void printRad();
		void printCoM();
		void printHeatMap();
		void printQuiver();
		void printDensity();
	void printColorwheel();
	void printContributions();

	//Keep track of OWN PERSONAL COPY of disks
		Disk* disks;
	
	//These are bins for the heatmap plot. The first records number of hits in a sector,
	//the second records ang_vels resulting from hits in a sector
		long ang_counts[HEATMAP_BINS];
		double ang_vels[HEATMAP_BINS];
	
	//Destructor
		~Stats(){
			if(disks) delete [] disks;
		}
	
	private:
		int num_of_disks;
		double boundpos[2];
		double boundvel[2];
	double contributions[3];
	//Colorwheel plot
	double colorwheel_forces[COLORWHEEL_BINS];
	int colorwheel_counters[COLORWHEEL_BINS];
	
	
	
	//Center of Mass
		double CoM[2];
	//Angle of M (from M-frame)
		double M_ang;
	
	//Density histogram bins
		int hist[DENSITY_BINS][DENSITY_BINS];
	
	//For quiver plots - first records number of velocities recorded,
	//second records sum (each component)
		int quiver_counts[QUIVER_BINS][QUIVER_BINS];
		double quiver_hist[QUIVER_BINS][QUIVER_BINS][2];
};






#endif
