
#ifndef  _ANIMATION_H_    /* only process this file once */
#define  _ANIMATION_H_

#define GLEW_STATIC
#define DELTA_T 1e-8
#include <iostream>
#include "circle.h"
#include <mutex>
#include <atomic>
#include <unistd.h>
#include "disk.h"
#include "cross.h"
#include "vec.h"
#include <assert.h>
#define M_FRAME 0 

class Animation{
	public:
	
		Animation(){  }
		Animation(int n){
			num_of_disks = n;
			disks = new Disk[n];
			notReady = true;
			drawing = false;
			boundpos[0] = 0;
			boundpos[1] = 0;
			total_time = 0;
			disks_buffer = new Disk[n];
			oa = 0;
		}
		void initialize();
		void setup();
		void compileShaders();
		void generateBuffers();
		void generateShapes();
		void setProjectionMatrices();
		void draw();
		void drawShapes();
		void setDisks(Disk* d, double* b, double* v, double inner_ang, double inner_vel);
		void moveDisks(double time);
		std::atomic<bool> notReady;
	std::atomic<bool> drawing;
	
	private:
	double total_time;
		GLuint s_VBO, s_VAO, s_EBO, shaderProgram, modelLoc, colorLoc, viewLoc;

	GLuint x_VBO, x_VAO, x_EBO;
		int width, height;
		GLFWwindow* window;
		Circle circle;
		Circle bound;
	Cross cross;
	Circle m_ball;
		Disk* disks;
	Disk* disks_buffer;
	double boundpos[2];
	double boundvel[2];
		int num_of_disks;
	double ia, iv, oa;
	
};
#endif
