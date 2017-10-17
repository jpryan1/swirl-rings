
#ifndef  _ANIMATION_H_    /* only process this file once */
#define  _ANIMATION_H_

#define GLEW_STATIC
#define DELTA_T 0.000003
#include <iostream>
#include "circle.h"
#include <mutex>
#include <atomic>
#include <unistd.h>
#include "disk.h"
#include "cross.h"
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
		}
		void initialize();
		void setup();
		void compileShaders();
		void generateBuffers();
		void generateShapes();
		void setProjectionMatrices();
		void draw();
		void drawShapes();
		void setDisks(Disk* d, double* b, double* v);
		void moveDisks(double time);
		std::atomic<bool> notReady;
	std::atomic<bool> drawing;
	
	private:
	double total_time;
		GLuint s_VBO, s_VAO, s_EBO, shaderProgram, modelLoc, colorLoc, viewLoc;
		GLuint b_VBO, b_VAO, b_EBO;
	GLuint m_VBO, m_VAO, m_EBO;
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
	
	
};
#endif
