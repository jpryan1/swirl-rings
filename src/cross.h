#ifndef  _CROSS_H_    /* only process this file once */
#define  _CROSS_H_
#define GLM_FORCE_RADIANS 1
#include <iostream>
#include <vector>
#include <cmath>
#include "glew.h"
#include <GLFW/glfw3.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#define PI 3.14159265359


class Cross{
	
public:
	std::vector<GLfloat> vertices;
	
	std::vector<GLushort> indices;
	

	static GLuint modelLoc;
	Cross(){}
	Cross(double r);
	void draw(GLfloat a, GLfloat b, GLfloat c, GLfloat ang);//, GLfloat c);
private:
	double radius;
};


#endif
