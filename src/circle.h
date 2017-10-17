#ifndef  _CIRCLE_H_    /* only process this file once */
#define  _CIRCLE_H_
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


class Circle{
	
public:
	std::vector<GLfloat> vertices;
	
	std::vector<GLushort> indices;
	

	static GLuint modelLoc;
	Circle(){}
	Circle(double r);
	void draw(GLfloat a, GLfloat b, GLfloat c);//, GLfloat c);
private:
	double radius;
};


#endif
