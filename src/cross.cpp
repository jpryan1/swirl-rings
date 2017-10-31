#include "cross.h"

GLuint Cross::modelLoc;

Cross::Cross(double radius){
	
	for(int i=-1; i<=1; i+=2){
		for(int j=-1; j<=1; j+=2){
			vertices.push_back(i* (radius));
			vertices.push_back(j* (radius/5));
			vertices.push_back(0);
		}
	}
	
	for(int i=-1; i<=1; i+=2){
		for(int j=-1; j<=1; j+=2){
			vertices.push_back(j* (radius/5));
			vertices.push_back(i* (radius));
			vertices.push_back(0);
		}
	}
	
	for(int i=0; i<2; i++){
		indices.push_back(0+i);
		indices.push_back(1+i);
		indices.push_back(2+i);
	}
	for(int i=0; i<2; i++){
		indices.push_back(4+i);
		indices.push_back(5+i);
		indices.push_back(6+i);
	}
	
	glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(GLfloat), &vertices[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*sizeof(GLushort), &indices[0], GL_STATIC_DRAW);

}


void Cross::draw(GLfloat a, GLfloat b, GLfloat c, GLfloat ang, GLfloat scale){
	
	glm::mat4 transform;
	transform = glm::translate(transform, glm::vec3(a,b,c));
	float angle =(180.0/PI)*ang;
	transform = glm::rotate(transform,ang, glm::vec3(0,0,1));
	transform = glm::scale(transform, glm::vec3(scale, scale, scale));
	
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(transform));
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_SHORT, (GLvoid*) 0);
}



//








