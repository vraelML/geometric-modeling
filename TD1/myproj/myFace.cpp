#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
	if (adjacent_halfedge == NULL ||
		adjacent_halfedge->next == NULL ||
		adjacent_halfedge->next->next == NULL)
		return;

	myPoint3D *p1 = adjacent_halfedge->source->point;
	myPoint3D *p2 = adjacent_halfedge->next->source->point;
	myPoint3D *p3 = adjacent_halfedge->next->next->source->point;

	normal->setNormal(p1, p2, p3);
}
