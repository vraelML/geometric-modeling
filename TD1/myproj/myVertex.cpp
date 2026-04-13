#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"
#include <set>

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
	normal->clear();
	if (originof == NULL) return;

	std::set<myFace *> visited_faces;
	myHalfedge *e = originof;

	while (e != NULL)
	{
		if (e->adjacent_face != NULL && visited_faces.insert(e->adjacent_face).second)
			*normal += *(e->adjacent_face->normal);

		if (e->twin == NULL || e->twin->next == originof) break;
		e = e->twin->next;
	}

	e = originof->prev;
	while (e != NULL && e->twin != NULL)
	{
		e = e->twin;
		if (e->adjacent_face != NULL && visited_faces.insert(e->adjacent_face).second)
			*normal += *(e->adjacent_face->normal);

		if (e->prev == NULL || e->prev == originof->prev) break;
		e = e->prev;
	}

	if (normal->length() > 0.0) normal->normalize();
}
