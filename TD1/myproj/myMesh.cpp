#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>
#include <GL/glew.h>
#include "myvector3d.h"

using namespace std;

namespace {
	vector<myVertex *> getFaceVertices(myFace *face)
	{
		vector<myVertex *> verts;
		if (face == NULL || face->adjacent_halfedge == NULL) return verts;

		myHalfedge *start = face->adjacent_halfedge;
		myHalfedge *e = start;
		do {
			if (e == NULL || e->source == NULL) {
				verts.clear();
				return verts;
			}
			verts.push_back(e->source);
			e = e->next;
		} while (e != NULL && e != start);

		if (e != start) verts.clear();
		return verts;
	}

	bool isPointInTriangle(double px, double py, double ax, double ay, double bx, double by, double cx, double cy)
	{
		double c1 = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
		double c2 = (cx - bx) * (py - by) - (cy - by) * (px - bx);
		double c3 = (ax - cx) * (py - cy) - (ay - cy) * (px - cx);
		return (c1 >= -1e-9 && c2 >= -1e-9 && c3 >= -1e-9) || (c1 <= 1e-9 && c2 <= 1e-9 && c3 <= 1e-9);
	}

	vector<vector<myVertex*>> earClipPolygon(const vector<myVertex*>& verts, myVector3D N)
	{
		vector<vector<myVertex*>> tris;
		int n = verts.size();
		if (n < 3) return tris;
		if (n == 3) {
			tris.push_back(verts);
			return tris;
		}

		myVector3D T;
		if (abs(N.dX) < abs(N.dY) && abs(N.dX) < abs(N.dZ)) {
			T = myVector3D(1.0, 0.0, 0.0);
		} else if (abs(N.dY) < abs(N.dZ)) {
			T = myVector3D(0.0, 1.0, 0.0);
		} else {
			T = myVector3D(0.0, 0.0, 1.0);
		}
		myVector3D U = T.crossproduct(N);
		U.normalize();
		myVector3D V = N.crossproduct(U);
		V.normalize();

		struct Point2D {
			double x, y;
			myVertex *orig;
		};
		vector<Point2D> pts(n);
		for (int i = 0; i < n; i++) {
			pts[i].x = verts[i]->point->X * U.dX + verts[i]->point->Y * U.dY + verts[i]->point->Z * U.dZ;
			pts[i].y = verts[i]->point->X * V.dX + verts[i]->point->Y * V.dY + verts[i]->point->Z * V.dZ;
			pts[i].orig = verts[i];
		}

		while (pts.size() > 3) {
			bool ear_found = false;
			int m = pts.size();
			for (int i = 0; i < m; i++) {
				int prev = (i - 1 + m) % m;
				int next = (i + 1) % m;

				double cross = (pts[i].x - pts[prev].x) * (pts[next].y - pts[i].y) - (pts[i].y - pts[prev].y) * (pts[next].x - pts[i].x);
				if (cross <= 1e-9) continue;

				bool has_point_inside = false;
				for (int j = 0; j < m; j++) {
					if (j == prev || j == i || j == next) continue;
					if (isPointInTriangle(pts[j].x, pts[j].y, pts[prev].x, pts[prev].y, pts[i].x, pts[i].y, pts[next].x, pts[next].y)) {
						has_point_inside = true;
						break;
					}
				}

				if (!has_point_inside) {
					vector<myVertex*> tri(3);
					tri[0] = pts[prev].orig;
					tri[1] = pts[i].orig;
					tri[2] = pts[next].orig;
					tris.push_back(tri);
					pts.erase(pts.begin() + i);
					ear_found = true;
					break;
				}
			}
			if (!ear_found) break;
		}

		if (pts.size() == 3) {
			vector<myVertex*> tri(3);
			tri[0] = pts[0].orig;
			tri[1] = pts[1].orig;
			tri[2] = pts[2].orig;
			tris.push_back(tri);
		} else {
			for (unsigned int j = 1; j + 1 < pts.size(); j++) {
				vector<myVertex*> tri(3);
				tri[0] = pts[0].orig;
				tri[1] = pts[j].orig;
				tri[2] = pts[j + 1].orig;
				tris.push_back(tri);
			}
		}
		return tris;
	}

	void rebuildMeshFromPolygons(myMesh *mesh, const vector<vector<myVertex *>> &polygons)
	{
		vector<myHalfedge *> old_halfedges;
		old_halfedges.swap(mesh->halfedges);
		vector<myFace *> old_faces;
		old_faces.swap(mesh->faces);

		for (unsigned int i = 0; i < old_halfedges.size(); i++) if (old_halfedges[i]) delete old_halfedges[i];
		for (unsigned int i = 0; i < old_faces.size(); i++) if (old_faces[i]) delete old_faces[i];

		map<myVertex *, int> vertex_to_id;
		for (unsigned int i = 0; i < mesh->vertices.size(); i++) {
			mesh->vertices[i]->originof = NULL;
			vertex_to_id[mesh->vertices[i]] = (int)i;
		}

		map<pair<int, int>, myHalfedge *> twin_map;
		for (unsigned int i = 0; i < polygons.size(); i++)
		{
			const vector<myVertex *> &poly = polygons[i];
			if (poly.size() < 3) continue;

			myFace *f = new myFace();
			mesh->faces.push_back(f);

			vector<myHalfedge *> face_edges(poly.size(), NULL);
			for (unsigned int j = 0; j < poly.size(); j++)
			{
				myHalfedge *h = new myHalfedge();
				h->source = poly[j];
				h->adjacent_face = f;
				if (h->source->originof == NULL) h->source->originof = h;

				mesh->halfedges.push_back(h);
				face_edges[j] = h;
			}

			for (unsigned int j = 0; j < face_edges.size(); j++)
			{
				unsigned int next_id = (j + 1) % face_edges.size();
				unsigned int prev_id = (j + face_edges.size() - 1) % face_edges.size();
				face_edges[j]->next = face_edges[next_id];
				face_edges[j]->prev = face_edges[prev_id];

				int src = vertex_to_id[poly[j]];
				int dst = vertex_to_id[poly[next_id]];
				pair<int, int> edge_key(src, dst);
				pair<int, int> opposite_key(dst, src);
				map<pair<int, int>, myHalfedge *>::iterator it = twin_map.find(opposite_key);
				if (it != twin_map.end()) {
					face_edges[j]->twin = it->second;
					it->second->twin = face_edges[j];
				}
				twin_map[edge_key] = face_edges[j];
			}

			f->adjacent_halfedge = face_edges[0];
		}
	}
}

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	vector<myHalfedge *>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL)
			break;
	}
	if (it != halfedges.end())
		cout << "Error! Not all edges have their twins!\n";
	else cout << "Each edge has a twin!\n";
}


bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v")
		{
			float x, y, z;
			myline >> x >> y >> z;

			myVertex *v = new myVertex();
			v->point = new myPoint3D(x, y, z);
			vertices.push_back(v);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			faceids.clear();
			while (myline >> u)
			{
				string vstr = u.substr(0, u.find("/"));
				if (vstr.empty()) continue;

				int obj_index = atoi(vstr.c_str());
				int vertex_id = obj_index > 0 ? obj_index - 1 : (int)vertices.size() + obj_index;
				if (vertex_id >= 0 && vertex_id < (int)vertices.size())
					faceids.push_back(vertex_id);
			}
			if (faceids.size() < 3) continue;

			myFace *f = new myFace();
			faces.push_back(f);

			vector<myHalfedge *> face_hedges(faceids.size(), NULL);
			for (unsigned int i = 0; i < faceids.size(); i++)
			{
				myHalfedge *h = new myHalfedge();
				h->source = vertices[faceids[i]];
				h->adjacent_face = f;
				if (h->source->originof == NULL) h->source->originof = h;

				halfedges.push_back(h);
				face_hedges[i] = h;
			}

			for (unsigned int i = 0; i < face_hedges.size(); i++)
			{
				unsigned int next_id = (i + 1) % face_hedges.size();
				unsigned int prev_id = (i + face_hedges.size() - 1) % face_hedges.size();
				face_hedges[i]->next = face_hedges[next_id];
				face_hedges[i]->prev = face_hedges[prev_id];

				pair<int, int> edge_key(faceids[i], faceids[next_id]);
				pair<int, int> opposite_key(faceids[next_id], faceids[i]);
				it = twin_map.find(opposite_key);
				if (it != twin_map.end()) {
					face_hedges[i]->twin = it->second;
					it->second->twin = face_hedges[i];
				}
				twin_map[edge_key] = face_hedges[i];
			}

			f->adjacent_halfedge = face_hedges[0];
		}
	}

	cout << "Vertices: " << vertices.size() << endl;
	cout << "Halfedges: " << halfedges.size() << endl;
	cout << "Faces: " << faces.size() << endl;
	normalize();

	return true;
}


void myMesh::computeNormals()
{
	for (unsigned int i = 0; i < faces.size(); i++)
		if (faces[i]) faces[i]->computeNormal();

	for (unsigned int i = 0; i < vertices.size(); i++)
		if (vertices[i]) vertices[i]->computeNormal();
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}


void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}


void myMesh::subdivisionCatmullClark()
{
	if (vertices.empty()) return;

	vector<myPoint3D*> face_points(faces.size());
	for (size_t i = 0; i < faces.size(); ++i) {
		myFace* f = faces[i];
		myHalfedge* start = f->adjacent_halfedge;
		myHalfedge* h = start;
		double sum_x = 0, sum_y = 0, sum_z = 0;
		int count = 0;
		do {
			sum_x += h->source->point->X;
			sum_y += h->source->point->Y;
			sum_z += h->source->point->Z;
			count++;
			h = h->next;
		} while (h != start);
		face_points[i] = new myPoint3D(sum_x / count, sum_y / count, sum_z / count);
	}

	for (size_t i = 0; i < halfedges.size(); ++i) {
		halfedges[i]->index = -1;
	}

	for (size_t i = 0; i < faces.size(); ++i) {
		faces[i]->index = i;
	}

	for (size_t i = 0; i < vertices.size(); ++i) {
		vertices[i]->index = i;
	}

	vector<myPoint3D*> edge_points;
	int edge_count = 0;
	for (size_t i = 0; i < halfedges.size(); ++i) {
		myHalfedge* h = halfedges[i];
		if (h->index != -1) continue;

		h->index = edge_count;
		if (h->twin) {
			h->twin->index = edge_count;
		}

		myPoint3D* ep = new myPoint3D();
		if (h->twin == NULL) {
			ep->X = 0.5 * (h->source->point->X + h->next->source->point->X);
			ep->Y = 0.5 * (h->source->point->Y + h->next->source->point->Y);
			ep->Z = 0.5 * (h->source->point->Z + h->next->source->point->Z);
		} else {
			myPoint3D* f1 = face_points[h->adjacent_face->index];
			myPoint3D* f2 = face_points[h->twin->adjacent_face->index];
			ep->X = 0.25 * (h->source->point->X + h->next->source->point->X + f1->X + f2->X);
			ep->Y = 0.25 * (h->source->point->Y + h->next->source->point->Y + f1->Y + f2->Y);
			ep->Z = 0.25 * (h->source->point->Z + h->next->source->point->Z + f1->Z + f2->Z);
		}
		edge_points.push_back(ep);
		edge_count++;
	}

	struct VertexInfo {
		int valence;
		double sum_face_x, sum_face_y, sum_face_z;
		double sum_edge_mid_x, sum_edge_mid_y, sum_edge_mid_z;
		bool is_boundary;
		vector<myVertex*> boundary_neighbors;
		VertexInfo() : valence(0), sum_face_x(0), sum_face_y(0), sum_face_z(0),
			sum_edge_mid_x(0), sum_edge_mid_y(0), sum_edge_mid_z(0), is_boundary(false) {}
	};
	vector<VertexInfo> v_info(vertices.size());

	for (size_t i = 0; i < faces.size(); ++i) {
		myFace* f = faces[i];
		myPoint3D* fp = face_points[i];
		myHalfedge* start = f->adjacent_halfedge;
		myHalfedge* h = start;
		do {
			int vi = h->source->index;
			v_info[vi].sum_face_x += fp->X;
			v_info[vi].sum_face_y += fp->Y;
			v_info[vi].sum_face_z += fp->Z;
			h = h->next;
		} while (h != start);
	}

	vector<bool> edge_processed(edge_count, false);
	for (size_t i = 0; i < halfedges.size(); ++i) {
		myHalfedge* h = halfedges[i];
		int idx = h->index;
		if (edge_processed[idx]) continue;
		edge_processed[idx] = true;

		myVertex* v1 = h->source;
		myVertex* v2 = h->next->source;

		v_info[v1->index].valence++;
		v_info[v2->index].valence++;

		double mid_x = 0.5 * (v1->point->X + v2->point->X);
		double mid_y = 0.5 * (v1->point->Y + v2->point->Y);
		double mid_z = 0.5 * (v1->point->Z + v2->point->Z);

		v_info[v1->index].sum_edge_mid_x += mid_x;
		v_info[v1->index].sum_edge_mid_y += mid_y;
		v_info[v1->index].sum_edge_mid_z += mid_z;

		v_info[v2->index].sum_edge_mid_x += mid_x;
		v_info[v2->index].sum_edge_mid_y += mid_y;
		v_info[v2->index].sum_edge_mid_z += mid_z;

		if (h->twin == NULL) {
			v_info[v1->index].is_boundary = true;
			v_info[v1->index].boundary_neighbors.push_back(v2);
			v_info[v2->index].is_boundary = true;
			v_info[v2->index].boundary_neighbors.push_back(v1);
		}
	}

	vector<myVertex*> new_face_vertices(faces.size());
	for (size_t i = 0; i < faces.size(); ++i) {
		myVertex* nv = new myVertex();
		nv->point = face_points[i];
		new_face_vertices[i] = nv;
	}

	vector<myVertex*> new_edge_vertices(edge_count);
	for (int i = 0; i < edge_count; ++i) {
		myVertex* nv = new myVertex();
		nv->point = edge_points[i];
		new_edge_vertices[i] = nv;
	}

	vector<myVertex*> new_vertex_vertices(vertices.size());
	for (size_t i = 0; i < vertices.size(); ++i) {
		myVertex* v = vertices[i];
		myVertex* nv = new myVertex();
		myPoint3D* vp = new myPoint3D();
		if (v_info[i].is_boundary && !v_info[i].boundary_neighbors.empty()) {
			double b_sum_x = 0, b_sum_y = 0, b_sum_z = 0;
			for (size_t j = 0; j < v_info[i].boundary_neighbors.size(); ++j) {
				b_sum_x += v_info[i].boundary_neighbors[j]->point->X;
				b_sum_y += v_info[i].boundary_neighbors[j]->point->Y;
				b_sum_z += v_info[i].boundary_neighbors[j]->point->Z;
			}
			double num_b = v_info[i].boundary_neighbors.size();
			vp->X = 0.75 * v->point->X + 0.25 * (b_sum_x / num_b);
			vp->Y = 0.75 * v->point->Y + 0.25 * (b_sum_y / num_b);
			vp->Z = 0.75 * v->point->Z + 0.25 * (b_sum_z / num_b);
		} else {
			double n = v_info[i].valence;
			if (n > 0) {
				double f_avg_x = v_info[i].sum_face_x / n;
				double f_avg_y = v_info[i].sum_face_y / n;
				double f_avg_z = v_info[i].sum_face_z / n;

				double e_avg_x = v_info[i].sum_edge_mid_x / n;
				double e_avg_y = v_info[i].sum_edge_mid_y / n;
				double e_avg_z = v_info[i].sum_edge_mid_z / n;

				vp->X = (f_avg_x + 2.0 * e_avg_x + (n - 3.0) * v->point->X) / n;
				vp->Y = (f_avg_y + 2.0 * e_avg_y + (n - 3.0) * v->point->Y) / n;
				vp->Z = (f_avg_z + 2.0 * e_avg_z + (n - 3.0) * v->point->Z) / n;
			} else {
				vp->X = v->point->X;
				vp->Y = v->point->Y;
				vp->Z = v->point->Z;
			}
		}
		nv->point = vp;
		new_vertex_vertices[i] = nv;
	}

	vector<vector<myVertex*>> new_polygons;
	for (size_t i = 0; i < faces.size(); ++i) {
		myFace* f = faces[i];
		vector<myHalfedge*> face_hedges;
		myHalfedge* start = f->adjacent_halfedge;
		myHalfedge* h = start;
		do {
			face_hedges.push_back(h);
			h = h->next;
		} while (h != start);

		int k = face_hedges.size();
		for (int j = 0; j < k; ++j) {
			myVertex* V_j = face_hedges[j]->source;
			myHalfedge* h_out = face_hedges[j];
			myHalfedge* h_in = face_hedges[(j - 1 + k) % k];

			vector<myVertex*> quad(4);
			quad[0] = new_vertex_vertices[V_j->index];
			quad[1] = new_edge_vertices[h_out->index];
			quad[2] = new_face_vertices[f->index];
			quad[3] = new_edge_vertices[h_in->index];
			new_polygons.push_back(quad);
		}
	}

	vector<myVertex*> all_new_vertices;
	all_new_vertices.reserve(new_vertex_vertices.size() + new_edge_vertices.size() + new_face_vertices.size());
	for (size_t i = 0; i < new_vertex_vertices.size(); ++i) {
		all_new_vertices.push_back(new_vertex_vertices[i]);
	}
	for (size_t i = 0; i < new_edge_vertices.size(); ++i) {
		all_new_vertices.push_back(new_edge_vertices[i]);
	}
	for (size_t i = 0; i < new_face_vertices.size(); ++i) {
		all_new_vertices.push_back(new_face_vertices[i]);
	}

	vector<myVertex*> old_vertices = vertices;
	vertices = all_new_vertices;

	rebuildMeshFromPolygons(this, new_polygons);

	for (size_t i = 0; i < old_vertices.size(); ++i) {
		if (old_vertices[i]) {
			if (old_vertices[i]->point) delete old_vertices[i]->point;
			delete old_vertices[i];
		}
	}

	computeNormals();
}


void myMesh::triangulate()
{
	vector<vector<myVertex *>> polygons;
	polygons.reserve(faces.size());

	for (unsigned int i = 0; i < faces.size(); i++)
	{
		vector<myVertex *> verts = getFaceVertices(faces[i]);
		if (verts.size() < 3) continue;

		if (verts.size() == 3) {
			polygons.push_back(verts);
			continue;
		}

		vector<vector<myVertex*>> tris = earClipPolygon(verts, *(faces[i]->normal));
		polygons.insert(polygons.end(), tris.begin(), tris.end());
	}

	rebuildMeshFromPolygons(this, polygons);
	computeNormals();
}

bool myMesh::triangulate(myFace *f)
{
	if (f == NULL) return false;

	vector<myVertex *> target = getFaceVertices(f);
	if (target.size() < 3) return false;
	if (target.size() == 3) return false;

	vector<vector<myVertex *>> polygons;
	polygons.reserve(faces.size() + target.size());

	for (unsigned int i = 0; i < faces.size(); i++)
	{
		myFace *curr = faces[i];
		vector<myVertex *> verts = getFaceVertices(curr);
		if (verts.size() < 3) continue;

		if (curr != f) {
			polygons.push_back(verts);
			continue;
		}

		vector<vector<myVertex*>> tris = earClipPolygon(verts, *(curr->normal));
		polygons.insert(polygons.end(), tris.begin(), tris.end());
	}

	rebuildMeshFromPolygons(this, polygons);
	computeNormals();
	return true;
}

void myMesh::simplify()
{
	if (vertices.size() < 4 || halfedges.empty()) return;

	size_t target_count = (size_t)(vertices.size() * 0.9);
	if (target_count >= vertices.size()) target_count = vertices.size() - 1;
	if (target_count < 4) target_count = 4;

	while (vertices.size() > target_count)
	{
		myHalfedge *best_edge = NULL;
		double min_dist = 1e30;

		for (unsigned int i = 0; i < halfedges.size(); i++) {
			myHalfedge *e = halfedges[i];
			if (e == NULL || e->source == NULL || e->twin == NULL || e->twin->source == NULL) continue;
			double d = e->source->point->dist(*(e->twin->source->point));
			if (d < min_dist) {
				min_dist = d;
				best_edge = e;
			}
		}

		if (best_edge == NULL) break;

		myVertex *v1 = best_edge->source;
		myVertex *v2 = best_edge->twin->source;

		v1->point->X = (v1->point->X + v2->point->X) / 2.0;
		v1->point->Y = (v1->point->Y + v2->point->Y) / 2.0;
		v1->point->Z = (v1->point->Z + v2->point->Z) / 2.0;

		vector<vector<myVertex *>> new_polygons;
		for (unsigned int i = 0; i < faces.size(); i++)
		{
			vector<myVertex *> face_verts = getFaceVertices(faces[i]);
			vector<myVertex *> rep_verts;
			for (unsigned int j = 0; j < face_verts.size(); j++)
			{
				myVertex *curr = (face_verts[j] == v2) ? v1 : face_verts[j];
				if (rep_verts.empty() || rep_verts.back() != curr) {
					rep_verts.push_back(curr);
				}
			}
			if (rep_verts.size() > 1 && rep_verts.front() == rep_verts.back()) {
				rep_verts.pop_back();
			}
			if (rep_verts.size() >= 3) {
				new_polygons.push_back(rep_verts);
			}
		}

		vertices.erase(std::remove(vertices.begin(), vertices.end(), v2), vertices.end());
		delete v2;

		rebuildMeshFromPolygons(this, new_polygons);
	}
}

