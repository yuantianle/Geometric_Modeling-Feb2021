/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__


#include "ply.h"
#include "icVector.H"
#include "glm/glm/glm.hpp"
#include "glm/glm/gtc/matrix_transform.hpp"
#include "glm/glm/gtc/type_ptr.hpp"
#include <vector>

const double EPS = 1.0e-6;
const double PI=3.1415926535898;

/* forward declarations */
class Triangle;
class Vertex;
class Edge;

class Corner {
public:
	double c;
	int index;

	Corner* o;  //opposite corner
	Corner* p;  //privious corner
	Corner* n;  //next corner
	Vertex* v; //corner's vertex
	Edge* e;
	Triangle* t;

//public:
//	Corner(double cc) {c = cc;}
};


class Vertex {
public:
  double x,y,z;
  int index;
  double r, g, b;   //heat color
  double tx, ty;    //texture coordinates  

  int ntris;
  Triangle **tris;
  int max_tris;

	icVector3 normal;
  void *other_props;

  double update_vert_pos[3];   //save the new vertex point location that need to update
  int edge_index;	//only for the vertex on the edge;

  double heatvalue;  //between 0-1
  double heatvaluesave; //for saving old temperature

  std::vector<double> neighborweight;

  double back_flag;  //for silhoutte

  double Amixed;
  double meancurvaturevalue[3]; // for mean curvature 
  double meancurvaturemagnitude;

  double gausscurvaturevalue; // for gauss curvature

  double principlecurvaturevector1[3];
  double principlecurvaturevector2[3];

  double TangentVector1[3];  // The two parpenticular vectors
  double TangentVector2[3];
  double Tensor[3];   // 0-l, 1-m, 2-n
  double TensorSave[3];

  glm::mat3 t1t2N;
  glm::mat3 Tglobal;
  glm::mat3 TglobalSave;

public:
  Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
};

class Edge {
public:
  double length;
  int index;
  Vertex* verts[2];
  int ntris;
  Triangle **tris;

  Vertex *new_vert;      //save the new vertex on the edge
  //Edge() { index = 0; ntris = 0; length = 0; tris = NULL; new_vert = NULL; }

  double zerovecpoint[3];
  double zeroflag;
};

class Triangle {
public:
  int index;
  int nverts;
  Vertex* verts[3];
  Edge* edges[3];

	//double angle[3];
  Corner* corners[3]; //the pointer from triangle to corners
	float area;

	icVector3 normal;
  void *other_props;
  double back_flag;
  double shadow_flag;
  //Triangle() { index = 0; nverts = 0; area = 0.0; other_props = NULL; normal = NULL; }
  

  double Tangent1[3];
  double Tangent2[3];
  double Tensor2D[3];
  glm::mat3 RotateMatrix;
  glm::mat3 Tensor3D;
  double principleTricurvaturevector1[3];
  double principleTricurvaturevector2[3];

  glm::vec3 CenterPoint;
  glm::vec3 StreamLinePoints[13];

};

//class CornerList {
//
//public:
//
//	int num_corners;
//	int max_corners;
//	Corner** corners;
//
//	CornerList(int max) {
//		num_corners = max;
//		corners = new Corner * [max_corners];
//		num_corners = 0;
//	}
//	void finalize() {
//		if (this != NULL) {
//			free(corners);
//			free(this);
//		}
//	}
//	void append(Corner* c)
//	{
//		int i;
//
//		/* first make sure there is enough room for new vertex */
//
//		if (num_corners >= max_corners) {
//			max_corners += 10;
//			Corner** tlist = new Corner * [max_corners];
//			for (i = 0; i < num_corners; i++)
//				tlist[i] = corners[i];
//			delete (corners);
//			corners = tlist;
//		}
//
//		/* add new vertex to list */
//
//		corners[num_corners] = c;
//		num_corners++;
//	}
//
//};

//class VertexList {
//
//public:
//
//  int num_verts;
//  int max_verts;
//  Vertex **verts;
//
//  VertexList (int max) {
//    max_verts = max;
//    verts = new Vertex *[max_verts];
//    num_verts = 0;
//  }
//	void finalize() {
//		if (this != NULL) {
//			free(verts);
//			free(this);
//		}
//	}
//	void append(Vertex *v)
//	{
//		int i;
//
//		/* first make sure there is enough room for new vertex */
//
//		if (num_verts >= max_verts) {
//			max_verts += 10;
//			Vertex **tlist = new Vertex *[max_verts];
//			for (i = 0; i < num_verts; i++)
//				tlist[i] = verts[i];
//			delete (verts);
//			verts = tlist;
//		}
//
//		/* add new vertex to list */
//
//		verts[num_verts] = v;
//		num_verts++;
//	}
//
//};
//
//class TriangleList {
//public:
//
//  int num_tris;
//  int max_tris;
//  Triangle **tris;
//
//  TriangleList (int max) {
//    max_tris = max;
//    tris = new Triangle *[max_tris];
//    num_tris = 0;
//  }
//	void append(Triangle *t)
//	{
//		int i;
//
//		/* first make sure there is enough room for new triangle */
//
//		if (num_tris >= max_tris) {
//			max_tris += 10;
//			Triangle **tlist = new Triangle *[max_tris];
//			for (i = 0; i < num_tris; i++)
//				tlist[i] = tris[i];
//			delete (tris);
//			tris = tlist;
//		}
//
//		/* add new triangle to list */
//
//		tris[num_tris] = t;
//		num_tris++;
//	}
//};
//
//class EdgeList {
//
//public:
//
//  int num_edges;
//  int max_edges;
//  Edge **edges;
//
//  EdgeList (int max) {
//    max_edges = max;
//    edges = new Edge *[max_edges];
//    num_edges = 0;
//  }
//	void append(Edge *e)
//	{
//		int i;
//
//		/* first make sure there is enough room for new edge */
//
//		if (num_edges >= max_edges) {
//			max_edges += 10;
//			Edge **tlist = new Edge *[max_edges];
//			for (i = 0; i < num_edges; i++)
//				tlist[i] = edges[i];
//			delete (edges);
//			edges = tlist;
//		}
//
//		/* add new edge to list */
//
//		edges[num_edges] = e;
//		num_edges++;
//	}
//};

class Polyhedron {
public:

	int index;
//------------------------------------
	Corner** clist;
	int ncorners;
	int max_corners;
//------------------------------------
  Triangle **tlist;  /* list of triangles */
  int ntris;
  int max_tris;

  Vertex **vlist;    /* list of vertices */
  int nverts;
  int max_verts;

  Edge **elist;      /* list of edges */
  int nedges;
  int max_edges;

	icVector3 center;
	double radius;
	double area;

	int seed;

  PlyOtherProp *vert_other,*face_other;

  unsigned int* tempA, * tempB, * tempC;

	void average_normals();

  void create_edge(Vertex *, Vertex *);
  void create_edges();
  //------------------------------------
  //---------Corner Structure-----------
  void create_corner(Vertex*, Vertex*, Vertex*);
  void create_corners();
  double calc_corner_angle(Vertex*, Vertex*, Vertex*);
  void setting_o_p_n();
  void corner_sort();
  void sort(unsigned int* A, unsigned int* B, unsigned int* C, unsigned int sid, unsigned int eid);
  //------------------------------------
  //----------Subdivision---------------
  void update_old_vert(Polyhedron* poly);
  void create_new_vert(Polyhedron* poly);
  void create_new_ivert(Polyhedron* poly);
  void regular_sub(Polyhedron* poly);
  void irregular_sub(Polyhedron* poly);
  //------------------------------------
  //-----------Smoothing----------------
  void update_old_vert_uniform(Polyhedron* poly);
  void update_old_vert_cord(Polyhedron* poly);
  void update_old_vert_meancurvature(Polyhedron* poly);
  void update_old_vert_meanvalue(Polyhedron* poly);

  void set_old_vert_cord_oldweight();       //weight not change
  void set_old_vert_meancurvature_oldweight();
  void set_old_vert_meanvalue_oldweight();
  void update_old_vert_all_old(Polyhedron* poly);

  void smooth1(Polyhedron* poly);
  void smooth(Polyhedron* poly);
  //------------------------------------
  //----------Function h----------------
  void setheat();
  int maxindex;// for the heat maxpoint
  int minindex;// for the heat minpoint
  void heat1(Polyhedron* poly, int maxindex, int minindex);  // Explicit
  void heat2(Polyhedron* poly, int maxindex, int minindex, int iterationtime);  // Implicit
  void heat3(Polyhedron* poly, int maxindex, int minindex);  // G-S
  void color_t_rainbow();
  int saddle(Vertex* vertex);
  void mcalculate();
  int Mf;
  //------------------------------------
  //----------Texturing----------------
  std::vector<int> texcoordindex;
  void settextcoord();
  void texmap(Polyhedron* poly);
  //------------------------------------
  //-----------Silhouette---------------
  void Silhouette1(Polyhedron* poly);   //for face method
  void Silhouette2(Polyhedron* poly);   //for vertex method
  //------------------------------------
  //-----------Curvature----------------
  void CalculateAmixed();
  void MeanCurvature();
  double maxmeancurv;
  double minmeancurv;
  void GaussCurvature();
  double maxgausscurv;
  double mingausscurv;
  void PrincipalCurvature();
  void CalculatePrincipalCurvature(int vindex);

  void GetTensor3D();
  void TensorSmoothing1();
  void GoBackTensor2D();
  void TensorSmoothing2(int iteratetime);

  void Hatching();
  void TriangleTensoring();
  void GoBackTriTensor2D();
  void CalculateTriPrincipalCurvature(int tindex);

  void SettingStreamlineFold();
  void SettingStreamlineUnfold();
  //------------------------------------


  int face_to_vertex_ref(Triangle *, Vertex *);
  void order_vertex_to_tri_ptrs(Vertex *);
  void vertex_to_tri_ptrs();
  Triangle *find_common_edge(Triangle *, Vertex *, Vertex *);  //for edges
  Triangle *other_triangle(Edge *, Triangle *);
	void calc_bounding_sphere();
	void calc_face_normals_and_area();
	void calc_edge_length();

	Polyhedron();
  Polyhedron(FILE *);
  void write_file(FILE *);

  void create_pointers();

	// initialization and finalization
	void initialize();
	void finalize();
};

#endif /* __LEARNPLY_H__ */

//if (longedgecount == 1)
//{
//	//push 2 new triangles
//	tlist[nowindex] = new Triangle;
//	tlist[nowindex]->nverts = 3;
//
//	Edge* longedge = new Edge;
//	int oppvertindex;
//	for (int j = 0; j < 3; j++)
//	{
//		if (poly->tlist[i]->edges[j]->new_vert != NULL)
//		{
//			longedge = poly->tlist[i]->edges[j];
//
//			for (int k = 0; k < 3; k++)
//			{
//				if (poly->tlist[i]->verts[k]->index != longedge->verts[0]->index && poly->tlist[i]->verts[k]->index != longedge->verts[1]->index)
//					oppvertindex = k;
//			}
//			break;
//		}
//	}
//
//	for (int p = 0; p < poly->max_verts; p++)
//		if (vlist[p]->index == poly->tlist[i]->verts[oppvertindex]->index)
//			tlist[nowindex]->verts[0] = vlist[p];
//
//	for (int p = 0; p < poly->max_verts; p++)
//		if (vlist[p]->index == poly->tlist[i]->verts[(oppvertindex + 1) % 3]->index)
//			tlist[nowindex]->verts[1] = vlist[p];
//
//	for (int p = poly->max_verts; p < max_verts; p++)
//		if (vlist[p]->edge_index == longedge->index)
//			tlist[nowindex]->verts[2] = vlist[p];
//	ntris += 1;
//
//
//	tlist[nowindex + 1] = new Triangle;
//	tlist[nowindex + 1]->nverts = 3;
//
//	for (int p = 0; p < poly->max_verts; p++)
//		if (vlist[p]->index == poly->tlist[i]->verts[oppvertindex]->index)
//			tlist[nowindex + 1]->verts[0] = vlist[p];
//
//	for (int p = poly->max_verts; p < max_verts; p++)
//		if (vlist[p]->edge_index == longedge->index)
//			tlist[nowindex + 1]->verts[1] = vlist[p];
//
//	for (int p = 0; p < poly->max_verts; p++)
//		if (vlist[p]->index == poly->tlist[i]->verts[(oppvertindex + 2) % 3]->index)
//			tlist[nowindex + 1]->verts[2] = vlist[p];
//	ntris += 1;
//
//	nowindex += (longedgecount + 1);
//}

