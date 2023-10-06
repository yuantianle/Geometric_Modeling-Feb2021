/*

Functions for learnply

Eugene Zhang, 2005
*/

/*

Edit by

Tianle Yuan, 2021
*/

/* Interface Guidance:
*  
1=> Original shape with tensor field
2=>	Regular subdivision
3=>	Irregular subdivision
4=> Smooth
5=> Heat diffusion
6=> Original shape with white color
7=> Shape with 3D checkerboard color 
or 3D checkerboard coordinates texture
9=> Silhouette
P=> Tensor Smoothing

s=> Swap smooth and flat
t=> Triangle mesh
*/

#include <stdlib.h>
#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <math.h>
#include "glut.h"
#include <string.h>
#include <fstream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "trackball.h"
#include "tmatrix.h"

#include "glm/glm/glm.hpp"
#include "glm/glm/gtc/matrix_transform.hpp"
#include "glm/glm/gtc/type_ptr.hpp"
#include <vector>

#include <iostream>
#include <cstdio>
#include <ctime>
#include <math.h> 
#include "bmptotexture.h"
#include "RIdeformation.h"


typedef Eigen::SparseMatrix<double> SpMat;  //declare a double sparse matrix type
typedef Eigen::Triplet<double> T; //triple:(row, column, value)
typedef glm::vec3 vec3;
static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;

int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float tridepth;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

int sflag = 1;    // for the smooth
int tflag = 1;    // for the triangle grids
int silhouetteflag = 1;  // for the silhouette grids

int maxflag = -1;

GLuint Tex0;
int textureflag = -1;
int heat2litertime = 0;

int chooseflag = -1; //for the anchor points


struct jitter_struct{
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125}, 
						  {0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375}, 
						  {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625}, 
						  {0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };

Polyhedron *poly;
Polyhedron * rpoly;
Polyhedron* irpoly;
Polyhedron* uniformsmoothpoly;
Polyhedron* cordsmoothpoly;

void init(void);
void inittexture(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char *argv[])
{
  char *progname;
  int num = 1;
	FILE *this_file;

  progname = argv[0];

	//this_file = fopen("../tempmodels/tetrahedron.ply", "r");
	//this_file = fopen("../tempmodels/octahedron.ply", "r");
    //this_file = fopen("../tempmodels/icosahedron.ply", "r");
	//this_file = fopen("../tempmodels/sphere.ply", "r");
	//this_file = fopen("../tempmodels/bunny.ply", "r");
	//this_file = fopen("../tempmodels/cow.ply", "r");
	//this_file = fopen("../tempmodels/feline.ply", "r");
	//this_file = fopen("../tempmodels/dragon.ply", "r");
	//this_file = fopen("../tempmodels/happy.ply", "r");
	//this_file = fopen("../tempmodels/torus.ply", "r");
	this_file = fopen("../tempmodels/cube.ply", "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );	

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	//poly->set_old_vert_cord_oldweight();
	poly->set_old_vert_meancurvature_oldweight();
	//poly->set_old_vert_meanvalue_oldweight();

	poly->setheat();
	poly->mcalculate();
	poly->color_t_rainbow();

	poly->settextcoord();

	poly->GaussCurvature();
	poly->MeanCurvature();
	poly->PrincipalCurvature();


	//poly->GetFrameForV();


	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition (20, 20);
	glutInitWindowSize (win_width, win_height); 
	glutCreateWindow ("Geometric Modeling");
	init ();
	glutKeyboardFunc (keyboard);
	glutDisplayFunc(display); 
	glutMotionFunc (motion);
	glutMouseFunc (mouse);

	inittexture();


	glutMainLoop(); 
	//---------------------------------------------------------------------




	poly->finalize();  // finalize everything

  return 0;    /* ANSI C requires main to return int. */
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0){
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0/3){
		col[0] = 1.0;
		col[1] = 1.0-percentage*3.0;
		col[2] = 1.0-percentage*3.0;
	}
	else if (percentage <= 2.0/3){
		col[0] = 1.0;
		col[1] = percentage*3.0-1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0/3){
		col[0] = 3.0-percentage*3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/******************************************************************************
Read in a polyhedron from a file.//--------------------------------------------Useful information of initializing------------
******************************************************************************/

Polyhedron::Polyhedron(FILE *file)
{
  int i,j;
  int elem_count;
  char *elem_name;

  /*** Read in the original PLY object ***/
  in_ply = read_ply (file);

  for (i = 0; i < in_ply->num_elem_types; i++) {

    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply (in_ply, i, &elem_count);

    if (equal_strings ("vertex", elem_name)) {

      /* create a vertex list to hold all the vertices */
      nverts = max_verts = elem_count;
      vlist = new Vertex *[nverts];

      /* set up for getting vertex elements */

      setup_property_ply (in_ply, &vert_props[0]);
      setup_property_ply (in_ply, &vert_props[1]);
      setup_property_ply (in_ply, &vert_props[2]);
      vert_other = get_other_properties_ply (in_ply, 
					     offsetof(Vertex_io,other_props));

      /* grab all the vertex elements */
      for (j = 0; j < nverts; j++) {
        Vertex_io vert;
		vert.other_props = NULL;
        get_element_ply (in_ply, (void *) &vert);

        /* copy info from the "vert" structure */
        vlist[j] = new Vertex (vert.x, vert.y, vert.z);
        vlist[j]->other_props = vert.other_props;
      }
    }
    else if (equal_strings ("face", elem_name)) {

      /* create a list to hold all the face elements */
      ntris = max_tris = elem_count;
      tlist = new Triangle *[ntris];

      /* set up for getting face elements */
      setup_property_ply (in_ply, &face_props[0]);
      face_other = get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

      /* grab all the face elements */
      for (j = 0; j < elem_count; j++) {
        Face_io face;
		face.other_props = NULL;
        get_element_ply (in_ply, (void *) &face);

        if (face.nverts != 3) {
          fprintf (stderr, "Face has %d vertices (should be three).\n",
                   face.nverts);
          exit (-1);
        }

        /* copy info from the "face" structure */
        tlist[j] = new Triangle;
        tlist[j]->nverts = 3;
        tlist[j]->verts[0] = (Vertex *) face.verts[0];
        tlist[j]->verts[1] = (Vertex *) face.verts[1];
        tlist[j]->verts[2] = (Vertex *) face.verts[2];
        tlist[j]->other_props = face.other_props;
      }
    }
    else
      get_other_element_ply (in_ply);
  }

  /* close the file */
  close_ply (in_ply);

  /* fix up vertex pointers in triangles */
  for (i = 0; i < ntris; i++) {
    tlist[i]->verts[0] = vlist[(int) tlist[i]->verts[0]];
    tlist[i]->verts[1] = vlist[(int) tlist[i]->verts[1]];
    tlist[i]->verts[2] = vlist[(int) tlist[i]->verts[2]];
  }

  /* get rid of triangles that use the same vertex more than once */

  for (i = ntris-1; i >= 0; i--) {

    Triangle *tri = tlist[i];
    Vertex *v0 = tri->verts[0];
    Vertex *v1 = tri->verts[1];
    Vertex *v2 = tri->verts[2];

    if (v0 == v1 || v1 == v2 || v2 == v0) {
		if (tlist[i])
		{
			free(tlist[i]);
			tlist[i] = NULL;
		}
      ntris--;
      tlist[i] = tlist[ntris];
    }
  }

  
}

/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/

void Polyhedron::write_file(FILE *file)
{
  int i;
  PlyFile *ply;
  char **elist;
  int num_elem_types;

  /*** Write out the transformed PLY object ***/

  elist = get_element_list_ply (in_ply, &num_elem_types);
  ply = write_ply (file, num_elem_types, elist, in_ply->file_type);

  /* describe what properties go into the vertex elements */

  describe_element_ply (ply, "vertex", nverts);
  describe_property_ply (ply, &vert_props[0]);
  describe_property_ply (ply, &vert_props[1]);
  describe_property_ply (ply, &vert_props[2]);
//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

  describe_element_ply (ply, "face", ntris);
  describe_property_ply (ply, &face_props[0]);

//  describe_other_properties_ply (ply, face_other,
//                                offsetof(Face_io,other_props));

//  describe_other_elements_ply (ply, in_ply->other_elems);

  copy_comments_ply (ply, in_ply);
	char mm[1024];
	sprintf(mm, "modified by learnply");
//  append_comment_ply (ply, "modified by simvizply %f");
	  append_comment_ply (ply, mm);
  copy_obj_info_ply (ply, in_ply);

  header_complete_ply (ply);

  /* set up and write the vertex elements */
  put_element_setup_ply (ply, "vertex");
  for (i = 0; i < nverts; i++) {
    Vertex_io vert;

    /* copy info to the "vert" structure */
    vert.x = vlist[i]->x;
    vert.y = vlist[i]->y;
    vert.z = vlist[i]->z;
    vert.other_props = vlist[i]->other_props;

    put_element_ply (ply, (void *) &vert);
  }

  /* index all the vertices */
  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  /* set up and write the face elements */
  put_element_setup_ply (ply, "face");

  Face_io face;
  face.verts = new int[3];
  
  for (i = 0; i < ntris; i++) {

    /* copy info to the "face" structure */
    face.nverts = 3;
    face.verts[0] = tlist[i]->verts[0]->index;
    face.verts[1] = tlist[i]->verts[1]->index;
    face.verts[2] = tlist[i]->verts[2]->index;
    face.other_props = tlist[i]->other_props;

    put_element_ply (ply, (void *) &face);
  }
  put_other_elements_ply (ply);

  close_ply (ply);
  free_ply (ply);
}

void Polyhedron::initialize(){
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	setting_o_p_n();
	seed = -1;
}

void Polyhedron::finalize(){
	int i;

	for (i=0; i<ntris; i++){
		if (tlist[i]->other_props)
		{
			free(tlist[i]->other_props);
			tlist[i]->other_props = NULL;
		}
		//-----------------------
		for (int j = 0; j < 3; j++)
		{
			if (tlist[i]->corners)
			{
				if (tlist[i]->corners[j])
				{
					free(tlist[i]->corners[j]);
					tlist[i]->corners[j] = NULL;
				}
			}
		}

		//-----------------------
		if (tlist[i])
		{
			free(tlist[i]);
			tlist[i] = NULL;
		}
	}
	for (i=0; i<nedges; i++) {
		if (elist[i]->tris)
		{
			free(elist[i]->tris);
			elist[i]->tris = NULL;
		}
		if (elist[i])
		{
			free(elist[i]);
			elist[i] = NULL;
		}
	}
	for (i=0; i<nverts; i++) {
		if (vlist[i]->tris)
		{
			free(vlist[i]->tris);
			vlist[i]->tris = NULL;
		}
		if (vlist[i]->other_props)
		{
			free(vlist[i]->other_props);
			vlist[i]->other_props = NULL;
		}
		if (vlist[i])
		{
			free(vlist[i]);
			vlist[i] = NULL;
		}
	}
	
	if (tlist)
	{
		free(tlist);
		tlist = NULL;
	}
	if (elist)
	{
		free(elist);
		elist = NULL;
	}
	if (vlist)
	{
		free(vlist);
		vlist = NULL;
	}
	//-----------------------
	if (clist)
	{
		free(clist);
		clist = NULL;
	}

	//-----------------------
	vert_other = NULL;
	face_other = NULL;
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face
******************************************************************************/

Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f2;
  Triangle *adjacent = NULL;

  /* look through all faces of the first vertex */

  for (i = 0; i < v1->ntris; i++) {
    f2 = v1->tris[i];
    if (f2 == f1)
      continue;
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < f2->nverts; j++) {

      /* look for a match */
      if (f2->verts[j] == v2) {

//#if 0
//	/* watch out for triple edges */
//
//        if (adjacent != NULL) {
//
//	  fprintf (stderr, "model has triple edges\n");
//
//	  fprintf (stderr, "face 1: ");
//	  for (k = 0; k < f1->nverts; k++)
//	    fprintf (stderr, "%d ", f1->iverts[k]);
//	  fprintf (stderr, "\nface 2: ");
//	  for (k = 0; k < f2->nverts; k++)
//	    fprintf (stderr, "%d ", f2->iverts[k]);
//	  fprintf (stderr, "\nface 3: ");
//	  for (k = 0; k < adjacent->nverts; k++)
//	    fprintf (stderr, "%d ", adjacent->iverts[k]);
//	  fprintf (stderr, "\n");
//
//	}

	/* if we've got a match, remember this face */
        adjacent = f2;
//#endif

//#if 1
	/* if we've got a match, return this face */
        return (f2);
//#endif

      }
    }
  }

  return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/

void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f;

  /* make sure there is enough room for a new edge */

  if (nedges >= max_edges) {

    max_edges += 100;
    Edge **list = new Edge *[max_edges];

    /* copy the old list to the new one */
	for (i = 0; i < nedges; i++)
	{
		list[i] = elist[i];
		if (elist[i])
		{
			free(elist[i]);
			elist[i] = NULL;
		}
	}
      

    /* replace list */
    //free (elist);
    elist = list;
  }

  /* create the edge */

  elist[nedges] = new Edge;
  Edge *e = elist[nedges];
  e->index = nedges;
  e->verts[0] = v1;
  e->verts[1] = v2;
  nedges++;

  /* count all triangles that will share the edge, and do this */
  /* by looking through all faces of the first vertex */

  e->ntris = 0;

  for (i = 0; i < v1->ntris; i++) {
    f = v1->tris[i];
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++) {
      /* look for a match */
      if (f->verts[j] == v2) {
        e->ntris++;
        break;
      }
    }
  }

  /* make room for the face pointers (at least two) */
  if (e->ntris < 2)
    e->tris = new Triangle *[2];
  else
    e->tris = new Triangle *[e->ntris];

  /* create pointers from edges to faces and vice-versa */

  e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

  for (i = 0; i < v1->ntris; i++) {

    f = v1->tris[i];

    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++)
      if (f->verts[j] == v2) {

        e->tris[e->ntris] = f;
        e->ntris++;

        if (f->verts[(j+1)%3] == v1)
          f->edges[j] = e;
        else if (f->verts[(j+2)%3] == v1)
          f->edges[(j+2)%3] = e;
        else {
          fprintf (stderr, "Non-recoverable inconsistancy in create_edge()\n");
          exit (-1);
        }

        break;  /* we'll only find one instance of v2 */
      }

  }
}

double Polyhedron::calc_corner_angle(Vertex* v1, Vertex* v2, Vertex* v3)
{
	double angle;
	Vertex* nowpoint = v1;
	Vertex* otherpoint1 = v2;
	Vertex* otherpoint2 = v3;

	vec3 vector1 = vec3(otherpoint1->x - nowpoint->x, otherpoint1->y - nowpoint->y, otherpoint1->z - nowpoint->z);
	vec3 vector2 = vec3(otherpoint2->x - nowpoint->x, otherpoint2->y - nowpoint->y, otherpoint2->z - nowpoint->z);
	//angle = 180.0*glm::acos(glm::dot(vector1, vector2)/(glm::length(vector1)* glm::length(vector2)))/3.1415926;
	angle = glm::acos(glm::dot(vector1, vector2) / (glm::length(vector1) * glm::length(vector2)));
	return angle;
}


/******************************************************************************
Create an corner.

Entry:
  v1,v2,v3 - 3 vertices of f1 that define corner
******************************************************************************/

void Polyhedron::create_corner(Vertex* v1, Vertex* v2, Vertex* v3)
{
	int i, j;
	Triangle* f;

	/* make sure there is enough room for a new corner */

	if (ncorners >= max_corners) {

		max_corners += 100;
		Corner** list = new Corner * [max_corners];

		/* copy the old list to the new one */
		for (i = 0; i < ncorners; i++)
		{
			list[i] = clist[i];
			if (clist[i])
			{
				free(clist[i]);
				clist[i] = NULL;
			}
		}
			

		/* replace list */
		//free(clist);
		clist = list;
	}

	/* create the corner */    //--------------------------------Place FIX all Corner MEMBERS---------------------------------

	clist[ncorners];// = new Corner;
	Corner* c = clist[ncorners];
	
	c->index = ncorners;
	c->c = calc_corner_angle(v1, v2, v3);
	//c->p = ; //should be later after we create all the corners
	//c->n = ;
	//c->o = ;//should be the last one based on max and min operation of c.p.v. c.n.v.
	c->v = v1;
	//c->e = ; //should do in the upper level;


	ncorners++;

}

void Polyhedron::setting_o_p_n()
{
	for (int i = 0; i < ntris; i++)
	{
		Triangle* f = tlist[i];
		for (int j = 0; j < 3; j++)
		{
			f->corners[j]->p = new Corner;
			f->corners[j]->p = f->corners[(j+2)%3];

			f->corners[j]->n = new Corner;
			f->corners[j]->n = f->corners[(j + 1) % 3];
		}
	}

	corner_sort();
}

/******************************************************************************
Create corners.//---------------------------------------------------------------Create corner stucture-------------------
******************************************************************************/
void Polyhedron::create_corners()
{
	int i, j;
	Triangle* f;
	Vertex* v1, * v2, * v3;
	double count = 0;

	/* count up how many edges we may require */
	count = ntris * 3.0;

	/*
	printf ("counted %f corners\n", count);
	*/

	/* create space for corner list */

	max_corners = (int)(count + 10);  /* leave some room for expansion */
	clist = new Corner * [max_corners];
	for (int i = 0; i < max_corners; i++)
	{
		clist[i] = new Corner;
		clist[i]->e = NULL;
		clist[i]->t = NULL;
		clist[i]->n = NULL; 
		clist[i]->p = NULL; 
		clist[i]->o = NULL;
	}


	ncorners = 0;

	/* zero out all the pointers from faces to corners */

	for (i = 0; i < ntris; i++)
		for (j = 0; j < 3; j++)
		{
			tlist[i]->corners[j] = NULL;
		}

	/* create all the corners by examining all the triangles */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < 3; j++) { //scan 3 angles of a triangle to create correspoding corners
			/* skip over corners that we've already created */
			if (f->corners[j])
				continue;
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			v3 = f->verts[(j + 2) % f->nverts];

			create_corner(v1, v2, v3);
			// make room for the face pointer 
			clist[ncorners-1]->e = f->edges[(j + 1) % f->nverts]; // set the edge of corner

			// make room for the face pointer 
			clist[ncorners-1]->t = tlist[i]; // set the triangle of corner
			
			f->corners[j] = clist[ncorners-1]; // set the pointer from face to corner
		}
	}
}

/******************************************************************************
Create edges.
******************************************************************************/

void Polyhedron::create_edges()
{
  int i,j;
  Triangle *f;
  Vertex *v1,*v2;
  double count = 0;

  /* count up how many edges we may require */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      Triangle *result = find_common_edge (f, v1, v2);
      if (result)
        count += 0.5;
      else
        count += 1;
    }
  }

  /*
  printf ("counted %f edges\n", count);
  */

  /* create space for edge list */

  max_edges = (int) (count + 10);  /* leave some room for expansion */
  elist = new Edge *[max_edges];
  nedges = 0;

  /* zero out all the pointers from faces to edges */

  for (i = 0; i < ntris; i++)
    for (j = 0; j < 3; j++)
      tlist[i]->edges[j] = NULL;

  /* create all the edges by examining all the triangles */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < 3; j++) {
      /* skip over edges that we've already created */
      if (f->edges[j])
        continue;
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      create_edge (v1, v2);
    }
  }
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_tri_ptrs()
{
  int i,j;
  Triangle *f;
  Vertex *v;

  /* zero the count of number of pointers to faces */

  for (i = 0; i < nverts; i++)
    vlist[i]->max_tris = 0;

  /* first just count all the face pointers needed for each vertex */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++)
      f->verts[j]->max_tris++;
  }

  /* allocate memory for face pointers of vertices */

  for (i = 0; i < nverts; i++) {
    vlist[i]->tris = (Triangle **)
		      malloc (sizeof (Triangle *) * vlist[i]->max_tris);
    vlist[i]->ntris = 0;
  }

  /* now actually create the face pointers */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v = f->verts[j];
      v->tris[v->ntris] = f;
      v->ntris++;
    }
  }
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
  /* search for any other triangle */

  for (int i = 0; i < edge->ntris; i++)
    if (edge->tris[i] != tri)
      return (edge->tris[i]);

  /* there is no such other triangle if we get here */
  return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered
******************************************************************************/

void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
  int i,j;
  Triangle *f;
  Triangle *fnext;
  int nf;
  int vindex;
  int boundary;
  int count;

  nf = v->ntris;
  f = v->tris[0];

  /* go backwards (clockwise) around faces that surround a vertex */
  /* to find out if we reach a boundary */

  boundary = 0;

  for (i = 1; i <= nf; i++) {

    /* find reference to v in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[j] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #1\n");
      exit (-1);
    }

    /* corresponding face is the previous one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* see if we've reached a boundary, and if so then place the */
    /* current face in the first position of the vertice's face list */

    if (fnext == NULL) {
      /* find reference to f in v */
      for (j = 0; j < v->ntris; j++)
        if (v->tris[j] == f) {
	  v->tris[j] = v->tris[0];
	  v->tris[0] = f;
	  break;
	}
      boundary = 1;
      break;
    }

    f = fnext;
  }

  /* now walk around the faces in the forward direction and place */
  /* them in order */

  f = v->tris[0];
  count = 0;

  for (i = 1; i < nf; i++) {

    /* find reference to vertex in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[(j+1) % f->nverts] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #2\n");
      exit (-1);
    }

    /* corresponding face is next one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* break out of loop if we've reached a boundary */
    count = i;
    if (fnext == NULL) {
      break;
    }

    /* swap the next face into its proper place in the face list */
    for (j = 0; j < v->ntris; j++)
      if (v->tris[j] == fnext) {
	v->tris[j] = v->tris[i];
	v->tris[i] = fnext;
	break;
      }

    f = fnext;
  }
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
  int j;
  int vindex = -1;

  for (j = 0; j < f->nverts; j++)
    if (f->verts[j] == v) {
      vindex = j;
      break;
    }

  return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.//---------------------------------------------Coner Pointer Location--------------
******************************************************************************/
// now the triangle list data and vetices data are all saved in the memberships.
void Polyhedron::create_pointers()
{
  int i;

  /* index the vertices and triangles */

  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  for (i = 0; i < ntris; i++) 
    tlist[i]->index = i;

  /* create pointers from vertices to triangles */
  vertex_to_tri_ptrs();

  /* make edges */
  create_edges();


  /* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
    order_vertex_to_tri_ptrs(vlist[i]);
		
	}
  /* index the edges */

  for (i = 0; i < nedges; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
    elist[i]->index = i;
	}

  create_corners();
}

void Polyhedron::calc_bounding_sphere()
{
  unsigned int i;
  icVector3 min, max;

  for (i=0; i<nverts; i++) {
    if (i==0)  {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
    }
    else {
			if (vlist[i]->x < min.entry[0])
			  min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
			  max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
			  min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
			  max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
			  min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
			  max.entry[2] = vlist[i]->z;
		}
  }
  center = (min + max) * 0.5;
  radius = length(center - min);
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i=0; i<nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1-v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
  Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i=0; i<ntris; i++){
		for (j=0; j<3; j++)
			length[j] = tlist[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2])/2.0;
		tlist[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		v2.set(vlist[tlist[i]->verts[1]->index]->x, vlist[tlist[i]->verts[1]->index]->y, vlist[tlist[i]->verts[1]->index]->z);
		v0.set(vlist[tlist[i]->verts[2]->index]->x, vlist[tlist[i]->verts[2]->index]->y, vlist[tlist[i]->verts[2]->index]->z);
		tlist[i]->normal = cross(v0-v1, v2-v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i=0; i<ntris; i++){
		icVector3 cent(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		signedvolume += dot(test-cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume<0) 
		orientation = 0;
	else {
		orientation = 1;
		for (i=0; i<ntris; i++)
			tlist[i]->normal *= -1.0;
	}
}


void Polyhedron::corner_sort()
{
	//std::vector <int> max_index;
	//std::vector <int> min_index;
	unsigned int * min_index, * max_index, * corner_index;
	min_index = (unsigned int*)malloc(sizeof(unsigned int) * (ncorners));
	max_index = (unsigned int*)malloc(sizeof(unsigned int) * (ncorners));
	corner_index = (unsigned int*)malloc(sizeof(unsigned int) * (ncorners));

	for (int i = 0; i < ncorners; i++)
	{
		Corner* corner = clist[i];
		int minindex = std::min(corner->p->v->index, corner->n->v->index);
		min_index[i] = minindex;
		int maxindex = std::max(corner->p->v->index, corner->n->v->index);
		max_index[i] = maxindex;
		corner_index[i] = i;
	}

	//for (int k = 0; k < ncorners; k++)
	//{
	//	for (int j = 0; j < ncorners; j++)
	//	{
	//		if (j != k)
	//		{ 
	//			if (max_index[k] == max_index[j] && min_index[k] == min_index[j])
	//			{
	//				clist[j]->o = clist[k];
	//				clist[k]->o = clist[j];
	//			}
	//		}
	//	}
	//}
	tempA = (unsigned int*)malloc(sizeof(unsigned int) * (ncorners + 1));
	tempB = (unsigned int*)malloc(sizeof(unsigned int) * (ncorners + 1));
	tempC = (unsigned int*)malloc(sizeof(unsigned int) * (ncorners + 1));
	sort(min_index, max_index, corner_index, 0 ,ncorners);
	free(tempA);
	tempA = NULL;
	free(tempB);
	tempB = NULL;
	free(tempC);
	tempC = NULL;
	
	for (int i = 0; i < ncorners-1; i += 2)
	{
		//std::cout << min_index[i] << " " << max_index[i] << " " << corner_index[i] <<std::endl;
		clist[corner_index[i]]->o = clist[corner_index[i+1]];
		clist[corner_index[i+1]]->o = clist[corner_index[i]];
	}



}

void Polyhedron::sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid){
  unsigned int i;
	unsigned int current1, current2, current0;

  if (sid>=eid)
		return;

	sort(A, B, C, sid, (sid+eid)/2);
	sort(A, B, C, (sid+eid)/2+1, eid);
	
	for (i=0; i<eid-sid+1; i++){
		tempA[i] = A[i+sid];
		tempB[i] = B[i+sid];
		tempC[i] = C[i+sid];
	}
	current1 = sid;
	current2 = (sid+eid)/2+1;
	current0 = sid;
	while ((current1<=(sid+eid)/2) && (current2<=eid)){
		if (tempA[current1-sid] < tempA[current2-sid]) {
			A[current0] = tempA[current1-sid];
			B[current0] = tempB[current1-sid];
			C[current0] = tempC[current1-sid];
			current1++;		
		}
		else if (tempA[current1-sid] > tempA[current2-sid]){
			A[current0] = tempA[current2-sid];
			B[current0] = tempB[current2-sid];
			C[current0] = tempC[current2-sid];
			current2++;		
		}
		else {
			if (tempB[current1-sid] < tempB[current2-sid]) {
				A[current0] = tempA[current1-sid];
				B[current0] = tempB[current1-sid];
				C[current0] = tempC[current1-sid];
				current1++;		
			} else {
				A[current0] = tempA[current2-sid];
				B[current0] = tempB[current2-sid];
				C[current0] = tempC[current2-sid];
				current2++;		
			}
		}
		current0++;
	}
	if (current1<=(sid+eid)/2){
		for (i=current1; i<=(sid+eid)/2; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}
	if (current2<=eid){
		for (i=current2; i<=eid; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}

}

void inittexture(void)
{
	int width, height;
	unsigned char* TextureArray0 = BmpToTexture("sample_bmp/221.bmp", &width, &height); //163 //checkerboard
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &Tex0);

	glBindTexture(GL_TEXTURE_2D, Tex0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexImage2D(GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureArray0);
}

void init(void) {
  /* select clearing color */ 

  glClearColor (0.0, 0.0, 0.0, 0.0);  // background
  glShadeModel (GL_FLAT);
  glPolygonMode(GL_FRONT, GL_FILL);

  glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
	// may need it
  glPixelStorei(GL_PACK_ALIGNMENT,1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0) 
		glFrontFace(GL_CW);
	else 
		glFrontFace(GL_CCW);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

  /* set escape key to exit */
  switch (key) {
    case 27:
			poly->finalize();  // finalize_everything
      exit(0);
      break;

		case '0':
			display_mode = 0;
			display();
			break;


		case '2':
			display_mode = 0;
			rpoly = new Polyhedron;
			rpoly->regular_sub(poly);
			poly = rpoly;
			poly->initialize(); // initialize everything

			poly->calc_bounding_sphere();
			poly->calc_face_normals_and_area();
			poly->average_normals();
			poly->setheat();
			poly->GaussCurvature();
			poly->MeanCurvature();
			poly->PrincipalCurvature();
			display();
			//poly->finalize();
			break;

		case '3':
			display_mode = 0;
			irpoly = new Polyhedron;
			irpoly->irregular_sub(poly);
			poly = irpoly;
			poly->initialize(); // initialize everything

			poly->calc_bounding_sphere();
			poly->calc_face_normals_and_area();
			poly->average_normals();
			poly->setheat();
			poly->GaussCurvature();
			poly->MeanCurvature();
			poly->PrincipalCurvature();
			display();
			//poly->finalize();
			break;

		case '4':
			display_mode = 0;
			

			//std::clock_t start;
			//double duration;
			//start = std::clock();
			for (int i = 0; i < 5; i++)
			{
				uniformsmoothpoly = new Polyhedron;
				uniformsmoothpoly->smooth(poly);

				//-------------------------for the unchange weight part--------------------------------
				for (int j = 0; j < poly->nverts; j++)
					uniformsmoothpoly->vlist[j]->neighborweight = poly->vlist[j]->neighborweight;
				//-------------------------------------------------------------------------------------
				poly->finalize();
				poly = uniformsmoothpoly;
				poly->initialize(); // initialize everything			
				poly->calc_bounding_sphere();
				poly->calc_face_normals_and_area();
				poly->average_normals();
				poly->setheat();
				poly->GaussCurvature();
				poly->MeanCurvature();
				poly->PrincipalCurvature();
			
			}
			//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			//std::cout << "\n*********************************************" << std::endl;
			//std::cout << "Time costs for 10 iterations: " << duration << "\n";
			//std::cout << "*********************************************" << std::endl;
			

		
			display();
			//
			
			break;

		case '5':
			display_mode = 5;

			std::clock_t start;
			double duration;
			start = std::clock();

			heat2litertime += 10000.;
			poly->heat2(poly, poly->maxindex, poly->minindex, heat2litertime);
			poly->color_t_rainbow();

			//for (int j = 0; j < 100; j++)
			//{
			//	//poly->heat1(poly, poly->maxindex, poly->minindex);
			//	poly->heat3(poly, poly->maxindex, poly->minindex);
			//	poly->color_t_rainbow();
			//}
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "\n*********************************************" << std::endl;
			std::cout << "Time costs for 1000 heat iterations: " << duration << "\n";
			std::cout << "*********************************************" << std::endl;

			display();
			break;

		case 'p':
			//heat2litertime += 10000.;
			//poly->TensorSmoothing2(heat2litertime);

			start = std::clock();

			poly->GetTensor3D();
			for (int j = 0; j < 200; j++)
			{
				//poly->TensorSmoothing1(poly, poly->maxindex, poly->minindex);
				poly->TensorSmoothing1();
			}
			poly->GoBackTensor2D();

			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << "\n*********************************************" << std::endl;
			std::cout << "Time costs for 5 iterations: " << duration << "\n";
			std::cout << "*********************************************" << std::endl;

			poly->SettingStreamlineFold();
			//poly->SettingStreamlineUnfold();

			display();
			break;

		case '6':
			display_mode = 6;
			display();
			break;

		case '7':
			display_mode = 7;
			display();
			break;

		case '8':
			display_mode = 8;
			textureflag = 1;
			for (int j = 0; j < 10; j++)
			{
				poly->texmap(poly);
			}
			display();
			break;

		case '9':
			display_mode = 9;
			silhouetteflag = -silhouetteflag;
			display();
			break;

		case '1':
			display_mode = 1;
			display();
			break;

		case 's':
			sflag = -sflag;
			display();
			break;

		case 't':
			tflag = -tflag;
			display();
			break;

		case 'x':
			switch(ACSIZE){
			case 1:
				ACSIZE = 16;
				break;

			case 16:
				ACSIZE = 1;
				break;

			default:
				ACSIZE = 1;
				break;
			}
			fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
			display();
			break;

		case '|':
			this_file = fopen("rotmat.txt", "w");
			for (i=0; i<4; i++) 
				fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
			fclose(this_file);
			break;

		case '^':
			this_file = fopen("rotmat.txt", "r");
			for (i=0; i<4; i++) 
				fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
			fclose(this_file);
			display();
			break;

	}
}

Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex *[max_verts];
	tlist = new Triangle *[max_tris];		
}


void multmatrix(const Matrix m)
{ 
  int i,j, index = 0;

  GLfloat mat[16];

  for ( i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mat[index++] = m[i][j];

  glMultMatrixf (mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


  glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix( rotmat );
	//
	glScalef(0.4, 0.4, 0.4);
	//glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	//glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
	
	
	//glRotatef(-90, 1, 0, 0);
	//glRotatef(-45, 0, 0, 1);
	
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch(mouse_mode){
	case -1:

		xsize = (float) win_width;
		ysize = (float) win_height;
	
		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat( rotmat, rvec );
		trackball( r, s_old, t_old, s, t );
		add_quats( r, rvec, rvec );
		quat_to_mat( rvec, rotmat );

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth=1.0e+20, current_depth;
	int seed_id=-1; 
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *) buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;
		
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	tridepth = smallest_depth;
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON) {
		switch(mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float) win_width;
				float ysize = (float) win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	} 
	else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
		  GLint hits;
		  GLint viewport[4];

		  glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
		  (void) glRenderMode(GL_SELECT);

		  glInitNames();
		  glPushName(0);

		  glMatrixMode(GL_PROJECTION);
	    glPushMatrix();
			glLoadIdentity();
/*  create 5x5 pixel picking region near cursor location */
	    gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3] - y),
                 1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix ();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
	    glPopMatrix();
		  glFlush();

	    hits = glRenderMode(GL_RENDER);
		  poly->seed = processHits(hits, selectBuf);
		  
//-------------------choose the index for the max/min temperature-----------------//
		if (textureflag == -1)
		{ 
		  if (maxflag == -1)
			  poly->maxindex = poly->seed;
		  else 
			  poly->minindex = poly->seed;
		  maxflag = -maxflag;
		}
		else
		{
			if (poly->texcoordindex.size() == 4) poly->texcoordindex.clear();
			poly->texcoordindex.push_back(poly->seed);
		}		
//--------------------------------------------------------------------------------//
		
//-------------------choose the index for the anchorpoints-----------------//
		if (chooseflag == -1)
		{
			poly->m_fixedindex = poly->tlist[poly->seed]->verts[0]->index;
			chooseflag = -chooseflag;
		}
		else
		{
			poly->m_editedintex = poly->tlist[poly->seed]->verts[0]->index;
			chooseflag = -chooseflag;
		}
//--------------------------------------------------------------------------------//
			display();
		}
	}

	else if (button == GLUT_RIGHT_BUTTON)
	{
		if (state == GLUT_DOWN) {
			//poly->RIDeforming(poly);
			//for (int i = 0; i < poly->nverts; i++)
			//{
			//	poly->vlist[i]->x = poly->vlist[i]->Newpoint.x;
			//	poly->vlist[i]->y = poly->vlist[i]->Newpoint.y;
			//	poly->vlist[i]->z = poly->vlist[i]->Newpoint.z;
			//}
			poly->m_iAroundNumber = 1;
			poly->m_iHeightNumber = 3;
			poly->InitLpls();
			poly->InitAxesAndMatrixSelf();
			poly->InitNewAxesMatrix();
			poly->CalculateNewAxes();
			poly->CalculateNewMatrixLpls();
			display();
		}
	}
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i=0; i<poly->ntris; i++) {
		Triangle *temp_t=poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};
		
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   
		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j=0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

void Polyhedron::set_old_vert_cord_oldweight()
{

	for (int i = 0; i < nverts; i++)
	{
		double n = vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (vlist[i] == vlist[i]->tris[j]->verts[k])
				{
					sum_eij += 1.0 / vlist[i]->tris[j]->edges[k]->length;
					break;
				}
			}
		}
		
		for (int j = 0; j < n; j++)
		{
			double wij = 0.0;
			for (int k = 0; k < 3; k++)
			{
				if (vlist[i] == vlist[i]->tris[j]->verts[k])
				{
					double eij = 1.0 / vlist[i]->tris[j]->edges[k]->length;
					wij = eij / sum_eij;
					vlist[i]->neighborweight.push_back(wij);
					break;
				}
			}
		}
	}
}


void Polyhedron::set_old_vert_meancurvature_oldweight()
{
	for (int i = 0; i < nverts; i++)
	{
		double dt = 0.1;
		double n =vlist[i]->ntris;
		double big = 0.5;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (vlist[i] == vlist[i]->tris[j]->verts[k])
				{
					for (int p = 0; p < 3; p++)
					{
						if (vlist[i]->tris[j]->corners[p]->e->index == vlist[i]->tris[j]->edges[k]->index)
						{
							c = vlist[i]->tris[j]->corners[p];
							break;
						}
					}

					if (sin(c->c) == 0 && sin(c->o->c) == 0)
						sum_eij += big;
					else if (sin(c->c) == 0 && sin(c->o->c) != 0)
						sum_eij += (big + cos(c->o->c) / sin(c->o->c)) / 2.0;
					else if (sin(c->c) != 0 && sin(c->o->c) == 0)
						sum_eij += (cos(c->c) / sin(c->c) + big) / 2.0;
					else
						sum_eij += (cos(c->c) / sin(c->c) + cos(c->o->c) / sin(c->o->c)) / 2.0;
					break;
				}
			}
		}

		for (int j = 0; j < n; j++)
		{
			Corner* c;
			double wij = 0.0;
			for (int k = 0; k < 3; k++)
			{
				if (vlist[i] == vlist[i]->tris[j]->verts[k])
				{
					for (int p = 0; p < 3; p++)
					{
						if (vlist[i]->tris[j]->corners[p]->e->index == vlist[i]->tris[j]->edges[k]->index)
						{
							c = vlist[i]->tris[j]->corners[p];
							break;
						}
					}

					double eij;
					if (sin(c->c) == 0 && sin(c->o->c) == 0)
						eij = big;
					else if (sin(c->c) == 0 && sin(c->o->c) != 0)
						eij = (big + cos(c->o->c) / sin(c->o->c)) / 2.0;
					else if (sin(c->c) != 0 && sin(c->o->c) == 0)
						eij = (cos(c->c) / sin(c->c) + big) / 2.0;
					else
						eij = (cos(c->c) / sin(c->c) + cos(c->o->c) / sin(c->o->c)) / 2.0;

					wij = eij / sum_eij;
					vlist[i]->neighborweight.push_back(wij);
					break;
				}
			}
		}
	}
}


void Polyhedron::set_old_vert_meanvalue_oldweight()
{
	for (int i = 0; i < nverts; i++)
	{
		double n = vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (vlist[i] == vlist[i]->tris[j]->verts[k])
					//if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					for (int p = 0; p < 3; p++)
					{
						if (vlist[i]->tris[j]->corners[p]->e->index == vlist[i]->tris[j]->edges[k]->index)
						{
							c = vlist[i]->tris[j]->corners[p];
							break;
						}
					}
					sum_eij += (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
					//c = poly->vlist[i]->tris[j]->corners[k];
					//sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;

				}
			}
		}

		for (int j = 0; j < n; j++)
		{
			Corner* c;
			double wij = 0.0;
			for (int k = 0; k < 3; k++)
			{
				if (vlist[i] == vlist[i]->tris[j]->verts[k])
					//if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					for (int p = 0; p < 3; p++)
					{
						if (vlist[i]->tris[j]->corners[p]->e->index == vlist[i]->tris[j]->edges[k]->index)
						{
							c = vlist[i]->tris[j]->corners[p];
							break;
						}
					}
					double eij = (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
					//c = poly->vlist[i]->tris[j]->corners[k];
					//double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
					wij = eij / sum_eij;
					//p_sumx += wij * (c->n->v->x - poly->vlist[i]->x);
					//p_sumy += wij * (c->n->v->y - poly->vlist[i]->y);
					//p_sumz += wij * (c->n->v->z - poly->vlist[i]->z);
					vlist[i]->neighborweight.push_back(wij);
				}
			}
		}
	}
}

void Polyhedron::update_old_vert_all_old(Polyhedron* poly)
{
	double dt = 0.9;
	for (int i = 0; i < poly->nverts; i++)
	{
		double n = poly->vlist[i]->ntris;
		double p_sumx = 0.0;
		double p_sumy = 0.0;
		double p_sumz = 0.0;
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
				{
					p_sumx += poly->vlist[i]->neighborweight[j] * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->x - poly->vlist[i]->x);
					p_sumy += poly->vlist[i]->neighborweight[j] * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->y - poly->vlist[i]->y);
					p_sumz += poly->vlist[i]->neighborweight[j] * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->z - poly->vlist[i]->z);
					break;
				}
			}
		}
		poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
		poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
		poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
	}
}
//===============================================================================================
//===============================================================================================

void cut_edge(Edge* e)
{
	Vertex* new_vert = e->verts[0];
}

void Polyhedron::update_old_vert(Polyhedron* poly)  // upgrade the location of the old vertices
{
	for (int i = 0; i < poly->nverts; i++)
	{
		double n = poly->vlist[i]->ntris;
		double beta = 1.0 / n * (0.625 - pow(0.375 + 0.25 * cos(2.0 * PI / n), 2.0));
		double p_sumx = 0;
		double p_sumy = 0;
		double p_sumz = 0;
		for (int j = 0; j < n; j++)
		{
			if (poly->vlist[i] != poly->vlist[i]->tris[j]->verts[0])
			{
				p_sumx += poly->vlist[i]->tris[j]->verts[0]->x;
				p_sumy += poly->vlist[i]->tris[j]->verts[0]->y;
				p_sumz += poly->vlist[i]->tris[j]->verts[0]->z;
			}
			if (poly->vlist[i] != poly->vlist[i]->tris[j]->verts[1])
			{
				p_sumx += poly->vlist[i]->tris[j]->verts[1]->x;
				p_sumy += poly->vlist[i]->tris[j]->verts[1]->y;
				p_sumz += poly->vlist[i]->tris[j]->verts[1]->z;
			}
			if (poly->vlist[i] != poly->vlist[i]->tris[j]->verts[2])
			{
				p_sumx += poly->vlist[i]->tris[j]->verts[2]->x;
				p_sumy += poly->vlist[i]->tris[j]->verts[2]->y;
				p_sumz += poly->vlist[i]->tris[j]->verts[2]->z;
			}
		}
		poly->vlist[i]->update_vert_pos[0] = (1.0 - n * beta) * poly->vlist[i]->x + beta * p_sumx;
		poly->vlist[i]->update_vert_pos[1] = (1.0 - n * beta) * poly->vlist[i]->y + beta * p_sumy;
		poly->vlist[i]->update_vert_pos[2] = (1.0 - n * beta) * poly->vlist[i]->z + beta * p_sumz;
	}

}


void Polyhedron::update_old_vert_uniform(Polyhedron* poly)
{

		for (int i = 0; i < poly->nverts; i++)
		{
			double dt = 0.9;
			double n = poly->vlist[i]->ntris;
			double eij = 1.0;
			double sum_eij = n;
			double wij = eij / sum_eij;
			double p_sumx = 0.0;
			double p_sumy = 0.0;
			double p_sumz = 0.0;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
					{
						p_sumx += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->x - poly->vlist[i]->x);
						p_sumy += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->y - poly->vlist[i]->y);
						p_sumz += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->z - poly->vlist[i]->z);
					}
				}
			}
			poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
			poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
			poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
		}

}


void Polyhedron::update_old_vert_cord(Polyhedron* poly)
{

		for (int i = 0; i < poly->nverts; i++)
		{
			double dt = 0.9;
			double n = poly->vlist[i]->ntris;
			double sum_eij = 0.0;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
					{
						sum_eij += 1.0 / poly->vlist[i]->tris[j]->edges[k]->length;
						break;
					}
				}
			}
			double p_sumx = 0.0;
			double p_sumy = 0.0;
			double p_sumz = 0.0;
			for (int j = 0; j < n; j++)
			{
				double wij = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
					{
						double eij = 1.0 / poly->vlist[i]->tris[j]->edges[k]->length;
						wij = eij / sum_eij;
						p_sumx += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->x - poly->vlist[i]->x);
						p_sumy += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->y - poly->vlist[i]->y);
						p_sumz += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->z - poly->vlist[i]->z);
						break;
					}
				}
			}
			poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
			poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
			poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
		}

}

void Polyhedron::update_old_vert_meancurvature(Polyhedron* poly)
{

		for (int i = 0; i < poly->nverts; i++)
		{
			double dt = 0.1;
			double n = poly->vlist[i]->ntris;
			double big = 1.;
			double sum_eij = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
					{
						c = poly->vlist[i]->tris[j]->corners[k];
						
						if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) < 0.000005)
							sum_eij += big;
						else if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) >= 0.000005)
							sum_eij += (big + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
						else if (sin(c->p->c) >= 0.000005 && sin(c->p->o->c) < 0.000005)
							sum_eij += (cos(c->p->c) / sin(c->p->c) + big) / 2.0;
						else 
							sum_eij += (cos(c->p->c) / sin(c->p->c) + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
						break;
					}
				}
			}
		
			double p_sumx = 0.0;
			double p_sumy = 0.0;
			double p_sumz = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				double wij = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
					{
						c = poly->vlist[i]->tris[j]->corners[k];
						double eij;
						if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) < 0.000005)
							eij = big;
						else if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) >= 0.000005)
							eij = (big + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
						else if (sin(c->p->c) >= 0.000005 && sin(c->p->o->c) < 0.000005)
							eij = (cos(c->p->c) / sin(c->p->c) + big) / 2.0;
						else
							eij = (cos(c->p->c) / sin(c->p->c) + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
						
						wij = eij / sum_eij;
						p_sumx += wij * (c->n->v->x - poly->vlist[i]->x);
						p_sumy += wij * (c->n->v->y - poly->vlist[i]->y);
						p_sumz += wij * (c->n->v->z - poly->vlist[i]->z);
						break;
					}
				}
			}
			poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
			poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
			poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
		}
	//for (int i = 0; i < poly->nverts; i++)
	//{
	//	double dt = 0.1;
	//	double n = poly->vlist[i]->ntris;
	//	double big = 0.5;
	//	double sum_eij = 0.0;
	//	for (int j = 0; j < n; j++)
	//	{
	//		Corner* c;
	//		for (int k = 0; k < 3; k++)
	//		{
	//			if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
	//			{
	//				for (int p = 0; p < 3; p++)
	//				{
	//					if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
	//					{
	//						c = poly->vlist[i]->tris[j]->corners[p];
	//						break;
	//					}
	//				}
	//
	//				if (sin(c->c) == 0 && sin(c->o->c) == 0)
	//					sum_eij += big;
	//				else if (sin(c->c) == 0 && sin(c->o->c) != 0)
	//					sum_eij += (big + cos(c->o->c) / sin(c->o->c)) / 2.0;
	//				else if (sin(c->c) != 0 && sin(c->o->c) == 0)
	//					sum_eij += (cos(c->c) / sin(c->c) + big) / 2.0;
	//				else
	//					sum_eij += (cos(c->c) / sin(c->c) + cos(c->o->c) / sin(c->o->c)) / 2.0;
	//				break;
	//			}
	//		}
	//	}
	//
	//	double p_sumx = 0.0;
	//	double p_sumy = 0.0;
	//	double p_sumz = 0.0;
	//	for (int j = 0; j < n; j++)
	//	{
	//		Corner* c;
	//		double wij = 0.0;
	//		for (int k = 0; k < 3; k++)
	//		{
	//			if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
	//			{
	//				for (int p = 0; p < 3; p++)
	//				{
	//					if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
	//					{
	//						c = poly->vlist[i]->tris[j]->corners[p];
	//						break;
	//					}
	//				}
	//
	//				double eij;
	//				if (sin(c->c) == 0 && sin(c->o->c) == 0)
	//					eij = big;
	//				else if (sin(c->c) == 0 && sin(c->o->c) != 0)
	//					eij = (big + cos(c->o->c) / sin(c->o->c)) / 2.0;
	//				else if (sin(c->c) != 0 && sin(c->o->c) == 0)
	//					eij = (cos(c->c) / sin(c->c) + big) / 2.0;
	//				else
	//					eij = (cos(c->c) / sin(c->c) + cos(c->o->c) / sin(c->o->c)) / 2.0;
	//
	//				wij = eij / sum_eij;
	//				p_sumx += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->x - poly->vlist[i]->x);
	//				p_sumy += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->y - poly->vlist[i]->y);
	//				p_sumz += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->z - poly->vlist[i]->z);
	//				break;
	//			}
	//		}
	//	}
	//	poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
	//	poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
	//	poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
	//}


}
void Polyhedron::update_old_vert_meanvalue(Polyhedron* poly)
{
	
		//for (int i = 0; i < poly->nverts; i++)
		//{
		//	double dt = 0.9;
		//	double n = poly->vlist[i]->ntris;
		//	double sum_eij = 0.0;
		//	for (int j = 0; j < n; j++)
		//	{
		//		Corner* c;
		//		for (int k = 0; k < 3; k++)
		//		{
		//			if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
		//			//if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
		//			{
		//				for (int p = 0; p < 3; p++)
		//				{
		//					if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
		//					{
		//						c = poly->vlist[i]->tris[j]->corners[p];
		//						break;
		//					}
		//				}
		//				sum_eij += (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
		//				//c = poly->vlist[i]->tris[j]->corners[k];
		//				//sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
		//				
		//			}
		//		}
		//	}
		//
		//	double p_sumx = 0.0;
		//	double p_sumy = 0.0;
		//	double p_sumz = 0.0;
		//	for (int j = 0; j < n; j++)
		//	{
		//		Corner* c;
		//		double wij = 0.0;
		//		for (int k = 0; k < 3; k++)
		//		{
		//			if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
		//			//if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
		//			{
		//				for (int p = 0; p < 3; p++)
		//				{
		//					if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
		//					{
		//						c = poly->vlist[i]->tris[j]->corners[p];
		//						break;
		//					}
		//				}
		//				double eij = (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
		//				//c = poly->vlist[i]->tris[j]->corners[k];
		//				//double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
		//				wij = eij / sum_eij;
		//				//p_sumx += wij * (c->n->v->x - poly->vlist[i]->x);
		//				//p_sumy += wij * (c->n->v->y - poly->vlist[i]->y);
		//				//p_sumz += wij * (c->n->v->z - poly->vlist[i]->z);
		//				p_sumx += wij * (poly->vlist[i]->tris[j]->verts[(k+1)%3]->x - poly->vlist[i]->x);
		//				p_sumy += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->y - poly->vlist[i]->y);
		//				p_sumz += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->z - poly->vlist[i]->z);
		//				
		//			}
		//		}
		//	}
		//	poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
		//	poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
		//	poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
		//}

	for (int i = 0; i < poly->nverts; i++)
	{
		double dt = 0.9;
		double n = poly->vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{

					c = poly->vlist[i]->tris[j]->corners[k];
					sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;

				}
			}
		}

		double p_sumx = 0.0;
		double p_sumy = 0.0;
		double p_sumz = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			double wij = 0.0;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{

					c = poly->vlist[i]->tris[j]->corners[k];
					double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
					wij = eij / sum_eij;
					p_sumx += wij * (c->n->v->x - poly->vlist[i]->x);
					p_sumy += wij * (c->n->v->y - poly->vlist[i]->y);
					p_sumz += wij * (c->n->v->z - poly->vlist[i]->z);
					
				}
			}
		}
		poly->vlist[i]->update_vert_pos[0] = poly->vlist[i]->x + dt * p_sumx;
		poly->vlist[i]->update_vert_pos[1] = poly->vlist[i]->y + dt * p_sumy;
		poly->vlist[i]->update_vert_pos[2] = poly->vlist[i]->z + dt * p_sumz;
	}


}

void Polyhedron::smooth1(Polyhedron* poly)
{
	std::clock_t start;
	double duration;
	start = std::clock();
	//update_old_vert_uniform(poly);
	//update_old_vert_cord(poly);
	update_old_vert_meancurvature(poly);    //1. create new points need old vert; 2. update old vert then
	//update_old_vert_meanvalue(poly);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "\n*********************************************" << std::endl;
	std::cout << "Time costs for 1 iteration: " << duration << "\n";
	std::cout << "*********************************************" << std::endl;
	Vertex** clean_vlist = new Vertex * [(int)(poly->max_verts)];
	for (int i = 0; i < poly->nverts; i++)
	{
		clean_vlist[i] = new Vertex(poly->vlist[i]->update_vert_pos[0], poly->vlist[i]->update_vert_pos[1], poly->vlist[i]->update_vert_pos[2]);
		clean_vlist[i]->index = poly->vlist[i]->index;
	}

	// create new vetex list
	for (int i = 0; i < poly->max_verts; i++)
	{
		poly->vlist[i]->x = clean_vlist[i]->x;
		poly->vlist[i]->y = clean_vlist[i]->y;
		poly->vlist[i]->z = clean_vlist[i]->z;
	}

}

void Polyhedron::smooth(Polyhedron* poly)
{
	max_tris = (int)(poly->max_tris);  /* leave some room for expansion */
	tlist = new Triangle * [max_tris];
	ntris = 0;//= new_poly->max_tris;

	max_verts = (int)(poly->max_verts);
	vlist = new Vertex * [max_verts];
	nverts = 0;// = new_poly->max_verts;

	//std::clock_t start;
	//double duration;
	//start = std::clock();
	
	//update_old_vert_uniform(poly);
	//update_old_vert_cord(poly);
	//update_old_vert_meancurvature(poly);    //1. create new points need old vert; 2. update old vert then
	update_old_vert_meanvalue(poly);
	//update_old_vert_all_old(poly);
	
	//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//std::cout << "\n*********************************************" << std::endl;
	//std::cout << "Time costs for 1 iteration: " << duration << "\n";
	//std::cout << "*********************************************" << std::endl;
	Vertex** clean_vlist = new Vertex * [(int)(poly->max_verts)];
	for (int i = 0; i < poly->nverts; i++)
	{
		clean_vlist[i] = new Vertex(poly->vlist[i]->update_vert_pos[0], poly->vlist[i]->update_vert_pos[1], poly->vlist[i]->update_vert_pos[2]);
		clean_vlist[i]->index = poly->vlist[i]->index;
		clean_vlist[i]->other_props = NULL;
	}

	// create new vetex list
	for (int i = 0; i < poly->max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
		vlist[i]->other_props = NULL;
		vlist[i] = clean_vlist[i];
		nverts += 1;
	}

	// create new triangle list
	for (int i = 0; i < poly->ntris; i++)    // i is the index of triangle
	{
		//Triangle* f = poly->tlist[i];

		//push new triangle with 3 central points
		tlist[i] = new Triangle;
		tlist[i]->nverts = 3;
		tlist[i]->other_props = NULL;


		//for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
		//{
		//	tlist[i * 4]->verts[j] = poly->tlist[i]->edges[j]->new_vert;
		//}
		//ntris += 1;
		bool find = false;
		for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
		{
			for (int p = 0; p < poly->max_verts; p++)
			{
				if (poly->tlist[i]->verts[j]->index == vlist[p]->index)
				{
					tlist[i]->verts[j] = vlist[p];
					break;
				}
			}
		}
		ntris += 1;

	}

	/* get rid of triangles that use the same vertex more than once */

	//for (int i = ntris - 1; i >= 0; i--) {
	//
	//	Triangle* tri = tlist[i];
	//	Vertex* v0 = tri->verts[0];
	//	Vertex* v1 = tri->verts[1];
	//	Vertex* v2 = tri->verts[2];
	//
	//	if (v0 == v1 || v1 == v2 || v2 == v0) {
	//		free(tlist[i]);
	//		ntris--;
	//		tlist[i] = tlist[ntris];
	//	}
	//}
}

//void Polyhedron::TensorSmoothing1()   //Explicit
//{
//
//	for (int i = 0; i < nverts; i++)
//	{
//		vlist[i]->TensorSave[0] = vlist[i]->Tensor[0];
//		vlist[i]->TensorSave[1] = vlist[i]->Tensor[1];
//		vlist[i]->TensorSave[2] = vlist[i]->Tensor[2];
//	}
//
//	//-----------------------Mean value coordinates-------------------------//
//	for (int i = 0; i < nverts; i++)
//	{
//
//		double dt = 1;
//		int n = vlist[i]->ntris;
//		double sum_eij = 0.0;
//		for (int j = 0; j < n; j++)
//		{
//			Corner* c;
//			for (int k = 0; k < 3; k++)
//			{
//				if (vlist[i] == vlist[i]->tris[j]->corners[k]->v)
//				{
//					c = vlist[i]->tris[j]->corners[k];
//					sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
//					break;
//				}
//			}
//		}
//
//		//if (vlist[i] != tlist[2]->verts[0] && vlist[i] != tlist[17]->verts[0])
//		//{
//			double p_sum0 = 0.0;
//			double p_sum1 = 0.0;
//			double p_sum2 = 0.0;
//
//			for (int j = 0; j < n; j++)
//			{
//				Corner* c;
//				double wij = 0.0;
//				for (int k = 0; k < 3; k++)
//				{
//					if (vlist[i] == vlist[i]->tris[j]->corners[k]->v)
//					{
//						c = vlist[i]->tris[j]->corners[k];
//						double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
//						wij = eij / sum_eij;//1. / vlist[i]->ntris;
//
//						glm::mat2 Mj = glm::mat2(c->n->v->TensorSave[0], c->n->v->TensorSave[1], c->n->v->TensorSave[1], c->n->v->TensorSave[2]);
//						glm::mat2 Mi = glm::mat2(vlist[i]->TensorSave[0], vlist[i]->TensorSave[1], vlist[i]->TensorSave[1], vlist[i]->TensorSave[2]);
//
//						Mj = Mi * Mj * inverse(Mi);
//
//						p_sum0 += wij * (Mj[0][0] - vlist[i]->TensorSave[0]);
//						p_sum1 += wij * (Mj[0][1] - vlist[i]->TensorSave[1]);
//						p_sum2 += wij * (Mj[1][1] - vlist[i]->TensorSave[2]);
//						break;
//					}
//				}
//			}
//			vlist[i]->Tensor[0] = vlist[i]->TensorSave[0] + dt * p_sum0;
//			vlist[i]->Tensor[1] = vlist[i]->TensorSave[1] + dt * p_sum1;
//			vlist[i]->Tensor[2] = vlist[i]->TensorSave[2] + dt * p_sum2;
//
//			CalculatePrincipalCurvature(i);
//		//}
//	}
//
//}

//void Polyhedron::TensorSmoothing1()   //Explicit - Global
//{
//	for (int i = 0; i < nverts; i++)
//	{
//		vlist[i]->TglobalSave = vlist[i]->Tglobal;
//	}
//	//-----------------------Mean value coordinates-------------------------//
//	for (int i = 0; i < nverts; i++)
//	{	
//
//		double dt = 1;
//		int n = vlist[i]->ntris;
//		double sum_eij = 0.0;
//		for (int j = 0; j < n; j++)
//		{
//			Corner* c;
//			for (int k = 0; k < 3; k++)
//			{
//				if (vlist[i] == vlist[i]->tris[j]->corners[k]->v)
//				{
//					c = vlist[i]->tris[j]->corners[k];
//					sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
//					break;
//				}
//			}
//		}
//
//		//if (vlist[i] != tlist[2]->verts[0] && vlist[i] != tlist[3]->verts[0])
//		//{
//			glm::mat3 p_sum = glm::mat3(0);
//
//			for (int j = 0; j < n; j++)
//			{
//				Corner* c;
//				double wij = 0.0;
//				for (int k = 0; k < 3; k++)
//				{
//					if (vlist[i] == vlist[i]->tris[j]->corners[k]->v)
//					{
//						c = vlist[i]->tris[j]->corners[k];
//						double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
//						wij = eij / sum_eij;//1. / vlist[i]->ntris;
//
//						p_sum += (float)wij * (c->n->v->TglobalSave - vlist[i]->TglobalSave);
//						break;
//					}
//				}
//			}
//			vlist[i]->Tglobal = vlist[i]->TglobalSave + (float)dt * p_sum;
//
//		//}
//	}
//
//}






//void Polyhedron::GetFrameForV()
//{
//	for (int i; i < nverts; i++)
//	{
//		Corner* c;
//		for (int k = 0; k < 3; k++)
//		{
//			if (vlist[i]== vlist[i]->tris[0]->corners[k]->v)
//			{
//				c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
//				break;
//			}
//		}
//		glm::vec3 nowpoint = glm::vec3(vlist[i]->x, vlist[i]->y, vlist[i]->z);
//		glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
//		glm::vec3 n = normalize(glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
//		glm::vec3 vij = (nextpoint - nowpoint);
//		glm::vec3 u1 = normalize(vij - dot(u1, n) * n);
//		glm::vec3 u2 = normalize(cross(u1, n));
//		vlist[i]->DiscreteFrame[0] = u1;
//		vlist[i]->DiscreteFrame[1] = n;
//		vlist[i]->DiscreteFrame[2] = u2;
//	}
//}


void Polyhedron::UpdateDiscreteFrame(Vertex* vi)
{

	Corner* c;
	for (int k = 0; k < 3; k++)
	{
		if (vi == vi->tris[0]->corners[k]->v)
		{
			c = vi->tris[0]->corners[k];   // record the first triangle's corner of Vi
			break;
		}
	}
	glm::vec3 nowpoint = glm::vec3(vi->x, vi->y, vi->z);
	glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
	glm::vec3 n = normalize(glm::vec3(vi->normal.x, vi->normal.y, vi->normal.z));
	glm::vec3 vij = (nextpoint - nowpoint);
	glm::vec3 u1 = normalize(vij - dot(u1, n) * n);
	glm::vec3 u2 = normalize(cross(u1, n));
	vi->DiscreteFrame[0] = u1;
	vi->DiscreteFrame[1] = n;
	vi->DiscreteFrame[2] = u2;

	Eigen::MatrixXd A(9* vi->ntris+3, 9 * vi->ntris);
	Eigen::VectorXd b(9* vi->ntris+3);
	b.setZero();

	for (int j = 0; j < vi->ntris; j++)
	{
		// Firstly, get the frame of vj.
		glm::vec3 nowpoint = glm::vec3(vi->x, vi->y, vi->z);
		glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		glm::vec3 nj = normalize(glm::vec3(c->n->v->normal.x, c->n->v->normal.y, c->n->v->normal.z));
		glm::vec3 uj1 = normalize(nowpoint - nextpoint);
		glm::vec3 uj2 = normalize(cross(uj1, nj));

		// Then use vi's frame to express vj's frame and get the coefficient
		float R00 = dot(uj1, u1); float R01 = dot(uj1, u2); float R02 = dot(uj1, n);
		float R10 = dot(uj2, u1); float R11 = dot(uj2, u2); float R12 = dot(uj2, n);
		float R20 = dot(nj, u1);  float R21 = dot(nj, u2);  float R22 = dot(nj, n);

		b[j * 3]     = uj1.x;  A(j * 3, 0)     = R00;  A(j * 3, 1)     = R01;  A(j * 3, 2)     = R02;
		b[j * 3 + 1] = uj1.y;  A(j * 3 + 1, 3) = R00;  A(j * 3 + 1, 4) = R01;  A(j * 3 + 1, 5) = R02;
		b[j * 3 + 2] = uj1.z;  A(j * 3 + 2, 6) = R00;  A(j * 3 + 2, 7) = R01;  A(j * 3 + 2, 8) = R02;

		b[j * 3 + 3] = uj2.x;  A(j * 3 + 3, 0) = R10;  A(j * 3 + 3, 1) = R11;  A(j * 3 + 3, 2) = R12;
		b[j * 3 + 4] = uj2.y;  A(j * 3 + 4, 3) = R10;  A(j * 3 + 4, 4) = R11;  A(j * 3 + 4, 5) = R12;
		b[j * 3 + 5] = uj2.z;  A(j * 3 + 5, 6) = R10;  A(j * 3 + 5, 7) = R11;  A(j * 3 + 5, 8) = R12;
		
		b[j * 3 + 6] = nj.x;   A(j * 3 + 6, 0) = R20;  A(j * 3 + 6, 1) = R21;  A(j * 3 + 6, 2) = R22;
		b[j * 3 + 7] = nj.y;   A(j * 3 + 7, 3) = R20;  A(j * 3 + 7, 4) = R21;  A(j * 3 + 7, 5) = R22;
		b[j * 3 + 8] = nj.z;   A(j * 3 + 8, 6) = R20;  A(j * 3 + 8, 7) = R21;  A(j * 3 + 8, 8) = R22;

		c = c->p->o->p;
	}

	glm::vec3 newpoint = glm::vec3(s_old, t_old, tridepth);    // The moved point coordinates
	glm::vec3 newvector = newpoint - nowpoint;
	float newb1 = dot(newvector, u1); float newb2 = dot(newvector, u2); float newn = dot(newvector, n);

	if (vi == tlist[seed]->verts[0] )
	{
		b[9 * vi->ntris] = newvector.x;     A(9 * vi->ntris, 0) = newb1;      A(9 * vi->ntris, 1) = newb2;      A(9 * vi->ntris, 2) = newn;
		b[9 * vi->ntris + 1] = newvector.y;   A(9 * vi->ntris + 1, 3) = newb1;  A(9 * vi->ntris + 1, 4) = newb2;  A(9 * vi->ntris + 1, 5) = newn;
		b[9 * vi->ntris + 2] = newvector.z;   A(9 * vi->ntris + 2, 6) = newb1;  A(9 * vi->ntris + 2, 7) = newb2;  A(9 * vi->ntris + 2, 8) = newn;
	}
	else
	{
		b[9 * vi->ntris] = 0;     A(9 * vi->ntris, 0) = 0;      A(9 * vi->ntris, 1) = 0;      A(9 * vi->ntris, 2) = 0;
		b[9 * vi->ntris + 1] = 0;   A(9 * vi->ntris + 1, 3) = 0;  A(9 * vi->ntris + 1, 4) = 0;  A(9 * vi->ntris + 1, 5) = 0;
		b[9 * vi->ntris + 2] = 0;   A(9 * vi->ntris + 2, 6) = 0;  A(9 * vi->ntris + 2, 7) = 0;  A(9 * vi->ntris + 2, 8) = 0;
	}

	Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
	vi->DiscreteFrame[0] = glm::vec3(x[0], x[3], x[6]);//glm::vec3(x[0], x[1], x[2]);
	vi->DiscreteFrame[1] = glm::vec3(x[1], x[4], x[7]);//glm::vec3(x[3], x[4], x[5]);
	vi->DiscreteFrame[2] = glm::vec3(x[2], x[5], x[8]);//glm::vec3(x[6], x[7], x[8]);
}


void Polyhedron::RIDeforming(Polyhedron* poly)
{
	for (int i = 0; i < poly->nverts; i++)
	{
		UpdateDiscreteFrame(poly->vlist[i]);
	}

	int neibornumb = 0;
	for (int i = 0; i < poly->nverts; i++)     // counts how many pairs of verteices are 1- neibor
	{
		neibornumb += poly->vlist[i]->ntris;
	}

	SpMat Ax(neibornumb + 1, poly->nverts);    
	Eigen::VectorXd bx(neibornumb + 1);
	bx.setZero();
	std::vector<T> tripletListx; 

	int rownumberflag = 0;
	for (int i = 0; i < poly->nverts; i++)
	{
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
			{
				c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
				break;
			}
		}

		glm::vec3 nowpoint = glm::vec3(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		glm::vec3 n = normalize(glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
		glm::vec3 vij = (nextpoint - nowpoint);
		glm::vec3 u1 = normalize(vij - dot(u1, n) * n);
		glm::vec3 u2 = normalize(cross(u1, n));
		
		for (int j = 0; j < poly->vlist[i]->ntris; j++)
		{
			glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
			glm::vec3 ij = (nextpoint - nowpoint);
			float ct1 = dot(ij, u1); float ct2 = dot(ij, u2); float cn = dot(ij, n);

			int indexi = poly->vlist[i]->index;
			int indexj = c->n->v->index;
			tripletListx.push_back(T(rownumberflag, indexi, -1));
			tripletListx.push_back(T(rownumberflag, indexj, 1));
			bx[rownumberflag] = ct1 * vlist[i]->DiscreteFrame[0].x + ct2 * vlist[i]->DiscreteFrame[1].x + cn * vlist[i]->DiscreteFrame[2].x;
			rownumberflag += 1;
		}

	}
	glm::vec3 editpoint = glm::vec3(s_old, t_old, tridepth);
	glm::vec3 editedpoint = glm::vec3(poly->tlist[seed]->verts[0]->x, poly->tlist[seed]->verts[0]->y, poly->tlist[seed]->verts[0]->z);
	glm::vec3 editvector = editpoint - editedpoint;

	tripletListx.push_back(T(neibornumb, poly->tlist[seed]->verts[0]->index, 1));
	bx[neibornumb] = editvector.x;
	Ax.setFromTriplets(tripletListx.begin(), tripletListx.end());

	SpMat Ax_transpose = Ax.transpose();
	SpMat AxAx = Ax_transpose * Ax;

	Eigen::SimplicialCholesky<SpMat>MatricesCholeskyX(AxAx);

	Eigen::VectorXd rightsideX = Ax_transpose * bx;
	Eigen::VectorXd Xs = MatricesCholeskyX.solve(rightsideX);


	//-------------------------------------------------------------------------------
	//---------------------------yyyyyyyyyyyyyyyyyyyyyyyyy---------------------------
	SpMat Ay(neibornumb + 1, poly->nverts);
	Eigen::VectorXd by(neibornumb + 1);
	by.setZero();
	std::vector<T> tripletListy;

	rownumberflag = 0;
	for (int i = 0; i < poly->nverts; i++)
	{
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
			{
				c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
				break;
			}
		}

		glm::vec3 nowpoint = glm::vec3(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		glm::vec3 n = normalize(glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
		glm::vec3 vij = (nextpoint - nowpoint);
		glm::vec3 u1 = normalize(vij - dot(u1, n) * n);
		glm::vec3 u2 = normalize(cross(u1, n));

		for (int j = 0; j < poly->vlist[i]->ntris; j++)
		{
			glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
			glm::vec3 ij = (nextpoint - nowpoint);
			float ct1 = dot(ij, u1); float ct2 = dot(ij, u2); float cn = dot(ij, n);

			int indexi = poly->vlist[i]->index;
			int indexj = c->n->v->index;
			tripletListy.push_back(T(rownumberflag, indexi, -1));
			tripletListy.push_back(T(rownumberflag, indexj, 1));
			by[rownumberflag] = ct1 * vlist[i]->DiscreteFrame[0].y + ct2 * vlist[i]->DiscreteFrame[1].y + cn * vlist[i]->DiscreteFrame[2].y;
			rownumberflag += 1;
		}

	}

	tripletListy.push_back(T(neibornumb, poly->tlist[seed]->verts[0]->index, 1));
	by[neibornumb] = editvector.y;
	Ay.setFromTriplets(tripletListy.begin(), tripletListy.end());

	SpMat Ay_transpose = Ay.transpose();
	SpMat AyAy = Ay_transpose * Ay;

	Eigen::SimplicialCholesky<SpMat>MatricesCholeskyY(AyAy);

	Eigen::VectorXd rightsideY = Ay_transpose * by;
	Eigen::VectorXd Ys = MatricesCholeskyY.solve(rightsideY);


	//-------------------------------------------------------------------------------
	//---------------------------zzzzzzzzzzzzzzzzzzzzzzzzz---------------------------
	SpMat Az(neibornumb + 1, poly->nverts);
	Eigen::VectorXd bz(neibornumb + 1);
	bz.setZero();
	std::vector<T> tripletListz;

	rownumberflag = 0;
	for (int i = 0; i < poly->nverts; i++)
	{
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
			{
				c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
				break;
			}
		}

		glm::vec3 nowpoint = glm::vec3(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		glm::vec3 n = normalize(glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
		glm::vec3 vij = (nextpoint - nowpoint);
		glm::vec3 u1 = normalize(vij - dot(u1, n) * n);
		glm::vec3 u2 = normalize(cross(u1, n));

		for (int j = 0; j < poly->vlist[i]->ntris; j++)
		{
			glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
			glm::vec3 ij = (nextpoint - nowpoint);
			float ct1 = dot(ij, u1); float ct2 = dot(ij, u2); float cn = dot(ij, n);

			int indexi = poly->vlist[i]->index;
			int indexj = c->n->v->index;
			tripletListz.push_back(T(rownumberflag, indexi, -1));
			tripletListz.push_back(T(rownumberflag, indexj, 1));
			bz[rownumberflag] = ct1 * vlist[i]->DiscreteFrame[0].z + ct2 * vlist[i]->DiscreteFrame[1].z + cn * vlist[i]->DiscreteFrame[2].z;
			rownumberflag += 1;
		}

	}

	tripletListz.push_back(T(neibornumb, poly->tlist[seed]->verts[0]->index, 1));
	bz[neibornumb] = editvector.z;
	Az.setFromTriplets(tripletListz.begin(), tripletListz.end());

	SpMat Az_transpose = Az.transpose();
	SpMat AzAz = Az_transpose * Az;

	Eigen::SimplicialCholesky<SpMat>MatricesCholeskyZ(AzAz);

	Eigen::VectorXd rightsideZ = Az_transpose * bz;
	Eigen::VectorXd Zs = MatricesCholeskyZ.solve(rightsideZ);


	for (int i = 0; i < poly->nverts; i++)
	{
		poly->vlist[i]->Newpoint.x = Xs[i];
		poly->vlist[i]->Newpoint.y = Ys[i];
		poly->vlist[i]->Newpoint.z = Zs[i];
	}





	//for (int i = 0; i < poly->nverts; i++)
	//{
	//	Corner* c;
	//	for (int k = 0; k < 3; k++)
	//	{
	//		if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
	//		{
	//			c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
	//			break;
	//		}
	//	}

		//Eigen::MatrixXd Ax(poly->vlist[i]->ntris + 1, 1);
		//Eigen::VectorXd bx(poly->vlist[i]->ntris + 1);
		//bx.setZero();
		//Eigen::MatrixXd Ay(poly->vlist[i]->ntris + 1, 1);
		//Eigen::VectorXd by(poly->vlist[i]->ntris + 1);
		//by.setZero();
		//Eigen::MatrixXd Az(poly->vlist[i]->ntris + 1, 1);
		//Eigen::VectorXd bz(poly->vlist[i]->ntris + 1);
		//bz.setZero();
		//
		//glm::vec3 nowpoint = glm::vec3(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		//glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		//glm::vec3 n = normalize(glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
		//glm::vec3 vij = (nextpoint - nowpoint);
		//glm::vec3 u1 = normalize(vij - dot(u1, n) * n);
		//glm::vec3 u2 = normalize(cross(u1, n));
		//vlist[i]->DiscreteFrame[0] = u1;
		//vlist[i]->DiscreteFrame[1] = n;
		//vlist[i]->DiscreteFrame[2] = u2;
		//
		//
		////---------------------------------------x-------------------------------------------
		//for (int j = 0; j < poly->vlist[i]->ntris; j++)
		//{
		//	glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		//	glm::vec3 ij = (nextpoint - nowpoint);
		//	// Then use vi's frame to express vj's frame and get the coefficient
		//	float ct1 = dot(ij, u1); float ct2 = dot(ij, u2); float cn = dot(ij, n);
		//	//int indexi = poly->vlist[i]->index;
		//	//int indexj = c->n->v->index;
		//	//tripletListx.push_back(T(indexi, indexi, -1));
		//	//tripletListx.push_back(T(indexi, indexj,1));
		//	Ax(j, 0) = -1;
		//	bx[j] = -(c->n->v->x) + ct1 * vlist[i]->DiscreteFrame[0].x + ct2 * vlist[i]->DiscreteFrame[1].x + cn * vlist[i]->DiscreteFrame[2].x;
		//	c = c->p->o->p;
		//}
		//
		//if (poly->vlist[i]->index = poly->tlist[seed]->verts[0]->index)
		//{
		//	Ax(poly->vlist[i]->ntris, 0) = -1;
		//	bx[poly->vlist[i]->ntris] = -s_old;
		//}
		//else
		//{
		//	Ax(poly->vlist[i]->ntris, 0) = 0;
		//	bx[poly->vlist[i]->ntris] = 0;
		//}
		//
		////---------------------------------------y-------------------------------------------
		//for (int j = 0; j < poly->vlist[i]->ntris; j++)
		//{
		//	glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		//	glm::vec3 ij = (nextpoint - nowpoint);
		//	// Then use vi's frame to express vj's frame and get the coefficient
		//	float ct1 = dot(ij, u1); float ct2 = dot(ij, u2); float cn = dot(ij, n);
		//	//int indexi = poly->vlist[i]->index;
		//	//int indexj = c->n->v->index;
		//	//tripletListx.push_back(T(indexi, indexi, -1));
		//	//tripletListx.push_back(T(indexi, indexj,1));
		//	Ay(j, 0) = -1;
		//	by[j] = -(c->n->v->y) + ct1 * vlist[i]->DiscreteFrame[0].y + ct2 * vlist[i]->DiscreteFrame[1].y + cn * vlist[i]->DiscreteFrame[2].y;
		//	c = c->p->o->p;
		//}
		//
		//if (poly->vlist[i]->index = poly->tlist[seed]->verts[0]->index)
		//{
		//	Ay(poly->vlist[i]->ntris, 0) = -1;
		//	by[poly->vlist[i]->ntris] = -t_old;
		//}
		//else
		//{
		//	Ay(poly->vlist[i]->ntris, 0) = 0;
		//	by[poly->vlist[i]->ntris] = 0;
		//}
		//
		////---------------------------------------z-------------------------------------------
		//for (int j = 0; j < poly->vlist[i]->ntris; j++)
		//{
		//	glm::vec3 nextpoint = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
		//	glm::vec3 ij = (nextpoint - nowpoint);
		//	// Then use vi's frame to express vj's frame and get the coefficient
		//	float ct1 = dot(ij, u1); float ct2 = dot(ij, u2); float cn = dot(ij, n);
		//	//int indexi = poly->vlist[i]->index;
		//	//int indexj = c->n->v->index;
		//	//tripletListx.push_back(T(indexi, indexi, -1));
		//	//tripletListx.push_back(T(indexi, indexj,1));
		//	Az(j, 0) = -1;
		//	bz[j] = -(c->n->v->z) + ct1 * vlist[i]->DiscreteFrame[0].z + ct2 * vlist[i]->DiscreteFrame[1].z + cn * vlist[i]->DiscreteFrame[2].z;
		//	c = c->p->o->p;
		//}
		//
		//if (poly->vlist[i]->index = poly->tlist[seed]->verts[0]->index)
		//{
		//	Az(poly->vlist[i]->ntris, 0) = -1;
		//	bz[poly->vlist[i]->ntris] = -tridepth;
		//}
		//else
		//{
		//	Az(poly->vlist[i]->ntris, 0) = 0;
		//	bz[poly->vlist[i]->ntris] = 0;
		//}
		//
		//Eigen::VectorXd ix = Ax.colPivHouseholderQr().solve(bx);
		//Eigen::VectorXd iy = Ay.colPivHouseholderQr().solve(by);
		//Eigen::VectorXd iz = Az.colPivHouseholderQr().solve(bz);
		//
		//poly->vlist[i]->Newpoint = glm::vec3(ix[0], iy[0], iz[0]);



	//}

}

void Polyhedron::Hatching()
{
	GetTensor3D();   // get the 3d tensor of each vertex in present status
	
	TriangleTensoring();    // triangle center tensor interpolation

	GoBackTriTensor2D();
}

void Polyhedron::TriangleTensoring()
{
	for (int i = 0; i < ntris; i++)
	{

		//----------------------for the rotate matrix---------------------------------
		glm::vec3 n = normalize(glm::vec3(tlist[i]->normal.x, tlist[i]->normal.y, tlist[i]->normal.z));
		glm::vec3 u1 = glm::vec3(tlist[i]->verts[1]->x - tlist[i]->verts[0]->x,
								 tlist[i]->verts[1]->y - tlist[i]->verts[0]->y,
								 tlist[i]->verts[1]->z - tlist[i]->verts[0]->z);
		u1 = normalize(u1);
		glm::vec3 u2 = normalize(cross(u1, n));

		tlist[i]->Tangent1[0] = u1.x;
		tlist[i]->Tangent1[1] = u1.y;
		tlist[i]->Tangent1[2] = u1.z;

		tlist[i]->Tangent2[0] = u2.x;
		tlist[i]->Tangent2[1] = u2.y;
		tlist[i]->Tangent2[2] = u2.z;

		tlist[i]->RotateMatrix = glm::mat3(tlist[i]->Tangent1[0], tlist[i]->Tangent2[0], tlist[i]->normal.x,
										   tlist[i]->Tangent1[1], tlist[i]->Tangent2[1], tlist[i]->normal.y,
										   tlist[i]->Tangent1[2], tlist[i]->Tangent2[2], tlist[i]->normal.z);
		//------------------------------------------------------------------------------

		tlist[i]->Tensor3D = (tlist[i]->verts[0]->Tglobal + tlist[i]->verts[1]->Tglobal + tlist[i]->verts[2]->Tglobal) / (float)3.;
	}
}

void Polyhedron::GoBackTriTensor2D()
{
	for (int i = 0; i < ntris; i++)
	{
		glm::mat3 Tensormat = (tlist[i]->RotateMatrix) * tlist[i]->Tensor3D * transpose(tlist[i]->RotateMatrix);

		tlist[i]->Tensor2D[0] = Tensormat[0][0];
		tlist[i]->Tensor2D[1] = Tensormat[0][1];
		tlist[i]->Tensor2D[2] = Tensormat[1][1];

		CalculateTriPrincipalCurvature(i);
	}
}

void Polyhedron::SettingStreamlineFold()   // get the stream line array; Fold one
{
	Hatching();   //Initialize the 3D principle directions of each triangle


	for (int i = 0; i < ntris; i++)
	{
		glm::vec3 tiv0 = glm::vec3(tlist[i]->verts[0]->x, tlist[i]->verts[0]->y, tlist[i]->verts[0]->z);
		glm::vec3 tiv1 = glm::vec3(tlist[i]->verts[1]->x, tlist[i]->verts[1]->y, tlist[i]->verts[1]->z);
		glm::vec3 tiv2 = glm::vec3(tlist[i]->verts[2]->x, tlist[i]->verts[2]->y, tlist[i]->verts[2]->z);
		
		glm::vec3 InitialPoint = (tiv0+ tiv1+ tiv2)/(float)3.;
		tlist[i]->CenterPoint = InitialPoint;    // for testing

		tlist[i]->StreamLinePoints[6] = InitialPoint;

		glm::vec3 Tmax = glm::vec3(tlist[i]->principleTricurvaturevector1[0],
								   tlist[i]->principleTricurvaturevector1[1],
								   tlist[i]->principleTricurvaturevector1[2]);

		glm::vec3 PointV0, PointV1, PointV2;
		//-------------edge 0 judgement-------------------
		PointV0 = tiv0 - InitialPoint;      //only for judgement
		PointV1 = tiv1 - InitialPoint;		
		glm::vec3 V0V1 = tiv1 - tiv0;
		
		//-------------edge 1 judgement-------------------
		PointV1 = tiv1 - InitialPoint;
		PointV2 = tiv2 - InitialPoint;
		glm::vec3 V1V2 = tiv2 - tiv1;

		//-------------edge 2 judgement-------------------
		PointV2 = tiv2 - InitialPoint;
		PointV0 = tiv0 - InitialPoint;
		glm::vec3 V2V0 = tiv0 - tiv2;


		int edgemaxflag;   //save the edge index
		int triangleflag;
		triangleflag = tlist[i]->index;

		if (dot(cross(PointV0, Tmax), cross(PointV1, Tmax)) < 0)  // go to edge0
		{
			edgemaxflag = tlist[i]->edges[0]->index;

			//glm::vec3 b = cross((InitialPoint - tiv0), V0V1);
			//glm::vec3 a = cross(V0V1, Tmax);
			//glm::vec3 p = b / a;

			glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv0), V0V1)) / (cross(V0V1, Tmax))) * Tmax;
			//if(length(P1)>20000)  P1 = InitialPoint + (length(cross((InitialPoint - tiv0), V0V1)) / length(cross(V0V1, Tmax))) * Tmax;
			tlist[i]->StreamLinePoints[7] = P1;
		}
		else if (dot(cross(PointV1, Tmax), cross(PointV2, Tmax)) < 0)  // go to edge1
		{
			edgemaxflag = tlist[i]->edges[1]->index;
			glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv1), V1V2)) / (cross(V1V2, Tmax))) * Tmax;
			tlist[i]->StreamLinePoints[7] = P1;
		}
		else if (dot(cross(PointV2, Tmax), cross(PointV0, Tmax)) < 0)     // go to edge2
		{
			edgemaxflag = tlist[i]->edges[2]->index;
			glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv2), V2V0)) / (cross(V2V0, Tmax))) * Tmax;
			tlist[i]->StreamLinePoints[7] = P1;
		}

		for (int j = 0; j < 5; j++)   // the last five triangles
		{
			if (elist[edgemaxflag]->tris[0] == tlist[triangleflag])   // if edge's triangle equals to the old triangle, update triangle flag
			{
				triangleflag = elist[edgemaxflag]->tris[1]->index;
			}
			else 
				triangleflag = elist[edgemaxflag]->tris[0]->index;
		
			int edgenow;
			for (int p = 0; p < 3; p++)
			{
				if (tlist[triangleflag]->edges[p] == elist[edgemaxflag]) edgenow = p;
			}
		
			glm::vec3 tiv0 = glm::vec3(tlist[triangleflag]->verts[0]->x, tlist[triangleflag]->verts[0]->y, tlist[triangleflag]->verts[0]->z);
			glm::vec3 tiv1 = glm::vec3(tlist[triangleflag]->verts[1]->x, tlist[triangleflag]->verts[1]->y, tlist[triangleflag]->verts[1]->z);
			glm::vec3 tiv2 = glm::vec3(tlist[triangleflag]->verts[2]->x, tlist[triangleflag]->verts[2]->y, tlist[triangleflag]->verts[2]->z);
		
			glm::vec3 InitialPoint = tlist[i]->StreamLinePoints[7+j];   //DONT CHANGE i
		
			glm::vec3 Tmax = glm::vec3(tlist[triangleflag]->principleTricurvaturevector1[0],
									   tlist[triangleflag]->principleTricurvaturevector1[1],
									   tlist[triangleflag]->principleTricurvaturevector1[2]);
		
			glm::vec3 PointV0, PointV1, PointV2;
		
			//-------------edge 0 judgement-------------------
			PointV0 = tiv0 - InitialPoint;      //only for judgement
			PointV1 = tiv1 - InitialPoint;
			glm::vec3 V0V1 = tiv1 - tiv0;
		
			//-------------edge 1 judgement-------------------
			PointV1 = tiv1 - InitialPoint;
			PointV2 = tiv2 - InitialPoint;
			glm::vec3 V1V2 = tiv2 - tiv1;
		
			//-------------edge 2 judgement-------------------
			PointV2 = tiv2 - InitialPoint;
			PointV0 = tiv0 - InitialPoint;
			glm::vec3 V2V0 = tiv0 - tiv2;
		
			if (edgenow != 0 && dot(cross(PointV0, Tmax), cross(PointV1, Tmax)) < 0)  // go to edge0
			{
				edgemaxflag = tlist[triangleflag]->edges[0]->index;
				glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv0), V0V1)) / (cross(V0V1, Tmax))) * Tmax;
				tlist[i]->StreamLinePoints[7+j+1] = P1;
			}
			else if (edgenow != 1 && dot(cross(PointV1, Tmax), cross(PointV2, Tmax)) < 0)  // go to edge1
			{
				edgemaxflag = tlist[triangleflag]->edges[1]->index;
				glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv1), V1V2)) / (cross(V1V2, Tmax))) * Tmax;
				tlist[i]->StreamLinePoints[7 + j + 1] = P1;
			}
			else if (edgenow != 2 && dot(cross(PointV2, Tmax), cross(PointV0, Tmax)) < 0)  // go to edge2
			{
				edgemaxflag = tlist[triangleflag]->edges[2]->index;
				glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv2), V2V0)) / (cross(V2V0, Tmax))) * Tmax;
				tlist[i]->StreamLinePoints[7 + j + 1] = P1;
			}
			
		}



		tiv0 = glm::vec3(tlist[i]->verts[0]->x, tlist[i]->verts[0]->y, tlist[i]->verts[0]->z);
		tiv1 = glm::vec3(tlist[i]->verts[1]->x, tlist[i]->verts[1]->y, tlist[i]->verts[1]->z);
		tiv2 = glm::vec3(tlist[i]->verts[2]->x, tlist[i]->verts[2]->y, tlist[i]->verts[2]->z);
		
		InitialPoint = (tiv0 + tiv1 + tiv2) / (float)3.;
		
		glm::vec3 Tmin = glm::vec3(poly->tlist[i]->principleTricurvaturevector2[0],
								   poly->tlist[i]->principleTricurvaturevector2[1],
								   poly->tlist[i]->principleTricurvaturevector2[2]);
		
		//-------------edge 0 judgement-------------------
		PointV0 = tiv0 - InitialPoint;      //only for judgement
		PointV1 = tiv1 - InitialPoint;
		V0V1 = tiv1 - tiv0;
		
		//-------------edge 1 judgement-------------------
		PointV1 = tiv1 - InitialPoint;
		PointV2 = tiv2 - InitialPoint;
		V1V2 = tiv2 - tiv1;
		
		//-------------edge 2 judgement-------------------
		PointV2 = tiv2 - InitialPoint;
		PointV0 = tiv0 - InitialPoint;
		V2V0 = tiv0 - tiv2;
		
		int edgeminflag;   //save the edge index
		int triangleflag1;
		triangleflag1 = tlist[i]->index;
		
		if (dot(cross(PointV0, Tmin), cross(PointV1, Tmin)) < 0)  // go to edge0
		{
			edgeminflag = tlist[i]->edges[0]->index;
			glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv0), V0V1) / cross(V0V1, Tmin)) * Tmin;
			tlist[i]->StreamLinePoints[5] = P1;
		}
		else if (dot(cross(PointV1, Tmin), cross(PointV2, Tmin)) < 0)  // go to edge1
		{
			edgeminflag = tlist[i]->edges[1]->index;
			glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv1), V1V2) / cross(V1V2, Tmin)) * Tmin;
			tlist[i]->StreamLinePoints[5] = P1;
		}
		else if (dot(cross(PointV2, Tmin), cross(PointV0, Tmin)) < 0)  // go to edge2
		{
			edgeminflag = tlist[i]->edges[2]->index;
			glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv2), V2V0) / cross(V2V0, Tmin)) * Tmin;
			tlist[i]->StreamLinePoints[5] = P1;
		}
		
		for (int j = 0; j < 5; j++)   // the last five triangles
		{
			if (elist[edgeminflag]->tris[0] == tlist[triangleflag1])   // if edge's triangle equals to the old triangle, update triangle flag
			{
				triangleflag1 = elist[edgeminflag]->tris[1]->index;
			}
			else
				triangleflag1 = elist[edgeminflag]->tris[0]->index;
		
			int edgenow;
			for (int p = 0; p < 3; p++)
			{
				if (tlist[triangleflag1]->edges[p] == elist[edgeminflag]) edgenow = p;
			}
		
			glm::vec3 tiv0 = glm::vec3(tlist[triangleflag1]->verts[0]->x, tlist[triangleflag1]->verts[0]->y, tlist[triangleflag1]->verts[0]->z);
			glm::vec3 tiv1 = glm::vec3(tlist[triangleflag1]->verts[1]->x, tlist[triangleflag1]->verts[1]->y, tlist[triangleflag1]->verts[1]->z);
			glm::vec3 tiv2 = glm::vec3(tlist[triangleflag1]->verts[2]->x, tlist[triangleflag1]->verts[2]->y, tlist[triangleflag1]->verts[2]->z);
		
			glm::vec3 InitialPoint = tlist[i]->StreamLinePoints[5 - j];   //DONT CHANGE i
		
			glm::vec3 Tmin = glm::vec3(tlist[triangleflag1]->principleTricurvaturevector2[0],
									   tlist[triangleflag1]->principleTricurvaturevector2[1],
									   tlist[triangleflag1]->principleTricurvaturevector2[2]);
		
			glm::vec3 PointV0, PointV1, PointV2;
		
			//-------------edge 0 judgement-------------------
			PointV0 = tiv0 - InitialPoint;      //only for judgement
			PointV1 = tiv1 - InitialPoint;
			glm::vec3 V0V1 = tiv1 - tiv0;
		
			//-------------edge 1 judgement-------------------
			PointV1 = tiv1 - InitialPoint;
			PointV2 = tiv2 - InitialPoint;
			glm::vec3 V1V2 = tiv2 - tiv1;
		
			//-------------edge 2 judgement-------------------
			PointV2 = tiv2 - InitialPoint;
			PointV0 = tiv0 - InitialPoint;
			glm::vec3 V2V0 = tiv0 - tiv2;
		
			if (edgenow != 0 && dot(cross(PointV0, Tmin), cross(PointV1, Tmin)) < 0)  // go to edge0
			{
				edgeminflag = tlist[triangleflag1]->edges[0]->index;
				glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv0), V0V1) / cross(V0V1, Tmin)) * Tmin;
				tlist[i]->StreamLinePoints[5 - j - 1] = P1;
			}
			else if (edgenow != 1 && dot(cross(PointV1, Tmin), cross(PointV2, Tmin)) < 0)  // go to edge1
			{
				edgeminflag = tlist[triangleflag1]->edges[1]->index;
				glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv1), V1V2) / cross(V1V2, Tmin)) * Tmin;
				tlist[i]->StreamLinePoints[5 - j - 1] = P1;
			}
			else if (edgenow != 2 && dot(cross(PointV2, Tmin), cross(PointV0, Tmin)) < 0)  // go to edge2
			{
				edgeminflag = tlist[triangleflag1]->edges[2]->index;
				glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv2), V2V0) / cross(V2V0, Tmin)) * Tmin;
				tlist[i]->StreamLinePoints[5 - j - 1] = P1;
			}
		
		}
	}
}

void Polyhedron::SettingStreamlineUnfold()  // get the stream line array; Unold one
{
	Hatching();   //Initialize the 3D principle directions of each triangle

	for (int i = 0; i < ntris; i++)
	{
		glm::vec3 tiv0 = glm::vec3(tlist[i]->verts[0]->x, tlist[i]->verts[0]->y, tlist[i]->verts[0]->z);
		glm::vec3 tiv1 = glm::vec3(tlist[i]->verts[1]->x, tlist[i]->verts[1]->y, tlist[i]->verts[1]->z);
		glm::vec3 tiv2 = glm::vec3(tlist[i]->verts[2]->x, tlist[i]->verts[2]->y, tlist[i]->verts[2]->z);

		glm::vec3 InitialPoint = (tiv0 + tiv1 + tiv2) / (float)3.;
		tlist[i]->CenterPoint = InitialPoint;    // for testing

		tlist[i]->StreamLinePoints[6] = InitialPoint;

		glm::vec3 Tmax = glm::vec3(tlist[i]->principleTricurvaturevector1[0],
			tlist[i]->principleTricurvaturevector1[1],
			tlist[i]->principleTricurvaturevector1[2]);

		glm::vec3 PointV0, PointV1, PointV2;
		//-------------edge 0 judgement-------------------
		PointV0 = tiv0 - InitialPoint;      //only for judgement
		PointV1 = tiv1 - InitialPoint;
		glm::vec3 V0V1 = tiv1 - tiv0;

		//-------------edge 1 judgement-------------------
		PointV1 = tiv1 - InitialPoint;
		PointV2 = tiv2 - InitialPoint;
		glm::vec3 V1V2 = tiv2 - tiv1;

		//-------------edge 2 judgement-------------------
		PointV2 = tiv2 - InitialPoint;
		PointV0 = tiv0 - InitialPoint;
		glm::vec3 V2V0 = tiv0 - tiv2;


		int edgemaxflag;   //save the edge index
		int triangleflag;
		triangleflag = tlist[i]->index;

		if (dot(cross(PointV0, Tmax), cross(PointV1, Tmax)) < 0)  // go to edge0
		{
			edgemaxflag = tlist[i]->edges[0]->index;

			//glm::vec3 b = cross((InitialPoint - tiv0), V0V1);
			//glm::vec3 a = cross(V0V1, Tmax);
			//glm::vec3 p = b / a;

			glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv0), V0V1)) / (cross(V0V1, Tmax))) * Tmax;
			//if(length(P1)>20000)  P1 = InitialPoint + (length(cross((InitialPoint - tiv0), V0V1)) / length(cross(V0V1, Tmax))) * Tmax;
			tlist[i]->StreamLinePoints[7] = P1;
		}
		else if (dot(cross(PointV1, Tmax), cross(PointV2, Tmax)) < 0)  // go to edge1
		{
			edgemaxflag = tlist[i]->edges[1]->index;
			glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv1), V1V2)) / (cross(V1V2, Tmax))) * Tmax;
			tlist[i]->StreamLinePoints[7] = P1;
		}
		else if (dot(cross(PointV2, Tmax), cross(PointV0, Tmax)) < 0)     // go to edge2
		{
			edgemaxflag = tlist[i]->edges[2]->index;
			glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv2), V2V0)) / (cross(V2V0, Tmax))) * Tmax;
			tlist[i]->StreamLinePoints[7] = P1;
		}

		glm::mat3 Told3D = tlist[i]->Tensor3D;

		for (int j = 0; j < 5; j++)   // the last five triangles
		{
			int triangleflagold;
			if (elist[edgemaxflag]->tris[0] == tlist[triangleflag])   // if edge's triangle equals to the old triangle, update triangle flag
			{
				triangleflag = elist[edgemaxflag]->tris[1]->index;
				triangleflagold = elist[edgemaxflag]->tris[0]->index;
			}
			else
			{
				triangleflag = elist[edgemaxflag]->tris[0]->index;
				triangleflagold = elist[edgemaxflag]->tris[1]->index;
			}
			int edgenow;
			for (int p = 0; p < 3; p++)
			{
				if (tlist[triangleflag]->edges[p] == elist[edgemaxflag]) edgenow = p;
			}

			glm::vec3 tiv0 = glm::vec3(tlist[triangleflag]->verts[0]->x, tlist[triangleflag]->verts[0]->y, tlist[triangleflag]->verts[0]->z);
			glm::vec3 tiv1 = glm::vec3(tlist[triangleflag]->verts[1]->x, tlist[triangleflag]->verts[1]->y, tlist[triangleflag]->verts[1]->z);
			glm::vec3 tiv2 = glm::vec3(tlist[triangleflag]->verts[2]->x, tlist[triangleflag]->verts[2]->y, tlist[triangleflag]->verts[2]->z);

			glm::vec3 InitialPoint = tlist[i]->StreamLinePoints[7 + j];   //DONT CHANGE i

//----------------------------------------------------------------------------------------------------------

			glm::mat3 TransformMat = tlist[triangleflag]->Tensor3D;
			glm::mat3 Tnow3D = transpose(TransformMat) * Told3D * (TransformMat);

			//---------------------back from 3D to 2D--------------------------
			glm::mat3 Tensormat = (tlist[triangleflag]->RotateMatrix) * Tnow3D * transpose(tlist[triangleflag]->RotateMatrix);
			//Told3D = Tensormat;



			//glm::mat2 Told2D = glm::mat2(tlist[triangleflagold]->Tensor2D[0], tlist[triangleflagold]->Tensor2D[1],
			//							 tlist[triangleflagold]->Tensor2D[1], tlist[triangleflagold]->Tensor2D[2]);
			//
			//float theta;
			//
			//glm::vec3 pdirection = normalize(glm::vec3(tlist[triangleflag]->principleTricurvaturevector1[0],
			//								 tlist[triangleflag]->principleTricurvaturevector1[1], 
			//								 tlist[triangleflag]->principleTricurvaturevector1[2]));
			//glm::vec3 oldpdirection = normalize(glm::vec3(tlist[triangleflagold]->principleTricurvaturevector1[0],
			//									tlist[triangleflagold]->principleTricurvaturevector1[1],
			//									tlist[triangleflagold]->principleTricurvaturevector1[2]));
			//
			////glm::vec3 pdirection = normalize(glm::vec3(tlist[triangleflag]->principleTricurvaturevector1[0],
			////			   							   tlist[triangleflag]->principleTricurvaturevector1[1], 
			////			   							   tlist[triangleflag]->principleTricurvaturevector1[2]));
			////
			////glm::vec3 oldpdirection = normalize(glm::vec3(tlist[i]->normal.x, tlist[i]->normal.y, tlist[i]->normal.z));
			//
			//
			//theta = acos(dot(pdirection, oldpdirection));
			//if (theta > PI / 2.) theta -= PI / 2.;
			//
			////glm::mat2 TransformMat = glm::mat2(cos(theta-PI/2.), -sin(theta - PI / 2.), sin(theta - PI / 2.), cos(theta - PI /2.));
			//glm::mat2 TransformMat = glm::mat2(cos(theta), -sin(theta), sin(theta), cos(theta));
			//glm::mat2 Tnow2D = transpose(TransformMat) * Told2D * (TransformMat);
			//
			//glm::mat2 Tensormat = Tnow2D;



			//-----------------Get the 3D Eig vectors------------------
			Eigen::MatrixXd M(2, 2);
			M(0, 0) = Tensormat[0][0];
			M(0, 1) = Tensormat[0][1];
			M(1, 0) = Tensormat[0][1];
			M(1, 1) = Tensormat[1][1];

			Eigen::EigenSolver<Eigen::MatrixXd> es(M);
			//Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);

			Eigen::VectorXcd v1 = es.eigenvectors().col(0);
			Eigen::VectorXcd v2 = es.eigenvectors().col(1);

			//std::cout << "The principle curvature vector is:\n"
			//	<< es.eigenvectors() << std::endl;
			//

			double x1 = v1[0].real();
			double y1 = v1[1].real();

			double x2 = v2[0].real();
			double y2 = v2[1].real();

			glm::vec3 u1 = glm::vec3(tlist[triangleflag]->Tangent1[0], tlist[triangleflag]->Tangent1[1], tlist[triangleflag]->Tangent1[2]);
			glm::vec3 u2 = glm::vec3(tlist[triangleflag]->Tangent2[0], tlist[triangleflag]->Tangent2[1], tlist[triangleflag]->Tangent2[2]);


			glm::vec3 principal1 = normalize((float)x1 * u1 + (float)y1 * u2);
			glm::vec3 principal2 = normalize((float)x2 * u1 + (float)y2 * u2);

			glm::vec3 Tmax;


			//--------------------------------------check new direction point to the inside of the new triangle----------------------------------------------------
			glm::vec3 normal = glm::vec3(tlist[triangleflag]->normal.x, tlist[triangleflag]->normal.y, tlist[triangleflag]->normal.z);
			glm::vec3 edge = glm::vec3(elist[edgemaxflag]->verts[1]->x - elist[edgemaxflag]->verts[0]->x,
				elist[edgemaxflag]->verts[1]->y - elist[edgemaxflag]->verts[0]->y,
				elist[edgemaxflag]->verts[1]->z - elist[edgemaxflag]->verts[0]->z);

			glm::vec3 edgenormal = cross(normal, edge);     // for check the direction of the new direction on a new plane wheather it point to the inside of the triangle

			glm::vec3 toppoint;
			for (int f = 0; f < 3; f++)
			{
				if (tlist[triangleflag]->verts[f] != elist[edgemaxflag]->verts[0] && tlist[triangleflag]->verts[f] != elist[edgemaxflag]->verts[1])
				{
					toppoint = glm::vec3(tlist[triangleflag]->verts[f]->x, tlist[triangleflag]->verts[f]->y, tlist[triangleflag]->verts[f]->z);
				}
			}
			glm::vec3 bottompoint = glm::vec3(elist[edgemaxflag]->verts[0]->x, elist[edgemaxflag]->verts[0]->y, elist[edgemaxflag]->verts[0]->z);

			if (dot((toppoint - bottompoint), edgenormal) < 0)
			{
				edgenormal = -edgenormal;
			}
			//-------------------------------------------------------------------------------------------------------------------------------------
			if (dot(principal1, u1) >= 0)
			{
				Tmax = glm::vec3(principal1.x, principal1.y, principal1.z);
			}
			else
			{
				Tmax = glm::vec3(principal2.x, principal2.y, principal2.z);
			}

			if (dot(Tmax, edgenormal) < 0)
			{
				Tmax = -Tmax;
			}

			//Tmax = glm::vec3(poly->tlist[triangleflag]->principleTricurvaturevector1[0],
			//poly->tlist[triangleflag]->principleTricurvaturevector1[1],
			//poly->tlist[triangleflag]->principleTricurvaturevector1[2]);

//------------------------------------------------------------------------------------------------------------
			glm::vec3 PointV0, PointV1, PointV2;

			//-------------edge 0 judgement-------------------
			PointV0 = tiv0 - InitialPoint;      //only for judgement
			PointV1 = tiv1 - InitialPoint;
			glm::vec3 V0V1 = tiv1 - tiv0;

			//-------------edge 1 judgement-------------------
			PointV1 = tiv1 - InitialPoint;
			PointV2 = tiv2 - InitialPoint;
			glm::vec3 V1V2 = tiv2 - tiv1;

			//-------------edge 2 judgement-------------------
			PointV2 = tiv2 - InitialPoint;
			PointV0 = tiv0 - InitialPoint;
			glm::vec3 V2V0 = tiv0 - tiv2;

			//int sss = j;
			//cross(PointV0, Tmax);
			//cross(PointV1, Tmax);


			if (edgenow != 0 && dot(cross(PointV0, Tmax), cross(PointV1, Tmax)) < 0)  // go to edge0
			{
				edgemaxflag = tlist[triangleflag]->edges[0]->index;
				glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv0), V0V1)) / (cross(V0V1, Tmax))) * Tmax;
				tlist[i]->StreamLinePoints[7 + j + 1] = P1;
			}
			else if (edgenow != 1 && dot(cross(PointV1, Tmax), cross(PointV2, Tmax)) < 0)  // go to edge1
			{
				edgemaxflag = tlist[triangleflag]->edges[1]->index;
				glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv1), V1V2)) / (cross(V1V2, Tmax))) * Tmax;
				tlist[i]->StreamLinePoints[7 + j + 1] = P1;
			}
			else if (edgenow != 2 && dot(cross(PointV2, Tmax), cross(PointV0, Tmax)) < 0)  // go to edge2
			{
				edgemaxflag = tlist[triangleflag]->edges[2]->index;
				glm::vec3 P1 = InitialPoint + ((cross((InitialPoint - tiv2), V2V0)) / (cross(V2V0, Tmax))) * Tmax;
				tlist[i]->StreamLinePoints[7 + j + 1] = P1;
			}
		}


		//--------------------------------------------------------------------------------------------
		//******************************************* MIN *********************************************
		//--------------------------------------------------------------------------------------------
		tiv0 = glm::vec3(tlist[i]->verts[0]->x, tlist[i]->verts[0]->y, tlist[i]->verts[0]->z);
		tiv1 = glm::vec3(tlist[i]->verts[1]->x, tlist[i]->verts[1]->y, tlist[i]->verts[1]->z);
		tiv2 = glm::vec3(tlist[i]->verts[2]->x, tlist[i]->verts[2]->y, tlist[i]->verts[2]->z);

		InitialPoint = (tiv0 + tiv1 + tiv2) / (float)3.;

		glm::vec3 Tmin = glm::vec3(poly->tlist[i]->principleTricurvaturevector2[0],
			poly->tlist[i]->principleTricurvaturevector2[1],
			poly->tlist[i]->principleTricurvaturevector2[2]);

		//-------------edge 0 judgement-------------------
		PointV0 = tiv0 - InitialPoint;      //only for judgement
		PointV1 = tiv1 - InitialPoint;
		V0V1 = tiv1 - tiv0;

		//-------------edge 1 judgement-------------------
		PointV1 = tiv1 - InitialPoint;
		PointV2 = tiv2 - InitialPoint;
		V1V2 = tiv2 - tiv1;

		//-------------edge 2 judgement-------------------
		PointV2 = tiv2 - InitialPoint;
		PointV0 = tiv0 - InitialPoint;
		V2V0 = tiv0 - tiv2;

		int edgeminflag;   //save the edge index
		int triangleflag1;
		triangleflag1 = tlist[i]->index;

		if (dot(cross(PointV0, Tmin), cross(PointV1, Tmin)) < 0)  // go to edge0
		{
			edgeminflag = tlist[i]->edges[0]->index;
			glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv0), V0V1) / cross(V0V1, Tmin)) * Tmin;
			tlist[i]->StreamLinePoints[5] = P1;
		}
		else if (dot(cross(PointV1, Tmin), cross(PointV2, Tmin)) < 0)  // go to edge1
		{
			edgeminflag = tlist[i]->edges[1]->index;
			glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv1), V1V2) / cross(V1V2, Tmin)) * Tmin;
			tlist[i]->StreamLinePoints[5] = P1;
		}
		else if (dot(cross(PointV2, Tmin), cross(PointV0, Tmin)) < 0)  // go to edge2
		{
			edgeminflag = tlist[i]->edges[2]->index;
			glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv2), V2V0) / cross(V2V0, Tmin)) * Tmin;
			tlist[i]->StreamLinePoints[5] = P1;
		}

		Told3D = tlist[i]->Tensor3D;

		for (int j = 0; j < 5; j++)   // the last five triangles
		{
			if (elist[edgeminflag]->tris[0] == tlist[triangleflag1])   // if edge's triangle equals to the old triangle, update triangle flag
			{
				triangleflag1 = elist[edgeminflag]->tris[1]->index;
			}
			else
				triangleflag1 = elist[edgeminflag]->tris[0]->index;

			int edgenow;
			for (int p = 0; p < 3; p++)
			{
				if (tlist[triangleflag1]->edges[p] == elist[edgeminflag]) edgenow = p;
			}

			glm::vec3 tiv0 = glm::vec3(tlist[triangleflag1]->verts[0]->x, tlist[triangleflag1]->verts[0]->y, tlist[triangleflag1]->verts[0]->z);
			glm::vec3 tiv1 = glm::vec3(tlist[triangleflag1]->verts[1]->x, tlist[triangleflag1]->verts[1]->y, tlist[triangleflag1]->verts[1]->z);
			glm::vec3 tiv2 = glm::vec3(tlist[triangleflag1]->verts[2]->x, tlist[triangleflag1]->verts[2]->y, tlist[triangleflag1]->verts[2]->z);

			glm::vec3 InitialPoint = tlist[i]->StreamLinePoints[5 - j];   //DONT CHANGE i

//----------------------------------------------------------------------------------------------------------
			
			glm::mat3 TransformMat = tlist[triangleflag1]->Tensor3D;
			glm::mat3 Tnow3D = transpose(TransformMat) * Told3D * TransformMat;

			//---------------------back to 2D--------------------------
			glm::mat3 Tensormat = (tlist[triangleflag1]->RotateMatrix) * Tnow3D * transpose(tlist[triangleflag1]->RotateMatrix);


			//-----------------Get the 3D Eig vectors------------------
			Eigen::MatrixXd M(2, 2);
			M(0, 0) = Tensormat[0][0];
			M(0, 1) = Tensormat[0][1];
			M(1, 0) = Tensormat[0][1];
			M(1, 1) = Tensormat[1][1];

			Eigen::EigenSolver<Eigen::MatrixXd> es(M);
			//Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);

			Eigen::VectorXcd v1 = es.eigenvectors().col(0);
			Eigen::VectorXcd v2 = es.eigenvectors().col(1);

			//std::cout << "The principle curvature vector is:\n"
			//	<< es.eigenvectors() << std::endl;
			//

			double x1 = v1[0].real();
			double y1 = v1[1].real();

			double x2 = v2[0].real();
			double y2 = v2[1].real();

			glm::vec3 u1 = glm::vec3(tlist[triangleflag1]->Tangent1[0], tlist[triangleflag1]->Tangent1[1], tlist[triangleflag1]->Tangent1[2]);
			glm::vec3 u2 = glm::vec3(tlist[triangleflag1]->Tangent2[0], tlist[triangleflag1]->Tangent2[1], tlist[triangleflag1]->Tangent2[2]);

			//double check = dot(glm::vec2(x1, y1), glm::vec2(x2, y2));

			glm::vec3 principal1 = normalize((float)x1 * u1 + (float)y1 * u2);
			glm::vec3 principal2 = normalize((float)x2 * u1 + (float)y2 * u2);


			glm::vec3 Tmin;
			//--------------------------------------check new direction point to the inside of the new triangle----------------------------------------------------
			glm::vec3 normal = glm::vec3(tlist[triangleflag1]->normal.x, tlist[triangleflag1]->normal.y, tlist[triangleflag1]->normal.z);
			glm::vec3 edge = glm::vec3(elist[edgeminflag]->verts[1]->x - elist[edgeminflag]->verts[0]->x,
				elist[edgeminflag]->verts[1]->y - elist[edgeminflag]->verts[0]->y,
				elist[edgeminflag]->verts[1]->z - elist[edgeminflag]->verts[0]->z);

			glm::vec3 edgenormal = cross(normal, edge);     // for check the direction of the new direction on a new plane wheather it point to the inside of the triangle

			glm::vec3 toppoint;
			for (int f = 0; f < 3; f++)
			{
				if (tlist[triangleflag1]->verts[f] != elist[edgeminflag]->verts[0] && tlist[triangleflag1]->verts[f] != elist[edgeminflag]->verts[1])
				{
					toppoint = glm::vec3(tlist[triangleflag1]->verts[f]->x, tlist[triangleflag1]->verts[f]->y, tlist[triangleflag1]->verts[f]->z);
				}
			}
			glm::vec3 bottompoint = glm::vec3(elist[edgeminflag]->verts[0]->x, elist[edgeminflag]->verts[0]->y, elist[edgeminflag]->verts[0]->z);

			if (dot((toppoint - bottompoint), edgenormal) < 0)
			{
				edgenormal = -edgenormal;
			}
			//-------------------------------------------------------------------------------------------------------------------------------------

			//double check = dot(u1, u2);

			if (dot(principal1, u1) >= 0)
			{
				Tmin = glm::vec3(principal2.x, principal2.y, principal2.z);
			}
			else
			{
				Tmin = glm::vec3(principal1.x, principal1.y, principal1.z);
			}

			if (dot(Tmin, edgenormal) < 0)
			{
				Tmin = -Tmin;
			}
			//------------------------------------------------------------------------------------------------------------

			glm::vec3 PointV0, PointV1, PointV2;

			//-------------edge 0 judgement-------------------
			PointV0 = tiv0 - InitialPoint;      //only for judgement
			PointV1 = tiv1 - InitialPoint;
			glm::vec3 V0V1 = tiv1 - tiv0;

			//-------------edge 1 judgement-------------------
			PointV1 = tiv1 - InitialPoint;
			PointV2 = tiv2 - InitialPoint;
			glm::vec3 V1V2 = tiv2 - tiv1;

			//-------------edge 2 judgement-------------------
			PointV2 = tiv2 - InitialPoint;
			PointV0 = tiv0 - InitialPoint;
			glm::vec3 V2V0 = tiv0 - tiv2;


			if (edgenow != 0 && dot(cross(PointV0, Tmin), cross(PointV1, Tmin)) < 0)  // go to edge0
			{
				edgeminflag = tlist[triangleflag1]->edges[0]->index;
				glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv0), V0V1) / cross(V0V1, Tmin)) * Tmin;
				tlist[i]->StreamLinePoints[5 - j - 1] = P1;
			}
			else if (edgenow != 1 && dot(cross(PointV1, Tmin), cross(PointV2, Tmin)) < 0)  // go to edge1
			{
				edgeminflag = tlist[triangleflag1]->edges[1]->index;
				glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv1), V1V2) / cross(V1V2, Tmin)) * Tmin;
				tlist[i]->StreamLinePoints[5 - j - 1] = P1;
			}
			else if (edgenow != 2 && dot(cross(PointV2, Tmin), cross(PointV0, Tmin)) < 0)  // go to edge2
			{
				edgeminflag = tlist[triangleflag1]->edges[2]->index;
				glm::vec3 P1 = InitialPoint + (cross((InitialPoint - tiv2), V2V0) / cross(V2V0, Tmin)) * Tmin;
				tlist[i]->StreamLinePoints[5 - j - 1] = P1;
			}

		}
	}
}



void Polyhedron::GetTensor3D()
{
	for (int i = 0; i < nverts; i++)
	{
		vlist[i]->t1t2N = glm::mat3(vlist[i]->TangentVector1[0], vlist[i]->TangentVector2[0], vlist[i]->normal.x,
									vlist[i]->TangentVector1[1], vlist[i]->TangentVector2[1], vlist[i]->normal.y,
									vlist[i]->TangentVector1[2], vlist[i]->TangentVector2[2], vlist[i]->normal.z);

		glm::mat3 lmn = glm::mat3(vlist[i]->Tensor[0], vlist[i]->Tensor[1], 0,
								  vlist[i]->Tensor[1], vlist[i]->Tensor[2], 0,
								  0,                   0,                   0);

		vlist[i]->Tglobal = transpose(vlist[i]->t1t2N) * lmn * (vlist[i]->t1t2N);
	}
}

void Polyhedron::TensorSmoothing1()   //Explicit - Global - easier corner
{
	for (int i = 0; i < nverts; i++)
	{
		vlist[i]->TglobalSave = vlist[i]->Tglobal;
	}
	//-----------------------Mean value coordinates  && Uniform-------------------------//
	for (int i = 0; i < nverts; i++)
	{
	
		double dt = 0.1;
		int n = vlist[i]->ntris;
		double sum_eij = 0.0;
	
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
			{
				c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
				break;
			}
		}
	
		for (int j = 0; j < n; j++)
		{
			sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
			c = c->p->o->p;      // corner to the next triangle
		}
	
	
		glm::mat3 p_sum = glm::mat3(0);
	
		for (int j = 0; j < n; j++)
		{
	
			double wij = 0.0;
	
			double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
			wij = eij / sum_eij; // 1. / (float)vlist[i]->ntris;
	
			p_sum += (float)wij * (c->n->v->TglobalSave - vlist[i]->TglobalSave);
			c = c->p->o->p;
	
		}
		
		vlist[i]->Tglobal = vlist[i]->TglobalSave + (float)dt * p_sum;
	}
	
	//---------------------------Cord----------------------------//
	//for (int i = 0; i < nverts; i++)
	//{
	//
	//	double dt = 0.9;
	//	int n = vlist[i]->ntris;
	//	double sum_eij = 0.0;
	//
	//	Corner* c;
	//	for (int k = 0; k < 3; k++)
	//	{
	//		if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
	//		{
	//			c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
	//			break;
	//		}
	//	}
	//
	//	for (int j = 0; j < n; j++)
	//	{
	//		sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
	//		sum_eij += 1.0 / c->p->e->length;
	//		c = c->p->o->p;      // corner to the next triangle
	//	}
	//
	//
	//	glm::mat3 p_sum = glm::mat3(0);
	//
	//	for (int j = 0; j < n; j++)
	//	{
	//
	//		double wij = 0.0;
	//
	//		double eij = 1.0 / c->p->e->length;
	//		wij = eij / sum_eij;
	//
	//		p_sum += (float)wij * (c->n->v->TglobalSave - vlist[i]->TglobalSave);
	//		c = c->p->o->p;
	//
	//	}
	//
	//	vlist[i]->Tglobal = vlist[i]->TglobalSave + (float)dt * p_sum;
	//}

	//-----------------------Mean curvature-------------------------//
	//for (int i = 0; i < nverts; i++)
	//{
	//
	//	double dt = 0.9;
	//	int n = vlist[i]->ntris;
	//	double sum_eij = 0.0;
	//	double big = 1.0;
	//
	//	Corner* c;
	//	for (int k = 0; k < 3; k++)
	//	{
	//		if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
	//		{
	//			c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
	//			break;
	//		}
	//	}
	//
	//	for (int j = 0; j < n; j++)
	//	{
	//		if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) < 0.000005)
	//			sum_eij += big;
	//		else if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) >= 0.000005)
	//			sum_eij += (big + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
	//		else if (sin(c->p->c) >= 0.000005 && sin(c->p->o->c) < 0.000005)
	//			sum_eij += (cos(c->p->c) / sin(c->p->c) + big) / 2.0;
	//		else
	//			sum_eij += (cos(c->p->c) / sin(c->p->c) + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
	//		c = c->p->o->p;      // corner to the next triangle
	//	}
	//
	//
	//	glm::mat3 p_sum = glm::mat3(0);
	//
	//	for (int j = 0; j < n; j++)
	//	{
	//
	//		double wij = 0.0;
	//
	//		double eij;
	//		if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) < 0.000005)
	//			eij = big;
	//		else if (sin(c->p->c) < 0.000005 && sin(c->p->o->c) >= 0.000005)
	//			eij = (big + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
	//		else if (sin(c->p->c) >= 0.000005 && sin(c->p->o->c) < 0.000005)
	//			eij = (cos(c->p->c) / sin(c->p->c) + big) / 2.0;
	//		else
	//			eij = (cos(c->p->c) / sin(c->p->c) + cos(c->p->o->c) / sin(c->p->o->c)) / 2.0;
	//		
	//		wij = eij / sum_eij;
	//
	//		p_sum += (float)wij * (c->n->v->TglobalSave - vlist[i]->TglobalSave);
	//		c = c->p->o->p;
	//
	//	}
	//
	//	vlist[i]->Tglobal = vlist[i]->TglobalSave + (float)dt * p_sum;
	//}



}

void Polyhedron::GoBackTensor2D()
{
	for (int i = 0; i < nverts; i++)
	{
		glm::mat3 Tensormat = (vlist[i]->t1t2N) * vlist[i]->Tglobal * transpose(vlist[i]->t1t2N);

		vlist[i]->Tensor[0] = Tensormat[0][0];
		vlist[i]->Tensor[1] = Tensormat[0][1];
		vlist[i]->Tensor[2] = Tensormat[1][1];

		CalculatePrincipalCurvature(i);
	}
}


void Polyhedron::TensorSmoothing2(int iterationtime)
{

	double dt = 1.;

	dt = dt * iterationtime;  // Set dt as the final time

	int m = poly->nverts;   // Number of unknown elements (equal to number of "vertices")

	//Assembly:
	std::vector<T> tripletList0;   // list of non-zeros coefficients the triple value
	std::vector<T> tripletList1;
	std::vector<T> tripletList2;
	Eigen::VectorXd b0(m);     // b vector at the right side of the equation
	Eigen::VectorXd b1(m);
	Eigen::VectorXd b2(m);

	//BuildProblem: push back coefficient

//----------------------------------------------------------------------------------------------------
	b0.setZero();
	b1.setZero();
	b2.setZero();

	for (int i = 0; i < m; i++)  // traverse row
	{
		int n = poly->vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					c = poly->vlist[i]->tris[j]->corners[k];
					sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
				}
			}
		}

		for (int j = 0; j < n; j++)   // traverse column
		{
			Corner* c;
			double wij = 0.0;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					c = poly->vlist[i]->tris[j]->corners[k];
					double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
					wij = eij / sum_eij;
					break;
				}
			}
			int a;
			//for (int q = 0; q < m; q++)
			//{
			//	if (poly->vlist[q] == c->p->v)
			//	{
			//		a = q;
			//		break;
			//	}
			//
			//}
			a = c->n->v->index;
			tripletList0.push_back(T(i, a, -dt * wij));
			tripletList1.push_back(T(i, a, -dt * wij));
			tripletList2.push_back(T(i, a, -dt * wij));
			b0[i] = poly->vlist[i]->Tensor[0];
			b1[i] = poly->vlist[i]->Tensor[1];
			b2[i] = poly->vlist[i]->Tensor[2];

			tripletList0.push_back(T(i, i, 1 + dt));
			tripletList1.push_back(T(i, i, 1 + dt));
			tripletList2.push_back(T(i, i, 1 + dt));
		}
	}


	//----------------------------------------------------------------------------------------------------
	SpMat A0(m, m);            // A matrix at the left side of the equation
	SpMat A1(m, m);
	SpMat A2(m, m);
	A0.setFromTriplets(tripletList0.begin(), tripletList0.end());  // sent the triple coefficient to matrix
	A1.setFromTriplets(tripletList1.begin(), tripletList1.end());
	A2.setFromTriplets(tripletList2.begin(), tripletList2.end());
	//Solve problem:

	Eigen::SparseLU<SpMat> solver0(A0);
	Eigen::SparseLU<SpMat> solver1(A1);
	Eigen::SparseLU<SpMat> solver2(A2);

	//Eigen::SimplicialCholesky<SpMat> solver(A);
	solver0.compute(A0);
	solver1.compute(A1);
	solver2.compute(A2);
	//if (solver.info() != Eigen::Success)
	//{
	//	std::cout << "Compute Matrix is error" << std::endl;
	//	return;
	//}

	Eigen::VectorXd x0 = solver0.solve(b0);
	Eigen::VectorXd x1 = solver1.solve(b1);
	Eigen::VectorXd x2 = solver2.solve(b2);

	//----do the saved value----//
	for (int i = 0; i < m; i++)
	{
			poly->vlist[i]->Tensor[0] = x0[i];
			poly->vlist[i]->Tensor[1] = x1[i];
			poly->vlist[i]->Tensor[2] = x2[i];
	}

}


void Polyhedron::setheat()
{
	maxindex = 3;
	minindex = 2;
	for (int i = 0; i < nverts; i++)
	{
		if (vlist[i] != tlist[maxindex]->verts[0] && vlist[i] != tlist[minindex]->verts[0])
		//if (vlist[i] != tlist[3]->verts[0] && vlist[i] != tlist[2]->verts[0])
		{
			vlist[i]->heatvalue = 0.0;
		}
	}
}

void Polyhedron::heat1(Polyhedron* poly, int maxindex, int minindex)
{

	//std::vector<int> array;
	
	poly->tlist[maxindex]->verts[0]->heatvalue = 1.0; //set the maxima
	poly->tlist[minindex]->verts[0]->heatvalue = 0.0; //set the minima


	for (int i = 0; i < poly->nverts; i++)
	{
		poly->vlist[i]->heatvaluesave = poly->vlist[i]->heatvalue;
	}

	for (int i = 0; i < poly->nverts; i++)
	{

			double dt = 1;
			int n = poly->vlist[i]->ntris;
			double sum_eij = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
					{
						c = poly->vlist[i]->tris[j]->corners[k];
						sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
						break;
					}
				}
			}
		if (poly->vlist[i] != poly->tlist[maxindex]->verts[0] && poly->vlist[i] != poly->tlist[minindex]->verts[0])
		{
			double p_sum = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				double wij = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
					{
						c = poly->vlist[i]->tris[j]->corners[k];
						double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
						wij = eij / sum_eij;
						p_sum += wij * (c->n->v->heatvaluesave - poly->vlist[i]->heatvaluesave);
						break;
					}
				}
			}
			poly->vlist[i]->heatvalue = poly->vlist[i]->heatvaluesave + dt * p_sum;
		}
	}

	

	poly->mcalculate();

}


void Polyhedron::heat2(Polyhedron* poly, int maxindex, int minindex, int iterationtime)
{

	poly->tlist[maxindex]->verts[0]->heatvalue = 1.0; //set the maxima
	poly->tlist[minindex]->verts[0]->heatvalue = 0.0; //set the minima

	double dt = 1.;

	dt = dt * iterationtime;  // Set dt as the final time

	int m = poly->nverts;   // Number of unknown elements (equal to number of "vertices")

	//Assembly:
	std::vector<T> tripletList;   // list of non-zeros coefficients the triple value
	Eigen::VectorXd b(m);     // b vector at the right side of the equation

	//BuildProblem: push back coefficient
	
//----------------------------------------------------------------------------------------------------
	b.setZero();

	for (int i = 0; i < m; i++)  // traverse row
	{ 
		int n = poly->vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					c = poly->vlist[i]->tris[j]->corners[k];
					sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
				}
			}
		}
		
		for (int j = 0; j < n; j++)   // traverse column
		{
			Corner* c;
			double wij = 0.0;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					c = poly->vlist[i]->tris[j]->corners[k];
					double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
					wij = eij / sum_eij;
					break;
				}
			}
			if (c->n->v != poly->tlist[maxindex]->verts[0] && c->n->v != poly->tlist[minindex]->verts[0])
			{
				int a;
				//for (int q = 0; q < m; q++)
				//{
				//	if (poly->vlist[q] == c->p->v)
				//	{
				//		a = q;
				//		break;
				//	}
				//
				//}
				a = c->n->v->index;
				tripletList.push_back(T(i, a, -dt * wij));
			}
			else if (c->n->v == poly->tlist[maxindex]->verts[0])
			{
				b[i] += dt * wij;
			}
			else b[i] += 0;

		}
		
		if (poly->vlist[i] == poly->tlist[maxindex]->verts[0])
		{
			b[i] += 1;
		}
		else b[i] += 0;
		tripletList.push_back(T(i, i, 1+dt));

		
	}


//----------------------------------------------------------------------------------------------------
	SpMat A(m, m);            // A matrix at the left side of the equation
	A.setFromTriplets(tripletList.begin(), tripletList.end());  // sent the triple coefficient to matrix

	//Solve problem:

	Eigen::SparseLU<SpMat> solver(A);
	//Eigen::SimplicialCholesky<SpMat> solver(A);
	solver.compute(A);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "Compute Matrix is error" << std::endl;
		return;
	}

	Eigen::VectorXd x = solver.solve(b);

	//----do the saved value----//
	for (int i = 0; i < m; i++)
	{
		if (poly->vlist[i] != poly->tlist[maxindex]->verts[0] && poly->vlist[i] != poly->tlist[minindex]->verts[0])
			poly->vlist[i]->heatvalue = x[i];
	}
	
	poly->mcalculate();
}

void Polyhedron::heat3(Polyhedron* poly, int maxindex, int minindex)
{

	//std::vector<int> array;

	poly->tlist[maxindex]->verts[0]->heatvalue = 1.0; //set the maxima
	poly->tlist[minindex]->verts[0]->heatvalue = 0.0; //set the minima


	for (int i = 0; i < poly->nverts; i++)
	{

		double dt = 1;
		int n = poly->vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
				{
					//for (int p = 0; p < 3; p++)
					//{
					//	if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
					//	{
					//		c = poly->vlist[i]->tris[j]->corners[p];
					//		break;
					//	}
					//}
					//sum_eij += (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
					c = poly->vlist[i]->tris[j]->corners[k];
					sum_eij += (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
				}
			}
		}
		if (poly->vlist[i] != poly->tlist[maxindex]->verts[0] && poly->vlist[i] != poly->tlist[minindex]->verts[0])
		{
			double p_sum = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				double wij = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->corners[k]->v)
					{
						//for (int p = 0; p < 3; p++)
						//{
						//	if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
						//	{
						//		c = poly->vlist[i]->tris[j]->corners[p];
						//		break;
						//	}
						//}
						//double eij = (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
						c = poly->vlist[i]->tris[j]->corners[k];
						double eij = (tan(c->c / 2.0) + tan(c->p->o->p->c / 2.0)) / 2.0;
						wij = eij / sum_eij;
						p_sum += wij * (c->n->v->heatvalue - poly->vlist[i]->heatvalue);
					}
				}
			}
			poly->vlist[i]->heatvalue = poly->vlist[i]->heatvalue + dt * p_sum;
		}
	}



	poly->mcalculate();

}

void Polyhedron::settextcoord()
{
	texcoordindex.push_back(0); //first high
	texcoordindex.push_back(1); //first low
	texcoordindex.push_back(2); //second high
	texcoordindex.push_back(3); //second low

	for (int i = 0; i < nverts; i++)
	{
		if (vlist[i] != tlist[texcoordindex[0]]->verts[0] && vlist[i] != tlist[texcoordindex[1]]->verts[0] )
			//if (vlist[i] != tlist[3]->verts[0] && vlist[i] != tlist[2]->verts[0])
		{
			vlist[i]->tx = 0.0;
		}
		if (vlist[i] != tlist[texcoordindex[2]]->verts[0] && vlist[i] != tlist[texcoordindex[3]]->verts[0])
		{
			vlist[i]->ty = 0.0;
		}
	}
}

void Polyhedron::texmap(Polyhedron* poly)
{
	int maxindex1 = poly->texcoordindex[0];
	int minindex1 = poly->texcoordindex[1];

	poly->tlist[maxindex1]->verts[0]->tx = 1.0; //set the maxima
	poly->tlist[minindex1]->verts[0]->tx = 0.0; //set the minima

	double dt = 0.9;

	for (int i = 0; i < poly->nverts; i++)
	{
		int n = poly->vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
				{
					for (int p = 0; p < 3; p++)
					{
						if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
						{
							c = poly->vlist[i]->tris[j]->corners[p];
							break;
						}
					}
					sum_eij += (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
				}
			}
		}
		if (poly->vlist[i] != poly->tlist[maxindex1]->verts[0] && poly->vlist[i] != poly->tlist[minindex1]->verts[0])
		{
			double p_sum = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				double wij = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
					{
						for (int p = 0; p < 3; p++)
						{
							if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
							{
								c = poly->vlist[i]->tris[j]->corners[p];
								break;
							}
						}
						double eij = (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
						wij =  eij / sum_eij;
						p_sum += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->tx - poly->vlist[i]->tx);
					}
				}
			}
			poly->vlist[i]->tx = poly->vlist[i]->tx + dt * p_sum;
		}
	}


	int maxindex2 = poly->texcoordindex[2];
	int minindex2= poly->texcoordindex[3];

	poly->tlist[maxindex2]->verts[0]->ty = 1.0; //set the maxima
	poly->tlist[minindex2]->verts[0]->ty = 0.0; //set the minima

	for (int i = 0; i < poly->nverts; i++)
	{
		int n = poly->vlist[i]->ntris;
		double sum_eij = 0.0;
		for (int j = 0; j < n; j++)
		{
			Corner* c;
			for (int k = 0; k < 3; k++)
			{
				if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
				{
					for (int p = 0; p < 3; p++)
					{
						if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
						{
							c = poly->vlist[i]->tris[j]->corners[p];
							break;
						}
					}
					sum_eij += (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
				}
			}
		}
		if (poly->vlist[i] != poly->tlist[maxindex2]->verts[0] && poly->vlist[i] != poly->tlist[minindex2]->verts[0])
		{
			double p_sum = 0.0;
			for (int j = 0; j < n; j++)
			{
				Corner* c;
				double wij = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (poly->vlist[i] == poly->vlist[i]->tris[j]->verts[k])
					{
						for (int p = 0; p < 3; p++)
						{
							if (poly->vlist[i]->tris[j]->corners[p]->e->index == poly->vlist[i]->tris[j]->edges[k]->index)
							{
								c = poly->vlist[i]->tris[j]->corners[p];
								break;
							}
						}
						double eij = (tan(c->n->c / 2.0) + tan(c->o->p->c / 2.0)) / 2.0;
						wij =  eij / sum_eij;
						p_sum += wij * (poly->vlist[i]->tris[j]->verts[(k + 1) % 3]->ty - poly->vlist[i]->ty);
					}
				}
			}
			poly->vlist[i]->ty = poly->vlist[i]->ty + dt * p_sum;
		}
	}

}


void Polyhedron::color_t_rainbow() 
{
	double maxt = 1.0, mint = 0.0;
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* v = poly->vlist[i];
		//cold
		if (v->heatvalue == 0.0f) {
			v->r = 1.0f;
			v->g = 0.0f;
			v->b = 1.0f;
		}
		else if (v->heatvalue <= .166 && v->heatvalue > 0.0f) {
			v->r = (.166 - v->heatvalue) / (.166f - 0.0f);
			v->g = 0.0f;
			v->b = 1.0f;
		}
		//cyan to blue
		else if (v->heatvalue <= .336 && v->heatvalue > .166) {
			v->r = 0.0;
			v->g = (v->heatvalue - .166) / (.336f - 0.166f);
			v->b = 1.0f;
		}
		//green to cyan
		else if (v->heatvalue <= .502 && v->heatvalue > .336) {
			v->r = 0.0;
			v->g = 1.0f;
			v->b = (.502 - v->heatvalue) / (.502f - 0.336f);
		}
		//yellow to green
		else if (v->heatvalue <= .668 && v->heatvalue > .502) {
			v->r = (v->heatvalue - .502) / (.668f - 0.502f);
			v->g = 1.0f;
			v->b = 0.0f;
		}
		//red to yellow?
		else if (v->heatvalue <= .834 && v->heatvalue > .668) {
			v->r = 1.0f;
			v->g = (.834 - v->heatvalue) / (.834f - 0.668f);
			v->b = 0.0;
		}
		//white to red
		else if (v->heatvalue < 1.0f && v->heatvalue > .834) {
			v->r = 1.0f;
			v->g = (v->heatvalue - .834f) / (1.0f - 0.834f);
			v->b = (v->heatvalue - .834f) / (1.0f - 0.834f);
		}
		else if (v->heatvalue == 1.0f) {
			v->r = 1.0f;
			v->g = 1.0f;
			v->b = 1.0f;
		}
	}
}

int Polyhedron::saddle(Vertex* vertex)
{
	std::vector<int> simbol;
	int pv = 0;   //p(v0)
	
	Corner* c;
	for (int j = 0; j < 3; j++)
	{
		if (vertex->tris[0]->corners[j]->v == vertex)
		{   
			c = vertex->tris[0]->corners[j]->n;
			for (int k = 0; k < vertex->ntris; k++)
			{
				if (c->v->heatvalue > vertex->heatvalue)
					simbol.push_back(1);
				else simbol.push_back(-1);
				c = c->n->o;
			}
			break;
		}
	}

	for (int i = 0; i < simbol.size()-1; i++)
	{
		if (simbol[i] > simbol[i + 1]) pv += 1;
	}

	if (pv == 0)  return pv;
	else return (pv - 1);
}

void Polyhedron::mcalculate()
{
	int maxima = 0;
	int minima = 0;
	int totalsad = 0;



	for (int i = 0; i < nverts; i++)
	{
		int sad = saddle(vlist[i]);
		int maxflag = 1;
		int minflag = 1;
		
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i]->tris[0]->corners[k]->v == vlist[i])
			{
				c = vlist[i]->tris[0]->corners[k]->n;
				break;
			}
		}

		//double i = c->n->v->heatvalue;
		//double i = c->n->v->heatvalue;
		//double i = c->n->v->heatvalue;
		//double i = c->n->v->heatvalue;


		for (int j = 0; j < vlist[i]->ntris; j++)
		{
			if (c->v->heatvalue > vlist[i]->heatvalue)
			{
				maxflag = 0;
				break;
			}
			c = c->n->o;
		}

		for (int j = 0; j < vlist[i]->ntris; j++)
		{
			if (c->v->heatvalue < vlist[i]->heatvalue)
			{
				minflag = 0;
				break;
			}
			c = c->n->o;
		}

		totalsad += sad;
		if (maxflag == 1) maxima += 1;
		if (minflag == 1) minima += 1;

		
	}

	int m = maxima + minima - totalsad;
	Mf = m;
}

void Polyhedron::create_new_vert(Polyhedron* poly)  //create new added vertices by travese all the edeges
{
	for (int i = 0; i < poly->nedges; i++)
	{
		Vertex* v = new Vertex(0,0,0);
		v->other_props = NULL;
		Corner* c;  // one of the opposite corner
		for (int j = 0; j < 3; j++)
		{
			if (poly->elist[i]->tris[0]->corners[j]->e == poly->elist[i])
			{
				c = poly->elist[i]->tris[0]->corners[j];
				break;
			}
		}

		v->x =0.5 * (c->n->v->x + c->p->v->x); //0.375 * (c->n->v->x + c->p->v->x) + 0.125 * (c->v->x + c->o->v->x);
		v->y =0.5 * (c->n->v->y + c->p->v->y); //0.375 * (c->n->v->y + c->p->v->y) + 0.125 * (c->v->y + c->o->v->y);
		v->z =0.5 * (c->n->v->z + c->p->v->z); //0.375 * (c->n->v->z + c->p->v->z) + 0.125 * (c->v->z + c->o->v->z);

		poly->elist[i]->new_vert = v;
	}
}

void Polyhedron::create_new_ivert(Polyhedron* poly)  //create new added vertices by travese all the edeges
{
	double average_elength = 0;
	for (int s = 0; s < poly->nedges; s++)
	{
		average_elength += poly->elist[s]->length;
	}
	average_elength /= (double)(poly->nedges);

	for (int i = 0; i < poly->nedges; i++)
	{
		Vertex* v = new Vertex(0, 0, 0);
		v->other_props = NULL;
		if (poly->elist[i]->length > average_elength)
		{ 
			Corner* c;  // one of the opposite corner
			for (int j = 0; j < 3; j++)
			{
				if (poly->elist[i]->tris[0]->corners[j]->e == poly->elist[i])
				{
					c = poly->elist[i]->tris[0]->corners[j];
					break;
				}
			}

			v->x = 0.5 * (c->n->v->x + c->p->v->x);//0.375 * (c->n->v->x + c->p->v->x) + 0.125 * (c->v->x + c->o->v->x);
			v->y = 0.5 * (c->n->v->y + c->p->v->y);//0.375 * (c->n->v->y + c->p->v->y) + 0.125 * (c->v->y + c->o->v->y);
			v->z = 0.5 * (c->n->v->z + c->p->v->z);//0.375 * (c->n->v->z + c->p->v->z) + 0.125 * (c->v->z + c->o->v->z);

			poly->elist[i]->new_vert = v;
		}
		else poly->elist[i]->new_vert = NULL;
	}
}

void Polyhedron::irregular_sub(Polyhedron* poly)
{
	create_new_ivert(poly);    //those two step cannot swap the sequence: 
	update_old_vert(poly);    //1. create new points need old vert; 2. update old vert then
	
	int edgevertcount = 0;
	for (int i = 0; i < poly->nedges; i++)
	{
		if (poly->elist[i]->new_vert != NULL)
			edgevertcount += 1;
	}

	max_tris = (int)(4 * poly->max_tris);  /* leave some room for expansion */
	tlist = new Triangle * [max_tris];
	ntris = 0;//= new_poly->max_tris;

	int biggest_number = (int)(poly->max_verts + poly->nedges);
	max_verts = (int)(poly->max_verts + edgevertcount);  //the whole points include NULL one
	vlist = new Vertex * [max_verts];
	nverts = 0;// = new_poly->max_verts;

	

	Vertex** clean_vlist = new Vertex * [(int)(poly->max_verts)];
	for (int i = 0; i < poly->nverts; i++)
	{
		//clean_vlist[i] = new Vertex(poly->vlist[i]->update_vert_pos[0], poly->vlist[i]->update_vert_pos[1], poly->vlist[i]->update_vert_pos[2]);
		clean_vlist[i] = new Vertex(poly->vlist[i]->x, poly->vlist[i]->y, poly->vlist[i]->z);
		clean_vlist[i]->other_props = NULL;
		clean_vlist[i]->edge_index = -1;
		clean_vlist[i]->index = poly->vlist[i]->index;
	}

	Vertex** clean_vlist1 = new Vertex * [(int)(biggest_number - poly->max_verts)];
	for (int i = poly->max_verts; i < biggest_number; i++)
	{
		if (poly->elist[i - poly->max_verts]->new_vert == NULL)
		{ 
			clean_vlist1[i - poly->max_verts] = NULL;
		}
		else
		{
			clean_vlist1[i - poly->max_verts] = new Vertex(poly->elist[i - poly->max_verts]->new_vert->x, poly->elist[i - poly->max_verts]->new_vert->y, poly->elist[i - poly->max_verts]->new_vert->z);
			clean_vlist1[i - poly->max_verts]->other_props = NULL;
			clean_vlist1[i - poly->max_verts]->edge_index = poly->elist[i - poly->max_verts]->index;
		}
		
	}

	// create new vetex list
	for (int i = 0; i < poly->max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
		vlist[i]->other_props = NULL;
		vlist[i] = clean_vlist[i];
		nverts += 1;
	}

	//int flag = poly->max_verts;
	//for (int i = poly->max_verts; i < max_verts; i++)
	//{
	//	if (clean_vlist1[flag - poly->max_verts] != NULL)
	//	{
	//		vlist[i] = new Vertex(0, 0, 0);
	//		vlist[i] = clean_vlist1[flag - poly->max_verts];//poly->elist[i- poly->max_verts]->new_vert;
	//		nverts += 1;
	//	}
	//	flag += 1;
	//}

	int i = poly->max_verts;
	for (int flag = poly->max_verts; flag < biggest_number; flag++)
	{
		if (clean_vlist1[flag - poly->max_verts] != NULL)
		{
			vlist[i] = new Vertex(0, 0, 0);
			vlist[i]->other_props = NULL;
			vlist[i] = clean_vlist1[flag - poly->max_verts];
			nverts += 1;
			i += 1;
		}
	}

	// create new triangle list
	int nowindex = 0;
	for (int i = 0; i < poly->ntris; i ++)    // i is the index of triangle
	{
		int longedgecount = 0;        //long edge count for each triangle
		for (int q = 0; q < 3; q++)
		{
			if (poly->tlist[i]->edges[q]->new_vert != NULL)
				longedgecount += 1;
		}

		//case1: when the number of long edges is 3
		if (longedgecount == 3)
		{
			
			//push new triangle with 3 central points
			tlist[nowindex] = new Triangle;
			tlist[nowindex]->nverts = 3;
			tlist[nowindex]->other_props = NULL;

			for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
			{
				for (int p = poly->max_verts; p < max_verts; p++)
				{
					if (poly->tlist[i]->edges[j]->index == vlist[p]->edge_index)
					{
						tlist[nowindex]->verts[j] = vlist[p];
						break;
					}
				}
			}
			ntris += 1;

			//push 3 new suround triangles
			for (int k = 0; k < 3; k++)         // k is the last 3 triangles
			{
				tlist[nowindex + k + 1] = new Triangle;
				tlist[nowindex + k + 1]->nverts = 3;
				tlist[nowindex + k + 1]->other_props = NULL;

				for (int p = poly->max_verts; p < max_verts; p++)
				{
					if (poly->tlist[i]->edges[k]->index == vlist[p]->edge_index)
					{
						tlist[nowindex + k + 1]->verts[0] = vlist[p];
						tlist[nowindex + k + 1]->verts[0]->ntris = 0;
					}
				}

				for (int p = poly->max_verts; p < max_verts; p++)
				{
					if (poly->tlist[i]->edges[(k + 2) % 3]->index == vlist[p]->edge_index)
					{
						tlist[nowindex + k + 1]->verts[1] = vlist[p];
						tlist[nowindex + k + 1]->verts[1]->ntris = 0;
					}
				}

				if (poly->tlist[i]->edges[k]->verts[0] == poly->tlist[i]->edges[(k + 2) % 3]->verts[0])
				{
					for (int p = 0; p < poly->max_verts; p++)
					{
						if (poly->tlist[i]->edges[k]->verts[0]->index == vlist[p]->index)
						{
							tlist[nowindex + k + 1]->verts[2] = vlist[p];
							break;
						}
					}

				}
				if (poly->tlist[i]->edges[k]->verts[0] == poly->tlist[i]->edges[(k + 2) % 3]->verts[1])
				{
					for (int p = 0; p < poly->max_verts; p++)
					{
						if (poly->tlist[i]->edges[k]->verts[0]->index == vlist[p]->index)
						{
							tlist[nowindex + k + 1]->verts[2] = vlist[p];
							break;
						}
					}
				}
				if (poly->tlist[i]->edges[k]->verts[1] == poly->tlist[i]->edges[(k + 2) % 3]->verts[0])
				{
					for (int p = 0; p < poly->max_verts; p++)
					{
						if (poly->tlist[i]->edges[k]->verts[1]->index == vlist[p]->index)
						{
							tlist[nowindex + k + 1]->verts[2] = vlist[p];
							break;
						}
					}
				}
				if (poly->tlist[i]->edges[k]->verts[1] == poly->tlist[i]->edges[(k + 2) % 3]->verts[1])
				{
					for (int p = 0; p < poly->max_verts; p++)
					{
						if (poly->tlist[i]->edges[k]->verts[1]->index == vlist[p]->index)
						{
							tlist[nowindex + k + 1]->verts[2] = vlist[p];
							break;
						}
					}
				}

				ntris += 1;
			}

			nowindex += (longedgecount+1);
		}
		//case2: when the number of long edges is 2
		if (longedgecount == 2)
		{
			//push new a new triangle (longcenter point, shortcenterpoint, reversepoint)
			tlist[nowindex] = new Triangle;
			tlist[nowindex]->nverts = 3;
			tlist[nowindex]->other_props = NULL;

			//for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
			//{
			//	for (int p = poly->max_verts; p < max_verts; p++)
			//	{
			//		if (poly->tlist[i]->edges[j]->index == vlist[p]->edge_index)
			//		{
			//			tlist[i * 3]->verts[j] = vlist[p];
			//			break;
			//		}
			//	}
			//}

			//scan the edge and get sequence
			Edge** nowedge = new Edge * [3]; //0--normal edge //1-- right long edge //2--left edge
			Edge* swapedge = new Edge;
			int commonindex;

			if (poly->tlist[i]->edges[0]->new_vert != NULL)
			{
				nowedge[1] = poly->tlist[i]->edges[0];
				if (poly->tlist[i]->edges[1]->new_vert != NULL)
				{ 
					nowedge[0] = poly->tlist[i]->edges[2];
					nowedge[2] = poly->tlist[i]->edges[1];
				}
				else
				{
					nowedge[0] = poly->tlist[i]->edges[1];
					nowedge[2] = poly->tlist[i]->edges[2];
				}
			}
			else
			{
				nowedge[0] = poly->tlist[i]->edges[0];
				nowedge[1] = poly->tlist[i]->edges[1];
				nowedge[2] = poly->tlist[i]->edges[2];
			}
			
			for (int m = 0; m < 2; m++)  // nowedge[1].verts
			{
				for (int n = 0; n < 2; n++)  // nowedge[2].verts
				{
					if (nowedge[1]->verts[m]->index == nowedge[2]->verts[n]->index)
						for (int k = 0; k < 3; k++)
						{
							if (poly->tlist[i]->verts[k]->index == nowedge[1]->verts[m]->index)
								commonindex = k;
						}		
				}
			}

			if (nowedge[1]->verts[0] != poly->tlist[i]->verts[(commonindex + 1) % 3]&& nowedge[1]->verts[1] != poly->tlist[i]->verts[(commonindex + 1) % 3])
			{
				swapedge = nowedge[1];
				nowedge[1] = nowedge[2];
				nowedge[2] = swapedge;
			}

			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[commonindex]->index)
					tlist[nowindex]->verts[0] = vlist[p];
			
			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == nowedge[1]->index)
					tlist[nowindex]->verts[1] = vlist[p];
			
			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == nowedge[2]->index)
					tlist[nowindex]->verts[2] = vlist[p];
			ntris += 1;


			tlist[nowindex +1] = new Triangle;
			tlist[nowindex + 1]->nverts = 3;
			tlist[nowindex+1]->other_props = NULL;
			
			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == nowedge[2]->index)
					tlist[nowindex + 1]->verts[0] = vlist[p];
			
			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == nowedge[1]->index)
					tlist[nowindex + 1]->verts[1] = vlist[p];
			
			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[(commonindex+2)%3]->index)
					tlist[nowindex + 1]->verts[2] = vlist[p];
		    ntris += 1;


			tlist[nowindex + 2] = new Triangle;
			tlist[nowindex + 2]->nverts = 3;
			tlist[nowindex+2]->other_props = NULL;
			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == nowedge[1]->index)
					tlist[nowindex + 2]->verts[0] = vlist[p];
			
			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[(commonindex + 1) % 3]->index)
					tlist[nowindex +2]->verts[1] = vlist[p];
			
			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[(commonindex + 2) % 3]->index)
					tlist[nowindex +2]->verts[2] = vlist[p];
			ntris += 1;
			nowindex += (longedgecount+1);

		}
		if (longedgecount == 1)
		{
			//push 2 new triangles
			tlist[nowindex] = new Triangle;
			tlist[nowindex]->nverts = 3;
			tlist[nowindex]->other_props = NULL;


			Edge* longedge = new Edge;
			int oppvertindex;
			for (int j = 0; j < 3; j++)
			{
				if (poly->tlist[i]->edges[j]->new_vert != NULL)
				{
					longedge = poly->tlist[i]->edges[j];

					for (int k = 0; k < 3; k++)
					{
						if (poly->tlist[i]->verts[k]->index != longedge->verts[0]->index && poly->tlist[i]->verts[k]->index != longedge->verts[1]->index)
							oppvertindex = k;
					}
					break;
				}
			}

			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[oppvertindex]->index)
					tlist[nowindex]->verts[0] = vlist[p];

			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[(oppvertindex + 1)%3]->index)
					tlist[nowindex]->verts[1] = vlist[p];

			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == longedge->index)
					tlist[nowindex]->verts[2] = vlist[p];
			ntris += 1;


			tlist[nowindex +1] = new Triangle;
			tlist[nowindex +1]->nverts = 3;
			tlist[nowindex + 1]->other_props = NULL;


			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[oppvertindex]->index)
					tlist[nowindex + 1]->verts[0] = vlist[p];

			for (int p = poly->max_verts; p < max_verts; p++)
				if (vlist[p]->edge_index == longedge->index)
					tlist[nowindex + 1]->verts[1] = vlist[p];

			for (int p = 0; p < poly->max_verts; p++)
				if (vlist[p]->index == poly->tlist[i]->verts[(oppvertindex + 2) % 3]->index)
					tlist[nowindex + 1]->verts[2] = vlist[p];
			ntris += 1;

			nowindex += (longedgecount+1);
		}
		if (longedgecount == 0)
		{
			//push original triangle with 3 central points
			tlist[nowindex] = new Triangle;
			tlist[nowindex]->nverts = 3;
			tlist[nowindex]->other_props = NULL;


			for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
			{
				for (int p = 0; p < poly->max_verts; p++)
				{
					if (poly->tlist[i]->verts[j]->index == vlist[p]->index)
					{
						tlist[nowindex]->verts[j] = vlist[p];
						break;
					}
				}
			}
			ntris += 1;
			nowindex += (longedgecount+1);
		}
	}


	///* get rid of triangles that use the same vertex more than once */
	//
	//for (int i = ntris - 1; i >= 0; i--) {
	//
	//	Triangle* tri = tlist[i];
	//	Vertex* v0 = tri->verts[0];
	//	Vertex* v1 = tri->verts[1];
	//	Vertex* v2 = tri->verts[2];
	//
	//	if (v0 == v1 || v1 == v2 || v2 == v0) {
	//		free(tlist[i]);
	//		ntris--;
	//		tlist[i] = tlist[ntris];
	//	}
	//}
}

void Polyhedron::regular_sub(Polyhedron* poly)
{	
	
	max_tris = (int)(4*poly->max_tris);  /* leave some room for expansion */
	tlist = new Triangle * [max_tris];
	ntris = 0;//= new_poly->max_tris;

	max_verts = (int)(poly->max_verts + poly->nedges);
	vlist = new Vertex * [max_verts];
	nverts = 0;// = new_poly->max_verts;

	create_new_vert(poly);    //those two step cannot swap the sequence: 
	update_old_vert(poly);    //1. create new points need old vert; 2. update old vert then
	
	Vertex ** clean_vlist = new Vertex * [(int)(poly->max_verts)];
	for (int i = 0; i < poly->nverts; i++)
	{
		//clean_vlist[i] = new Vertex(poly->vlist[i]->update_vert_pos[0], poly->vlist[i]->update_vert_pos[1], poly->vlist[i]->update_vert_pos[2]);
		clean_vlist[i] = new Vertex(poly->vlist[i]->x, poly->vlist[i]->y, poly->vlist[i]->z);
		clean_vlist[i]->other_props = NULL;
		clean_vlist[i]->edge_index = -1;
		clean_vlist[i]->index = poly->vlist[i]->index;
	}

	Vertex** clean_vlist1 = new Vertex * [(int)(max_verts- poly->max_verts)];
	for (int i = poly->max_verts; i < max_verts; i++)
	{
		clean_vlist1[i - poly->max_verts] = new Vertex(poly->elist[i - poly->max_verts]->new_vert->x, poly->elist[i - poly->max_verts]->new_vert->y, poly->elist[i - poly->max_verts]->new_vert->z);
		clean_vlist1[i - poly->max_verts]->other_props = NULL;
		clean_vlist1[i - poly->max_verts]->edge_index = poly->elist[i - poly->max_verts]->index;
	}

	// create new vetex list
	for (int i = 0; i < poly->max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
		vlist[i]->other_props = NULL;
		vlist[i] = clean_vlist[i];
		nverts += 1;
	}
	for (int i = poly->max_verts; i < max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
		vlist[i]->other_props = NULL;
		vlist[i] = clean_vlist1[i - poly->max_verts];//poly->elist[i- poly->max_verts]->new_vert;
		nverts += 1;
	}

	// create new triangle list
	for (int i = 0; i < poly->ntris; i++)    // i is the index of triangle
	{
		//Triangle* f = poly->tlist[i];

		//push new triangle with 3 central points
		tlist[i * 4] = new Triangle;
		tlist[i * 4]->nverts = 3;
		tlist[i * 4]->other_props = NULL;


		//for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
		//{
		//	tlist[i * 4]->verts[j] = poly->tlist[i]->edges[j]->new_vert;
		//}
		//ntris += 1;

		for (int j = 0; j < 3; j++)         // j is the index of the vertex of the center triangle   
		{
			for (int p = poly->max_verts; p < max_verts; p++)
			{
				if (poly->tlist[i]->edges[j]->index == vlist[p]->edge_index)
				{
					tlist[i * 4]->verts[j] = vlist[p];
					break;
				}
			}
		}
		ntris += 1;
		
		//push 3 new suround triangles
		for (int k = 0; k < 3; k++)         // k is the last 3 triangles
		{
			tlist[i * 4 + k + 1] = new Triangle;
			tlist[i * 4 + k + 1]->nverts = 3;
			tlist[i * 4 + k + 1]->other_props = NULL;


			for (int p = poly->max_verts; p < max_verts; p++)
			{
				if (poly->tlist[i]->edges[k]->index == vlist[p]->edge_index)
				{
					tlist[i * 4 + k + 1]->verts[0] = vlist[p];
					tlist[i * 4 + k + 1]->verts[0]->ntris = 0;
				}
			}

			for (int p = poly->max_verts; p < max_verts; p++)
			{
				if (poly->tlist[i]->edges[(k + 2) % 3]->index == vlist[p]->edge_index)
				{
					tlist[i * 4 + k + 1]->verts[1] = vlist[p];
					tlist[i * 4 + k + 1]->verts[1]->ntris = 0;
				}
			}

			if (poly->tlist[i]->edges[k]->verts[0] == poly->tlist[i]->edges[(k + 2) % 3]->verts[0])
			{
				for (int p = 0; p < poly->max_verts; p++)
				{
					if (poly->tlist[i]->edges[k]->verts[0]->index == vlist[p]->index)
					{
						tlist[i * 4 + k + 1]->verts[2] = vlist[p];
						break;
					}
				}
				
			}
			if (poly->tlist[i]->edges[k]->verts[0] == poly->tlist[i]->edges[(k + 2) % 3]->verts[1])
			{
				for (int p = 0; p < poly->max_verts; p++)
				{
					if (poly->tlist[i]->edges[k]->verts[0]->index == vlist[p]->index)
					{
						tlist[i * 4 + k + 1]->verts[2] = vlist[p];
						break;
					}
				}
			}
			if (poly->tlist[i]->edges[k]->verts[1] == poly->tlist[i]->edges[(k + 2) % 3]->verts[0])
			{
				for (int p = 0; p < poly->max_verts; p++)
				{
					if (poly->tlist[i]->edges[k]->verts[1]->index == vlist[p]->index)
					{
						tlist[i * 4 + k + 1]->verts[2] = vlist[p];
						break;
					}
				}
			}
			if (poly->tlist[i]->edges[k]->verts[1] == poly->tlist[i]->edges[(k + 2) % 3]->verts[1])
			{
				for (int p = 0; p < poly->max_verts; p++)
				{
					if (poly->tlist[i]->edges[k]->verts[1]->index == vlist[p]->index)
					{
						tlist[i * 4 + k + 1]->verts[2] = vlist[p];
						break;
					}
				}
			}
			//tlist[i * 4 + k + 1]->verts[2]->ntris = 0;

			//tlist[i * 4 + k + 1] = new Triangle;
			//tlist[i * 4 + k + 1]->nverts = 3;
			//tlist[i * 4 + k + 1]->verts[0] = poly->tlist[i]->edges[k]->new_vert;
			//tlist[i * 4 + k + 1]->verts[0]->ntris = 0;
			//
			//tlist[i * 4 + k + 1]->verts[1] = poly->tlist[i]->edges[(k + 2) % 3]->new_vert;
			//tlist[i * 4 + k + 1]->verts[1]->ntris = 0;
			//
			//if (poly->tlist[i]->edges[k]->verts[0] == poly->tlist[i]->edges[(k + 2) % 3]->verts[0])
			//	tlist[i * 4 + k + 1]->verts[2] = new Vertex(poly->tlist[i]->edges[k]->verts[0]->x, poly->tlist[i]->edges[k]->verts[0]->y, poly->tlist[i]->edges[k]->verts[0]->z);//poly->tlist[i]->edges[k]->verts[0];
			//if (poly->tlist[i]->edges[k]->verts[0] == poly->tlist[i]->edges[(k + 2) % 3]->verts[1])
			//	tlist[i * 4 + k + 1]->verts[2] = new Vertex(poly->tlist[i]->edges[k]->verts[0]->x, poly->tlist[i]->edges[k]->verts[0]->y, poly->tlist[i]->edges[k]->verts[0]->z);//poly->tlist[i]->edges[k]->verts[0];
			//if (poly->tlist[i]->edges[k]->verts[1] == poly->tlist[i]->edges[(k + 2) % 3]->verts[0])
			//	tlist[i * 4 + k + 1]->verts[2] = new Vertex(poly->tlist[i]->edges[k]->verts[1]->x, poly->tlist[i]->edges[k]->verts[1]->y, poly->tlist[i]->edges[k]->verts[1]->z);//poly->tlist[i]->edges[k]->verts[1];
			//if (poly->tlist[i]->edges[k]->verts[1] == poly->tlist[i]->edges[(k + 2) % 3]->verts[1])
			//	tlist[i * 4 + k + 1]->verts[2] = new Vertex(poly->tlist[i]->edges[k]->verts[1]->x, poly->tlist[i]->edges[k]->verts[1]->y, poly->tlist[i]->edges[k]->verts[1]->z);//poly->tlist[i]->edges[k]->verts[1];
			//tlist[i * 4 + k + 1]->verts[2]->ntris = 0;

			////Triangle* t = tlist[i * 4 + k + 1];  //for debug point
			////Vertex* v = new Vertex(poly->tlist[i]->edges[k]->verts[1]->x, poly->tlist[i]->edges[k]->verts[1]->y, poly->tlist[i]->edges[k]->verts[1]->z);
			////tlist[i * 4 + k + 1]->verts[1] = v;
			
			ntris += 1;
		}
	}

	/* get rid of triangles that use the same vertex more than once */

	for (int i = ntris - 1; i >= 0; i--) {

		Triangle* tri = tlist[i];
		Vertex* v0 = tri->verts[0];
		Vertex* v1 = tri->verts[1];
		Vertex* v2 = tri->verts[2];

		if (v0 == v1 || v1 == v2 || v2 == v0) {
			if (tlist[i])
			{
				free(tlist[i]);
				tlist[i] = NULL;
			}
			ntris--;
			tlist[i] = tlist[ntris];
		}
	}
	//for (int i = 0; i < new_poly->ntris; i++)
	//{
	//	std::cout << "triangle" << i << std::endl;
	//	for (int j = 0; j < 3; j++)
	//	{
	//		std::cout << new_poly->tlist[i]->verts[j]->x << "," << new_poly->tlist[i]->verts[j]->y << "," << new_poly->tlist[i]->verts[j]->z << std::endl;
	//	}
	//}
}

void Polyhedron::Silhouette1(Polyhedron* poly)
{		
	glm::mat4x4 Q;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			Q[i][j] = rotmat[i][j];

	glm::mat4x4 QINVE = glm::inverse(Q);
	glm::vec3 eye_position = vec3(0.0, 0.0, -10);
	glm::vec4 eye_pos = glm::vec4(eye_position, 1.0);
	eye_position = vec3((QINVE * eye_pos).x, (QINVE * eye_pos).y, (QINVE * eye_pos).z);

	for (int i = 0; i < poly->ntris; i++)
	{
		double meanvertx = (poly->tlist[i]->verts[0]->x+ poly->tlist[i]->verts[1]->x+ poly->tlist[i]->verts[2]->x)/3.0;
		double meanverty = (poly->tlist[i]->verts[0]->y + poly->tlist[i]->verts[1]->y + poly->tlist[i]->verts[2]->y) / 3.0;
		double meanvertz = (poly->tlist[i]->verts[0]->z + poly->tlist[i]->verts[1]->z + poly->tlist[i]->verts[2]->z) / 3.0;

		glm::vec3 facecenter = glm::vec3(meanvertx, meanverty, meanvertz);
		//glm::vec3 facecenter1 = glm::vec3(poly->tlist[i]->verts[0]->x, poly->tlist[i]->verts[0]->y, poly->tlist[i]->verts[0]->z);
		//glm::vec3 facecenter2 = glm::vec3(poly->tlist[i]->verts[1]->x, poly->tlist[i]->verts[1]->y, poly->tlist[i]->verts[1]->z);
		//glm::vec3 facecenter3 = glm::vec3(poly->tlist[i]->verts[2]->x, poly->tlist[i]->verts[2]->y, poly->tlist[i]->verts[2]->z);
		
		glm::vec3 ray = normalize(eye_position - facecenter);
		//glm::vec3 ray1 = normalize(eye_position - facecenter1);
		//glm::vec3 ray2 = normalize(eye_position - facecenter2);
		//glm::vec3 ray3 = normalize(eye_position - facecenter3);

		glm::vec3 normal = normalize(glm::vec3(poly->tlist[i]->normal.x, poly->tlist[i]->normal.y, poly->tlist[i]->normal.z));

		if (glm::dot(normal,ray) < 0.) poly->tlist[i]->back_flag = -1;
		else poly->tlist[i]->back_flag = 1;
		 
		//--------------for fun---------------//
		glm::vec3 light = vec3(5.5, 0.0, 0.0);
		glm::vec4 light_pos = glm::vec4(light, 1.0);
		light_pos = QINVE * light_pos;
		light = vec3(light_pos.x, light_pos.y, light_pos.z);
		//if (glm::dot(normal, normalize(light - facecenter)) < 0.) poly->tlist[i]->shadow_flag = -1;
		 poly->tlist[i]->shadow_flag = glm::dot(normal, normalize(light - facecenter));
		//------------------------------------------

	}
	//glBegin(GL_LINES);
	//glVertex3f(0, 0, 0);
	//glVertex3f(eye_position.x, eye_position.y, eye_position.z);
	//glEnd();

}

void Polyhedron::Silhouette2(Polyhedron* poly)
{

	glm::mat4x4 Q;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			Q[i][j] = rotmat[i][j];

	glm::mat4x4 QINVE = glm::inverse(Q);
	glm::vec3 eye_position = vec3(0.0, 0.0, -10.0);
	glm::vec4 eye_pos = glm::vec4(eye_position, 1.0);
	eye_position = vec3((QINVE * eye_pos).x, (QINVE * eye_pos).y, (QINVE * eye_pos).z);

	for (int i = 0; i < poly->nverts; i++)
	{
		glm::vec3 vertlocation = glm::vec3(poly->vlist[i]->x, poly->vlist[i]->y, poly->vlist[i]->z);
		glm::vec3 ray = eye_position - vertlocation;

		glm::vec3 normal = glm::vec3(poly->vlist[i]->normal.entry[0], poly->vlist[i]->normal.entry[1], poly->vlist[i]->normal.entry[2]);

		//if (glm::dot(ray, normal) < 0) poly->vlist[i]->back_flag = -1;
		//else poly->vlist[i]->back_flag = 1;

		//--------------threshold---------------
		poly->vlist[i]->back_flag = glm::dot(ray, normal);
		//--------------------------------------
	}

	for (int i = 0; i < poly->nedges; i++)
	{
		poly->elist[i]->zeroflag = 0;
	}
	//--------------new--------------------
	for (int i = 0; i < poly->nedges; i++)
	{
		if (poly->elist[i]->verts[0]->back_flag * poly->elist[i]->verts[1]->back_flag < 0)
		{
			double w0 = abs(poly->elist[i]->verts[0]->back_flag) / (abs(poly->elist[i]->verts[0]->back_flag) + abs(poly->elist[i]->verts[1]->back_flag));
			double w1 = abs(poly->elist[i]->verts[1]->back_flag) / (abs(poly->elist[i]->verts[0]->back_flag) + abs(poly->elist[i]->verts[1]->back_flag));
			poly->elist[i]->zerovecpoint[0] = w0 * poly->elist[i]->verts[0]->x + w1 * poly->elist[i]->verts[1]->x;
			poly->elist[i]->zerovecpoint[1] = w0 * poly->elist[i]->verts[0]->y + w1 * poly->elist[i]->verts[1]->y;
			poly->elist[i]->zerovecpoint[2] = w0 * poly->elist[i]->verts[0]->z + w1 * poly->elist[i]->verts[1]->z;
			poly->elist[i]->zeroflag = -1;
		}
	}
}

void Polyhedron::CalculateAmixed()
{
	for (int i = 0; i < nverts; i++)
	{
		vlist[i]->Amixed = 0;
	}

	for (int i = 0; i < nverts; i++)
	{
		Corner *c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i]->tris[0]->corners[k]->v == vlist[i])
			{
				c = vlist[i]->tris[0]->corners[k];
				break;
			}
		}

		for (int j = 0; j <vlist[i]->ntris; j++)
		{
			if (c->c <= PI / 2. && c->n->c <= PI / 2. && c->n->n->c <= PI / 2.)    //non-obtuse
			{
				double distance = pow((vlist[i]->x - c->n->v->x),2)+ pow((vlist[i]->y - c->n->v->y), 2) + pow((vlist[i]->z - c->n->v->z), 2);
				vlist[i]->Amixed += 0.125 * (cos(c->p->c) / sin(c->p->c) + cos(c->p->o->c) / sin(c->p->o->c)) * distance;
			}
			else if (c->c > PI / 2)
			{
				double p = 0.5*(c->e->length+ c->n->e->length + c->n->n->e->length);  //Halen formula
				double S = sqrt(p * (p - c->e->length) * (p - c->n->e->length) * (p - c->n->n->e->length));
				vlist[i]->Amixed += S / 2.;
			}
			else
			{
				double p = 0.5 * (c->e->length + c->n->e->length + c->n->n->e->length);  //Halen formula
				double S = sqrt(p * (p - c->e->length) * (p - c->n->e->length) * (p - c->n->n->e->length));
				vlist[i]->Amixed += S / 4.;
			}		
			c = c->p->o->p;
		}
	}
}

void Polyhedron::MeanCurvature()
{
	maxmeancurv = -10.;
	minmeancurv = 10;

	CalculateAmixed();
	for (int i = 0; i < nverts; i++)
	{
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i]->tris[0]->corners[k]->v == vlist[i])
			{
				c = vlist[i]->tris[0]->corners[k];
				break;
			}
		}
		vlist[i]->meancurvaturevalue[0] = 0;
		vlist[i]->meancurvaturevalue[1] = 0;
		vlist[i]->meancurvaturevalue[2] = 0;
		for (int j = 0; j < vlist[i]->ntris; j++)
		{
			double angle;

			if (sin(c->p->c) == 0. && sin(c->p->o->c) == 0.)
				angle = 2;
			else if (sin(c->p->c) == 0. && sin(c->p->o->c) != 0.)
				angle = (1 + cos(c->p->o->c) / sin(c->p->o->c));
			else if (sin(c->p->c) != 0. && sin(c->p->o->c) == 0. )
				angle = (cos(c->p->c) / sin(c->p->c) + 1);
			else
				angle = (cos(c->p->c) / sin(c->p->c) + cos(c->p->o->c) / sin(c->p->o->c));

			vlist[i]->meancurvaturevalue[0] += (1. / (2. * vlist[i]->Amixed))*angle* (vlist[i]->x - c->n->v->x);
			vlist[i]->meancurvaturevalue[1] += (1. / (2. * vlist[i]->Amixed)) * angle * (vlist[i]->y - c->n->v->y);
			vlist[i]->meancurvaturevalue[2] += (1. / (2. * vlist[i]->Amixed)) * angle * (vlist[i]->z - c->n->v->z);
			c = c->p->o->p;
		}
		//vlist[i]->meancurvaturemagnitude = sqrt(pow(vlist[i]->meancurvaturevalue[0], 2) + pow(vlist[i]->meancurvaturevalue[1], 2) + pow(vlist[i]->meancurvaturevalue[2], 2)) / 2;
		vlist[i]->meancurvaturemagnitude = 0.5 * dot(glm::vec3(vlist[i]->meancurvaturevalue[0], vlist[i]->meancurvaturevalue[1], vlist[i]->meancurvaturevalue[2]), 
													 glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
		if (vlist[i]->meancurvaturemagnitude > maxmeancurv) maxmeancurv = vlist[i]->meancurvaturemagnitude;
		if (vlist[i]->meancurvaturemagnitude < minmeancurv) minmeancurv = vlist[i]->meancurvaturemagnitude;
	}
}

void Polyhedron::GaussCurvature()
{
	maxgausscurv = -10.;
	mingausscurv = 10.;
	CalculateAmixed();
	for (int i = 0; i < nverts; i++)
	{
		double totalangle = 0;
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i]->tris[0]->corners[k]->v ==vlist[i])
			{
				c =vlist[i]->tris[0]->corners[k];
				break;
			}
		}
		for (int j = 0; j <vlist[i]->ntris; j++)
		{
			totalangle += c->c;
			c = c->p->o->p;
		}
		vlist[i]->gausscurvaturevalue = (2.0 * PI - totalangle) / vlist[i]->Amixed;
		if (vlist[i]->gausscurvaturevalue > maxgausscurv) maxgausscurv = vlist[i]->gausscurvaturevalue;
		if (vlist[i]->gausscurvaturevalue < mingausscurv) mingausscurv = vlist[i]->gausscurvaturevalue;
	}
}

void Polyhedron::PrincipalCurvature()
{
	CalculateAmixed();
	for (int i = 0; i < nverts; i++)
	{
		double totalangle = 0;
		Corner* c;
		for (int k = 0; k < 3; k++)
		{
			if (vlist[i]->tris[0]->corners[k]->v == vlist[i])
			{
				c = vlist[i]->tris[0]->corners[k];
				break;
			}
		}

		Eigen::MatrixXd A(vlist[i]->ntris+1, 3);
		Eigen::VectorXd b(vlist[i]->ntris + 1);     // b vector at the right side of the equation
		b.setZero();

		glm::vec3 n = normalize(glm::vec3(vlist[i]->normal.x, vlist[i]->normal.y, vlist[i]->normal.z));
		//glm::vec3 n = normalize(glm::vec3(vlist[i]->meancurvaturevalue[0], vlist[i]->meancurvaturevalue[1], vlist[i]->meancurvaturevalue[2]));
		//can n be the mean curvature normal or vertex normal?

		glm::vec3 u1;
		glm::vec3 u2;

		glm::vec3 xjxi = glm::vec3((c->n->v->x - vlist[i]->x), (c->n->v->y - vlist[i]->y), (c->n->v->z - vlist[i]->z));
		u1 = xjxi - dot(xjxi, n) * n;

		//glm::vec3 check = ( u1 / length(u1));

		u1 = normalize(u1 / length(u1));
		u2 = normalize(cross(u1,n)); //normalize(cross(n,u1));

		//double check = dot(u1, n);

		vlist[i]->TangentVector1[0] = u1.x;
		vlist[i]->TangentVector1[1] = u1.y;
		vlist[i]->TangentVector1[2] = u1.z;

		vlist[i]->TangentVector2[0] = u2.x;
		vlist[i]->TangentVector2[1] = u2.y;
		vlist[i]->TangentVector2[2] = u2.z;


		for (int j = 0; j < vlist[i]->ntris; j++)
		{			
			glm::vec3 d1; 
			glm::vec3 d2;
			
			glm::vec3 xjmxi = glm::vec3((c->n->v->x - vlist[i]->x), (c->n->v->y - vlist[i]->y), (c->n->v->z - vlist[i]->z));
			d1 = xjmxi - dot(xjmxi, n) * n;
			d1 = normalize(d1 / length(d1));

			//std::cout << "The is:\n"
			//	<< d1.x<< ","<< d1.y << ","<<d1.z << "," << d2.x << "," << d2.y << "," << d2.z << std::endl;

			double kijN = 2. * dot((-xjmxi), n) / pow(length(-xjmxi),2);
			double bj = kijN;

			double m = dot(d1, u1);   //d1 project on u1 multiple length(u1) = 1    => d1's project length on u1
			double n = dot(d1, u2);

			double mm = m*m;
			double mn = m*n;
			double nn = n*n;

			//Assembly:

			//BuildProblem: push back coefficient
			b[j] = bj;
			
			A(j, 0) = mm;
			A(j, 1) = 2. * mn;
			A(j, 2) = nn;
			
			c = c->p->o->p;
		}

		double kH = vlist[i]->meancurvaturemagnitude;
		double bnp1 = 2. * kH;
		b[vlist[i]->ntris] = bnp1;
		A(vlist[i]->ntris, 0) = 1.;
		A(vlist[i]->ntris, 1) = 1.;
		A(vlist[i]->ntris, 2) = 0.;
		//std::cout << "Here is the matrix A:\n" << A << std::endl;
		//std::cout << "Here is the right hand side b:\n" << b << std::endl;

		Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
		//std::cout << "The solution using the QR decomposition is:\n"
		//	<< A.colPivHouseholderQr().solve(b) << std::endl;

		vlist[i]->Tensor[0] = x[0];
		vlist[i]->Tensor[1] = x[1];
		vlist[i]->Tensor[2] = x[2];

		//-----------------Principle vector getting----------------------
		//Eigen::EigenSolver<Eigen::MatrixXd> es(M);
		//Eigen::VectorXcd v1 = es.eigenvectors().col(0);
		//Eigen::VectorXcd v2 = es.eigenvectors().col(1);
		//
		////std::cout << "The principle curvature vector is:\n"
		////	<< es.eigenvectors() << std::endl;
		////
		//
		//double x1 = v1[0].real();
		//double y1 = v1[1].real();
		//
		//double x2 = v2[0].real();
		//double y2 = v2[1].real();
		//
		//
		//
		//glm::vec3 principal1 = normalize((float)x1 * u1 + (float)y1 * u2);
		//glm::vec3 principal2 = normalize((float)x2 * u1 + (float)y2 * u2);
		//
		//vlist[i]->principlecurvaturevector1[0] = principal1.x;
		//vlist[i]->principlecurvaturevector1[1] = principal1.y;
		//vlist[i]->principlecurvaturevector1[2] = principal1.z;
		//
		//vlist[i]->principlecurvaturevector2[0] = principal2.x;
		//vlist[i]->principlecurvaturevector2[1] = principal2.y;
		//vlist[i]->principlecurvaturevector2[2] = principal2.z;

		CalculatePrincipalCurvature(i);

	}
}


void Polyhedron::CalculateTriPrincipalCurvature(int tindex)
{
	Eigen::MatrixXd M(2, 2);
	M(0, 0) = tlist[tindex]->Tensor2D[0];
	M(0, 1) = tlist[tindex]->Tensor2D[1];
	M(1, 0) = tlist[tindex]->Tensor2D[1];
	M(1, 1) = tlist[tindex]->Tensor2D[2];

	Eigen::EigenSolver<Eigen::MatrixXd> es(M);
	//Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);

	Eigen::VectorXcd v1 = es.eigenvectors().col(0);
	Eigen::VectorXcd v2 = es.eigenvectors().col(1);

	//std::cout << "The principle curvature vector is:\n"
	//	<< es.eigenvectors() << std::endl;
	//

	double x1 = v1[0].real();
	double y1 = v1[1].real();

	double x2 = v2[0].real();
	double y2 = v2[1].real();

	glm::vec3 u1 = glm::vec3(tlist[tindex]->Tangent1[0], tlist[tindex]->Tangent1[1], tlist[tindex]->Tangent1[2]);
	glm::vec3 u2 = glm::vec3(tlist[tindex]->Tangent2[0], tlist[tindex]->Tangent2[1], tlist[tindex]->Tangent2[2]);

	//double check = dot(glm::vec2(x1, y1), glm::vec2(x2, y2));

	glm::vec3 principal1 = normalize((float)x1 * u1 + (float)y1 * u2);
	glm::vec3 principal2 = normalize((float)x2 * u1 + (float)y2 * u2);

	//double check = dot(u1, u2);

	if (dot(principal1, u1) >= 0)
	{
		tlist[tindex]->principleTricurvaturevector1[0] = principal1.x;
		tlist[tindex]->principleTricurvaturevector1[1] = principal1.y;
		tlist[tindex]->principleTricurvaturevector1[2] = principal1.z;
		  						
		tlist[tindex]->principleTricurvaturevector2[0] = principal2.x;
		tlist[tindex]->principleTricurvaturevector2[1] = principal2.y;
		tlist[tindex]->principleTricurvaturevector2[2] = principal2.z;
	}							
	else						
	{									
		tlist[tindex]->principleTricurvaturevector2[0] = principal1.x;
		tlist[tindex]->principleTricurvaturevector2[1] = principal1.y;
		tlist[tindex]->principleTricurvaturevector2[2] = principal1.z;
					
		tlist[tindex]->principleTricurvaturevector1[0] = principal2.x;
		tlist[tindex]->principleTricurvaturevector1[1] = principal2.y;
		tlist[tindex]->principleTricurvaturevector1[2] = principal2.z;

	}
}
void Polyhedron::CalculatePrincipalCurvature(int vindex)
{
	Eigen::MatrixXd M(2, 2); 
	M(0,0) = vlist[vindex]->Tensor[0]; 
	M(0,1) = vlist[vindex]->Tensor[1];
	M(1,0) = vlist[vindex]->Tensor[1];
	M(1,1) = vlist[vindex]->Tensor[2];

	Eigen::EigenSolver<Eigen::MatrixXd> es(M);
	//Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);

	Eigen::VectorXcd v1 = es.eigenvectors().col(0);
	Eigen::VectorXcd v2 = es.eigenvectors().col(1);

	//std::cout << "The principle curvature vector is:\n"
	//	<< es.eigenvectors() << std::endl;
	//

	double x1 = v1[0].real();
	double y1 = v1[1].real();
					 
	double x2 = v2[0].real();
	double y2 = v2[1].real();

	glm::vec3 u1 = glm::vec3(vlist[vindex]->TangentVector1[0], vlist[vindex]->TangentVector1[1], vlist[vindex]->TangentVector1[2]);
	glm::vec3 u2 = glm::vec3(vlist[vindex]->TangentVector2[0], vlist[vindex]->TangentVector2[1], vlist[vindex]->TangentVector2[2]);

	//double check = dot(glm::vec2(x1, y1), glm::vec2(x2, y2));

	glm::vec3 principal1 = normalize((float)x1 * u1 + (float)y1 * u2);
	glm::vec3 principal2 = normalize((float)x2 * u1 + (float)y2 * u2);

	//double check = dot(u1, u2);

	if (dot(principal1, u1) >= 0)
	{
		vlist[vindex]->principlecurvaturevector1[0] = principal1.x;
		vlist[vindex]->principlecurvaturevector1[1] = principal1.y;
		vlist[vindex]->principlecurvaturevector1[2] = principal1.z;

		vlist[vindex]->principlecurvaturevector2[0] = principal2.x;
		vlist[vindex]->principlecurvaturevector2[1] = principal2.y;
		vlist[vindex]->principlecurvaturevector2[2] = principal2.z;
	}
	else
	{
	
		vlist[vindex]->principlecurvaturevector2[0] = principal1.x;
		vlist[vindex]->principlecurvaturevector2[1] = principal1.y;
		vlist[vindex]->principlecurvaturevector2[2] = principal1.z;

		vlist[vindex]->principlecurvaturevector1[0] = principal2.x;
		vlist[vindex]->principlecurvaturevector1[1] = principal2.y;
		vlist[vindex]->principlecurvaturevector1[2] = principal2.z;
	
	}

}	

void colorMap(double colorvalue, double &r, double &g, double &b)
{
	const double rone = 0.8;
	const double gone = 1.0;
	const double bone = 1.0;
	double x = colorvalue;
	x = (colorvalue < 0 ? 0 : (x > 1 ? 1 : x));
	if (x < 1. / 8.) 
	{
		r = 0; 
		g = 0; 
		b = bone * (0.5 + (x) / (1. / 8.) * 0.5);
	}
	else if (x < 3. / 8.)
	{
		r = 0; 
		g = gone * (x - 1. / 8.) / (3. / 8. - 1. / 8.); 
		b = bone;
	}
	else if (x < 5. / 8.)
	{
		r = rone * (x - 3. / 8.) / (5. / 8. - 3. / 8.);
		g = gone;
		b = gone * (x - 3. / 8.) / (5. / 8. - 3. / 8.);
	}
	else if (x < 7. / 8.)
	{
		r = rone;
		g = gone - (x - 5. / 8.) / (7. / 8. - 5. / 8.);
		b = 0;
	}
	else 
	{
		r = rone * (x - 7. / 8.) / (1. - 7. / 8.)*0.5;
		g = 0;
		b = 0;
	}

}

//void colorMap(double colorvalue, double& r, double& g, double& b)
//{
//	double maxt = 1.0, mint = 0.0;
//
//		//cold
//		if (colorvalue == 0.0f) {
//			r = 1.0f;
//			g = 0.0f;
//			b = 1.0f;
//		}
//		else if (colorvalue <= .166 && colorvalue > 0.0f) {
//			r = (.166 - colorvalue) / (.166f - 0.0f);
//			g = 0.0f;
//			b = 1.0f;
//		}
//		//cyan to blue
//		else if (colorvalue <= .336 && colorvalue > .166) {
//			r = 0.0;
//			g = (colorvalue - .166) / (.336f - 0.166f);
//			b = 1.0f;
//		}
//		//green to cyan
//		else if (colorvalue <= .502 && colorvalue > .336) {
//			r = 0.0;
//			g = 1.0f;
//			b = (.502 - colorvalue) / (.502f - 0.336f);
//		}
//		//yellow to green
//		else if (colorvalue <= .668 && colorvalue > .502) {
//			r = (colorvalue - .502) / (.668f - 0.502f);
//			g = 1.0f;
//			b = 0.0f;
//		}
//		//red to yellow?
//		else if (colorvalue <= .834 && colorvalue > .668) {
//			r = 1.0f;
//			g = (.834 - colorvalue) / (.834f - 0.668f);
//			b = 0.0;
//		}
//		//white to red
//		else if (colorvalue < 1.0f && colorvalue > .834) {
//			r = 1.0f;
//			g = (colorvalue - .834f) / (1.0f - 0.834f);
//			b = (colorvalue - .834f) / (1.0f - 0.834f);
//		}
//		else if (colorvalue == 1.0f) {
//			r = 1.0f;
//			g = 1.0f;
//			b = 1.0f;
//		}
//	
//}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];
	int count_valence = 0;
	double count_angle = 0;

  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (sflag != 1) glShadeModel(GL_SMOOTH);
	else glShadeModel(GL_FLAT);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1); 
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if (tflag != 1)
	{ 
		for (i = 0; i < this_poly->nedges; i++)
		{
			glBegin(GL_LINES);
			glColor4f(0, 0, 0,0.5);
			glVertex3d(this_poly->elist[i]->verts[0]->x, this_poly->elist[i]->verts[0]->y, this_poly->elist[i]->verts[0]->z);
			glColor4f(0, 0, 0, 0.5);
			glVertex3d(this_poly->elist[i]->verts[1]->x, this_poly->elist[i]->verts[1]->y, this_poly->elist[i]->verts[1]->z);
			glEnd();
		}
	}


	if (silhouetteflag != 1)
	{
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(1.4);
		poly->Silhouette1(this_poly);
		for (i = 0; i < this_poly->nedges; i++)
		{
			if (this_poly->elist[i]->tris[0]->back_flag * this_poly->elist[i]->tris[1]->back_flag < 0)
			{
				glBegin(GL_LINES);
				glColor3f(0, 0, 0);
				glVertex3d(this_poly->elist[i]->verts[0]->x, this_poly->elist[i]->verts[0]->y, this_poly->elist[i]->verts[0]->z);
				glColor3f(0, 0, 0);
				glVertex3d(this_poly->elist[i]->verts[1]->x, this_poly->elist[i]->verts[1]->y, this_poly->elist[i]->verts[1]->z);
				glEnd();
			}
		}
		//--------------for fun (shadows)-----------------//
		//for (i = 0; i < this_poly->ntris; i++)
		//{
		//	if (this_poly->tlist[i]->shadow_flag <= 1.0 && this_poly->tlist[i]->shadow_flag > 0.3)
		//	{
		//		glm::vec3 normal = glm::vec3(this_poly->tlist[i]->normal.x, this_poly->tlist[i]->normal.y, this_poly->tlist[i]->normal.z)/(float)200;
		//		glBegin(GL_TRIANGLES);
		//		glColor3f(1,1,1);
		//		glVertex3d(this_poly->tlist[i]->verts[0]->x+normal.x, this_poly->tlist[i]->verts[0]->y + normal.y, this_poly->tlist[i]->verts[0]->z + normal.z);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[1]->x + normal.x, this_poly->tlist[i]->verts[1]->y + normal.y, this_poly->tlist[i]->verts[1]->z + normal.z);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[2]->x + normal.x, this_poly->tlist[i]->verts[2]->y + normal.y, this_poly->tlist[i]->verts[2]->z + normal.z);
		//		glEnd();
		//		
		//		normal = glm::vec3(this_poly->tlist[i]->normal.x, this_poly->tlist[i]->normal.y, this_poly->tlist[i]->normal.z) / (float)80;
		//		glBegin(GL_TRIANGLES);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[0]->x + normal.x, this_poly->tlist[i]->verts[0]->y + normal.y, this_poly->tlist[i]->verts[0]->z + normal.z);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[1]->x + normal.x, this_poly->tlist[i]->verts[1]->y + normal.y, this_poly->tlist[i]->verts[1]->z + normal.z);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[2]->x + normal.x, this_poly->tlist[i]->verts[2]->y + normal.y, this_poly->tlist[i]->verts[2]->z + normal.z);
		//		glEnd();
		//
		//		glBegin(GL_TRIANGLES);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[0]->x, this_poly->tlist[i]->verts[0]->y, this_poly->tlist[i]->verts[0]->z);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[1]->x, this_poly->tlist[i]->verts[1]->y, this_poly->tlist[i]->verts[1]->z);
		//		glColor3f(1, 1, 1);
		//		glVertex3d(this_poly->tlist[i]->verts[2]->x, this_poly->tlist[i]->verts[2]->y, this_poly->tlist[i]->verts[2]->z);
		//		glEnd();
		//	}
		//}
		////----------------------------------------//



		//poly->Silhouette2(this_poly);
		//for (i = 0; i < this_poly->nedges; i++)
		//{
		//	//if (this_poly->elist[i]->verts[0]->back_flag * this_poly->elist[i]->verts[1]->back_flag < 0)
		//	if (abs(this_poly->elist[i]->verts[0]->back_flag + this_poly->elist[i]->verts[1]->back_flag) < 0.5)
		//	{
		//		glBegin(GL_LINES);
		//		glColor3f(1, 0, 0);
		//		glVertex3d(this_poly->elist[i]->verts[0]->x, this_poly->elist[i]->verts[0]->y, this_poly->elist[i]->verts[0]->z);
		//		glColor3f(1, 0, 0);
		//		glVertex3d(this_poly->elist[i]->verts[1]->x, this_poly->elist[i]->verts[1]->y, this_poly->elist[i]->verts[1]->z);
		//		glEnd();
		//	}
		//
		//	
		//}
		//-------------------------------new------------------------------------------
		//poly->Silhouette2(this_poly);
		//glBegin(GL_LINES);
		//for (i = 0; i < this_poly->ntris; i++)
		//{
		//	for (int j = 0; j < 3; j++)
		//	{
		//		if (this_poly->tlist[i]->edges[j]->zeroflag == -1)
		//		{
		//			glColor3f(0, 0, 0);
		//			glVertex3d(this_poly->tlist[i]->edges[j]->zerovecpoint[0], this_poly->tlist[i]->edges[j]->zerovecpoint[1], this_poly->tlist[i]->edges[j]->zerovecpoint[2]);
		//		}
		//	}
		//}
		//glEnd();
		//----------------------------------------------------------------------------
		glDisable(GL_LINE_SMOOTH);
	}


	//------------------tensor field---------------------------//
	//glEnable(GL_LINE_SMOOTH);
	//glLineWidth(1.4);
	//for (int q = 0; q < this_poly->nverts; q++)
	//{	
	//	glm::vec3 v = glm::vec3(this_poly->vlist[q]->x, this_poly->vlist[q]->y, this_poly->vlist[q]->z);
	//
	//	//glBegin(GL_LINES);
	//	//glColor3f(1.0, 0.0, 0.0);
	//	//glm::vec3 p1 = glm::vec3(this_poly->vlist[q]->principlecurvaturevector1[0], this_poly->vlist[q]->principlecurvaturevector1[1], this_poly->vlist[q]->principlecurvaturevector1[2]);
	//	//p1 = p1 / (float)30.;//
	//	//glVertex3d(0.5*p1.x + v.x, 0.5 * p1.y + v.y, 0.5 * p1.z + v.z);
	//	//glVertex3d(-0.5 * p1.x + v.x, -0.5 * p1.y + v.y, -0.5 * p1.z + v.z);
	//	////glVertex3d(p1.x + v.x, p1.y + v.y, p1.z + v.z);
	//	////glVertex3d(v.x, v.y, v.z);
	//	//glEnd();
	//	
	//	glBegin(GL_LINES);
	//	glColor3f(0.0, 0.0, 1.0);
	//	glm::vec3 p2 = glm::vec3(this_poly->vlist[q]->principlecurvaturevector2[0], this_poly->vlist[q]->principlecurvaturevector2[1], this_poly->vlist[q]->principlecurvaturevector2[2]);
	//	p2 = p2 / (float)30.;
	//	glVertex3d(0.5 * p2.x + v.x, 0.5 * p2.y + v.y, 0.5 * p2.z + v.z);
	//	glVertex3d(-0.5 * p2.x + v.x, -0.5 * p2.y + v.y, -0.5 * p2.z + v.z);
	//	//glVertex3d(p2.x + v.x, p2.y + v.y, p2.z + v.z);
	//	//glVertex3d(v.x, v.y, v.z);
	//	glEnd();
	//
	//	//---------------------------------Mean curvature visualization-----------------------------------
	//	//glBegin(GL_LINES);
	//	//glColor3f(0.0, 0.0, 1.0);
	//	//glm::vec3 p2 = glm::vec3(this_poly->vlist[q]->meancurvaturevalue[0], this_poly->vlist[q]->meancurvaturevalue[1], this_poly->vlist[q]->meancurvaturevalue[2]);
	//	//p2 = p2 / (float)40.;
	//	////glVertex3d(0.5 * p2.x + v.x, 0.5 * p2.y + v.y, 0.5 * p2.z + v.z);
	//	////glVertex3d(-0.5 * p2.x + v.x, -0.5 * p2.y + v.y, -0.5 * p2.z + v.z);
	//	//glVertex3d(p2.x + v.x, p2.y + v.y, p2.z + v.z);
	//	//glVertex3d(v.x, v.y, v.z);
	//	//glEnd();
	//}
	//glDisable(GL_LINE_SMOOTH);

	//---------------------------------For triangle tensor test ---------------------------------------------
	for (int i = 0; i < this_poly->ntris; i++)
	{
		glm::vec3 cp = this_poly->tlist[i]->CenterPoint;

		//-------------test tensor-------------------
		//glBegin(GL_LINES);
		//glColor3f(1.0, 0.0, 0.0);
		//glm::vec3 p1 = glm::vec3(this_poly->tlist[i]->principleTricurvaturevector1[0], this_poly->tlist[i]->principleTricurvaturevector1[1], this_poly->tlist[i]->principleTricurvaturevector1[2]);
		//p1 = p1 / (float)30.;//
		//glVertex3d(p1.x + cp.x, p1.y + cp.y, p1.z + cp.z);
		//glVertex3d(cp.x, cp.y, cp.z);
		//glEnd();
		//
		//glBegin(GL_LINES);
		//glColor3f(0.0, 0.0, 1.0);
		//glm::vec3 p2 = glm::vec3(this_poly->tlist[i]->principleTricurvaturevector2[0], this_poly->tlist[i]->principleTricurvaturevector2[1], this_poly->tlist[i]->principleTricurvaturevector2[2]);
		//p2 = p2 / (float)30.;
		//glVertex3d(p2.x + cp.x, p2.y + cp.y, p2.z + cp.z);
		//glVertex3d(cp.x, cp.y, cp.z);
		//glEnd();
		//--------------test StreamLine-----------------
		if (this_poly->tlist[i]->shadow_flag <= 0.2)
		{	
			glBegin(GL_LINES);
			for (int q = 6; q < 12; q++)
			{
				glm::vec3 p1 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q].x, this_poly->tlist[i]->StreamLinePoints[q].y, this_poly->tlist[i]->StreamLinePoints[q].z);
				glm::vec3 p2 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q + 1].x, this_poly->tlist[i]->StreamLinePoints[q + 1].y, this_poly->tlist[i]->StreamLinePoints[q + 1].z);
			
				glColor3f(0.0, 0., 0.0);
				if (length(p1 - p2) < 0.09)
				{
					glVertex3d(p1.x, p1.y, p1.z);
					glVertex3d(p2.x, p2.y, p2.z);
				}
			}
			glEnd();
			
			//glBegin(GL_LINES);
			//for (int q = 0; q < 6; q++)
			//{
			//	glm::vec3 p1 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q].x, this_poly->tlist[i]->StreamLinePoints[q].y, this_poly->tlist[i]->StreamLinePoints[q].z);
			//	glm::vec3 p2 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q + 1].x, this_poly->tlist[i]->StreamLinePoints[q + 1].y, this_poly->tlist[i]->StreamLinePoints[q + 1].z);
			//
			//	glColor3f(0.0, 0.0, 0.0);
			//	if (length(p1 - p2) < 0.09)
			//	{
			//		glVertex3d(p1.x, p1.y, p1.z);
			//		glVertex3d(p2.x, p2.y, p2.z);
			//	}
			//}
			//glEnd();
		}
		if (this_poly->tlist[i]->shadow_flag <= -0.5)
		{
			//glBegin(GL_LINES);
			//for (int q = 6; q < 12; q++)
			//{
			//	glm::vec3 p1 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q].x, this_poly->tlist[i]->StreamLinePoints[q].y, this_poly->tlist[i]->StreamLinePoints[q].z);
			//	glm::vec3 p2 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q + 1].x, this_poly->tlist[i]->StreamLinePoints[q + 1].y, this_poly->tlist[i]->StreamLinePoints[q + 1].z);
			//
			//	glColor3f(0.0, 0., 0.0);
			//	if (length(p1 - p2) < 0.09)
			//	{
			//		glVertex3d(p1.x, p1.y, p1.z);
			//		glVertex3d(p2.x, p2.y, p2.z);
			//	}
			//}
			//glEnd();

			glBegin(GL_LINES);
			for (int q = 0; q < 6; q++)
			{
				glm::vec3 p1 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q].x, this_poly->tlist[i]->StreamLinePoints[q].y, this_poly->tlist[i]->StreamLinePoints[q].z);
				glm::vec3 p2 = glm::vec3(this_poly->tlist[i]->StreamLinePoints[q + 1].x, this_poly->tlist[i]->StreamLinePoints[q + 1].y, this_poly->tlist[i]->StreamLinePoints[q + 1].z);
			
				glColor3f(0.0, 0.0, 0.0);
				if (length(p1 - p2) < 0.08)
				{
					glVertex3d(p1.x, p1.y, p1.z);
					glVertex3d(p2.x, p2.y, p2.z);
				}
			}
			glEnd();
		}
		//glm::vec3 p = glm::vec3(this_poly->tlist[i]->StreamLinePoints[7].x, this_poly->tlist[i]->StreamLinePoints[7].y, this_poly->tlist[i]->StreamLinePoints[7].z);
		//glVertex3d(p.x, p.y, p.z);
		//glVertex3d(cp.x, cp.y, cp.z);
	
	
	}


	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
      glLoadName(i+1);
	
		Triangle *temp_t=this_poly->tlist[i];
	
		switch (display_mode) {
		case 0:
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			} else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			
			glBegin(GL_TRIANGLES);//glBegin(GL_POLYGON);//
			for (j=0; j<3; j++) {
	
				Vertex *temp_v = temp_t->verts[j];
			//glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
			
				if (i == this_poly->seed)
				{ 
					glColor3f(0.0, 0.0, 1.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				else
				{ 
					//---------------------------------for deficit vetex' triangle showing----------------------------------
					//int count = temp_v->ntris;
					//if (count == 1) { glColor3f(1.0, 0.0, 1.0);}// if the valence is small than 6 then the point get red
					//else if (count == 2) { glColor3f(1.0, 0.0, 1.0); }
					//else if (count == 3) { glColor3f(1.0, 0.0, 0.0); }
					//else if (count == 4) { glColor3f(1.0, 1.0, 0.0); }
					//else if (count == 5) { glColor3f(0.0, 1.0, 0.0); }
					//else if (count == 7) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 8) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 9) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 10) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 11) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 12) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 13) { glColor3f(0.0, 0.0, 1.0); }
					//else if (count == 14) { glColor3f(0.0, 0.0, 1.0); }
	
					//glColor4f(0.7, 0.7, 0.0, 0.4);
					glColor4f(0.5, 0.8, 1, 0.4);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}	
			}
			glEnd();

			
			if (this_poly->m_fixedindex >= 0 && this_poly->m_editedintex >= 0)
			{
				glBegin(GL_TRIANGLES);
				for (j = 0; j < this_poly->vlist[this_poly->m_editedintex]->ntris; j++) {
					Triangle* temp_t1 = this_poly->vlist[this_poly->m_editedintex]->tris[j];
					glColor4f(0.0, 0.0, 1.0, 0.8);
					glVertex3d(temp_t1->verts[0]->x, temp_t1->verts[0]->y, temp_t1->verts[0]->z);
					glVertex3d(temp_t1->verts[1]->x, temp_t1->verts[1]->y, temp_t1->verts[1]->z);
					glVertex3d(temp_t1->verts[2]->x, temp_t1->verts[2]->y, temp_t1->verts[2]->z);
				}
				glEnd();

				glBegin(GL_TRIANGLES);
				for (j = 0; j < this_poly->vlist[this_poly->m_fixedindex]->ntris; j++) {
					Triangle* temp_t1 = this_poly->vlist[this_poly->m_fixedindex]->tris[j];
					glColor4f(1.0, 0.0, 0.0, 0.8);
					glVertex3d(temp_t1->verts[0]->x, temp_t1->verts[0]->y, temp_t1->verts[0]->z);
					glVertex3d(temp_t1->verts[1]->x, temp_t1->verts[1]->y, temp_t1->verts[1]->z);
					glVertex3d(temp_t1->verts[2]->x, temp_t1->verts[2]->y, temp_t1->verts[2]->z);
				}
				glEnd();
			}
			break;
		//---------------------------------for checkerboard texture showing------------------------------------
		case 7:
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, Tex0);
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++) {
				Vertex* temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				float L = 0.5;//0.25;//0.1;////1;
				float R, G, B;
				int xf = floor(temp_v->x / L);
				int yf = floor(temp_v->y / L);
				int zf = floor(temp_v->z / L);
				if ((xf % 2) == 0) R = 1.0; else R = 0.0;
				if ((yf % 2) == 0) G = 1.0; else G = 0.0;
				if ((zf % 2) == 0) B = 1.0; else B = 0.0;
				//glColor3f(R, G, B);


				float loops = 5;//1;
				//------------------------------------Tianle TEXTURE------------------------------------
				if ((xf % 2) != 0)
				{
					glTexCoord2f(temp_v->y* loops, temp_v->z * loops);
					if ((zf % 2) != 0) glTexCoord2f(temp_v->y * G * loops, temp_v->z * loops);
					if ((yf % 2) != 0) glTexCoord2f(temp_v->y * loops, temp_v->z * B * loops);
					if ((yf % 2) != 0&& (zf % 2) != 0) glTexCoord2f(temp_v->y * loops, temp_v->z * loops);
				}
				else
				{
					glTexCoord2f(temp_v->x * loops, temp_v->z * loops);
					if ((zf % 2) == 0 && (yf % 2) != 0) glTexCoord2f(temp_v->x * R * loops, temp_v->z * loops);
					if ((yf % 2) == 0 && (zf % 2) != 0) glTexCoord2f(temp_v->x * R * loops, temp_v->z * loops);
					if ((yf % 2) == 0 && (zf % 2) == 0) glTexCoord2f(temp_v->x * loops, temp_v->z * loops);
				}
				//--------------------------------------------------------------------------------------
				glColor3f(1, 1, 1);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			glDisable(GL_TEXTURE_2D);
			break;
		//---------------------------------for heat equation showing------------------------------------
		case 5:
			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; j++) {
				Vertex* temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				//if (temp_v->heatvalue >= 0.9)
				//	glColor3f(1.0, 1.0, 1.0);
				//if (temp_v->heatvalue >= 0.91 && temp_v->heatvalue < 0.95)
				//	glColor3f(1.0, 0.6, 0.6);
				//if (temp_v->heatvalue >= 0.84 && temp_v->heatvalue < 0.91)
				//	glColor3f(1.0, 0.3, 0.3);
				//if (temp_v->heatvalue >= 0.7 && temp_v->heatvalue < 0.77)
				//	glColor3f(1.0, 0.0, 0.0);
				//if (temp_v->heatvalue >= 0.63 && temp_v->heatvalue < 0.7)
				//	glColor3f(1.0, 0.5, 0.0);
				//if (temp_v->heatvalue >= 0.56 && temp_v->heatvalue < 0.63)
				//	glColor3f(1.0, 1.0, 0.0);
				//if (temp_v->heatvalue >= 0.49 && temp_v->heatvalue < 0.56)
				//	glColor3f(0.5, 1.0, 0.0);
				//if (temp_v->heatvalue >= 0.42 && temp_v->heatvalue < 0.49)
				//	glColor3f(0.0, 1.0, 0.0);
				//if (temp_v->heatvalue >= 0.35 && temp_v->heatvalue < 0.42)
				//	glColor3f(0.0, 1.0, 0.5);
				//if (temp_v->heatvalue >= 0.28 && temp_v->heatvalue < 0.35)
				//	glColor3f(0.0, 1.0, 1.0);
				//if (temp_v->heatvalue >= 0.21 && temp_v->heatvalue < 0.28)
				//	glColor3f(0.0, 0.5, 1.0);
				//if (temp_v->heatvalue >= 0.14 && temp_v->heatvalue < 0.21)
				//	glColor3f(0.0, 0.0, 1.0);
				//if (temp_v->heatvalue >= 0.07 && temp_v->heatvalue < 0.14)
				//	glColor3f(0.5, 0.0, 1.0);
				//if (temp_v->heatvalue < 0.07)
				//	glColor3f(0.5, 0.0, 0.5);

				glColor3f(temp_v->r, temp_v->g, temp_v->b);

				//int x = temp_v->heatvalue * 256. * 256. * 256.;
				//int r = x % 256;
				//int g = ((x - r) / 256 % 256);
				//int b = ((x - r - g * 256) / (256 * 256) % 256);
				//float rf = (float)r / 256.;
				//float gf = (float)g / 256.;
				//float bf = (float)b / 256.;
				//glColor3f(rf, gf, bf);
				//std::cout << rf << "," << gf << "," << bf << "," << std::endl;
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();

			break;

		case 8:

			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, Tex0);
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++) {
				Vertex* temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				if (i == this_poly->seed)
				{
					glColor3f(0.0, 0.0, 1.0);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
				glColor3f(0.8, 0.8, 0.8);
				glTexCoord2f(temp_v->tx*20, temp_v->ty*20);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			glDisable(GL_TEXTURE_2D);

			//glEnable(GL_TEXTURE_2D);
			//glBindTexture(GL_TEXTURE_2D, Tex0);
			//glBegin(GL_TRIANGLES);
			//glColor3d(1, 0, 0);
			//glTexCoord2f(0., 0.);
			//glVertex3d(1, 1, 0);
			//glTexCoord2f(1., 0.);
			//glVertex3d(1, -1, 0);
			//glTexCoord2f(1., 1.);
			//glVertex3d(-1, -1, 0);
			//glTexCoord2f(0., 1.);
			//glVertex3d(-1, 1, 0);
			//glEnd();
			//glDisable(GL_TEXTURE_2D);
			break;
	
		//case 9:
		//	glBegin(GL_POLYGON);
		//	for (j=0; j<3; j++) {
		//		mat_diffuse[0] = 1.0;
		//		mat_diffuse[1] = 0.0;
		//		mat_diffuse[2] = 0.0;
		//		mat_diffuse[3] = 1.0;
		//
		//		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		//
		//		Vertex *temp_v = temp_t->verts[j];
		//		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		//
		//		glColor3f(1.0, 0.0, 0.0);
		//		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		//	}
		//	glEnd();
		//	break;

		case 9:
			glDisable(GL_LIGHTING);
			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; j++) {
				Vertex* temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(0.9, 0.9, 0.9);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 1:
			
			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; j++) {
				Vertex* temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		//--------------------method 1-------------------------
				if (temp_v->gausscurvaturevalue >= 5)
					glColor3f(1.0, 0.0, 0.0);
				if (temp_v->gausscurvaturevalue >= 3 && temp_v->gausscurvaturevalue < 5)
					glColor3f(1.0, 0.5, 0.0);
				if (temp_v->gausscurvaturevalue >= 1 && temp_v->gausscurvaturevalue < 3)
					glColor3f(1.0, 1.0, 0.0);
				if (temp_v->gausscurvaturevalue > 0 && temp_v->gausscurvaturevalue < 1)
					glColor3f(0.5, 1.0, 0.0);
				if (temp_v->gausscurvaturevalue == 0)
					glColor3f(0.0, 1.0, 0.0);
				if (temp_v->gausscurvaturevalue >= -1 && temp_v->gausscurvaturevalue < 0)
					glColor3f(0.0, 1.0, 0.5);
				if (temp_v->gausscurvaturevalue >= -3 && temp_v->gausscurvaturevalue < -1)
					glColor3f(0.0, 1.0, 1.0);
				if (temp_v->gausscurvaturevalue >= -5 && temp_v->gausscurvaturevalue < -3)
					glColor3f(0.0, 0.5, 1.0);
				if (temp_v->gausscurvaturevalue < -5)
					glColor3f(0.0, 0.0, 1.0);
			
			//	//double x = temp_v->gausscurvaturevalue - this_poly->mingausscurv / (this_poly->maxgausscurv - this_poly->mingausscurv); //map the value to [0,1]
			//	//if (x >= 0.875)
			//	//	glColor3f(1.0, 0.0, 0.0);
			//	//if (x >= 0.75 && x < 0.875)
			//	//	glColor3f(1.0, 0.5, 0.0);
			//	//if (x >= 0.625 && x < 0.75)
			//	//	glColor3f(1.0, 1.0, 0.0);
			//	//if (x > 0.5 && x < 0.625)
			//	//	glColor3f(0.5, 1.0, 0.0);
			//	//if (x == 0.5)
			//	//	glColor3f(0.0, 1.0, 0.0);
			//	//if (x >= 0.375 && x < 0.5)
			//	//	glColor3f(0.0, 1.0, 0.5);
			//	//if (x >= 0.25 && x < 0.375)
			//	//	glColor3f(0.0, 1.0, 1.0);
			//	//if (x >= 0.125 && x < 0.25)
			//	//	glColor3f(0.0, 0.5, 1.0);
			//	//if (x < 0.125)
			//	//	glColor3f(0.0, 0.0, 1.0);
		//--//-------------------method 2--------------------------
			//	//double x = temp_v->gausscurvaturevalue - this_poly->mingausscurv / (this_poly->maxgausscurv - this_poly->mingausscurv); //map the value to [0,1]
			//	//double r, g, b;
			//	//colorMap(x, r, g, b);
			//	//glColor3f(r, g, b);
		//--//-----------------------------------------------------
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();

			//glBegin(GL_TRIANGLES);
			//for (j = 0; j < 3; j++) {
			//	Vertex* temp_v = temp_t->verts[j];
			//	glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
			//	
			//	double value = temp_v->meancurvaturemagnitude;//sqrt(pow(temp_v->meancurvaturevalue[0], 2)+ pow(temp_v->meancurvaturevalue[1], 2)+ pow(temp_v->meancurvaturevalue[2], 2)) / 2;
			//	
			//	if (value >= 5)
			//		glColor3f(1.0, 0.0, 0.0);
			//	if (value >= 3 && value < 5)
			//		glColor3f(1.0, 0.5, 0.0);
			//	if (value >= 1 && value < 3)
			//		glColor3f(1.0, 1.0, 0.0);
			//	if (value > 0 && value < 1)
			//		glColor3f(0.5, 1.0, 0.0);
			//	if (value == 0)
			//		glColor3f(0.0, 1.0, 0.0);
			//	if (value >= -1 && value < 0)
			//		glColor3f(0.0, 1.0, 0.5);
			//	if (value >= -3 && value < -1)
			//		glColor3f(0.0, 1.0, 1.0);
			//	if (value >= -5 && value < -3)
			//		glColor3f(0.0, 0.5, 1.0);
			//	if (value < -5)
			//		glColor3f(0.0, 0.0, 1.0);
			//
		//--//-------------------method 2--------------------------
			////double x = temp_v->meancurvaturemagnitude - this_poly->minmeancurv / (this_poly->maxmeancurv - this_poly->minmeancurv); //map the value to [0,1]
			////double r, g, b;
			////colorMap(x, r, g, b);
			////glColor3f(r, g, b);
		//--//-----------------------------------------------------
			//
			//	glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			//}
			//glEnd();

			break;
		}
	}

	//--------------------------------------------------for deficit points-------------------------------------------------------
	//glBegin(GL_POINTS);
	//for (i = 0; i < this_poly->nverts; i++) {
	//	int count = this_poly->vlist[i]->ntris;
	//	if ( count == 1) { glColor3f(1.0, 0.0, 1.0); count_valence += 5;}// if the valence is small than 6 then the point get red
	//	else if (count == 2) { glColor3f(1.0, 0.0, 1.0); count_valence += 4;}
	//	else if (count == 3) { glColor3f(1.0, 0.0, 0.0); count_valence += 3;}
	//	else if (count == 4) { glColor3f(1.0, 1.0, 0.0); count_valence += 2;}
	//	else if (count == 5) { glColor3f(0.0, 1.0, 0.0); count_valence += 1;}
	//	else if (count == 7) { glColor3f(0.0, 0.0, 1.0); count_valence -= 1; }
	//	else if (count == 8) { glColor3f(0.0, 0.0, 1.0); count_valence -= 2; }
	//	else if (count == 9) { glColor3f(0.0, 0.0, 1.0); count_valence -= 3; }
	//	else if (count == 10) { glColor3f(0.0, 0.0, 1.0); count_valence -= 4; }
	//	else if (count == 11) { glColor3f(0.0, 0.0, 1.0); count_valence -= 5; }
	//	else if (count == 12) { glColor3f(0.0, 0.0, 1.0); count_valence -= 6; }
	//	else if (count == 13) { glColor3f(0.0, 0.0, 1.0); count_valence -= 7; }
	//	else if (count == 14) { glColor3f(0.0, 0.0, 1.0); count_valence -= 8; }
	//	glVertex3d(this_poly->vlist[i]->x, this_poly->vlist[i]->y, this_poly->vlist[i]->z);
	//}
	//glEnd();

	//----------------------------------------------------for angle deficit------------------------------------------------------
	//glBegin(GL_POINTS);
	//glPointSize(20);
	//for (i = 0; i < this_poly->nverts; i++) {
	//	double sumangle = 0;
	//	for (j = 0; j < this_poly->vlist[i]->ntris; j++) {
	//		for (int k = 0; k < 3; k++) {
	//			if (this_poly->vlist[i]->tris[j]->verts[k] == this_poly->vlist[i])
	//			{
	//				double angle;
	//				Vertex* nowpoint = this_poly->vlist[i];
	//				std::vector <Vertex*> otherpoints;
	//				if (k != 0) otherpoints.push_back(this_poly->vlist[i]->tris[j]->verts[0]);
	//				if (k != 1) otherpoints.push_back(this_poly->vlist[i]->tris[j]->verts[1]);
	//				if (k != 2) otherpoints.push_back(this_poly->vlist[i]->tris[j]->verts[2]);
	//				vec3 vector1 = vec3(otherpoints[0]->x - nowpoint->x, otherpoints[0]->y - nowpoint->y, otherpoints[0]->z - nowpoint->z);
	//				vec3 vector2 = vec3(otherpoints[1]->x - nowpoint->x, otherpoints[1]->y - nowpoint->y, otherpoints[1]->z - nowpoint->z);
	//				angle = 180.0*glm::acos(glm::dot(vector1, vector2)/(glm::length(vector1)* glm::length(vector2)))/3.1415926;
	//				sumangle += angle;
	//				break;
	//			}
	//		}
	//	}
	//	int count = 360 - sumangle;
	//	if (count > 0)
	//		glColor3f(0.0, 0.0, 1.0);
	//	else
	//		glColor3f(1.0, 0.0, 0.0);
	//	glVertex3d(this_poly->vlist[i]->x, this_poly->vlist[i]->y, this_poly->vlist[i]->z);
	//	
	//	count_angle += (360.0 - sumangle);
	//}
	//glEnd();


	//Vertex* v1 = this_poly->tlist[208]->corners[0]->v;
	//Vertex* v2 = this_poly->tlist[208]->corners[1]->v;
	//Vertex* v3 = this_poly->tlist[208]->corners[2]->v;
	//glBegin(GL_TRIANGLES);
	//glColor3f(1.0, 0.0, 0.0);
	//glVertex3d(v1->x, v1->y, v1->z);
	//glColor3f(0.0, 1.0, 0.0);
	//glVertex3d(v2->x, v2->y, v2->z);
	//glColor3f(0.0, 0.0, 1.0);
	//glVertex3d(v3->x, v3->y, v3->z);
	//glEnd();


	//std::cout << "\nTest 1" << std::endl;
	//Corner* c = this_poly->tlist[2]->corners[2];
	//std::cout << "corner.c ="<< c->c / PI * 180.0 << std::endl;
	//for (int i = 0; i < this_poly->tlist[2]->corners[2]->v->ntris ; i++)
	//{
	//	c = c->p->o->p;
	//	std::cout << "neighbor"<< i << " .c ="<< c->c/PI*180.0 << std::endl;
	//}
	//
	//std::cout << "\nTest 2"<< std::endl;
	//std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c / PI * 180.0 << std::endl;
	//std::cout << "c.n.n.n.c =" << this_poly->tlist[2]->corners[2]->n->n->n->c / PI * 180.0 << std::endl;
	//
	//std::cout << "\nTest 3" << std::endl;
	//std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c / PI * 180.0 << std::endl;
	//std::cout << "c.p.p.p.c =" << this_poly->tlist[2]->corners[2]->p->p->p->c / PI * 180.0 << std::endl;
	//
	//std::cout << "\nTest 4" << std::endl;
	//std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c / PI * 180.0 << std::endl;
	//std::cout << "c.n.p.c =" << this_poly->tlist[2]->corners[2]->n->p->c / PI * 180.0 << std::endl;
	//std::cout << "c.p.n.c =" << this_poly->tlist[2]->corners[2]->p->n->c / PI * 180.0 << std::endl;
	//
	//std::cout << "\nTest 5" << std::endl;
	//std::cout << "corner.e.length =" << this_poly->tlist[2]->corners[2]->e->length << std::endl;
	//std::cout << "corner.o.e.length =" << this_poly->tlist[2]->corners[2]->o->e->length << std::endl;
	//
	//std::cout << "\nTest 6" << std::endl;
	//std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c / PI * 180.0 << std::endl;
	//std::cout << "corner.o.o.c =" << this_poly->tlist[2]->corners[2]->o->o->c / PI * 180.0 << std::endl;
	//
	//
	//std::cout<<"\nTotal Vertices valence deficit: "<<count_valence<<std::endl;
	//std::cout << "Total Angless valence deficit: " << count_angle << std::endl;
	//std::cout << "Total number of Edges: " << this_poly->nedges << std::endl;
	//std::cout << "Total number of Vertices: " << this_poly->nverts << std::endl;
	//std::cout << "Total number of Triangles: " << this_poly->ntris << std::endl;
	//std::cout << "Total number of Corners: " << this_poly->ncorners << std::endl;
	//


	//std::cout << "\nMax + Min - Saddle: " << this_poly->Mf << std::endl;
	//if (this_poly->m_editedintex >= 0 && this_poly->m_fixedindex >= 0)
	//{
	//	std::cout << "blue: " << this_poly->vlist[this_poly->m_editedintex]->x << ", " << this_poly->vlist[this_poly->m_editedintex]->y << ", " << this_poly->vlist[this_poly->m_editedintex]->z << std::endl;
	//	std::cout << "red: " << this_poly->vlist[this_poly->m_fixedindex]->x << ", " << this_poly->vlist[this_poly->m_fixedindex]->y << ", " << this_poly->vlist[this_poly->m_fixedindex]->z << std::endl;
	//}
	std::cout << "random: " << this_poly->vlist[0]->x << ", " << this_poly->vlist[0]->y << ", " << this_poly->vlist[0]->z << std::endl;
	//std::cout << "seed: " << this_poly->seed<< std::endl;
	//std::cout << "anchor red index: " << this_poly->m_fixedindex << std::endl;
	//std::cout << "anchor blue index: " << this_poly->m_editedintex << std::endl;


	//std::cout << "\n vert0:"<< this_poly->tlist[0]->verts[0]->x << ", " << this_poly->tlist[0]->verts[0]->y << ", " << this_poly->tlist[0]->verts[0]->z << std::endl;
	//std::cout << "\n vert1:" << this_poly->tlist[0]->verts[1]->x << ", " << this_poly->tlist[0]->verts[1]->y << ", " << this_poly->tlist[0]->verts[1]->z << std::endl;
	//std::cout << "\n vert2:" << this_poly->tlist[0]->verts[2]->x << ", " << this_poly->tlist[0]->verts[2]->y << ", " << this_poly->tlist[0]->verts[2]->z << std::endl;
	//
	//std::cout << "\n ct0:" << this_poly->tlist[0]->corners[0]->v->x << ", " << this_poly->tlist[0]->corners[0]->v->y << ", " << this_poly->tlist[0]->corners[0]->v->z << std::endl;
	////std::cout << "\n ct0:" << this_poly->tlist[0]->corners[1]->v->x << ", " << this_poly->tlist[0]->corners[1]->v->y << ", " << this_poly->tlist[0]->corners[1]->v->z << std::endl;
	////std::cout << "\n ct0:" << this_poly->tlist[0]->corners[2]->v->x << ", " << this_poly->tlist[0]->corners[2]->v->y << ", " << this_poly->tlist[0]->corners[2]->v->z << std::endl;
	//std::cout << "\n ct0:" << this_poly->tlist[0]->corners[0]->n->v->x << ", " << this_poly->tlist[0]->corners[0]->n->v->y << ", " << this_poly->tlist[0]->corners[0]->n->v->z << std::endl;
	//std::cout << "\n ct0:" << this_poly->tlist[0]->corners[0]->n->n->v->x << ", " << this_poly->tlist[0]->corners[0]->n->n->v->y << ", " << this_poly->tlist[0]->corners[0]->n->n->v->z << std::endl;
	//std::cout << "\n normal:" << this_poly->tlist[0]->normal.entry[0] << ", " << this_poly->tlist[0]->normal.entry[1] << ", " << this_poly->tlist[0]->normal.entry[2] << std::endl;

	//for (int i = 0; i < this_poly->nverts; i++)
	//std::cout << "\n---------------------------------------------------------------------" << std::endl;
	//for (int i = 0; i < this_poly->nverts; i++)
	//	std::cout << "A vertex "<< i <<" coordinates in Feline model: " << this_poly->vlist[i]->x<<"," << this_poly->vlist[i]->y << "," << this_poly->vlist[i]->z << std::endl;
	//std::cout << "---------------------------------------------------------------------" << std::endl;
	//	std::cout << this_poly->vlist[i]->tx <<","<< this_poly->vlist[i]->ty << std::endl;
		//std::cout << this_poly->vlist[i]->heatvalue << std::endl;

}

void display(void)
{
  GLint viewport[4];
  int jitter;

  //glClearColor (0.6, 0.6, 0.6, 1.0);  // background for rendering color coding and lighting
  glClearColor(1., 1., 1., 1.0);

  glGetIntegerv (GL_VIEWPORT, viewport);
 
  glClear(GL_ACCUM_BUFFER_BIT);
  for (jitter = 0; jitter < ACSIZE; jitter++) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);//set_view(GL_RENDER, new_poly);
    glPushMatrix ();
		switch(ACSIZE){
		case 1:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;

		case 16:
			glTranslatef (ji16[jitter].x*2.0/viewport[2], ji16[jitter].y*2.0/viewport[3], 0.0);
			break;

		default:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);//set_scene(GL_RENDER, new_poly);
		display_shape(GL_RENDER, poly);//display_shape(GL_RENDER, new_poly);
    glPopMatrix ();
    //glAccum(GL_ACCUM, 1.0/ACSIZE);
  }
  //glAccum (GL_RETURN, 1.0);
  glFlush();
  glutSwapBuffers();
 	glFinish();
}

void Polyhedron::average_normals()
{
	int i, j;

	for (i=0; i<nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j=0; j<vlist[i]->ntris; j++) 
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		normalize(vlist[i]->normal);
	}
}


