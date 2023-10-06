/*

Functions for learnply

Eugene Zhang, 2005
*/

#include <stdlib.h>
#include <stdio.h>
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
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

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

void init(void);
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

	this_file = fopen("../tempmodels/bunny.ply", "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );	

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	//npoly = new Polyhedron;
	//npoly->regular_sub(poly);
	//npoly->initialize();
	//npoly->calc_bounding_sphere();
	////npoly->calc_face_normals_and_area();
	////npoly->average_normals();
	//
	//new_poly = new Polyhedron;
	//new_poly->regular_sub(npoly);
	//new_poly->initialize();
	//new_poly->calc_bounding_sphere();
	////new_poly->calc_face_normals_and_area();
	////new_poly->average_normals();

	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition (20, 20);
	glutInitWindowSize (win_width, win_height); 
	glutCreateWindow ("Geometric Modeling");
	init ();
	glutKeyboardFunc (keyboard);
	glutDisplayFunc(display); 
	glutMotionFunc (motion);
	glutMouseFunc (mouse);
	glutMainLoop(); 

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
      free (tlist[i]);
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
		free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i=0; i<nedges; i++) {
		free(elist[i]->tris);
		free(elist[i]);
	}
	for (i=0; i<nverts; i++) {
		free(vlist[i]->tris);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	free(tlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
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
      list[i] = elist[i];

    /* replace list */
    free (elist);
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
	angle = 180.0*glm::acos(glm::dot(vector1, vector2)/(glm::length(vector1)* glm::length(vector2)))/3.1415926;
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
			list[i] = clist[i];

		/* replace list */
		free(clist);
		clist = list;
	}

	/* create the corner */    //--------------------------------Place FIX all Corner MEMBERS---------------------------------

	clist[ncorners] = new Corner;
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
			f->corners[j]->p = f->corners[(j+1)%3];

			f->corners[j]->n = new Corner;
			f->corners[j]->n = f->corners[(j + 2) % 3];
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
	ncorners = 0;

	/* zero out all the pointers from faces to corners */

	for (i = 0; i < ntris; i++)
		for (j = 0; j < 3; j++)
			tlist[i]->corners[j] = NULL;

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
			clist[ncorners-1]->e = new Edge;
			clist[ncorners-1]->e = f->edges[(j + 1) % f->nverts]; // set the edge of corner

			// make room for the face pointer 
			clist[ncorners-1]->t = new Triangle;
			clist[ncorners-1]->t = tlist[i]; // set the triangle of corner
			
			f->corners[j] = new Corner;
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
	free(tempB);
	free(tempC);
	
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

		case '1':
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
			display();
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
			display();
			break;

		case '4':
			display_mode = 4;
			display();
			break;

		case '5':
			display_mode = 5;
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
			display();
			break;

		case '9':
			display_mode = 9;
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

	glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
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
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
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
	} else if (button == GLUT_MIDDLE_BUTTON) {
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
void Polyhedron::create_new_vert(Polyhedron* poly)  //create new added vertices by travese all the edeges
{
	for (int i = 0; i < poly->nedges; i++)
	{
		Vertex* v = new Vertex(0,0,0);
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
			clean_vlist1[i - poly->max_verts]->edge_index = poly->elist[i - poly->max_verts]->index;
		}
		
	}

	// create new vetex list
	for (int i = 0; i < poly->max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
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
		clean_vlist[i]->edge_index = -1;
		clean_vlist[i]->index = poly->vlist[i]->index;
	}

	Vertex** clean_vlist1 = new Vertex * [(int)(max_verts- poly->max_verts)];
	for (int i = poly->max_verts; i < max_verts; i++)
	{
		clean_vlist1[i - poly->max_verts] = new Vertex(poly->elist[i - poly->max_verts]->new_vert->x, poly->elist[i - poly->max_verts]->new_vert->y, poly->elist[i - poly->max_verts]->new_vert->z);
		clean_vlist1[i - poly->max_verts]->edge_index = poly->elist[i - poly->max_verts]->index;
	}

	// create new vetex list
	for (int i = 0; i < poly->max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
		vlist[i] = clean_vlist[i];
		nverts += 1;
	}
	for (int i = poly->max_verts; i < max_verts; i++)
	{
		vlist[i] = new Vertex(0, 0, 0);
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
			free(tlist[i]);
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
	glShadeModel(GL_SMOOTH);
	//glShadeModel(GL_FLAT);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1); 
	
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
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i == this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
				{ 
					//---------------------------------for checkerboard texture showing------------------------------------
					//float L = 0.25; 
					//float R, G, B;
					//int xf = floor(temp_v->x / L);
					//if ((xf % 2) == 0) R = 1.0; else R = 0.0;
					//int yf = floor(temp_v->y / L);
					//if ((yf % 2) == 0) G = 1.0; else G = 0.0;
					//int zf = floor(temp_v->z / L);
					//if ((zf % 2) == 0) B = 1.0; else B = 0.0;
					//glColor3f(R,G,B); 
	
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
					glColor3f(0.5, 0.5, 0.5);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
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
	
		case 10:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
		
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
	
				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
	}

	//--------------------------------------------------for deficit points-------------------------------------------------------
	glBegin(GL_POINTS);
	for (i = 0; i < this_poly->nverts; i++) {
		int count = this_poly->vlist[i]->ntris;
		if ( count == 1) { glColor3f(1.0, 0.0, 1.0); count_valence += 5;}// if the valence is small than 6 then the point get red
		else if (count == 2) { glColor3f(1.0, 0.0, 1.0); count_valence += 4;}
		else if (count == 3) { glColor3f(1.0, 0.0, 0.0); count_valence += 3;}
		else if (count == 4) { glColor3f(1.0, 1.0, 0.0); count_valence += 2;}
		else if (count == 5) { glColor3f(0.0, 1.0, 0.0); count_valence += 1;}
		else if (count == 7) { glColor3f(0.0, 0.0, 1.0); count_valence -= 1; }
		else if (count == 8) { glColor3f(0.0, 0.0, 1.0); count_valence -= 2; }
		else if (count == 9) { glColor3f(0.0, 0.0, 1.0); count_valence -= 3; }
		else if (count == 10) { glColor3f(0.0, 0.0, 1.0); count_valence -= 4; }
		else if (count == 11) { glColor3f(0.0, 0.0, 1.0); count_valence -= 5; }
		else if (count == 12) { glColor3f(0.0, 0.0, 1.0); count_valence -= 6; }
		else if (count == 13) { glColor3f(0.0, 0.0, 1.0); count_valence -= 7; }
		else if (count == 14) { glColor3f(0.0, 0.0, 1.0); count_valence -= 8; }
		glVertex3d(this_poly->vlist[i]->x, this_poly->vlist[i]->y, this_poly->vlist[i]->z);
	}
	glEnd();

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


	std::cout << "\nTest 1" << std::endl;
	Corner* c = this_poly->tlist[2]->corners[2];
	std::cout << "corner.c ="<< c->c << std::endl;
	for (int i = 0; i < this_poly->tlist[2]->corners[2]->v->ntris ; i++)
	{
		c = c->p->o->p;
		std::cout << "neighbor"<< i << " .c ="<< c->c << std::endl;
	}
	
	std::cout << "\nTest 2"<< std::endl;
	std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c << std::endl;
	std::cout << "c.n.n.n.c =" << this_poly->tlist[2]->corners[2]->n->n->n->c<< std::endl;
	
	std::cout << "\nTest 3" << std::endl;
	std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c << std::endl;
	std::cout << "c.p.p.p.c =" << this_poly->tlist[2]->corners[2]->p->p->p->c << std::endl;
	
	std::cout << "\nTest 4" << std::endl;
	std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c << std::endl;
	std::cout << "c.n.p.c =" << this_poly->tlist[2]->corners[2]->n->p->c << std::endl;
	std::cout << "c.p.n.c =" << this_poly->tlist[2]->corners[2]->p->n->c << std::endl;
	
	std::cout << "\nTest 5" << std::endl;
	std::cout << "corner.e.length =" << this_poly->tlist[2]->corners[2]->e->length << std::endl;
	std::cout << "corner.o.e.length =" << this_poly->tlist[2]->corners[2]->o->e->length << std::endl;
	
	std::cout << "\nTest 6" << std::endl;
	std::cout << "corner.c =" << this_poly->tlist[2]->corners[2]->c << std::endl;
	std::cout << "corner.o.o.c =" << this_poly->tlist[2]->corners[2]->o->o->c << std::endl;
	
	
	std::cout<<"\nTotal Vertices valence deficit: "<<count_valence<<std::endl;
	std::cout << "Total Angless valence deficit: " << count_angle << std::endl;
	std::cout << "Total number of Edges: " << this_poly->nedges << std::endl;
	std::cout << "Total number of Vertices: " << this_poly->nverts << std::endl;
	std::cout << "Total number of Triangles: " << this_poly->ntris << std::endl;
	std::cout << "Total number of Corners: " << this_poly->ncorners << std::endl;


}

void display(void)
{
  GLint viewport[4];
  int jitter;

  glClearColor (0.6, 0.6, 0.6, 1.0);  // background for rendering color coding and lighting
 
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
    glAccum(GL_ACCUM, 1.0/ACSIZE);
  }
  glAccum (GL_RETURN, 1.0);
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


