#include "../obj3d.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

struct edge{
	int face;
	int a;
	int b;
	int mid_point;
	bool sliced;
};

struct face{
	int edges[3];
};

struct edges_info{
	edge *edges;
	face *faces;
};


/*
mesh3D diamond{
	{{0,1,0}, {-1,0,1}, {1,0,1}, {1,0,-1}, {-1,0,-1}, {0, -1, 0}},
	{{0,1,0}, {-1,0,1}, {1,0,1}, {1,0,-1}, {-1,0,-1}, {0, -1, 0}},
	0,
	0,
	{{0,1,0}, {-1,0,1}, {1,0,1}, {1,0,-1}, {-1,0,-1}, {0, -1, 0}},
	0,
	0,
	{0,0,0},
	1,
	6,
	6,
	8,
	{1,2,0, 2,3,0, 3,4,0, 4,1,0, 1,2,5, 2,3,5, 3,4,5, 4,1,5},
	{0,0,0,0}
};
*/

edges_info make_edge(mesh3D *obj){

	edge *edges = (edge*)malloc(sizeof(edge)*obj->nb_f);
	face *faces = (face*)malloc(sizeof(face)*obj->nb_f/3);
	int edge_pointer = 0;

	for (int f = 0; f < obj->nb_f; f += 3){
		for (int i = 0; i<3; i++){ 
			bool ok = true;
			int a = obj->face[f + i];
			int b = obj->face[f + ((i+1)%3)];

			for (int e=0; e < edge_pointer; e++){
				if ((edges[e].a == a && edges[e].b == b) || (edges[e].a == b && edges[e].b == a)){
					faces[f/3].edges[i] = e;
					ok = false;
					break;
				}
			}
			if (ok){
				faces[f/3].edges[i] = edge_pointer;
				edges[edge_pointer].face = f;
				edges[edge_pointer].a = a;
				edges[edge_pointer].b = b;
				edges[edge_pointer].sliced = false;
				edge_pointer++;
			}
		}
	}
	printf("nb edges %d \n", edge_pointer);
	//for (int i= 0; i< edge_pointer; i++)printf("a:%d  b:%d face:%d\n",edges[i].a, edges[i].b, edges[i].face);
	for (int i= 0; i< obj->nb_f/3; i++)printf(":%d  :%d :%d\n",faces[i].edges[0], faces[i].edges[1], faces[i].edges[2]);
	realloc(edges, sizeof(edge)*edge_pointer);
	edges_info out;
	out.edges = edges;
	out.faces = faces;

	return out;
}


void tessel(mesh3D *obj){
	int nb_f = obj->nb_f;
	int nb = obj->nb;
	printf("mmmmm \n");

	int *tp_face = (int*)malloc(sizeof(int)*nb_f*4);
	vector *tp_vertics = (vector*)malloc(sizeof(vector)*nb*4);
	vector *tp_normal = (vector*)malloc(sizeof(vector)*nb*4);
	printf("done malloc \n");
	//memcpy (tp_face, obj->face, sizeof(int)*nb_f);
	memcpy(tp_vertics, obj->vertics, sizeof(vector)*nb);
	memcpy(tp_normal, obj->normals, sizeof(vector)*nb);
	//obj->face = tp_face;
	obj->vertics = tp_vertics;
	obj->normals = tp_normal;
	obj->vertics_trans = (vector*)malloc(sizeof(vector)*nb*3);

	edges_info edgeF = make_edge(obj);

	printf("done realloc \n");
	int vertex_pointer = nb;
	int face_pointer = 0;

	for (int face = 0; face < obj->nb_f; face += 3){
		int new_vertics[3];

		for (int i=0; i<3; i++){
			edge* current_edge = edgeF.edges + edgeF.faces[face/3].edges[i];

			if (current_edge -> sliced){
				printf("  edge %d, midpoint %d\n", edgeF.faces[face/3].edges[i], current_edge -> mid_point);
				new_vertics[i] = current_edge -> mid_point;
			}else{
				new_vertics[i] = vertex_pointer;
		 		obj->vertics[new_vertics[i]] = add(obj->vertics + obj->face[face + i], obj->vertics + obj->face[face + ((i + 1)%3)]);
				normalize(obj->vertics + new_vertics[i] );

		 		obj->normals[new_vertics[i]] = add(obj->normals + obj->face[face + i], obj->normals + obj->face[face + ((i + 1)%3)]);
				normalize(obj->normals + new_vertics[i] );
				current_edge -> mid_point = vertex_pointer;
				current_edge -> sliced = true;
				vertex_pointer++;
			}
		}
		tp_face[face_pointer++] = obj->face[face];
		tp_face[face_pointer++] = new_vertics[0];
		tp_face[face_pointer++] = new_vertics[2];

		tp_face[face_pointer++] = new_vertics[0];
		tp_face[face_pointer++] = obj->face[face + 1];
		tp_face[face_pointer++] = new_vertics[1];

		tp_face[face_pointer++] = new_vertics[2];
		tp_face[face_pointer++] = new_vertics[1];
		tp_face[face_pointer++] = obj->face[face + 2];

		tp_face[face_pointer++] = new_vertics[0];
		tp_face[face_pointer++] = new_vertics[1];
		tp_face[face_pointer++] = new_vertics[2];

	}
	printf("nb_point %d nb_face %d \n", vertex_pointer, face_pointer);
	obj->nb = vertex_pointer;
	obj->nb_f = face_pointer;
	obj->face = tp_face;
	for (int i = 0; i < obj->nb_f; i += 3){printf("face: %d %d %d  %d\n", obj->face[i], obj->face[i+1], obj->face[i+2], (i/3)%3);}
}
