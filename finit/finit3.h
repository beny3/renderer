#include <iostream>
#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "math.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define DT 1

//using namespace Eigen;

struct constraint{
	vector p1;
	vector p2;

	vector n;
	matrix3 material_base;

	vector2D rest1;
	vector2D rest2;
	vector2D rest3;

	//Matrix2f base;

};

//a bosser
struct edge_info{
	vector2D base_vector1;
	vector2D base_vector2;
	int e1;
	int e2;
	int f1;
	int f2;
	int vert1;
	int vert2;
	vector2D rest1;
	vector2D rest2;
	vector2D rest3;

	vector2D rest4;
	vector2D rest5;
	vector2D rest6;
};

struct edge_constraints{
	edge_info *edge_infos;
	int nb;
};

void make_edge_bases(vector *n1_p, vector*n2_p, vector *local_x, vector *local_z, vector *edge, vector *v1, vector *v2){

	vector n1 = cross(*v1, *edge);
	vector n2 = cross(*edge, *v2);
	normalize(&n1);
	normalize(&n2);

	*local_z = add(&n1, &n2);
	normalize(local_z);
	*local_x =  cross(*local_z, *edge);
	*n1_p = cross(n1, *edge);
	*n2_p = cross(*edge, n2);
	normalize(n1_p);
	normalize(n2_p);
	normalize(local_x);
}


edge_constraints find_edge(mesh3D *obj, int *nb_constraints){
	int nb = obj->nb_f/3;
	edge_constraints constraints;
	constraints.edge_infos = new edge_info[obj->nb_f];
	edge_info *edge_infos = constraints.edge_infos;
	vector *v = obj->vertics;
	int pointer = 0;

	char *marker = new char[nb*nb]();

	for (int faceA = 0; faceA < obj->nb_f; faceA += 3){

		for (int edge = 0; edge < 3; edge++){

			int edge0 = obj->face[faceA + edge];
			int edge1 = obj->face[faceA + (edge + 1)%3];

			for (int face = 0; face < obj->nb_f; face += 3){
				for (int i = 0; i < 3; i++){
					if ( face!= faceA && !marker[nb*faceA/3 + face/3] && !marker[nb*face/3 + faceA/3] &&
							((obj->face[face + i] == edge0 && obj->face[face + ((i+1)%3)] == edge1) ||
							 (obj->face[face + i] == edge1 && obj->face[face + ((i+1)%3)] == edge0))){

						marker[nb*faceA/3 + face/3] = 1;
						marker[nb*face/3 + faceA/3] = 1;

						edge_infos[pointer].e1 = edge0;
						edge_infos[pointer].e2 = edge1;

						edge_infos[pointer].f1 = faceA/3;
						edge_infos[pointer].f2 = face/3;

						int vert1 = obj->face[faceA + (edge + 2)%3];
						int vert2 = obj->face[face + (i + 2)%3];

						edge_infos[pointer].vert1 = vert1;
						edge_infos[pointer].vert2 = vert2;

						vector center = add(v + obj->face[faceA], v + obj->face[faceA + 1]);
						acc(&center, v + obj->face[faceA + 2]);
						scalar( &center, 1./3);

						vector cote1 = minus(v + obj->face[faceA]    , &center);
						vector cote2 = minus(v + obj->face[faceA + 1], &center);
						vector cote3 = minus(v + obj->face[faceA + 2], &center);

						vector p1 = minus(&obj->vertics[edge0], &obj->vertics[edge1]);
						normalize(&p1);
						vector n1_p, n2_p, local_x, local_z;

						vector vec_vert1 = minus(v + vert1, v + edge0);
						vector vec_vert2 = minus(v + vert2, v + edge0);

						make_edge_bases(&n1_p, &n2_p, &local_x, &local_z, &p1, &vec_vert1, &vec_vert2);
						//printVect(local_z);
						edge_infos[pointer].base_vector1 = vec2D(dot(&n1_p, &local_x), dot(&n1_p, &local_z));
						edge_infos[pointer].base_vector2 = vec2D(dot(&n2_p, &local_x), dot(&n2_p, &local_z)); 
						
						edge_infos[pointer].rest1 = vec2D(dot(&cote1, &p1), dot(&cote1, &n1_p));
						edge_infos[pointer].rest2 = vec2D(dot(&cote2, &p1), dot(&cote2, &n1_p));
						edge_infos[pointer].rest3 = vec2D(dot(&cote3, &p1), dot(&cote3, &n1_p));

						center = add(v + obj->face[face], v + obj->face[face + 1]);
						acc(&center, v + obj->face[face + 2]);
						scalar( &center, 1./3);

						cote1 = minus(v + obj->face[face]    , &center);
						cote2 = minus(v + obj->face[face + 1], &center);
						cote3 = minus(v + obj->face[face + 2], &center);

						edge_infos[pointer].rest4 = vec2D(dot(&cote1, &p1), dot(&cote1, &n2_p));
						edge_infos[pointer].rest5 = vec2D(dot(&cote2, &p1), dot(&cote2, &n2_p));
						edge_infos[pointer].rest6 = vec2D(dot(&cote3, &p1), dot(&cote3, &n2_p));

						nb_constraints[obj->face[faceA    ]]++;
						nb_constraints[obj->face[faceA + 1]]++;
						nb_constraints[obj->face[faceA + 2]]++;

						nb_constraints[obj->face[face    ]]++;
						nb_constraints[obj->face[face + 1]]++;
						nb_constraints[obj->face[face + 2]]++;
						pointer++;
					}
				}
			}
		}
	}
	constraints.nb = pointer;
	return constraints;
}


constraint *create_constraint(mesh3D *obj, int *nb_constraints){
	constraint *constraints= new constraint[obj->nb_f];
	int *face = obj->face;
	vector *v = obj->vertics;

	//memset(nb_constraints, 0, sizeof(int)*obj->nb);

	int c = 0;
	for (int i=0; i< obj->nb_f; i+=3){
		vector center = add(v + face[i + 1],v + face[i]);
		acc( &center, v + face[i + 2]);
		scalar( &center, 1./3);

		vector edge1 = minus(v + face[i + 1], &center);
		vector edge2 = minus(v + face[i + 2], &center);
		vector edge3 = minus(v + face[i], &center);

		vector n = cross( edge1, edge2);
		normalize(&n);

		constraints[c].p1 = copy(&edge1);
		normalize(&constraints[c].p1);
		constraints[c].p2 = cross(constraints[c].p1, n);
		normalize(&constraints[c].p2);

		//make_plane(edge1, edge2, constraints[c].p1, constraints[c].p2);

		constraints[c].rest1 = vec2D(dot(&edge1, &constraints[c].p1), dot(&edge1, &constraints[c].p2));
		constraints[c].rest2 = vec2D(dot(&edge2, &constraints[c].p1), dot(&edge2, &constraints[c].p2));
		constraints[c].rest3 = vec2D(dot(&edge3, &constraints[c].p1), dot(&edge3, &constraints[c].p2));

		//Matrix2f G = Map<Matrix2f >((float*)(&constraints[c].rest1));
		//constraints[c].base = G.inverse();

		//vector test = linear(&constraints[c].p1, &constraints[c].p2, constraints[c].rest1.x, constraints[c].rest1.y);
		//acc(&test, v + face[i]);
		//printVect(test);
		//printVect(v[face[i + 1]]);

		float mat[9]={ edge1.x, edge2.x, n.x, 
					edge1.y, edge2.y, n.y, 
					edge1.z, edge2.z, n.z};

		//print_mat9(mat);
		constraints[c].material_base = gauss(mat);
		nb_constraints[face[i    ]]++;
		nb_constraints[face[i + 1]]++;
		nb_constraints[face[i + 2]]++;
		c++;
	}
	return constraints;
}


float set_momentum(vector *velocities, vector *pos, int n){
	for (int i=0; i<n; i++){
		velocities[i].x = pos[i].x + velocities[i].x*DT;
		velocities[i].y = pos[i].y + velocities[i].x*DT;
		velocities[i].z = pos[i].z + velocities[i].z*DT;
	}
}


void get_momentum(vector *new_p, vector *vel, int n){
	for (int i=0; i<n; i++){
		sub_overWrite_left(new_p + i, vel + i);
	}
}


void project_constraint(mesh3D *obj, constraint *constraints, vector *projected_points){

	int *face = obj->face;
	vector *v = obj->vertics;
	int c = 0;
	for (int i=0; i< obj->nb_f; i+=3){
		//printVect(v[face[i+1]]);

		vector center = add(v + face[i + 1],v + face[i]);
		acc( &center, v + face[i + 2]);
		scalar( &center, 1./3);

		vector edge1 = minus(v + face[i + 1], &center);
		vector edge2 = minus(v + face[i + 2], &center);
		vector edge3 = minus(v + face[i], &center);
		vector p1, p2;

		vector n = cross( edge1, edge2);
		normalize(&n);
		p1 = copy(&edge1);
		normalize(&p1);
		p2 = cross(p1, n);
		normalize(&p2);
		
		/*
		Matrix2f X;
		X << dot(&edge1, &p1), dot(&edge2, &p1),
			dot(&edge1, &p2), dot(&edge2, &p2);

		//Matrix2f G = Map<Matrix2f >((float*)(&constraints[c].rest1));		//std::cout << G << "\n";

		JacobiSVD<Matrix2f> svd(X*constraints[c].base, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Matrix2f X_proj = svd.matrixU()*svd.matrixV().transpose();
		
		projected_points[face[i + 1]].x += v[face[i]].x + p1.x*X_proj(0,0) + p2.x*X_proj(1,0);
		projected_points[face[i + 1]].y += v[face[i]].y + p1.y*X_proj(0,0) + p2.y*X_proj(1,0);
		projected_points[face[i + 1]].z += v[face[i]].z + p1.z*X_proj(0,0) + p2.z*X_proj(1,0);

		projected_points[face[i + 2]].x += v[face[i]].x + p1.x*X_proj(0,1) + p2.x*X_proj(1,1);
		projected_points[face[i + 2]].y += v[face[i]].y + p1.y*X_proj(0,1) + p2.y*X_proj(1,1);
		projected_points[face[i + 2]].z += v[face[i]].z + p1.z*X_proj(0,1) + p2.z*X_proj(1,1);
		*/

		/*
		vector test = linear(&p1, &p2,constraints[c].rest3.x, constraints[c].rest3.y);
		acc(&test, &center);
		printVect(test);
		//printVect(edge1);
		printVect(v[face[i]]);
		printf("\n");
		*/

		projected_points[face[i + 1]].x += center.x + p1.x*constraints[c].rest1.x + p2.x*constraints[c].rest1.y;
		projected_points[face[i + 1]].y += center.y + p1.y*constraints[c].rest1.x + p2.y*constraints[c].rest1.y;
		projected_points[face[i + 1]].z += center.z + p1.z*constraints[c].rest1.x + p2.z*constraints[c].rest1.y;

		projected_points[face[i + 2]].x += center.x + p1.x*constraints[c].rest2.x + p2.x*constraints[c].rest2.y;
		projected_points[face[i + 2]].y += center.y + p1.y*constraints[c].rest2.x + p2.y*constraints[c].rest2.y;
		projected_points[face[i + 2]].z += center.z + p1.z*constraints[c].rest2.x + p2.z*constraints[c].rest2.y;

		projected_points[face[i]].x += center.x + p1.x*constraints[c].rest3.x + p2.x*constraints[c].rest3.y;
		projected_points[face[i]].y += center.y + p1.y*constraints[c].rest3.x + p2.y*constraints[c].rest3.y;
		projected_points[face[i]].z += center.z + p1.z*constraints[c].rest3.x + p2.z*constraints[c].rest3.y;

		c++;
	}
}

void project_edges(mesh3D *obj, edge_constraints constraints, vector *projected_points ){

	edge_info *edge_infos = constraints.edge_infos;
	vector *v = obj->vertics;
	int *debug = new int[obj->nb];

	for (int i=0; i< constraints.nb; i++){

		int faceA = edge_infos[i].f1*3;
		int face = edge_infos[i].f2*3;

		int edge0 = edge_infos[i].e1;
		int edge1 = edge_infos[i].e2;

		int vert1 = edge_infos[i].vert1;
		int vert2 = edge_infos[i].vert2;

		vector center = add(v + obj->face[faceA], v + obj->face[faceA + 1]);
		acc(&center, v + obj->face[faceA + 2]);
		scalar( &center, 1./3);

		vector p1 = minus(&obj->vertics[edge0], &obj->vertics[edge1]);
		normalize(&p1);
		vector n1_p, n2_p, local_x, local_z;

		vector vec_vert1 = minus(v + vert1, v + edge0);
		vector vec_vert2 = minus(v + vert2, v + edge0);

		make_edge_bases(&n1_p, &n2_p, &local_x, &local_z, &p1, &vec_vert1, &vec_vert2);

		n1_p = linear(&local_x, &local_z, edge_infos[i].base_vector1.x, edge_infos[i].base_vector1.y);
		n2_p = linear(&local_x, &local_z, edge_infos[i].base_vector2.x, edge_infos[i].base_vector2.y);
		//printVect(local_z);
		//printVect(local_x);
		//edge_infos[pointer].rest1 = vec2D(dot(&cote1, &p1), dot(&cote1, &n1_p));

		vector new_1 = linear(&p1, &n1_p, edge_infos[i].rest1.x, edge_infos[i].rest1.y);
		vector new_2 = linear(&p1, &n1_p, edge_infos[i].rest2.x, edge_infos[i].rest2.y);
		vector new_3 = linear(&p1, &n1_p, edge_infos[i].rest3.x, edge_infos[i].rest3.y);
		acc(&new_1, &center);
		acc(&new_2, &center);
		acc(&new_3, &center);
		acc(&projected_points[obj->face[faceA    ]], &new_1);
		acc(&projected_points[obj->face[faceA + 1]], &new_2);
		acc(&projected_points[obj->face[faceA + 2]], &new_3);

		//printVect( v[obj->face[faceA + 2]]);
		//printVect( new_3);

		center = add(v + obj->face[face], v + obj->face[face + 1]);
		acc(&center, v + obj->face[face + 2]);
		scalar( &center, 1./3);

		new_1 = linear(&p1, &n2_p, edge_infos[i].rest4.x, edge_infos[i].rest4.y);
		new_2 = linear(&p1, &n2_p, edge_infos[i].rest5.x, edge_infos[i].rest5.y);
		new_3 = linear(&p1, &n2_p, edge_infos[i].rest6.x, edge_infos[i].rest6.y);
		acc(&new_1, &center);
		acc(&new_2, &center);
		acc(&new_3, &center);

		//printVect( v[obj->face[face + 2]]);
		//printVect( new_3 );

		acc(&projected_points[obj->face[face    ]], &new_1);
		acc(&projected_points[obj->face[face + 1]], &new_2);
		acc(&projected_points[obj->face[face + 2]], &new_3);
	}
}


void global_solve(mesh3D *obj, vector *projected_points, int *nb_constraints, vector *velocities){
	//printf("%d \n", nb_constraints[2]);
	vector *v = obj->vertics;
	float m = 6/obj->nb;

	float alpha = 1.-1/24.;

	for (int i=0; i< obj->nb; i++){
		v[i] = velocities[i];
		scalar(v + i, alpha*m);
		acc(v + i, projected_points + i, (1-alpha));
		scalar(v + i, 1./(alpha*m + (1-alpha)*nb_constraints[i]));
	} 
}
