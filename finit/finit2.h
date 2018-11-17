#include "Eigen/Dense"
#include <iostream>
#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "math.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define DT 0.1

using namespace Eigen;

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

struct edge_info{
	vector n1;
	vector n2;
	int f1;
	int f2;
	int e1;
	int e2;
	int vert1;
	int vert2;
	float face_base[9];
};

struct edge_constraints{
	edge_info *edge_infos;
	int nb;
};

void make_edge_base(float mat[9], vector *edge1, vector *edge2, vector *v1, vector *v2){

	vector *baseF = (vector*)mat;
	baseF[0] = minus(edge1, edge2);
	normalize(&baseF[0]);
	vector n1 = cross(*v1, baseF[0]);
	vector n2 = cross(baseF[0], *v2);

	normalize(&n1);
	normalize(&n2);

	baseF[1] = add(&n1, &n2);
	normalize(&baseF[1]);
	baseF[2] = cross(baseF[0], baseF[1]);

}


edge_constraints find_edge(mesh3D *obj, int *nb_constraints){
	int nb = obj->nb_f/3;
	edge_constraints constraints;
	constraints.edge_infos = new edge_info[obj->nb_f];
	edge_info *edge_infos = constraints.edge_infos;

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

						edge_infos[pointer].f1 = faceA/3;
						edge_infos[pointer].f2 = face/3;

						edge_infos[pointer].e1 = edge0;
						edge_infos[pointer].e2 = edge1;

						//printf("i %d :", pointer);
						//printVect(obj->vertics[edge1]);

						edge_infos[pointer].vert1 = obj->face[faceA + (edge + 2)%3];
						edge_infos[pointer].vert2 = obj->face[face + (i + 2)%3];
						vector v1 = minus(&obj->vertics[edge_infos[pointer].vert1], &obj->vertics[edge0]);
						vector v2 = minus(&obj->vertics[edge_infos[pointer].vert2], &obj->vertics[edge0]);

						make_edge_base(edge_infos[pointer].face_base,&obj->vertics[edge0], &obj->vertics[edge1],  &v1, &v2);
						//matrix3 bt = transpose9(edge_infos[pointer].face_base);
						//matrix3 ide = mult_m3(bt.data, edge_infos[pointer].face_base);
						print_mat9(edge_infos[pointer].face_base);

						edge_infos[pointer].n1 = mult_vm9( edge_infos[pointer].face_base, obj->normals_f[faceA/3]);
						edge_infos[pointer].n2 = mult_vm9( edge_infos[pointer].face_base, obj->normals_f[face/3]);
						//printf("%d %d | %d %d \n", faceA/3+1, face/3+1 , edge_infos[pointer].vert1 + 1 , edge_infos[pointer].vert2 + 1);

						nb_constraints[edge_infos[pointer].vert1] += 1;
						nb_constraints[edge_infos[pointer].vert2] += 1;
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
	//pas optimal
	float tmp_matrix[9];
	matrix3 local_base;	
	edge_info *edge_infos = constraints.edge_infos;

	for (int i=0; i< constraints.nb; i++){
		vector e_to_v1 = minus(&obj->vertics[edge_infos[i].vert1], &obj->vertics[edge_infos[i].e1]);
		vector e_to_v2 = minus(&obj->vertics[edge_infos[i].vert2], &obj->vertics[edge_infos[i].e1]);

		make_edge_base(tmp_matrix, &obj->vertics[edge_infos[i].e1], &obj->vertics[edge_infos[i].e2], &e_to_v1, &e_to_v2);
		local_base = transpose9(tmp_matrix);
		vector n1 = mult_vm9(local_base.data , edge_infos[i].n1);
		vector n2 = mult_vm9(local_base.data , edge_infos[i].n2);

		//matrix3 ide = mult_m3(local_base.data, tmp_matrix);
		//printVect(obj->vertics[edge_infos[i].e1]);
		//print_mat9(tmp_matrix);
		//print_mat9(edge_infos[i].face_base);
		//print_mat9(ide.data);

		//printf("edge %d \n", edge_infos[i].e2);
		//printVect(obj->normals_f[edge_infos[i].f1]);
		//printVect(n1);

		//printf("dot prod %f \n", dot(&n1, &e_to_v));
		scalar(&n1, 0.1*dot(&n1, &e_to_v1));
		vector proj = minus( &e_to_v1, &n1);
		vector new_vertics = add(&obj->vertics[edge_infos[i].e1], &proj);
		acc(&projected_points[edge_infos[i].vert1], &new_vertics);
		//printf("%d  \n", edge_infos[i].vert1);
		//printVect(n1);
		//printVect(new_vertics);
		//printVect(obj->vertics[edge_infos[i].vert1]);

		scalar(&n2,  0.1*dot(&n2, &e_to_v2));
		proj = minus( &e_to_v2, &n2);
		new_vertics = add(&obj->vertics[edge_infos[i].e1], &proj);
		acc(&projected_points[edge_infos[i].vert2], &new_vertics);
		//printf("%d  \n", edge_infos[i].vert2);
		//printVect(n2);
		//printVect(new_vertics);
		//printVect(obj->vertics[edge_infos[i].vert2]);

	}
}


void global_solve(mesh3D *obj, vector *projected_points, int *nb_constraints, vector *velocities){

	vector *v = obj->vertics;
	float m = 8/obj->nb;

	for (int i=0; i< obj->nb; i++){
		v[i] = velocities[i];
		scalar(v + i, m);
		acc(v + i, projected_points + i);
		scalar(v + i, 1./(m + nb_constraints[i]));
	} 
}
