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
	char dir;

	//Matrix2f base;

};


constraint *create_constraint(mesh3D *obj, int *nb_constraints){
	constraint *constraints= new constraint[obj->nb_f];
	int *face = obj->face;
	vector *v = obj->vertics;

	memset(nb_constraints, 0, sizeof(int)*obj->nb);

	int c = 0;
	for (int i=0; i< obj->nb_f; i+=3){
		vector center = add(v + face[i + 1],v + face[i]);
		acc( &center, v + face[i + 2]);
		scalar( &center, 1./3);

		vector edges[3];

		edges[0] = minus(v + face[i], &center);
		edges[1] = minus(v + face[i + 1], &center);
		edges[2] = minus(v + face[i + 2], &center);

		vector n = cross( edges[1], edges[2]);
		normalize(&n);
		int dir = 1;
		/*
		float max;
		vector axe = vec(0,1,0);
		if (dot(&edges[0], &axe) > dot(&edges[1], &axe)){
			dir = 1; 
			max = dot(&edges[1], &axe);
		}
		if (max > dot(&edges[2], &axe))dir = 2;
		*/

		constraints[c].p1 = copy(&edges[dir]);
		normalize(&constraints[c].p1);
		constraints[c].p2 = cross(constraints[c].p1, n);
		normalize(&constraints[c].p2);

		constraints[c].dir = dir;
		//make_plane(edges[1], edges[2], constraints[c].p1, constraints[c].p2);

		constraints[c].rest1 = vec2D(dot(&edges[0], &constraints[c].p1), dot(&edges[0], &constraints[c].p2));
		constraints[c].rest2 = vec2D(dot(&edges[1], &constraints[c].p1), dot(&edges[1], &constraints[c].p2));
		constraints[c].rest3 = vec2D(dot(&edges[2], &constraints[c].p1), dot(&edges[2], &constraints[c].p2));

		//Matrix2f G = Map<Matrix2f >((float*)(&constraints[c].rest1));
		//constraints[c].base = G.inverse();

		//vector test = linear(&constraints[c].p1, &constraints[c].p2, constraints[c].rest1.x, constraints[c].rest1.y);
		//acc(&test, v + face[i]);
		//printVect(test);
		//printVect(v[face[i + 1]]);

		float mat[9]={ edges[1].x, edges[2].x, n.x, 
					edges[1].y, edges[2].y, n.y, 
					edges[1].z, edges[2].z, n.z};

		//print_mat9(mat);
		constraints[c].material_base = gauss(mat);
		nb_constraints[face[i    ]]++;
		nb_constraints[face[i + 1]]++;
		nb_constraints[face[i + 2]]++;
		c++;
	}
	return constraints;
}

void momentum(vector *vertics, int n, vector *v){
	for (int i=0; i<n; i++){
		acc( vertics + i, v);
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
	for (int i=0; i< obj->nb; i++){
		projected_points[i].x = 0;
		projected_points[i].y = 0;
		projected_points[i].z = 0;
	}

	int c = 0;
	for (int i=0; i< obj->nb_f; i+=3){
		//printVect(v[face[i+1]]);

		vector center = add(v + face[i + 1],v + face[i]);
		acc( &center, v + face[i + 2]);
		scalar( &center, 1./3);

		vector edges[3];

		edges[0] = minus(v + face[i    ], &center);
		edges[1] = minus(v + face[i + 1], &center);
		edges[2] = minus(v + face[i + 2], &center);
		vector p1, p2;

		vector n = cross(edges[1], edges[2]);
		normalize(&n);
		p1 = copy(&edges[constraints[c].dir]);
		normalize(&p1);
		p2 = cross(p1, n);
		normalize(&p2);
		
		/*
		Matrix2f X;
		X << dot(&edges[1], &p1), dot(&edges[2], &p1),
			dot(&edges[1], &p2), dot(&edges[2], &p2);

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
		//printVect(edges[1]);
		printVect(v[face[i]]);
		printf("\n");
		*/

		projected_points[face[i    ]].x += center.x + p1.x*constraints[c].rest1.x + p2.x*constraints[c].rest1.y;
		projected_points[face[i    ]].y += center.y + p1.y*constraints[c].rest1.x + p2.y*constraints[c].rest1.y;
		projected_points[face[i    ]].z += center.z + p1.z*constraints[c].rest1.x + p2.z*constraints[c].rest1.y;

		projected_points[face[i + 1]].x += center.x + p1.x*constraints[c].rest2.x + p2.x*constraints[c].rest2.y;
		projected_points[face[i + 1]].y += center.y + p1.y*constraints[c].rest2.x + p2.y*constraints[c].rest2.y;
		projected_points[face[i + 1]].z += center.z + p1.z*constraints[c].rest2.x + p2.z*constraints[c].rest2.y;

		projected_points[face[i + 2]].x += center.x + p1.x*constraints[c].rest3.x + p2.x*constraints[c].rest3.y;
		projected_points[face[i + 2]].y += center.y + p1.y*constraints[c].rest3.x + p2.y*constraints[c].rest3.y;
		projected_points[face[i + 2]].z += center.z + p1.z*constraints[c].rest3.x + p2.z*constraints[c].rest3.y;

		c++;
	}
}

void global_solve(mesh3D *obj, vector *projected_points, int *nb_constraints, vector *velocities){

	vector *v = obj->vertics;

	float m = 0.1;
	for (int i=0; i< obj->nb; i++){
		acc(v + i, velocities + i, DT);
		scalar(v + i, m);
		acc(v + i, projected_points + i);
		scalar(v + i, 1./(m + nb_constraints[i]));
	} 



}
