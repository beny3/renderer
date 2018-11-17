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

using namespace Eigen;

struct constraint{
	vector p1;
	vector p2;

	vector n;
	matrix3 material_base;

	vector2D rest1;
	vector2D rest2;

	//Matrix2f base;

};


void make_plane(vector edge1, vector &edge2, vector &p1, vector &p2){
		p1 = edge1;
		//printVect(edge1);
		//printVect(p1);
		normalize(&p1);
		float f = dot(&edge2, &p1);

		p2.x = edge2.x - p1.x*f;
		p2.y = edge2.y - p1.y*f;
		p2.z = edge2.z - p1.z*f;


		//printVect(p1);
		//printVect(p2);
		//printf("%f \n", dot(&p1, &p2));
		normalize(&p2);


}


constraint *create_constraint(mesh3D *obj, int *nb_constraints){
	constraint *constraints= new constraint[obj->nb_f];
	int *face = obj->face;
	vector *v = obj->vertics;

	memset(nb_constraints, 0, sizeof(int)*obj->nb);

	int c = 0;
	for (int i=0; i< obj->nb_f; i+=3){
		vector edge1 = minus(v + face[i + 1],v + face[i]);
		vector edge2 = minus(v + face[i + 2],v + face[i]);

		vector n = cross( edge1, edge2);
		normalize(&n);
		constraints[c].p1 = copy(&edge1);
		normalize(&constraints[c].p1);
		constraints[c].p2 = cross(constraints[c].p1, n);
		normalize(&constraints[c].p2);

		//make_plane(edge1, edge2, constraints[c].p1, constraints[c].p2);

		constraints[c].rest1 = vec2D(dot(&edge1, &constraints[c].p1), dot(&edge1, &constraints[c].p2));
		constraints[c].rest2 = vec2D(dot(&edge2, &constraints[c].p1), dot(&edge2, &constraints[c].p2));

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

void momentum(vector *vertics, int n, vector *v){
	for (int i=0; i<n; i++){
		acc( vertics + i, v);
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

		vector edge1 = minus(v + face[i + 1],v + face[i]);
		vector edge2 = minus(v + face[i + 2],v + face[i]);
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

		
		projected_points[face[i + 1]].x += v[face[i]].x + p1.x*constraints[c].rest1.x + p2.x*constraints[c].rest1.y;
		projected_points[face[i + 1]].y += v[face[i]].y + p1.y*constraints[c].rest1.x + p2.y*constraints[c].rest1.y;
		projected_points[face[i + 1]].z += v[face[i]].z + p1.z*constraints[c].rest1.x + p2.z*constraints[c].rest1.y;

		projected_points[face[i + 2]].x += v[face[i]].x + p1.x*constraints[c].rest2.x + p2.x*constraints[c].rest2.y;
		projected_points[face[i + 2]].y += v[face[i]].y + p1.y*constraints[c].rest2.x + p2.y*constraints[c].rest2.y;
		projected_points[face[i + 2]].z += v[face[i]].z + p1.z*constraints[c].rest2.x + p2.z*constraints[c].rest2.y;
		
		acc(projected_points + face[i], v + face[i]);
		c++;
	}
}

void global_solve(mesh3D *obj, vector *projected_points, int *nb_constraints){
	vector *v = obj->vertics;
	for (int i=0; i< obj->nb; i++){
		v[i] = projected_points[i];
		scalar(v + i, 1./nb_constraints[i]);
	} 

}
