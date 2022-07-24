
#include "../obj3d.h"
#include "math.h"
/* include some silly stuff */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define SIMPLEX_PUSH_AND_TEST(i)\
	simplex->data[i] = support( &simplex->D, A, na, B, nb);\
	if ( dot(&simplex->data[i], &simplex->D) < 0 ){\
			simplex->s_last = -1;\
			return 0;\
	}

struct tuple{
	vector3D a;
	vector3D b;
	vector3D c;
};

struct Simplex{
	vector3D D;
	vector3D data[4];
	int s_last = -1;
};

float distance(vector3D *A, vector3D *D, vector3D *Dn){
	vector3D AXD = cross(*A, *Dn);
	return dot(D, A);// + (1./sqrt(dot(&AXD, &AXD)));
};

vector3D support2(vector3D *D, vector3D *A, int na, vector3D *B, int nb){
	float max = -100000000;
	vector3D Dn = *D;
	normalize(&Dn);
	vector3D *AB = new vector3D[nb*na];
	vector3D *arMax = AB;

	for (int a=0; a<na; ++a){
		for (int b=0; b<nb; ++b){
			AB[a*nb + b]=minus(&A[a], &B[b]);
			float s = distance(&AB[a*nb + b] , D, &Dn);
			if (s > max){
				max = s;
				arMax = &AB[a*nb + b];
			}
		}
	}

	return *arMax;
}

vector3D support(vector3D *D, vector3D *A, int na, vector3D *B, int nb){
	vector3D D_ = inv(D);
	vector3D Dn = *D;
	normalize(&Dn);
	vector3D Dn_ = inv(&Dn);

	float s  = distance(A, D, &Dn);
	float max = s;

	vector3D *argmaxA = A; 
	vector3D *argmaxB = B;

	for (int i=1; i<na; i++){
		s = distance(&A[i], D, &Dn);

		if( s > max){
			max = s;
			argmaxA = &A[i];
		}
	}

	s = distance(B, &D_, &Dn_);
	max = s;

	for (int i=1; i<nb; i++){
		s = distance(&B[i], &D_, &Dn_);

		if( s > max){
			max = s;
			argmaxB = &B[i];
		}
	}

	return minus(argmaxA, argmaxB);
}

void do_simplex(Simplex *simplexS, int *sign_n){

	int s_last = simplexS->s_last;
	vector3D *simplex = simplexS->data;

	if(s_last == -1){
		//simplexS->D = vec(0,1,0);
		simplexS->s_last = 0;

	}else if(s_last == 0){
		simplexS->D = inv(&simplex[0]);
		simplexS->s_last = 1;

	}else if(s_last == 1){
	
		vector3D ab = minus( &simplex[0], &simplex[1]);
		vector3D ao = inv(&simplex[1]);

		if ( dot(&ab, &ao) >= 0 ){
			vector3D abXao = cross(ab, ao);
			simplexS->D = cross(abXao, ab);
			simplexS->s_last = 2;
		}else{
			printf("of_line \n");
			simplexS->D = ao;
			simplexS->s_last = 1;
		}

	}else if (s_last == 2){

		//cas triangle
		vector3D ao = inv(&simplex[2]);
		vector3D ab = minus(&simplex[1], &simplex[2]);
		vector3D ac = minus(&simplex[0], &simplex[2]);

		vector3D n = cross(ab, ac);
		vector3D nXac = cross(n, ac);
		vector3D abXn = cross(ab, n);

		//outside plan aco
		if (dot(&nXac, &ao) > 0){
			if ( dot(&ac, &ao) > 0 ){
				vector3D acXao = cross(ac, ao);
				simplexS->D = cross(acXao, ac);
				simplex[1]=simplex[2];
				simplexS->s_last = 2;		
			}else{
				simplexS->D = ao;
				simplex[0]=simplex[2];
				simplexS->s_last = 1;
			}
		//outside plan abo
		}else if(dot(&abXn, &ao) > 0 ){

		 	if ( dot(&ab, &ao) > 0 ){
				printf("out_abo\n");
				vector3D abXao = cross(ab, ao);
				simplexS->D = cross(abXao, ab);
				simplex[0]=simplex[1];
				simplex[1]=simplex[2];
				simplexS->s_last = 2;			
			}else{
				printf("out_a \n");
				simplexS->D = ao;
				simplex[0]=simplex[2];
				simplexS->s_last = 1;
			}
		//inside the 2 plan
		}else{
			if ( dot(&n, &ao) > 0 ){
				//printf("in_top\n");
				simplexS->D = n;
				*sign_n = 1;
			}else{
				//printf("in_bottom\n");
				simplexS->D = inv(&n);
				*sign_n = -1;
			}
			simplexS->s_last = 3;
		}
		
	}else{ //s_last == 3
		
		vector3D bc = minus(&simplex[1], &simplex[0]);
		vector3D ab = minus(&simplex[2], &simplex[1]);
		vector3D ac = minus(&simplex[0], &simplex[2]);

		vector3D top_to_base = minus(&simplex[1], &simplex[3]);
		vector3D tpX = cross(bc, top_to_base, *sign_n);
		
		printf("orto %f %f \n", tpX.x*bc.x + tpX.y*bc.y + tpX.z*bc.z, dot(&tpX, &top_to_base));
 
		if ( dot(&simplex[3], &tpX) < 0 ){
			
			simplex[0]=simplex[1];
			simplex[1]=simplex[2];
			simplex[2]=simplex[3];
			simplexS->s_last = 2;
			simplexS->D = tpX;
			return;
		}
		//out a b top
		top_to_base = minus(&simplex[2], &simplex[3]);
		tpX = cross(ab, top_to_base, *sign_n);

		if ( dot( &tpX, &simplex[3]) < 0 ){
			simplex[2] = simplex[3];
			simplexS->s_last = 2;
			simplexS->D = tpX;
			return;
		}
		//out a c top
		top_to_base = minus(&simplex[0], &simplex[3]);
		tpX = cross(ac, top_to_base, *sign_n);

		if ( dot( &tpX, &simplex[3]) < 0 ){
			simplex[1]=simplex[2];
			simplex[2]=simplex[3];
			simplexS->s_last = 2;
			simplexS->D = tpX;
			return;
		}
		simplexS->s_last=4;
	}

}


int gjk(vector3D *A, int na, vector3D *B, int nb, Simplex *simplex, vector3D *AB){
	
	vector3D D_prev = simplex->D;
	/*
	vector3D invS;
	vector3D ab;
	vector3D ac;
	vector3D n;
	vector3D ao;
	*/

	for (int a=0; a<na; ++a){
		for (int b=0; b<nb; ++b){
			AB[a*nb + b]=minus(&A[a], &B[b]);
		}
	}
	//simplex->s_last=-1;
	int sign = 1;
	int i = 0;
	do_simplex(simplex, &sign);

	while( simplex->s_last < 4 && i < 50){
		SIMPLEX_PUSH_AND_TEST(simplex->s_last);
		do_simplex(simplex, &sign);
		printf(" i%d last%d \n", i, simplex->s_last);	
		i++;
	}

	//SIMPLEX_PUSH_AND_TEST(3);
	//simplex->D = D_prev;
	if (simplex->s_last < 4 ) return 0;
	simplex->s_last=-1;
	

	return 1;
}



vector3D colid2(physic  *pyA, physic *pyB, colition_hull *A, colition_hull *B,vector3D *prev){

	vector3D pA = add(&pyA->position, &pyA->velocity);
	vector3D pB = add(&pyB->position, &pyB->velocity);

	
	vector3D out = vec(0,0,0);
	vector3D *normalB = B->mesh.normals;

	vector3D col = minus( &pB, &pA);
	vector3D tpa, tpb, tp;

	// vector3D vert_A_max =0; 
	// vert_B_min = 0;
	float max = -100000000;
	float min = 100000000;
	int count = 0;
	

	for (int i=0; i < A->mesh.nb ; i++){
		for (int j=0; j < B->mesh.nb ; j++){
		
			tp = minus(A->mesh.vertics_trans + i, B->mesh.vertics_trans + j);
			float maybe = dot(&tp, &tp);
			if (maybe < 200){
				acc(&out, A->mesh.vertics_trans + i);
				acc(&out, B->mesh.vertics_trans + j);
				count+=2;
				/*
				min = maybe;
				vert_B_min = j;
				vert_A_max = i;
				*/
			}
		}
	}
	if(count > 0) scalar(&out, 1./count);
	printf("count %d\n", count);
	//printVect(out);

	//vector3D out = linear((A->mesh.vertics_trans + vert_A_max), (B->mesh.vertics_trans + vert_B_min), 0.5, 0.5);
	//vector3D dist = minus(A->mesh.vertics_trans + vert_A_max, B->mesh.vertics_trans + vert_B_min);

	if (dot(&pyA->velocity, &col) > 0 && count > 0){
		printf("				BAAANNNNNNNNNG\n");
	
		vector3D r =  minus(&out, &pyA->position);
		vector3D force = cross(r,  pyA->moment);
	
		acc(&force, &pyA->velocity);
			
		vector3D torque = cross(force, r);
		vector3D delta_velocity = r;
		normalize(&delta_velocity);
		scalar(&delta_velocity, dot(&force, &delta_velocity));
		normalize(&col);//a effacer peut etre
		
		acc(&pyA->velocity , &delta_velocity, -1);

		pyA->angle = 0;

		acc(&pyA->moment, &torque, -0.0002);
		
	}else{
		
	}
	
	return out;

}

char visited[9604];

vector3D colid_fast(physic  *pyA, physic *pyB, colition_hull *A, colition_hull *B,vector3D *prev){

	
	vector3D pA = add(&pyA->position, &pyA->velocity);
	vector3D pB = add(&pyB->position, &pyB->velocity);
	
	//char *visited=new char[A->mesh.nb*B->mesh.nb];
	//memset(visited, 0, A->mesh.nb*B->mesh.nb);

	vector3D col = minus( &pB, &pA);
	vector3D out = vec(0,0,0);

	if (B->mesh.radius + A->mesh.radius < dot(&col, &col)){
		return out;
	}

	vector3D delta_v = minus( &pyA->velocity, &pyB->velocity);

	float min = 10000000;	
	int a, b, currentA=0, currentB=0;
	vector3D dist;

	int count = 0;
	int step = 0;
	int found = 1;

	for (int i=0; i < A->mesh.nb ; i++){
		for (int j=0; j < B->mesh.nb ; j++){
		
			vector3D tp = minus(A->mesh.vertics_trans + i, B->mesh.vertics_trans + j);
			float maybe = dot(&tp, &tp);
			if (maybe < 200){
				acc(&out, A->mesh.vertics_trans + i);
				acc(&out, B->mesh.vertics_trans + j);
				count+=2;
				/*
				min = maybe;
				vert_B_min = j;
				vert_A_max = i;
				*/
			}
		}
	}



	/*
	while (found){
	
		a = currentA*A->faces_per_vertex;
		b = currentB*B->faces_per_vertex;
		visited[currentA*B->mesh.nb + currentB] = 1;
		int begingB = b;
		found = 0;
		
		while (A->vect_graph_table[a] != -1 && a < A->faces_per_vertex*A->mesh.nb){
			b = begingB;
			while (B->vect_graph_table[b] != -1  && b < B->faces_per_vertex*B->mesh.nb){
					step++;
					dist = minus(B->mesh.vertics_trans + B->vect_graph_table[b], A->mesh.vertics_trans + A->vect_graph_table[a]);
					float maybe = dot(&dist, &dist);

					if ( (maybe < min || maybe < 200) && !visited[A->vect_graph_table[a]*B->mesh.nb + B->vect_graph_table[b]]){
			
						min = maybe;
						currentA = A->vect_graph_table[a];
						currentB = B->vect_graph_table[b];
						found = 1;
					}
					//printf("dist %f min %f\n", maybe, min);

					if ( maybe < 200  && !visited[A->vect_graph_table[a]*B->mesh.nb + B->vect_graph_table[b]]){
						acc( &out, A->mesh.vertics_trans + A->vect_graph_table[a]);
						acc( &out, B->mesh.vertics_trans + B->vect_graph_table[b]);
						count +=2;
								
					}	

				b++;		
			}			
			a++;		
		}
		step++;
	}
	*/
	A->current_vert = currentA;
	B->current_vert = currentB;

	//printf("step %d  count%d total%d\n", step, count, B->mesh.nb*A->mesh.nb);
	//delete[] visited;
	if(count > 0) scalar(&out, 1./count);//

	if (dot(&delta_v, &col)> 0 && count){ //dot(&delta_v, &col)>0 &&;
		
		//printf("				BAAANNNNNNNNNG\n");
		//torque a->b
		//out.z = 0;

		vector3D rA =  minus(&out, &pyA->position);
		vector3D angular_velA = cross(rA,  pyA->moment);
		
		
		//torque b->a
		vector3D rB =  minus(&out, &pyB->position);
		vector3D angular_velB = cross(rB, pyB->moment);
		
		vector3D angular_vel = minus(&angular_velA, &angular_velB);

		acc(&angular_vel, &delta_v);

		vector3D torqueA = cross(angular_vel, rA);
		vector3D torqueB = cross(angular_vel, rB);

		normalize(&rA);
		scalar(&rA, dot(&angular_vel, &rA));
		acc(&pyA->velocity , &rA,  -0.8);

		normalize(&rB);
		scalar(&rB, dot(&angular_vel, &rB));

		acc(&pyB->velocity , &rB,   0.8);

		//pyA->angle = 0;

		acc(&pyA->moment, &torqueA, -0.0002);	
		acc(&pyB->moment, &torqueB, 0.0002);	
		
	}else{
		
	}
	/*
	tuple tu;
	tu.a = *(A->mesh.vertics_trans + A->vect_graph_table[a-1]);
	tu.b = *(B->mesh.vertics_trans + B->vect_graph_table[b-1]);
	tu.c = out;
	*/
	return  out;

}


vector3D colid_box(physic  *pyA, physic *pyB, bounding_box A, bounding_box B){
	//append in the local reference of A
	vector3D edge[3];
	vector3D corner, current_edge;
	vector3D out = vec(0,0,0);
	int count = 0;

	float max = -1;
	
	vector3D zero = mult_vm(pyB->mat.data, B.zero);
	zero = inv_vm(pyA->mat.data, zero );
	vector3D zero_limit;
	zero_limit.y = A.zero.y;
	zero_limit.x = A.zero.x;
	zero_limit.z = A.zero.z;
	
	edge[0]=vec( B.limits[0]*pyB->mat.data[1], B.limits[0]*pyB->mat.data[5], B.limits[0]*pyB->mat.data[9]);
	edge[1]=vec( B.limits[1]*pyB->mat.data[0], B.limits[1]*pyB->mat.data[4], B.limits[1]*pyB->mat.data[8]);
	edge[2]=vec( B.limits[2]*pyB->mat.data[2], B.limits[2]*pyB->mat.data[6], B.limits[2]*pyB->mat.data[10]);
	//printf("init \n");
	//printVect(zero);

	for (int i=0; i<3; i++){
		edge[i]  = inv_vm(pyA->mat.data, edge[i]);
		//printVect(edge[i]);
	}
	//printf("corner \n");
	for (int i=0; i<4; i++){
		if (i < 3){
			corner = add(&zero, edge + i);
			acc(&corner, edge + (i + 1)%3);
			//printVect(corner);
		}else{
			corner = zero;
			//printVect(corner);
		}
		//printf("edges \n");
		for (int j=0; j<3; j++){
			if (j==i || j==(i+1)%3){
				current_edge = inv(&edge[j]);
			}else{
				current_edge = edge[j];	
			}
			//if (dot(&current_edge, &pyB->velocity)> 0){
			//printVect(current_edge);
			for(int dim=0; dim<3; dim++){
				float lambda = (A.limits[dim] + zero_limit.el[dim] - corner.el[dim])/(current_edge .el[dim]);
				//printf("lambda %f, limitA %f \n",lambda, A.limits[dim] + zero_limit.el[dim]);
				if (lambda >= 0 && lambda <= 1){
					vector3D colid = add(&corner, &current_edge, lambda);
					int next_dim = (dim + 1)%3;
					if (colid.el[next_dim] > zero_limit.el[next_dim] && colid.el[next_dim] < A.limits[next_dim] + zero_limit.el[next_dim]){
						next_dim = (dim + 2)%3;
						if (colid.el[next_dim] > zero_limit.el[next_dim] && colid.el[next_dim] < A.limits[next_dim] + zero_limit.el[next_dim]){
							acc(&out, &colid);
							count++;
						}
					}
				}

				lambda = (zero_limit.el[dim] - corner.el[dim])/(current_edge.el[dim]);
				if (lambda >= 0 && lambda <= 1){
					vector3D colid = add(&corner, &current_edge, lambda);
					int next_dim = (dim + 1)%3;
					if (colid.el[next_dim] > zero_limit.el[next_dim] && colid.el[next_dim] < A.limits[next_dim] + zero_limit.el[next_dim]){
						next_dim = (dim + 2)%3;
						if (colid.el[next_dim] > zero_limit.el[next_dim] && colid.el[next_dim] < A.limits[next_dim] + zero_limit.el[next_dim]){
							acc(&out, &colid);
							count++;
						}
					}						
				}
				//}
			}
		}
	}
	scalar(&out, 1./count);
	return out;
}



