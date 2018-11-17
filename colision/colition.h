
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
	vector a;
	vector b;
	vector c;
};

struct Simplex{
	vector D;
	vector data[4];
	int s_last = -1;
};

float distance(vector *A, vector *D, vector *Dn){
	vector AXD = cross(*A, *Dn);
	return dot(D, A);// + (1./sqrt(dot(&AXD, &AXD)));
};

vector support2(vector *D, vector *A, int na, vector *B, int nb){
	float max = -100000000;
	vector Dn = *D;
	normalize(&Dn);
	vector *AB = new vector[nb*na];
	vector *arMax = AB;

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

vector support(vector *D, vector *A, int na, vector *B, int nb){
	vector D_ = inv(D);
	vector Dn = *D;
	normalize(&Dn);
	vector Dn_ = inv(&Dn);

	float s  = distance(A, D, &Dn);
	float max = s;

	vector *argmaxA = A; 
	vector *argmaxB = B;

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
	vector *simplex = simplexS->data;

	if(s_last == -1){
		//simplexS->D = vec(0,1,0);
		simplexS->s_last = 0;

	}else if(s_last == 0){
		simplexS->D = inv(&simplex[0]);
		simplexS->s_last = 1;

	}else if(s_last == 1){
	
		vector ab = minus( &simplex[0], &simplex[1]);
		vector ao = inv(&simplex[1]);

		if ( dot(&ab, &ao) >= 0 ){
			vector abXao = cross(ab, ao);
			simplexS->D = cross(abXao, ab);
			simplexS->s_last = 2;
		}else{
			printf("of_line \n");
			simplexS->D = ao;
			simplexS->s_last = 1;
		}

	}else if (s_last == 2){

		//cas triangle
		vector ao = inv(&simplex[2]);
		vector ab = minus(&simplex[1], &simplex[2]);
		vector ac = minus(&simplex[0], &simplex[2]);

		vector n = cross(ab, ac);
		vector nXac = cross(n, ac);
		vector abXn = cross(ab, n);

		//outside plan aco
		if (dot(&nXac, &ao) > 0){
			if ( dot(&ac, &ao) > 0 ){
				vector acXao = cross(ac, ao);
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
				vector abXao = cross(ab, ao);
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
		
		vector bc = minus(&simplex[1], &simplex[0]);
		vector ab = minus(&simplex[2], &simplex[1]);
		vector ac = minus(&simplex[0], &simplex[2]);

		vector top_to_base = minus(&simplex[1], &simplex[3]);
		vector tpX = cross(bc, top_to_base, *sign_n);
		
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


int gjk(vector *A, int na, vector *B, int nb, Simplex *simplex, vector *AB){
	
	vector D_prev = simplex->D;
	/*
	vector invS;
	vector ab;
	vector ac;
	vector n;
	vector ao;
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



vector colid2(physic  *pyA, physic *pyB, colition_hull *A, colition_hull *B,vector *prev){

	vector pA = add(&pyA->position, &pyA->velocity);
	vector pB = add(&pyB->position, &pyB->velocity);

	
	vector out = vec(0,0,0);
	vector *normalB = B->mesh.normals;

	vector col = minus( &pB, &pA);
	vector tpa, tpb, tp;

	// vector vert_A_max =0; 
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

	//vector out = linear((A->mesh.vertics_trans + vert_A_max), (B->mesh.vertics_trans + vert_B_min), 0.5, 0.5);
	//vector dist = minus(A->mesh.vertics_trans + vert_A_max, B->mesh.vertics_trans + vert_B_min);

	if (dot(&pyA->velocity, &col) > 0 && count > 0){
		printf("				BAAANNNNNNNNNG\n");
	
		vector r =  minus(&out, &pyA->position);
		vector force = cross(r,  pyA->moment);
	
		acc(&force, &pyA->velocity);
			
		vector torque = cross(force, r);
		vector delta_velocity = r;
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

vector colid_fast(physic  *pyA, physic *pyB, colition_hull *A, colition_hull *B,vector *prev){

	
	vector pA = add(&pyA->position, &pyA->velocity);
	vector pB = add(&pyB->position, &pyB->velocity);
	
	//char *visited=new char[A->mesh.nb*B->mesh.nb];
	//memset(visited, 0, A->mesh.nb*B->mesh.nb);

	vector col = minus( &pB, &pA);
	vector out = vec(0,0,0);

	if (B->mesh.radius + A->mesh.radius < dot(&col, &col)){
		return out;
	}

	vector delta_v = minus( &pyA->velocity, &pyB->velocity);

	float min = 10000000;	
	int a, b, currentA=0, currentB=0;
	vector dist;

	int count = 0;
	int step = 0;
	int found = 1;

	for (int i=0; i < A->mesh.nb ; i++){
		for (int j=0; j < B->mesh.nb ; j++){
		
			vector tp = minus(A->mesh.vertics_trans + i, B->mesh.vertics_trans + j);
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

		vector rA =  minus(&out, &pyA->position);
		vector angular_velA = cross(rA,  pyA->moment);
		
		
		//torque b->a
		vector rB =  minus(&out, &pyB->position);
		vector angular_velB = cross(rB, pyB->moment);
		
		vector angular_vel = minus(&angular_velA, &angular_velB);

		acc(&angular_vel, &delta_v);

		vector torqueA = cross(angular_vel, rA);
		vector torqueB = cross(angular_vel, rB);

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


vector colid_box(physic  *pyA, physic *pyB, bounding_box A, bounding_box B){
	//append in the local reference of A
	vector edge[3];
	vector corner, current_edge;
	vector out = vec(0,0,0);
	int count = 0;

	float max = -1;
	
	vector zero = mult_vm(pyB->mat.data, B.zero);
	zero = inv_vm(pyA->mat.data, zero );
	vector zero_limit;
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
					vector colid = add(&corner, &current_edge, lambda);
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
					vector colid = add(&corner, &current_edge, lambda);
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



