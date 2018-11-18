#ifndef KDTREE
#define KDTREE
#include <stdlib.h>
#include <stdio.h>
#include "../vector.h"

struct node{
	char plan;
	float limit;
	int point = -1;
	node *left = NULL;
	node *right = NULL;
};
//beuurk
struct vector_id : vector {
	vector_id(void){
	}
	vector_id(vector &in){
		this->x = in.x;
		this->y = in.y;
		this->z = in.z;
	}
	int id;
};

template<int DIM>
int compare(const void * a, const void * b)
{
  if ( ((vector_id*)a)->el[DIM] >((vector_id*)b)->el[DIM] ) return 1;
  if ( ((vector_id*)a)->el[DIM] ==((vector_id*)b)->el[DIM] ) return 0; 
  if ( ((vector_id*)a)->el[DIM] <((vector_id*)b)->el[DIM] ) return -1; 
}




node pool[4000];
node *pool_pointer = pool;

	void build_rec_KD(vector_id *points, node *root,  int nb){

	char plan = root->plan;
	int next_plan = (plan+1)%3;
	
	if(plan == 0){
		qsort(points, nb, sizeof(vector_id), compare<0>);
	}else if (plan == 1){
		qsort(points, nb, sizeof(vector_id), compare<1>);
	}else{
		qsort(points, nb, sizeof(vector_id), compare<2>);
	}
	
	int mid_point = nb/2;
	root->limit = points[mid_point].el[plan];
	root->point = points[mid_point].id;
	
	int nb_right =  nb - mid_point - (nb%2);
	int nb_left =  mid_point - (1 - nb%2);
	root->left = NULL;
	root->right = NULL;
	
	if (nb_left > 0){
		node *left = pool_pointer++;
		root->left = left;
		left->plan = next_plan;
		build_rec_KD(points, left, nb_left);
	}
	
    if (nb_right > 0){
		node *right = pool_pointer++;
		root->right = right;
		right->plan = next_plan;
		build_rec_KD(points + mid_point + nb%2, right,  nb_right);
	}
}

node *build_KD(vector *points, int nb){
	pool_pointer = pool;
	node *root = pool_pointer++;
	float limit = 0;

	vector_id *points_id = new vector_id[nb];
	for(int i=0; i<nb; i++){
		points_id[i]=points[i];
		points_id[i].id = i;
	}
	root->plan = 0;
	
	build_rec_KD(points_id, root, nb);
	delete[] points_id;
	return root;
}

void find(node *root, vector *vertics, int *verif){
	if (root->left){
		find(root->left, vertics, verif);
	}
	if (root->point > -1){
		printVect(vertics[root->point]);
		 //printf("leaf %d \n", root->point);
		verif[root->point]++;
	}

	if (root->right){
		find(root->right, vertics, verif);
	}
}

void line_query(node *root, vector *vertics, float min, float max, char plan, int *count, int depth){
	printf("call %d depth %d \n", (*count)++, depth);
	if (root->plan != plan){
		if (root->left){
		 	line_query(root->left, vertics, min, max, plan, count, depth + 1);
		}
		if (root->right){
		 	line_query(root->right, vertics, min, max, plan, count, depth + 1);
		}
	}else{
		if (root->point > -1 && vertics[root->point].x >= min && vertics[root->point].x <= max) {
			printf("x %f nb call %d\n", vertics[root->point].x, *count); 
			
		}else{
			if (min <= root->limit && root->left){
				 line_query(root->left, vertics, min, max, plan, count, depth + 1);
			}
			if (max >= root->limit && root->right){
				 line_query(root->right, vertics, min, max, plan, count, depth + 1);
			}
		}
	}
}

void linear_box_query(vector *vertics, vector2D limits[3], int nb){
	printVect(limits[0]);
	printVect(limits[1]);
	printVect(limits[2]);
	for (int i=0; i< nb; i++){
		bool inside_box = true;
		
		for (int plan = 0; plan < 3; plan++){
			inside_box = inside_box && (vertics[i].el[plan] >=  limits[plan].x && vertics[i].el[plan] <=  limits[plan].y);
			//printf("vector el %f \n", vertics[i].el[plan]);
		}
		if (inside_box){
			printVect(vertics[i]);
		}
	}
}

void box_query(node *root, vector *vertics, vector2D limits[3], int *count, int depth){

	//printf("call %d depth %d \n", (*count)++, depth);
	
	if (root->point > -1 && 
	limits[0].x <= vertics[root->point].x && 
	limits[0].y >= vertics[root->point].x
	&&
	limits[1].x <= vertics[root->point].y && 
	limits[1].y >= vertics[root->point].y
	&&
	limits[2].x <= vertics[root->point].z && 
	limits[2].y >= vertics[root->point].z){
		
		printVect(vertics[root->point]);
		//printf("x %f nb call %d\n", vertics[root->point].x, *count); 
		
	}else{
		
		if (limits[root->plan].x <= root->limit && root->left){
			 box_query(root->left, vertics, limits, count, depth + 1);
		}
		if (limits[root->plan].y >= root->limit && root->right){
			 box_query(root->right, vertics, limits, count, depth + 1);
		}
	}	
}

void nearest_linear(vector *vertics, vector *point, int nb){
	float min = 1000;
	
	for (int i = 0; i< nb; i++){
		vector dist = minus(vertics + i, point);
		float n = normalize(&dist);
		if(n < min) {
			min = n;
			printf("min %f \n", min);
		}
		
	}
	
}

void nearest(node *root, vector *vertics, vector *point, int *count, int depth, float *dist_min){
	//printf("call %d depth %d \n", (*count)++, depth);
	
	vector dist = minus(point, vertics + root->point);
	float maybe_min = normalize(&dist);
	if (maybe_min < *dist_min){
		*dist_min = maybe_min;
	 	//printf("min %f point %d \n", *dist_min, root->point);	
	}
	if ( point->el[root->plan] > root->limit){
		if (point->el[root->plan] + *dist_min >= root->limit && root->right){
			 nearest(root->right, vertics, point, count, depth + 1, dist_min);
		}
		if (point->el[root->plan] - *dist_min <= root->limit && root->left){
			 nearest(root->left, vertics, point, count, depth + 1, dist_min);
		}
	}else{
		if (point->el[root->plan] - *dist_min <= root->limit && root->left){
			 nearest(root->left, vertics, point, count, depth + 1, dist_min);
		}
		if (point->el[root->plan] + *dist_min >= root->limit && root->right){
			 nearest(root->right, vertics, point, count, depth + 1, dist_min);
		}
	}
		
}


#endif 
