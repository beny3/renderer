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

void build_rec_KD(vector *points, node *root, node **fathers_point,  int nb){
	node *left = new node();
	node *right = new node();
	root->left = left;
	root->right = right;
	char plan = root->plan;
	float limit = root->limit;
	int nb_left = 0;
	int nb_right = 0;
	float medianLeft = 0;
	float medianRight = 0;
	int next_plan = (plan+1)%3;

	for (int i=0; i< nb; i++){
		if (fathers_point[i] == root){
			if (points[i].el[plan] < limit){
				medianLeft += points[i].el[next_plan];
				fathers_point[i] = left;
				left->point = i;
				nb_left++;
			}else{
				medianRight += points[i].el[next_plan];
				fathers_point[i] = right;
				right->point = i;
				nb_right++;
			}
		}
	}

	right->limit = medianRight/nb_right;
	right->plan = next_plan;
	left->limit = medianLeft/nb_left;
	left->plan = next_plan;

	if (nb_right > 1){
		right->point = -1;
		build_rec_KD(points, right, fathers_point,  nb);
	}

	if (nb_left > 1){
		left->point = -1;
		build_rec_KD(points, left, fathers_point,  nb);
	}
}

node *build_KD(vector *points, int nb){
	node *root = new node();
	float limit = 0;

	node **fathers_point = new node*[nb];
	for (int i=0; i< nb; i++){
		fathers_point[i] = root;
		limit+= points[i].x;
	}
	limit/=nb;
	root->limit = limit;
	root->plan = 0;
	build_rec_KD(points, root, fathers_point,  nb);
	return root;
}

void erase(node *root){
	if (root->left){
		erase(root->left);
	}

	if (root->right){
		erase(root->right);
	}
	delete[] root;
}

void find(node *root, int *verif){
	if (root->left){
		find(root->left, verif);
	}
	if (root->point > -1){
		 printf("leaf %d \n", root->point);
		verif[root->point]++;
	}

	if (root->right){
		find(root->right, verif);
	}
}

void line_query(node *root, vector *vertics, float line, char plan, int *count){

	//printf("call %d\n", (*count)++);
	if (root->plan != plan){
		if (root->left){
		 	line_query(root->left, vertics, line, plan, count);
		}
		if (root->right){
		 	line_query(root->right, vertics, line, plan, count);
		}
	}else{
		 if ( line < root->limit ){
			if (root->left){
			 	line_query(root->left, vertics, line, plan, count);
			}
		}else{
			if (root->left){
			 	line_query(root->left, vertics, line, plan, count);
			}
			if (root->right){
			 	line_query(root->right, vertics, line, plan, count);
			}
			if (root->point > -1) {
			//printf("x %f \n", vertics[root->point].x); 
			*(count)++;
			}
		}
	}
}
#endif 
