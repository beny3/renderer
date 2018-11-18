#ifndef KDTREE
#define KDTREE
#include <stdlib.h>
#include <stdio.h>
#include "../vector.h"

struct node{
	char plan;
	float limit;
	int id;
	int point = -1;
	node *left = NULL;
	node *right = NULL;
};

node pool[4000];
node *pool_pointer = pool;

void build_rec_KD(vector *points, node *root, node **fathers_point,  int nb){

	int stk_ptr = 0;
	node *stk[200];
	stk[stk_ptr++]=root;

	while(stk_ptr){
		root = stk[--stk_ptr];

		node *left = pool_pointer++;
		node *right = pool_pointer++;
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
			stk[stk_ptr++] = right;
			//build_rec_KD(points, right, fathers_point,  nb);
		}

		if (nb_left > 1){
			left->point = -1;
			stk[stk_ptr++] = left;
			//build_rec_KD(points, left, fathers_point,  nb);
		}
	}
}

node *build_KD(vector *points, int nb){
	for (int i=0; i<4000; i++){
		pool[i].id = i;
	}

	pool_pointer = pool;
	node *root = pool_pointer++;
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
	printf("call %d depth %d \n", (*count)++, depth);

	float min = *dist_min;

	node *stack[200];
	float test_stack[200];

	int stk_pointer = 0;
	test_stack[stk_pointer] = -1;
	stack[stk_pointer] = root;
	stk_pointer++;

	float left_condition  = point->el[root->plan] - root->limit;
	float right_condition = root->limit - point->el[root->plan];

	int k = 0;

	int max_stack = 0;

	while(stk_pointer){
		k++;
		do{
			root = stack[--stk_pointer];
		}while(test_stack[stk_pointer] > min && stk_pointer);

		///if(stk_pointer > max_stack)max_stack = stk_pointer;

		if (root->point > -1){
			vector dist = minus(point, vertics + root->point);
			float maybe_min = normalize(&dist);
			if(maybe_min < min){
				min = maybe_min;
				printf("min %f point %d  call %d\n", min, root->point, k);
			}
		
		}else{
			if ( point->el[root->plan] > root->limit){

				if (left_condition <= min && root->left){
					stack[stk_pointer] = root->left;
					test_stack[stk_pointer] = left_condition;
					stk_pointer++;
				}
				if (right_condition <= min && root->right){
					stack[stk_pointer] = root->right;
					test_stack[stk_pointer] = right_condition;
					stk_pointer++;
				}

			}else{
				if (right_condition <= min && root->right){
					stack[stk_pointer] = root->right;
					test_stack[stk_pointer] = right_condition;
					stk_pointer++;
				}
				if (left_condition <= min && root->left){
					stack[stk_pointer] = root->left;
					test_stack[stk_pointer] = left_condition;
					stk_pointer++;
				}
			}
		}
	}
	//printf("stack_height %d \n", max_stack);
}


#endif 
