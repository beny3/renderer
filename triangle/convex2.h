#include "../obj3d.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void make_convex(mesh3D *obj, char *been_there, char *which_shape, int index_first_face, char convex_shape_index, int maxNb){

	if (index_first_face >= obj->nb_f){
		printf("starting index out of range \n");
		return ;
	}

	int* current_face = obj->face + index_first_face;
	vector *current_vertex = obj->vertics + *current_face;
	
	int border_face[1000];
	int border_face_index = 0;

	int nb_f = obj->nb_f/3;

	int *convex_index = new int[nb_f];
	vector *convex_vertics = new vector[nb_f];
	vector *convex_normal = new vector[nb_f];

	int convex_nb = 1;
	int head = 1;
	int tail = 0;
	convex_index[0] = index_first_face;
	convex_vertics[0] = obj->vertics[*current_face];
	convex_normal[0] = obj->normals_f[0];

	memset(been_there, 0, nb_f);

	//printf("%d %d \n", index_first_face, convex_shape_index);
	which_shape[index_first_face/3] = convex_shape_index;
	been_there[index_first_face/3] = 1;

	int list_edge_index0[] = {0, 1, 2};
	int list_edge_index1[] = {1, 2, 0};

	int last_face;
	float lambda;

	while(head != tail && convex_nb <= maxNb){

		for (int edge = 0; edge < 3; edge++){

			int edge0 = current_face[edge];
			int edge1 = current_face[(edge + 1)%3];

			//find the neigbouring points of current face
			for (int face = 0; face < obj->nb_f; face += 3){
				if (!been_there[face/3] && !which_shape[face/3] && obj->face + face != current_face){

					for (int i = 0; i < 3; i++){

						if ((obj->face[face + i] == edge0 && obj->face[face + ((i+1)%3)] == edge1) ||
						    (obj->face[face + i] == edge1 && obj->face[face + ((i+1)%3)] == edge0)  ){

							current_vertex = obj->vertics + obj->face[face + ((i+2)%3)];
							been_there[face/3] = 1;
							vector r;

							//find if we add the current point, the shape would stil be convex
							char is_convex = 0;
							int vertex = 0;
							for (vertex = 0; vertex < convex_nb; vertex++){
								r = minus(current_vertex ,convex_vertics + vertex);
								is_convex += (dot(&r,  convex_normal + vertex) > 0.0001);

								if (is_convex){
									border_face[ border_face_index++ ] =  face;
									//printf("colide %d lambda %f\n", face, lambda);
									break;
								}
							}
							if(is_convex){
								vertex = convex_nb;
							}else{
								vertex = 0;
							}

							for (vertex; vertex < convex_nb; vertex++){
								r = minus(convex_vertics + vertex, current_vertex);
								is_convex += (dot(&r,  obj->normals_f + face/3) > 0.0001);

								if (is_convex){
									border_face[ border_face_index++ ] =  face;
									//printf("colide %d lambda %f\n", face, lambda);
									break;
								}
							}


							if(!is_convex){
								//printf("push_face %d \n", face);
								last_face = face;
								which_shape[face/3] = convex_shape_index;
								convex_index[head] = face;
								head = (head + 1)%nb_f;

								convex_vertics[convex_nb] = *current_vertex;
								convex_normal[convex_nb] = obj->normals_f[face/3];
								convex_nb++;
							}
							break;
						}
					} 
				}
			}
		}
		//take the next current face from the queue
		current_face = obj->face + convex_index[tail];
		current_vertex = obj->vertics + *current_face;
		tail = (tail+1)%nb_f;
	}
	//printf("last face %d \n", last_face*3);
	//which_shape[last_face] = 0;
	delete[] convex_index;

	for(int i=0; i<nb_f; i++){
		printf("%d %d \n", 3*i, which_shape[i]);
	}
	printf(" \n");

	for(int i=0; i<border_face_index; i++){
		printf("%d %d \n", i, border_face[i]);
	}
	printf(" \n");

	printf("nb_border %d \n", border_face_index);

	printf("nb convex %d %d\n",convex_nb, border_face[ border_face_index-1]);

	if (border_face_index){
		//make_convex(obj, been_there, which_shape, border_face[ border_face_index-1]/3, convex_shape_index + 1);
	}
}
