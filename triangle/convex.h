#include "../obj3d.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void make_convex(mesh3D *obj, char *been_there, char *which_shape, int index_first_face, char convex_shape_index){

	if (index_first_face >= obj->nb_f){
		printf("starting index out of range \n");
		return ;
	}

	int* current_face = obj->face + index_first_face;
	vector *current_vertex = obj->vertics + *current_face;
	
	int border_face[1000];
	int border_face_index = 0;
	int *convex_index = new int[obj->nb_f/3];
	vector *convex_vertics = new vector[obj->nb_f/3];
	vector *convex_normal = new vector[obj->nb_f/3];

	int convex_nb = 1;
	int stack_pointer = 1;
	convex_index[0] = index_first_face;
	convex_vertics[0] = obj->vertics[*current_face];
	convex_normal[0] = obj->normals_f[0];

	memset(been_there, 0, obj->nb_f/3);

	printf("%d %d \n", index_first_face, convex_shape_index);
	which_shape[index_first_face] = convex_shape_index;
	been_there[index_first_face] = 1;

	int list_edge_index0[] = {0, 1, 2};
	int list_edge_index1[] = {1, 2, 0};

	int last_face;
	float lambda;

	while(stack_pointer){

		for (int edge = 0; edge < 3; edge++){

			int edge0 = current_face[edge];
			int edge1 = current_face[(edge + 1)%3];

			for (int face = 0; face < obj->nb_f; face += 3){

				if (!been_there[face/3] && !which_shape[face/3] && obj->face + face != current_face){

					for (int i = 0; i < 3; i++){

						if ((obj->face[face + i] == edge0 && obj->face[face + ((i+1)%3)] == edge1) ||
						    (obj->face[face + i] == edge1 && obj->face[face + ((i+1)%3)] == edge0)  ){

							current_vertex = obj->vertics + obj->face[face + ((i+2)%3)];

							been_there[face/3] = 1;
							vector r;
							char is_convex = 0;

							for (int vertex = 0; vertex < convex_nb; vertex++){
								vector r = minus(convex_vertics + vertex, current_vertex);

								for (int conv_face = 0; conv_face < convex_nb; conv_face++){
									if (conv_face != vertex){
										vector *test_vert = convex_vertics + conv_face;
										//printVect(*test_vert);
										is_convex += plane_intersection(test_vert, convex_normal + conv_face, current_vertex, &r, &lambda);

										if (is_convex){
											border_face[ border_face_index++ ] =  face;
											printf("colide %d lambda %f\n", face, lambda);
											break;
										}
									}
								}
								if (is_convex){
									break;
								}
							}

							if(!is_convex){
								which_shape[face/3] = convex_shape_index;
								convex_index[stack_pointer++] = face;
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
		stack_pointer--;
		current_face = obj->face + convex_index[stack_pointer];
		current_vertex = obj->vertics + *current_face;
	}
	printf("last face %d \n", last_face*3);
	//which_shape[last_face] = 0;
	delete[] convex_index;

	for(int i=0; i<obj->nb_f/3; i++){
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
