/* include the X library headers */
#include "../x11_wrapper.h"

#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "./triangle.h"
#include "./octree.h"
#include "math.h"
/* include some silly stuff */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void draw_from_list(std::vector<vector3D>& list, mesh3D& dot, bitmap &buffer){
    float mat[16];
	for (auto& v: list){
		rotate4(0, 0, 0, mat);
		scale(mat, 10);
		make_mov4(mat, v);
		mult_vm(mat, dot.vertics , dot.vertics_trans,  dot.nb);
		draw_mesh4(dot, buffer, mat, dot.which_shape);
	}
}


int main () {

	int total = 10000;

	int iter = 0;
	mesh3D obj  = make_triangle_obj("tri");
	mesh3D dot  = make_diamond("dot");
	mesh3D dot2 = make_diamond("dot2");
	//read_mesh3D(&obj, "test.obj");
	//readBMP(&obj.texture, "../toyplane.bmp");

	bitmap buffer(W, H);
	
	init_x(&buffer);

	float mat[16];
	float mat2[16];
	float mat3[16];

	while(iter < total) {

		wait_for_key();
		reset(&buffer);
		float angle = 0.2;//iter*0.03;
		/////////////
		rotate4(0, 0, angle, mat);
		scale(mat, 200);
		//mat[7] = iter;
		////////////
		rotate4(0, 0, 0, mat2);
		scale(mat2, 20);
		////////////
		rotate4(0, 10, 0, mat3);
		scale(mat3, 40);
		
		mult_vm(mat, obj.vertics, obj.vertics_trans, obj.nb);
		triangle tri = make_triangle(&obj, 0);
		vector3D point = vec(300, 100,0);
		vector3D projVec = proj(&tri, &point);
		vector3D proj = triangle_in_box(&tri, &projVec, &point);
		make_mov4(mat2, proj);
		make_mov4(mat3, projVec);
	
		mult_vm(mat3, dot2.vertics, dot2.vertics_trans, dot2.nb);
		mult_vm(mat2, dot.vertics , dot.vertics_trans,  dot.nb);
		draw_mesh4(dot2,buffer, mat3, dot2.which_shape);
		draw_mesh4(dot, buffer, mat2, dot.which_shape);
		draw_mesh4(obj, buffer, mat, obj.which_shape);
		clock_t t1 = clock();
		nb_call = 0;
		Octree octree = make_octree(&obj);
		clock_t t2 = clock();
		printf(" temps octree %f nbcall %d size %u %u \n", (double)(t2 - t1)/CLOCKS_PER_SEC, nb_call, sizeof(Node)*octree.nodes.size(), sizeof(Node));
		std::vector<vector3D> out = dfsOctree(&octree);	
		draw_from_list(out, dot, buffer);
		printf("nb %u \n", out.size());
		putImage(&buffer);
		iter++;
	}

	
	return 0;
}

