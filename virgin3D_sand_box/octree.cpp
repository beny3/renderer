/* include the X library headers */
#include "../x11_wrapper.h"

#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "math.h"
/* include some silly stuff */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main () {

	int total = 10000;

	int iter = 0;
	mesh3D obj = make_triangle("tri");
	//read_mesh3D(&obj, "test.obj");
	//readBMP(&obj.texture, "../toyplane.bmp");

	bitmap buffer(W, H);
	
	init_x(&buffer);

	float mat[16];

	while(iter < total) {

		wait_for_key();
		reset(&buffer);
		float angle = iter*0.03;
		rotate4(0, 0, angle, mat);
		scale(mat, 200);
		mult_vm(mat, obj.vertics, obj.vertics_trans, obj.nb);
        //print_mat(mat);
		draw_mesh4(obj, buffer, mat, obj.which_shape);
		putImage(&buffer);
		iter++;
	}

	
	return 0;
}

