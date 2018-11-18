/* include the X library headers */
//#include "Eigen/Dense"

#include "kdtree2.h"
#include "finit3.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <time.h>

#define H 800
#define W 800
#define COLOR(r, g, b) ((b) | ((g)<<8) | ((r)<<16))
#define MAP_S H*W*3


/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;

void draw(int *bitmap){
	for(int i=0; i<W*H*3; i+=3){
		XSetForeground(dis,gc, COLOR(bitmap[i],bitmap[i+1],bitmap[i+2]));
		XDrawPoint(dis, win, gc, (i/3)%W, (i/3)/W);
	}
}

void draw(vector *v, int *f, int nb_f){

	for (int i=0; i<nb_f; i+=3){
		XDrawLine(dis, win, gc, v[f[i]].x , v[f[i]].y, v[f[i+1]].x, v[f[i+1]].y);
		XDrawLine(dis, win, gc, v[f[i+1]].x , v[f[i+1]].y, v[f[i+2]].x, v[f[i+2]].y);
		XDrawLine(dis, win, gc, v[f[i+2]].x , v[f[i+2]].y, v[f[i]].x, v[f[i]].y);
	}
}

void gradient(unsigned  char *bitmap, int iter, int h, int w){
	for(int x=0; x<h; x++){
		for(int y=0; y<w; y++){
			bitmap[4*(y*w + x)]   = (x/2 + iter)%255;
			bitmap[4*(y*w + x)+1] = 0;
			bitmap[4*(y*w + x)+2] = (y/2 + iter)%255;
		}
	}
}

void init_x() {
/* get the colors black and white (see section for details) */        
	unsigned long black,white;

	dis=XOpenDisplay((char *)0);
   	screen=DefaultScreen(dis);
	black=BlackPixel(dis,screen),
	white=WhitePixel(dis, screen);
   	win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0, H, W, 5,black, white);
	XSetStandardProperties(dis,win,"Howdy","Hi",None,NULL,0,NULL);
	XSelectInput(dis, win, ExposureMask|ButtonPressMask|KeyPressMask);
     gc=XCreateGC(dis, win, 0,0);        
	XSetBackground(dis,gc,white);
	XSetForeground(dis,gc,black);
	XClearWindow(dis, win);
	XMapRaised(dis, win);
};

void close_x() {
	XFreeGC(dis, gc);
	XDestroyWindow(dis,win);
	XCloseDisplay(dis);	
	exit(1);				
};

void redraw() {
	XClearWindow(dis, win);
};

//using namespace Eigen;
int main () {
	XEvent event;		/* the XEvent declaration !!! */
	KeySym key;		/* a dealie-bob to handle KeyPress Events */	

	int total = 5000;
	init_x();
	int iter = 0;
	mesh3D obj;
	mesh3D ico;
	read_mesh3D(&obj, "test.obj");
	read_mesh3D(&ico, "tetra_long_super_fin.obj");

	char *which_shape = new char[ico.nb_f/3];
	memset(which_shape, 0, ico.nb_f/3);
	//for(int i=0; i<ico.nb_f; i+=3){printf("face %d %d %d \n",ico.face[i],ico.face[i+1], ico.face[i+2]);}
	//readBMP(&obj.texture, "../toyplane.bmp");

	bitmap buffer;
	buffer.data = (unsigned char*)malloc(H*W*4);
	buffer.zbuffer = (float*)malloc(H*W*sizeof(float));
	buffer.w = W;
	buffer.h = H;
	XImage *image = XCreateImage(dis, DefaultVisual(dis, DefaultScreen(dis)), DefaultDepth(dis, DefaultScreen(dis)), ZPixmap, 0, (char*)buffer.data, buffer.w, buffer.h, 8, 0);

	float mat[16];
	memset(mat,0,16*sizeof(float));
	//mat[0]=1;
	//mat[4]=1;

	vector drawing[3];
	drawing[0]= vec(10,10, 0);
	drawing[1]= vec(100,100, 0);
	drawing[2]= vec(0,100, 0);
	vector n = vec(0, 0.7071, 0.7071);

	int faces[]={0, 1, 2};

	//inverse 2d test
	vector u = vec(1, 0, 0);
	vector v = vec(0.1,-1.8, 4);
	vector p = linear(&u, &v, 3.1, 1.75);
	vector2D uv = solve_uv(&u, &v, &p);
	printVect( p );
	printVect( uv );

	//test svd
	vector2D x = vec2D(1,2);
	vector2D y = vec2D(4,-1);
	printf("svd test \n");
	svd(&x, &y);

	//intersect test 
	vector o = vec(0,50,-50);
	vector r = vec(0, 0, 25);
	float lambda;

	printf("%d \n", plane_intersection(drawing, &n, &o, &r, &lambda));

	//test KD tree
	int *verif = new int[ico.nb]();
	int count=0;
	//for (float l=-0.3; l < 0.4; l+= 0.001){
		float dist = 2;
		vector test_point = vec(0.001158, -0.430398, -0.53);

		clock_t start = clock();

		node *root = build_KD(ico.vertics, ico.nb);
		//box_query(root, ico.vertics, (vector2D*)limit, &count, 0);
		//printf("\n*************\n linear search \n");
		//linear_box_query(ico.vertics, (vector2D*)limit, ico.nb);
		nearest(root, ico.vertics, &test_point, &count, 0, &dist);

		clock_t end = clock();

		printf("temps kd %f %d\n", (float)(end -start)/CLOCKS_PER_SEC, count);
		printf("linear nearest search \n");
		nearest_linear( ico.vertics, &test_point, ico.nb);
		//erase(root); 
	//}

	//find(root, verif);
	//for(int i=0; i<ico.nb; i++)printf("%d \n", verif[i]);
	//test gaussbeny@beny:~/renderer/finit$ beny@beny:~/renderer/finit$ 
	float mat9[9] = {10,2,3,40,15,6,7,8,19};
	///gauss(mat9);

	//test constraint
	float a = 0.1;

	mat[0] = (cos(a));
	mat[1] = (sin(a));
		
	mat[10] = 1;

	mat[4] = (-sin(a));
	mat[5] = (cos(a));

	int *nb_constraints = new int[ico.nb]();
	vector *next_velocities = new vector[ico.nb];
	vector *velocities = new vector[ico.nb];
	vector *momentum = new vector[ico.nb];
	set_zero(next_velocities, ico.nb);
	set_zero(velocities, ico.nb);

	edge_constraints edge_constraint = find_edge(&ico, nb_constraints);
	//constraint *constraints = create_constraint( &ico, nb_constraints);

	vector target1 = ico.vertics[129];
	vector target2 = ico.vertics[6];
	//target1.x += 0.5;
	target2.x += 0.1;
	//mult_vm(mat, ico.vertics, ico.nb);
	vector *projected_points = new vector[ico.nb];

	/*
	for (int i = 0; i < 2; i++){ 
 		project_constraint(&ico, constraints, projected_points);
		global_solve(&ico, projected_points, nb_constraints);
	}
	*/

	for (int i=0; i<ico.nb; i++){
		ico.vertics_trans[i].x = ico.vertics[i].x*200;
		ico.vertics_trans[i].y = ico.vertics[i].y*200;
		ico.vertics_trans[i].z = ico.vertics[i].z*200;
	}

	draw_mesh4(ico, buffer, mat, which_shape);
	/*
	Matrix3f A;
	A << 10, 2, 3, 4, 15, 6, 7, 8, 19;
	std::cout << A << "\n";
	*/
	int debug = 0;
	make_scale(mat, 400);
	vector dir = vec(200, -400, 0);
	make_mov4(mat, dir);

	while(iter < total) { //total 
		
		//XNextEvent(dis, &event);
		/*
		if (event.type == KeyPress){

			if (event.xkey.keycode == 86){
			//
			};
			if (event.xkey.keycode == 82){
			//
			};

		} 
		*/
		//redraw();
		//printVect(velocities, ico.nb);
		//target.x = startx + 0.3*sin(iter*0.01);
		reset(&buffer);
		node *root = build_KD(ico.vertics, ico.nb);
		start = clock();
		dist = 2;
		count = 0;
		nearest(root, ico.vertics, &test_point, &count, 0, &dist);
		end = clock();
		printf("temps kd %f %d\n", (float)(end -start)/CLOCKS_PER_SEC, count);
		start = clock();
		nearest_linear( ico.vertics, &test_point, ico.nb);
		end = clock();
		printf("temps linear %f %d\n", (float)(end -start)/CLOCKS_PER_SEC, count);
		if(iter > 100){
			clock_t start = clock();

			copy(ico.vertics, next_velocities, ico.nb);
			vector g = vec(0,9,0);
			acc_one_to_many(ico.vertics, &g, 0.0004, ico.nb);
			acc(ico.vertics, velocities, 0.999, ico.nb);
			copy(ico.vertics, momentum, ico.nb);

			for (int i = 0; i < 80; i++){ 
				//printf("compute solver\n");
				set_zero(projected_points, ico.nb);
				//project_constraint(&ico, constraints, projected_points);
				project_edges(&ico, edge_constraint, projected_points );
				global_solve(&ico, projected_points, nb_constraints, momentum);
				
				
				//ico.vertics[129] = target1;
				//ico.vertics[129].x = target1.x + sin((iter-5)*0.03);
				ico.vertics[6] = target2;
			}
			clock_t end = clock();
			double seconds = (double)(end-start)/ CLOCKS_PER_SEC;
			//printf("temps %f \n", seconds);
			make_normal(&ico);
			//next_velocities[1].x = 1/DT; //4*(startx + sin(iter*0.003) - ico.vertics[1].x)/DT;


		//printVect(ico.vertics, ico.nb);
		//next_velocities[3].y = ((starty + cos(iter*0.003) - ico.vertics[3].x))/DT;
		//printf("vel %f \n", next_velocities[3].x );

		get_momentum(ico.vertics, next_velocities, ico.nb);
		swap(&velocities, &next_velocities);
		}

		mult_vm(mat, ico.vertics, ico.vertics_trans, ico.nb);
		draw_mesh4(ico, buffer, mat, which_shape);

		XPutImage(dis, win, gc, image, 0, 0, 0, 0, buffer.w, buffer.h);
		//draw(drawing, faces, 1);
		usleep(1500);
		XSync( dis, 0);
		//printf("iter %d \n", iter);
		iter++;
	}

	
	return 0;
}


