/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "math.h"
/* include some silly stuff */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "convex2.h"
#include "tesselisation.h"

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



int main () {
	XEvent event;		/* the XEvent declaration !!! */
	KeySym key;		/* a dealie-bob to handle KeyPress Events */	

	int total = 10000;
	init_x();
	int iter = 0;
	mesh3D obj;
	mesh3D ico;
	read_mesh3D(&obj, "test.obj");
	read_mesh3D(&ico, "sphere.obj");
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

	//intersect test 
	vector o = vec(0,50,-50);
	vector r = vec(0, 0, 25);
	float lambda;

	printf("%d \n", plane_intersection(drawing, &n, &o, &r, &lambda));

	//make_convex test
	char *been_there = new char[obj.nb_f/3];
	char *which_shape = new char[obj.nb_f/3];
	memset(which_shape, 0, obj.nb_f/3);
	make_convex(&obj,been_there, which_shape , 0, 1, 100);

	//tessel test
	mesh3D sphere = make_diamond("sphere");
	tessel(&sphere);
	tessel(&sphere);
	tessel(&sphere);
	//make_normal(sphere);

	char *which_shape2 = new char[sphere.nb_f/3];
	
	//edge *edges = make_edge(sphere);

	int debug = 0;

	while(iter < total) {
		XNextEvent(dis, &event);

		if (event.type == KeyPress){

			if (event.xkey.keycode == 86){
				printf("clik");
				memset(which_shape, 0, obj.nb_f/3);
				make_convex(&obj,been_there, which_shape , 0, 1, ++debug);
			};
			if (event.xkey.keycode == 82){
				memset(which_shape, 0, obj.nb_f/3);
				make_convex(&obj,been_there, which_shape , 0, 1, --debug);
			};


		} 


		redraw();
		
		reset(&buffer);
		float a = iter*0.03;

		mat[0] = (cos(a)*200);
		mat[2] = (sin(a)*200);
		
		mat[5] = 200;

		mat[8] = (-sin(a)*200);
		mat[10] = (cos(a)*200);
		//print_mat(mat);

		//mult_vm(mat, obj.vertics, obj.vertics_trans, obj.nb);
		//draw_mesh4(obj, buffer, mat, which_shape);

		//mult_vm(mat, ico.vertics, ico.vertics_trans, ico.nb);
		//draw_mesh4(ico, buffer, mat, which_shape2);

		mult_vm(mat, sphere.vertics, sphere.vertics_trans, sphere.nb);
		draw_mesh4(sphere, buffer, mat, which_shape2);

		XPutImage(dis, win, gc, image, 0, 0, 0, 0, buffer.w, buffer.h);
		//draw(drawing, faces, 1);
		usleep(20000);
		XSync( dis, 0);
		iter++;
	}

	
	return 0;
}


