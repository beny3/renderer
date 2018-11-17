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
	read_mesh3D(&obj, "test.obj");
	//readBMP(&obj.texture, "../toyplane.bmp");

	bitmap buffer;
	buffer.data = (unsigned char*)malloc(H*W*4);
	buffer.zbuffer = (float*)malloc(H*W*sizeof(float));
	buffer.w = W;
	buffer.h = H;
	XImage *image = XCreateImage(dis, DefaultVisual(dis, DefaultScreen(dis)), DefaultDepth(dis, DefaultScreen(dis)), ZPixmap, 0, (char*)buffer.data, buffer.w, buffer.h, 8, 0);

	float mat[9];
	float mat4[16];

	memset(mat, 0, 9*sizeof(float));
	memset(mat4, 0, 16*sizeof(float));
	mat4[10]=1;
	mat4[15]=1;

	//test inverse
	float angle = PI/7;

	mat4[0]=(cos(angle));
	mat4[4]=(sin(angle));

	mat4[1]=(-sin(angle));
	mat4[5]=(cos(angle));
	vector dir = vec(100, 100, 5.3);

	vector u = vec(4,1,3);
	normalize(&u);
/*
	make_mov4(mat4, dir);
	make_rot4_u(mat4, &u, angle);
	print_mat(mat4);
	matrix inv_mat = inv_m(mat4);
	print_mat(inv_mat.data);
	matrix ide = mult_m(mat4 , inv_mat.data);
	print_mat(ide.data);
*/
	char *which_shape = new char[obj.nb_f/3];
	memset(which_shape, 0, obj.nb_f/3);

	while(iter < total) {
		XNextEvent(dis, &event);

		if (event.type == KeyPress){

			if (event.xkey.keycode == 86){
				
			};
			if (event.xkey.keycode == 82){

			};
		} 

		redraw();
		reset(&buffer);
		angle = iter*0.03;

		mat[0]=(cos(angle)*70);
		mat[5]=(sin(angle)*70);
		
		mat[5]=70;

		mat[2]=(-sin(angle)*70);
		mat[10]=(cos(angle)*70);

		matrix model = make_rot4( 0, PI/4, 0);
		dir = vec(0.2, 0.2, 0);
		make_mov4(model.data , dir);

		mult_vm(mat, obj.vertics, obj.vertics_trans, obj.nb);
		//make_ide4(mat4);
		print_mat(model.data);
		draw_cube(model.data, buffer);
		//draw_mesh4(obj, buffer, mat, which_shape);
		XPutImage(dis, win, gc, image, 0, 0, 0, 0, buffer.w, buffer.h);
		//draw(drawing, faces, 1);
		usleep(20000);
		XSync( dis, 0);
		iter++;
	}

	
	return 0;
}

