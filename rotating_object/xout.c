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

	
	init_x();
	int iter = 0;
	mesh3D obj;
	read_mesh3D(&obj, "toyplane.obj");
	readBMP(&obj.texture, "toyplane.bmp");
	

	bitmap buffer;
	buffer.data = (unsigned char*)malloc(H*W*4);
	buffer.zbuffer = (float*)malloc(H*W*sizeof(float));
	buffer.w = W;
	buffer.h = H;

	reset(&buffer);

	vector2D face[3];
	face[0].x = 20;
	face[0].y = 20;

	face[1].x = 20;
	face[1].y = 300;

	face[2].x = 400;
	face[2].y = 300;
	float mat[9];
	memset(mat,0,9*sizeof(float));
	mat[0]=1;
	mat[4]=1;
	
	

	//draw_mesh(obj, buffer, mat);
	XImage *image = XCreateImage(dis, DefaultVisual(dis, DefaultScreen(dis)), DefaultDepth(dis, DefaultScreen(dis)), ZPixmap, 0, (char*)buffer.data, buffer.w, buffer.h, 8, 0);
	double seconds = 0;
	/* look for events forever... */
	int total = 500;
	clock_t start = clock();

	while(iter < total) {		
		/* get the next event and stuff it into our event variable.
		   Note:  only events we set the mask for are detected!
		*/
	
		/* the window was exposed redraw it! */
		
		
		//XNextEvent(dis, &event);
		//gradient(buffer.data, iter, buffer.h, buffer.w);
		//usleep(4000);

		reset(&buffer);
		float a = iter*0.03;
		float c = 1.8;

		mat[0]=(cos(a)*0.2);
		mat[6]=(sin(a)*0.2);
		
		mat[4]=1*0.2;

		mat[2]=(-sin(a)*0.2);
		mat[8]=(cos(a)*0.2);

		matrix mat2 = make_rot4( a, PI/3, 0);
		matrix mat3 = mat2;	
		matrix mat4 = mat2;
		scale(mat2.data, 20);
		vector pos = vec(200, 0, 0);
		make_mov4(mat4.data, pos);

		mult_vm(mat2.data, obj.vertics, obj.vertics_trans, obj.nb);
		draw_mesh4(obj, buffer, mat3.data);

		//mult_vm(mat4.data, obj.vertics, obj.vertics_trans, obj.nb);
		//draw_mesh4(obj, buffer, mat3.data);

		XPutImage(dis, win, gc, image, 0, 0, 0, 0, buffer.w, buffer.h);
		//redraw();
		/*
		if (event.type==Expose && event.xexpose.count==0) {
			
			redraw();
		}
		
		if (event.type==ButtonPress) {
		// tell where the mouse Button was Pressed 
			int x=event.xbutton.x,
			    y=event.xbutton.y;

			XSetForeground(dis,gc, COLOR(255,0,255));
			XPutImage(dis, win, gc, image, 0, 0, 0, 0, w, h);
			//XFillRectangle(dis, win, gc, x, y, 20, 20);
			//XDrawPoint(dis, win, gc, x, y);
		}
		*/
		iter++;
	}

	clock_t end = clock();
	seconds = (double)(end - start) / CLOCKS_PER_SEC;
	printf("temps %f \n", seconds/total);
	
	return 0;
}


