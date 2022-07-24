/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "vector.h"
#include "bitmap.h"
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
XImage* image;
XEvent event;		/* the XEvent declaration !!! */
KeySym key;		/* a dealie-bob to handle KeyPress Events */	

void draw(int *bitmap){
	for(int i=0; i<W*H*3; i+=3){
		XSetForeground(dis,gc, COLOR(bitmap[i],bitmap[i+1],bitmap[i+2]));
		XDrawPoint(dis, win, gc, (i/3)%W, (i/3)/W);
	}
}

void draw(vector3D *v, int *f, int nb_f){

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

void init_x(bitmap* buffer) {
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
	image = XCreateImage(dis, DefaultVisual(dis, DefaultScreen(dis)), DefaultDepth(dis, DefaultScreen(dis)), ZPixmap, 0, (char*)buffer->data, buffer->w, buffer->h, 8, 0);
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

void  putImage(bitmap *buffer){
	XPutImage(dis, win, gc, image, 0, 0, 0, 0, buffer->w, buffer->h);
}

void wait_for_key(){
    XNextEvent(dis, &event);

	if (event.type == KeyPress){

		if (event.xkey.keycode == 86){
			
		};
		if (event.xkey.keycode == 82){

		};
	} 
}
