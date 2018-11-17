/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "math.h"
#include "../anim.h"
/* include some silly stuff */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define H 800
#define W 800
#define COLOR(r, g, b) ((b) | ((g)<<8) | ((r)<<16))
#define MAP_S H*W*3


/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;


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
	int total = 60;
	int start = 0;
	vector p1 = vec(100,500,0);
	vector p2 = vec(300,400,0);
	vector p3 = vec(500,400,0);
	vector p4 = vec(700,500,0);
	
	sequence seq;
	sequence seq2;
	seq.nb = 4;
	seq2.nb = 4;
	keyFrame frames[4];
	keyFrame frames2[4];

	seq.frames = frames;
	seq2.frames = frames2;

	frames[0].p = p1;
	frames[0].temps = 0;
	frames[1].p = p2;
	frames[1].temps = 100;
	frames[2].p = p3;
	frames[2].temps = 200;
	frames[3].p = p4;
	frames[3].temps = 300;

	frames2[0].p = p1;
	frames2[0].temps = 0;
	frames2[1].p = p2;
	frames2[1].temps = 100;
	frames2[2].p = p3;
	frames2[2].temps = 200;
	frames2[3].p = p4;
	frames2[3].temps = 300;

	float frame_list[12] = {100,500,0, 300,400,0, 500,400,0, 700,500,0};
	sequence seq3 = make_sequence((vector*)frame_list, 4, 20);
	make_tangents( &seq3 );

	printf("%f %f %f \n", frames[0].t.x, frames[0].t.y, frames[0].t.z);
	printf("%f %f %f \n", frames[1].t.x, frames[1].t.y, frames[1].t.z);
	printf("%f %f %f \n", frames[2].t.x, frames[2].t.y, frames[2].t.z);
	vector p;
	vector pp;

	while(1) {

		if (start == 0)XNextEvent(dis, &event);
	
		if (event.type==ButtonPress && start == 0) {
			start = 1;
			printf("start %d \n", start);
		}
		if (start){;
			usleep(20000);
			//printf("iter %d \n", iter);
			redraw();
			//interpolate(p,  &seq, iter);
			interpolate(&pp, &seq3, iter);
			//XDrawRectangle(dis, win, gc, p.x , p.y, 10, 10);
			XDrawRectangle(dis, win, gc, pp.x , pp.y, 20, 20);
			XSync(dis, 0);
			iter++;
		}
	}
	printf("%f %f %f \n", p.x, p.y, p.z);
	
	return 0;
}


