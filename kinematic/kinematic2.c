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
#define GETX(n) 2*(n)
#define GETY(n) 2*(n) + 1
#define N 2

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
			bitmap[4*(y*w + x)]=(x/2 + iter)%255;
			bitmap[4*(y*w + x)+1]=0;
			bitmap[4*(y*w + x)+2]=(y/2 + iter)%255;
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
	char text[255];
	
	init_x();
	int x=0;
	int y=0;
	int n = N;
	float boneA[]={0,100,100,100,200,100};
	vector2D *bones = (vector2D*)boneA;
	vector2D norm[N+1];
	int first = 0;

	for (int i=1; i<=n; i++){

		XSetForeground(dis,gc, 0xFF0000);
		XDrawLine(dis, win, gc, (int)bones[i].x, (int)bones[i].y, (int)bones[i-1].x, (int)bones[i-1].y);
	}	
	/* look for events forever... */
	while(1) {		
		/* get the next event and stuff it into our event variable.
		   Note:  only events we set the mask for are detected!
		*/
		XNextEvent(dis, &event);
	
	
		if (event.type==ButtonPress) {
			//redraw();
			if (first){
				for (int converge=0; converge<5; converge++){
				bones[n].x= event.xbutton.x;
				bones[n].y= event.xbutton.y;
				
				
				for (int i=n; i>1; i--){
					
					norm[i].x =  (bones[i].y - bones[i].y);
					norm[i].y =  (bones[i-1].x - bones[i].x);
					
					vector2D bone;
					
					bone.x = bones[i-1].x  - bones[i].x;
					bone.y = bones[i-1].y  - bones[i].y;
					float proj_n = dot2D( &norm[i], &bone);
					if(proj_n < 0){
						bone.x -= 2*proj_n*norm[i].x;
						bone.y -= 2*proj_n*norm[i].y;
					}

					normalize2D(&bone);
					scalar2D(&bone,100);
					bones[i-1]=add2D(&bones[i], &bone);

				}
				for (int i=0; i<n; i++){
						
					vector2D bone;
					bone.x = bones[i+1].x - bones[i].x;
					bone.y = bones[i+1].y - bones[i].y;
					normalize2D(&bone);
					scalar2D(&bone,100);
					bones[i+1]=add2D(&bones[i], &bone);
				}
				}	
			}
			first = 1;
			redraw();
			for (int i=1; i<=n; i++){
				

				XSetForeground(dis,gc, 0x0000FF);
				XDrawLine(dis, win, gc, bones[i].x, bones[i].y, bones[i].x + norm[i].x, bones[i].y + norm[i].y);

				XSetForeground(dis,gc, 0xFF0000);
				XDrawLine(dis, win, gc, (int)bones[i].x, (int)bones[i].y, (int)bones[i-1].x, (int)bones[i-1].y);
			}				
		}
	}


	
	return 0;
}


