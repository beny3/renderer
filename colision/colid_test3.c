/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "colition.h"
#include <time.h>
#include "../3d.h"
#define H 800
#define W 800
#define COLOR(r, g, b) ((b) | ((g)<<8) | ((r)<<16))
#define MAP_S H*W*3
#define NBS 4
#define NB 16

/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;

void draw(mesh3D *obj){
	vector *v = obj->vertics_trans;
	int *f = obj->face;
	int nb_f = obj->nb_f;
	double d = 500;

	for (int i=0; i<nb_f; i+=3){
		double z1=v[f[i]].z+100, z2=v[f[i+1]].z+100, z3=v[f[i+2]].z+100;
		XDrawLine(dis, win, gc, v[f[i]].x*d/(d + z1)   + H/2, v[f[i]].y*d/(d + z1)  + W/2, v[f[i+1]].x*d/(d + z2) + W/2, v[f[i+1]].y*d/(d + z2) + H/2);
		XDrawLine(dis, win, gc, v[f[i+1]].x*d/(d + z2) + H/2, v[f[i+1]].y*d/(d + z2)+ W/2, v[f[i+2]].x*d/(d + z3) + W/2, v[f[i+2]].y*d/(d + z3) + H/2);
		XDrawLine(dis, win, gc, v[f[i+2]].x*d/(d + z3) + H/2, v[f[i+2]].y*d/(d + z3)+ W/2, v[f[i]].x*d/(d + z1)   + W/2, v[f[i]].y*d/(d + z1)   + H/2);
	}
}

void drawNormal(mesh3D *obj, int shift){
	vector *v = obj->vertics_trans;
	vector *normals = obj->normals;
	int nb = obj->nb;
	XSetForeground(dis,gc, 0x4444FF);
	for (int i=0; i<nb; i++){
		XDrawLine(dis, win, gc, v[i].x + shift, v[i].y + shift, v[i].x + normals[i].x*15 + shift, v[i].y + normals[i].y*15 + shift);

	}
}

void draw(vector *v, int nb, int shift){

	for (int i=0; i<nb; ++i){
		XFillRectangle(dis, win, gc, v[i].x + shift, v[i].y + shift, 5, 5);
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
	mesh3D rect[NB];
	physic phy[NB];
	colition_hull hull[NB];

	bitmap buffer;
	buffer.data = (unsigned char*)malloc(H*W*4);
	buffer.zbuffer = (float*)malloc(H*W*sizeof(float));
	buffer.w = W;
	buffer.h = H;

	reset(&buffer);

	for (int i=0; i<NB; i++){
	
		read_mesh3D(&rect[i], "cube3.obj");
		phy[i].position = vec( (i%NBS - NBS/2)*W/NBS, (i/NBS - NBS/2)*H/NBS  , 0);
		phy[i].velocity = vec(0, 0,0);
		phy[i].moment = vec(0, 0, 0);
		//phy[i].moment = vec(0, PI/16*((i+1)%2), PI/16*(i%2));
		phy[i].angle = 0;
		hull[i]= make_graph(rect[i]);

	}
	phy[0].velocity = vec(10, 13,0);
	
	//printf("colid_test %d \n", gjk(rectA.vertics, rectA.nb, rectB.vertics, rectB.nb) );

	/* look for events forever... */
	matrix m2 = make_rot4(0.,0.,PI/2);
	matrix m22 = make_rot4(0.,PI/2,0.);
	matrix m = make_rot4(1,2,3);
	matrix eye;
	make_ide4(eye.data);
	float phi = 0;
	float rho = PI/2;
	vector c = vec(0,0,0);

	XImage *image = XCreateImage(dis, DefaultVisual(dis, DefaultScreen(dis)), DefaultDepth(dis, DefaultScreen(dis)), ZPixmap, 0, (char*)buffer.data, buffer.w, buffer.h, 8, 0);
	int frame=0;
	
	clock_t start = clock();
	while(frame<200) {
	
		//XNextEvent(dis, &event);

		for (int i=0; i<NB; i++){
			euler(&phy[i], 1);
			mult_vm(phy[i].mat.data, rect[i].vertics, rect[i].vertics_trans, rect[i].nb);
	
			XSetForeground(dis,gc, 0xFF0000);
			draw_mesh4(hull[i].mesh, buffer, phy[i].mat.data);

			if (phy[i].position.x < -W/2 || phy[i].position.x > W/2 )phy[i].velocity.x *= -1;
			if (phy[i].position.y < -H/2 || phy[i].position.y > H/2 )phy[i].velocity.y *= -1;
			if (phy[i].position.z < -W/8 || phy[i].position.z > W/8 )phy[i].velocity.z *= -1;
				
		}
		
		for (int i=0; i<NB; i++){
			for (int j=i+1; j<NB; j++){
				vector c = colid_fast(&phy[i], &phy[j], &hull[i], &hull[j], &c);
				
			}
		}		
		//printf("temps %f \n", seconds);

		XPutImage(dis, win, gc, image, 0, 0, 0, 0, buffer.w, buffer.h);	
		reset(&buffer);
		frame++;
	}
	clock_t end = clock();
	double seconds = (double)(end - start) / CLOCKS_PER_SEC;
	printf("temps %f \n", seconds);

	
	return 0;
}


