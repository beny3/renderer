/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "colition.h"
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


void draw(mesh3D *obj, int shift){
	vector *v = obj->vertics_trans;
	int *f = obj->face;
	int nb_f = obj->nb_f;

	for (int i=0; i<nb_f; i+=3){
		XDrawLine(dis, win, gc, v[f[i]].x + shift, v[f[i]].y + shift, v[f[i+1]].x + shift, v[f[i+1]].y + shift);
		XDrawLine(dis, win, gc, v[f[i+1]].x + shift, v[f[i+1]].y + shift, v[f[i+2]].x + shift, v[f[i+2]].y + shift);
		XDrawLine(dis, win, gc, v[f[i+2]].x + shift, v[f[i+2]].y + shift, v[f[i]].x + shift, v[f[i]].y + shift);
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
	
	mesh3D rectA;// = make_rect(-10, -10, -10, 40, 40, 40);
	//rectA.vertics[3].x +=20;
	read_mesh3D(&rectA, "cube3.obj");

	physic phyA;
	phyA.position = vec(0,0,0);
	phyA.velocity = vec(0,0,0);
	phyA.moment = vec(0,0,0);
	physic phyB;
	phyB.position = vec(100,70,0);
	phyB.velocity = vec(-1, 0, 0);
	phyB.moment = vec(0, 1, 0);
	phyB.angle = 0;


	mesh3D rectB;// = make_rect(-10, -10, -10, 50, 50, 50);
	read_mesh3D(&rectB, "cube3.obj");

	colition_hull hullB = make_graph(rectB);
	colition_hull hullA = make_graph(rectA);
	
	print_hull(&hullB);

	vector D = vec(0,1,0);

	vector AB[64];
	vector AB_t[64];
	vector normals[3]; 

	Simplex simplex;

	mesh3D tetra;
	tetra.vertics = simplex.data;
	tetra.vertics_trans = new vector[4];
	tetra.nb = 4;
	tetra.nb_f = 12;
	tetra.face = new int[12];
	int *f = tetra.face;
	f[0]=0; f[1]=1; f[2]=2;
	f[3]=0; f[4]=1; f[5]=3;
	f[6]=0; f[7]=2; f[8]=3;
	f[9]=1; f[10]=2; f[11]=3;

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
	
	while(1) {		
	
		XNextEvent(dis, &event);

		if (event.type == KeyPress){
			if (event.xkey.keycode == 10){
				rho+=PI/256.;
				m = make_rot4(0.,rho,0.);
			};
			if (event.xkey.keycode == 11){
				rho-=PI/256.;
				m = make_rot4(0.,rho,0.);
			};
			if (event.xkey.keycode == 86){
				phi+=PI/256.;
				m2 = make_rot4(0.,0., phi);
			};
			if (event.xkey.keycode == 82){
				phi-=PI/256.;
				m2 = make_rot4(0.,0., phi);
			};

		} 
	
		if (event.type == KeyPress || event.type==ButtonPress) {
			redraw();
			
			//print_mat(m2.data);
			
			vector pos = vec(event.xbutton.x-200, event.xbutton.y-200, 0);
			
			euler(&phyB, 1);
			euler(&phyA, 1);

			make_mov4(m.data, phyB.position);
			vector pos2 = vec(400,   0, 0);
			vector pos3 = vec(  0, 400, 0);
			make_mov4(m2.data, pos2);
			make_mov4(m22.data, pos3);

			matrix m3 = mult_m(m2.data, phyB.mat.data);
			matrix m4 = mult_m(m22.data, phyB.mat.data);

			mult_vm(phyB.mat.data, rectB.vertics, rectB.vertics_trans, rectB.nb);
			mult_vm(phyA.mat.data, rectA.vertics, rectA.vertics_trans, rectA.nb);

			XSetForeground(dis,gc, 0xFF0000);
			draw(&rectB, 200);
			draw(&rectA, 200);
			drawNormal(&rectA, 200);

			clock_t start = clock();
			vector c = vec(0, 0, 0);
			//c = colid_fast(&phyA, &phyB, &hullA, &hullB, &c);
			c = colid_fast(&phyB, &phyA, &hullB, &hullA, &c);
			clock_t end = clock();

			double seconds = (double)(end - start) / CLOCKS_PER_SEC;
			printf("temps %f \n", seconds);
			
			
			XFillRectangle(dis, win, gc, c.x +200, c.y+200, 10, 10);
			//XFillRectangle(dis, win, gc, ab.b.x +200, ab.b.y+200, 10, 10);
		
			mult_vm(m3.data, rectB.vertics, rectB.vertics_trans, rectB.nb);
			mult_vm(m2.data, rectA.vertics, rectA.vertics_trans, rectA.nb);
			draw(&rectB, 200);
			draw(&rectA, 200);
	


			mult_vm(m4.data, rectB.vertics, rectB.vertics_trans, rectB.nb);
			mult_vm(m22.data, rectA.vertics, rectA.vertics_trans, rectA.nb);

			vector c1 = mult_vm(m22.data, c);
			draw(&rectB, 200);
			draw(&rectA, 200);
			XFillRectangle(dis, win, gc, c1.x +200, c1.y+200, 10, 10);
			//XFillRectangle(dis, win, gc, b1.x +200, b1.y+200, 10, 10);
			//printf("colid_test %d \n",  c);
	
		}
	}

	
	
	return 0;
}


