/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "colition.h"

#define H 800
#define W 800
#define COLOR(r, g, b) ((b) | ((g)<<8) | ((r)<<16))
#define MAP_S H*W*3


/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;


void draw(vector *v, int *f, int nb_f, int shift){

	for (int i=0; i<nb_f; i+=3){

		XDrawLine(dis, win, gc, v[f[i]].x + shift, v[f[i]].y + shift, v[f[i+1]].x + shift, v[f[i+1]].y + shift);
		XDrawLine(dis, win, gc, v[f[i+1]].x + shift, v[f[i+1]].y + shift, v[f[i+2]].x + shift, v[f[i+2]].y + shift);
		XDrawLine(dis, win, gc, v[f[i+2]].x + shift, v[f[i+2]].y + shift, v[f[i]].x + shift, v[f[i]].y + shift);
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
	
	mesh3D rectA = make_rect(-10, -10, -10, 40, 40, 40);
	//rectA.vertics[3].x +=20;
	
	mesh3D rectB = make_rect(-10, -10, -10, 50, 50, 50);

	colition_hull hullB = make_graph(rectB);
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
	matrix m2 = make_rot4(0,0,PI/2);
	matrix m22 = make_rot4(0,PI/2,0);
	matrix m = make_rot4(1,2,3);
	float phi = 0;
	float rho = PI/2;

	while(1) {		
	
		XNextEvent(dis, &event);

		if (event.type == KeyPress){
			if (event.xkey.keycode == 10){
				rho+=PI/256.;
				m = make_rot4(0,rho,0);
			};
			if (event.xkey.keycode == 11){
				rho-=PI/256.;
				m = make_rot4(0,rho,0);
			};
			if (event.xkey.keycode == 86){
				phi+=PI/256.;
				m2 = make_rot4(0,0, phi);
			};
			if (event.xkey.keycode == 82){
				phi-=PI/256.;
				m2 = make_rot4(0,0, phi);
			};

		} 
	
		if (event.type == KeyPress || event.type==ButtonPress) {
			redraw();
			
			//print_mat(m2.data);

			vector pos = vec(event.xbutton.x-200, event.xbutton.y-200, 0);
			make_mov4(m.data, pos);
			vector pos2 = vec(400,   0, 0);
			vector pos3 = vec(  0, 400, 0);
			make_mov4(m2.data, pos2);
			make_mov4(m22.data, pos3);

			matrix m3 = mult_m(m2.data, m.data);
			matrix m4 = mult_m(m22.data, m.data);

			mult_vm(m.data, rectB.vertics, rectB.vertics_trans, rectB.nb);

			XSetForeground(dis,gc, 0xFF0000);
			draw(rectB.vertics_trans, rectB.face, rectB.nb_f, 200);
			draw(rectA.vertics, rectA.face, rectA.nb_f, 200);
			
			int c = gjk(rectB.vertics_trans, rectB.nb, rectA.vertics , rectA.nb, &simplex, AB);
		
			mult_vm(m3.data, rectB.vertics, rectB.vertics_trans, rectB.nb);
			mult_vm(m2.data, rectA.vertics, rectA.vertics_trans, rectA.nb);
			draw(rectB.vertics_trans, rectB.face, rectB.nb_f, 200);
			draw(rectA.vertics_trans, rectA.face, rectA.nb_f, 200);


			mult_vm(m4.data, rectB.vertics, rectB.vertics_trans, rectB.nb);
			mult_vm(m22.data, rectA.vertics, rectA.vertics_trans, rectA.nb);
			draw(rectB.vertics_trans, rectB.face, rectB.nb_f, 200);
			draw(rectA.vertics_trans, rectA.face, rectA.nb_f, 200);

			XSetForeground(dis,gc, 0x888888);
			draw(AB, 64, 200);
			mult_vm(m2.data,AB, AB_t, 64);
			draw(AB_t, 64, 200);
			mult_vm(m22.data,AB, AB_t, 64);
			draw(AB_t, 64, 200);


			XFillRectangle(dis, win, gc, 200,  200, 3, 3);
			XFillRectangle(dis, win, gc, 600,  200, 3, 3);

			if (c>0){
				XSetForeground(dis,gc, 0x0000FF);
				mult_vm(m2.data, tetra.vertics, tetra.vertics_trans, tetra.nb);
				draw(tetra.vertics, tetra.face, tetra.nb_f, 200);
				draw(tetra.vertics_trans, tetra.face, tetra.nb_f, 200);
				XSetForeground(dis,gc, 0x4444FF);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[0].x +200,  tetra.vertics_trans[0].y+200, 10, 10);
				XFillRectangle(dis, win, gc, tetra.vertics[0].x +200,  tetra.vertics[0].y+200, 10, 10);
				XSetForeground(dis,gc, 0x44FF44);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[1].x +200,  tetra.vertics_trans[1].y+200, 10, 10);
				XFillRectangle(dis, win, gc, tetra.vertics[1].x +200,  tetra.vertics[1].y+200, 10, 10);
				XSetForeground(dis,gc, 0xFF4444);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[2].x +200,  tetra.vertics_trans[2].y+200, 10, 10);
				XFillRectangle(dis, win, gc, tetra.vertics[2].x +200,  tetra.vertics[2].y+200, 10, 10);
				XSetForeground(dis,gc, 0x444444);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[3].x +200,  tetra.vertics_trans[3].y+200, 10, 10);
				XFillRectangle(dis, win, gc, tetra.vertics[3].x +200,  tetra.vertics[3].y+200, 10, 10);
			
				vector centre = add(&tetra.vertics[0], &tetra.vertics[1]);
				acc(&centre, &tetra.vertics[3]);
				scalar(&centre, 0.33);
			
				XSetForeground(dis,gc, 0xF200FF);
				XDrawLine(dis, win, gc, centre.x + 200, 
								    centre.y + 200, 
								    centre.x + normals[0].x + 200, 
								    centre.y + normals[1].y + 200);

				centre = mult_vm(m2.data, centre);
				normals[0] = mult_vm(m2.data, normals[0]);

				XDrawLine(dis, win, gc, centre.x + 200, 
								    centre.y + 200, 
								    centre.x + normals[0].z + 200, 
								    centre.y + normals[0].y + 200);

				mult_vm(m22.data, tetra.vertics, tetra.vertics_trans, tetra.nb);

				draw(tetra.vertics_trans, tetra.face, tetra.nb_f, 200);
				XSetForeground(dis,gc, 0x4444FF);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[0].x +200,  tetra.vertics_trans[0].y+200, 10, 10);
				XSetForeground(dis,gc, 0x44FF44);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[1].x +200,  tetra.vertics_trans[1].y+200, 10, 10);
				XSetForeground(dis,gc, 0xFF4444);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[2].x +200,  tetra.vertics_trans[2].y+200, 10, 10);
				XSetForeground(dis,gc, 0x444444);
				XFillRectangle(dis, win, gc, tetra.vertics_trans[3].x +200,  tetra.vertics_trans[3].y+200, 10, 10);
			}

			printf("colid_test %d \n",  c);
	
		}
	}


	
	return 0;
}


