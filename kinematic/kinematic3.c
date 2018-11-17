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
#define PI (3.14159265)
#define NB_ENTITIES 4



typedef struct transform_data{
	int idx=0;
	vector positions[NB_ENTITIES];
	vector rotations[NB_ENTITIES];
	matrix mats[NB_ENTITIES];
	

	void resetEntity( ){
		idx = 0;
	}

	int make_entity(float x, float y, float z, float rx, float ry, float rz){
		positions[idx]=vec(x, y, z);
		rotations[idx]=vec(rx, ry, rz);
		this->idx++%NB_ENTITIES;
		return idx-1;
	}

	void compute_mat( ){
		for (int i=0; i<idx+1; i++){
			vector *c = &rotations[i];
			mats[i] = make_rot4(c->x, c->y, c->z);
		}
		for (int i=0; i<idx+1; i++){
			make_mov4(mats[i].data, positions[i]);
		}
	}
	
} tranform_data;

typedef struct chaine{
	chaine(vector *debut, int n){
		this->debut = debut;
		this->n = n;
	}
	vector *debut;
	int n;
	vector target;
} chaine;

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


void inversKin(vector2D *bones, vector2D *norm, int n, vector2D &target){
	bones[n].x= target.x;
	bones[n].y= target.y;
	vector2D total = minus2D(&bones[n], &bones[0]);
	float length = dot2D(&total, &total);		

	for (int i=1; i<=n; i++){
		norm[i].x =  (bones[i].y - bones[i-1].y);
		norm[i].y =  (bones[i-1].x - bones[i].x);
		normalize2D(&norm[i]);
		scalar2D(&norm[i], length/128);
	}
				
	for (int converge=0; converge<10; converge++){
		bones[n].x= target.x;
		bones[n].y= target.y;
		
		for (int i=n; i>1; i--){
			vector2D bone;				
			bone.x = (bones[i-1].x + (converge==0? norm[i-1].x:0) ) - bones[i].x;
			bone.y = (bones[i-1].y + (converge==0? norm[i-1].y:0) ) - bones[i].y;
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

void inversKin_step1(chaine &ch){
	vector t1 = minus( &ch.target, &ch.debut[0]);
	//t1.x = 0;
	normalize(&t1);
	vector t2;
	if (t1.y!=0 || t1.y!= 0){
		t2 = vec(0, t1.z, -t1.y);
	}else{
		t2 = vec(0, 1, 0);
	}
	normalize(&t2);
	vector t3 = cross(t2, t1);

	float mat[16]={t3.x, t1.x, t2.x,0, 
				t3.y, t1.y, t2.y,0, 
				t3.z, t1.z, t2.z,0, 
				0   , 0   , 0   ,1};

	
	//printf("orto? %f %f %f \n", dot(&t1, &t2), dot(&t2, &t3), dot(&t1, &t3) );
	 print_mat(mat);
	//print_mat(mult_m(mat, transpose(mat).data).data);
	vector tp = minus(&ch.target, &ch.debut[0]);
	///vector target_p = inv_vm(mat, tp);
	vector2D target_p;

	float normy = 1./sqrt(mat[5]*mat[5] + mat[9]*mat[9]);
	float normz = 1./sqrt(mat[6]*mat[6] + mat[10]*mat[10]);

	target_p.y = tp.y*mat[5]*normy + tp.z*mat[9]*normy + ch.debut[0].y;
	target_p.x = ch.target.x;
	
	float phi=0.5;
	float cosP = cos(phi);
	float sinP = sin(phi);
	

	vector2D bones[N+1];
	vector2D norm[N+1];

	for (int i=0; i<= ch.n; i++){

		//bones[i].x = (ch.debut[i].x-ch.debut[0].x)*cosP - (ch.debut[i].y-ch.debut[0].y)*mat[6]*normz*sinP - (ch.debut[i].z-ch.debut[0].z)*mat[10]*normz*sinP + ch.debut[0].x;
		bones[i].x = ch.debut[i].x;
		bones[i].y = (ch.debut[i].y-ch.debut[0].y)*mat[5]*normy + (ch.debut[i].z-ch.debut[0].z)*mat[9]*normy + ch.debut[0].y;	
	}

	inversKin(bones, norm, ch.n, target_p);

	for (int i=1; i<= ch.n; i++){
		
		ch.debut[i].x =  bones[i].x;
		ch.debut[i].y = (bones[i].y - bones[0].y)*mat[5]*normy  +  ch.debut[0].y;
		ch.debut[i].z = (bones[i].y - bones[0].y)*mat[9]*normy  +  ch.debut[0].z;
		
		/*
		ch.debut[i].x = (bones[i].x - bones[0].x)*cosP   +  ch.debut[0].x;
		ch.debut[i].y = (bones[i].x - bones[0].x)*mat[6]*normz*sinP  + (bones[i].y - bones[0].y)*mat[5]*normy  +  ch.debut[0].y;
		ch.debut[i].z = (bones[i].x - bones[0].x)*mat[10]*normz*sinP + (bones[i].y - bones[0].y)*mat[9]*normy  +  ch.debut[0].z;
		*/
 
	}
	
	for (int i=1; i<=ch.n; i++){
				
		XSetForeground(dis,gc, 0x0000FF);
		//XDrawLine(dis, win, gc, bones[i].x, bones[i].y, bones[i].x + norm[i].x, bones[i].y + norm[i].y);

		XSetForeground(dis,gc, 0xFF0000);
		XDrawLine(dis, win, gc, (int)ch.debut[i].x, (int)ch.debut[i].y, (int)ch.debut[i-1].x, (int)ch.debut[i-1].y);

		XDrawLine(dis, win, gc, (int)ch.debut[i].z, (int)ch.debut[i].y, (int)ch.debut[i-1].z, (int)ch.debut[i-1].y);

		XSetForeground(dis,gc, 0x0000FF);
		XDrawLine(dis, win, gc, (int)bones[i].x, (int)bones[i].y, (int)bones[i-1].x, (int)bones[i-1].y);
	}

	XSetForeground(dis,gc, 0x880033);

	XDrawLine(dis, win, gc, 0, 200, 600, 200);		

	XDrawLine(dis, win, gc, (int)ch.target.x, (int)ch.target.y, ch.target.x + t1.x*10, ch.target.y + t1.y*10);
	XDrawLine(dis, win, gc, (int)ch.target.x, (int)ch.target.y, ch.target.x + t3.x*10, ch.target.y + t3.y*10);
	XDrawLine(dis, win, gc, (int)ch.target.x, (int)ch.target.y, ch.target.x + t2.x*10, ch.target.y + t2.y*10); 
	
	XDrawLine(dis, win, gc, (int)ch.target.z, (int)ch.target.y, ch.target.z + t1.z*10, ch.target.y + t1.y*10);
	XDrawLine(dis, win, gc, (int)ch.target.z, (int)ch.target.y, ch.target.z + t3.z*10, ch.target.y + t3.y*10);
	XDrawLine(dis, win, gc, (int)ch.target.z, (int)ch.target.y, ch.target.z + t2.z*10, ch.target.y + t2.y*10); 

	XSetForeground(dis,gc, 0x0000FF);
	XDrawLine(dis, win, gc, target_p.x, (int)target_p.y, target_p.x + t1.x*10, target_p.y + t1.y*10);
	XDrawLine(dis, win, gc, target_p.x, (int)target_p.y, target_p.x + t3.x*10, target_p.y + t3.y*10);
	XDrawLine(dis, win, gc, target_p.x, (int)target_p.y, target_p.x + t2.x*10, target_p.y + t2.y*10);
	
	t1=inv_vm(mat, t1);
	t2=inv_vm(mat, t2);
	t3=inv_vm(mat, t3);

	XDrawLine(dis, win, gc, 10, (int)target_p.y, 10 + t1.z*10, target_p.y + t1.y*10);
	XDrawLine(dis, win, gc, 10, (int)target_p.y, 10 + t3.z*10, target_p.y + t3.y*10);
	XDrawLine(dis, win, gc, 10, (int)target_p.y, 10 + t2.z*10, target_p.y + t2.y*10); 
	
}


void test_vect(){

	float mat2[16];
	matrix mat1 = make_rot4(0, PI/6, PI/6);
	float *mat = mat1.data; 
	matrix mat1_inv = transpose(mat);

	vector a = vec(1,2,3);
	vector b = vec(0,1,3);

	make_ide4(mat2);
	make_mov4(mat2, a);
	print_mat(mat2);
	print_mat(mat);
	print_mat(mat1_inv.data);
	print_mat(mult_m(mat, transpose(mat).data).data);

	matrix c = mult_m(mat2,mat);
	print_mat(c.data);
	vector r= mult_vm(mat2, a);
	vector r2 = mult_vm(mat, b);
	vector r3 = mult_vm(c.data, a);
	vector r4 = inv_vm(c.data, r3);
	c = transpose(mat2);
	print_mat(c.data);
	c = transpose(mat);
	print_mat(c.data);
	r2= mult_vm(mat, r2);
	printf("%f %f %f\n",r.x,r.y,r.z);
	printf("%f %f %f\n",r2.x,r2.y,r2.z);
	printf("%f %f %f\n",r3.x,r3.y,r3.z);
	printf("%f %f %f\n",r4.x,r4.y,r4.z);
}



int main () {
	XEvent event;		/* the XEvent declaration !!! */
	KeySym key;		/* a dealie-bob to handle KeyPress Events */	
	char text[255];
	
	init_x();
	int x=0;
	int y=0;
	int n = N;
	//init bones
	float bone3DA[] = {300,100,0 ,300,200,0 ,300,400,0};
	vector *bones3D = (vector*)bone3DA;
	chaine leg(bones3D, 2);
	leg.target = vec(100,100,100);

	float boneA[]={300,100,300,200,300,400};
	vector2D *bones = (vector2D*)boneA;
	vector2D norm[N+1];


	transform_data entities_transform;
	entities_transform.make_entity(100,0,0 ,0,PI/6, PI/6);
	entities_transform.compute_mat();
	print_mat(entities_transform.mats[0].data);


	int first = 0;
	test_vect();	

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
			redraw();
			if (first){
				 leg.target.x = event.xbutton.x;
				 leg.target.y = event.xbutton.y;

				 inversKin_step1(leg);
				 vector2D target = vec2D(event.xbutton.x, event.xbutton.y);
				 inversKin(bones, norm, n, target);
			}

			first = 1;
			
			/*
			for (int i=1; i<=n; i++){
				
				XSetForeground(dis,gc, 0x0000FF);
				XDrawLine(dis, win, gc, bones[i].x, bones[i].y, bones[i].x + norm[i].x, bones[i].y + norm[i].y);

				XSetForeground(dis,gc, 0xFF0000);
				XDrawLine(dis, win, gc, (int)bones[i].x, (int)bones[i].y, (int)bones[i-1].x, (int)bones[i-1].y);
			}
			*/				
		}
	}


	
	return 0;
}


