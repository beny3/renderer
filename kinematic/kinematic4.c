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
		this->rest = new vector[n];
		for (int i=0; i<=n; i++){
			this->rest[i]=debut[i];
		}
		this->n = n;
		mats = new float[n*16];
	}
	float *mats;
	vector *debut;
	vector *rest;
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

void draw(vector *v, int *f, int nb_f){

	for (int i=0; i<nb_f; i+=3){
		XDrawLine(dis, win, gc, v[f[i]].x , v[f[i]].y, v[f[i+1]].x, v[f[i+1]].y);
		XDrawLine(dis, win, gc, v[f[i+1]].x , v[f[i+1]].y, v[f[i+2]].x, v[f[i+2]].y);
		XDrawLine(dis, win, gc, v[f[i+2]].x , v[f[i+2]].y, v[f[i]].x, v[f[i]].y);
	}
}

mesh3D make_rect(float x, float y, float z, float sx, float sy, float sz){
	mesh3D rect;
	rect.nb = 8;
	rect.nb_f = 12*3;
	rect.vertics = new vector[8];
	rect.vertics_trans = new vector[8];
	rect.face = new int[12*3];
	vector *v = rect.vertics;
	int *f = rect.face;
	v[0] = vec(x     , y     , z);
	v[1] = vec(x + sx, y     , z);
	v[2] = vec(x + sx, y + sy, z);
	v[3] = vec(x     , y + sy, z);

	v[4] = vec(x     , y     , z + sz);
	v[5] = vec(x + sx, y     , z + sz);
	v[6] = vec(x + sx, y + sy, z + sz);
	v[7] = vec(x     , y + sy, z + sz);

	f[0]=0;f[1]=1;f[2]=2;
	f[3]=0;f[4]=2;f[5]=3;

	f[6]=4;f[7]=5;f[8]=6;
	f[9]=4;f[10]=6;f[11]=7;

	f[12]=3;f[13]=2;f[14]=6;
	f[15]=3;f[16]=6;f[17]=7;

	f[18]=0;f[19]=1;f[20]=5;
	f[21]=0;f[22]=5;f[23]=4;

	f[24]=1;f[25]=5;f[26]=6;
	f[27]=1;f[28]=6;f[29]=2;

	f[30]=0;f[31]=4;f[32]=7;
	f[33]=0;f[34]=7;f[35]=3;

	return rect;
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


void inversKin(vector2D *bones, int n, vector2D &target){
	bones[n].x= target.x;
	bones[n].y= target.y;
		

	for (int converge=0; converge<5; converge++){
		bones[n].x= target.x;
		bones[n].y= target.y;
		
		for (int i=n; i>1; i--){
			vector2D bone;				
			bone.x = (bones[i-1].x + (converge==0? 1:0) ) - bones[i].x;
			bone.y =  bones[i-1].y  - bones[i].y;
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

	float phi=-0.1;
	float cosP = cos(phi);
	float sinP = sin(phi);

	vector t1 = minus( &ch.target, &ch.debut[0]);
	normalize(&t1);
	vector t2;
	if (t1.y!=0 || t1.y!= 0){
		t2 = vec(0, t1.z, -t1.y);
	}else{
		t2 = vec(0, 1, 0);
	}
	normalize(&t2);
	vector t3 = cross(t2, t1);
	vector normal = vec(t2.x*cosP - t3.x*sinP, t2.y*cosP - t3.y*sinP, t2.z*cosP - t3.z*sinP); 

	float mat[16]={t3.x*cosP + t2.x*sinP, t1.x, normal.x, ch.debut[0].x, 
				t3.y*cosP + t2.y*sinP, t1.y, normal.y, ch.debut[0].y, 
				t3.z*cosP + t2.z*sinP, t1.z, normal.z, ch.debut[0].z, 
				0   , 0   , 0   ,1};


	t1 = minus( &ch.debut[ch.n], &ch.debut[0]);
	normalize(&t1);
	t2;
	if (t1.y!=0 || t1.y!= 0){
		t2 = vec(0, t1.z, -t1.y);
	}else{
		t2 = vec(0, 1, 0);
	}
	normalize(&t2);
	t3 = cross(t2, t1);



	float mat2[16]={t3.x, t1.x, t2.x, ch.debut[0].x, 
				 t3.y, t1.y, t2.y, ch.debut[0].y, 
				 t3.z, t1.z, t2.z, ch.debut[0].z, 
				0   , 0   , 0   ,1};

	
	vector2D target_p = inv_vm_p(mat, ch.target);
	vector2D bones[N+1];

	for (int i=0; i<= ch.n; i++){
		bones[i] = inv_vm_p(mat2, ch.debut[i]);	
	}

	/*
	for (int i=1; i<=ch.n; i++){
				
		XSetForeground(dis,gc, 0x808080);
		XDrawLine(dis, win, gc, (int)bones[i].x + 100, (int)bones[i].y, (int)bones[i-1].x + 100, (int)bones[i-1].y);
	}
	*/
	inversKin(bones, ch.n, target_p);

	for (int i=1; i<= ch.n; i++){	
		ch.debut[i]=mult_vm(mat, bones[i]);
	}

	for (int i=0; i<ch.n; i++){

		// rest rotat matrix
		/*
		t1 = minus(&ch.rest[i+1], &ch.rest[i]);
		normalize(&t1);
		if (t1.y!=0 || t1.y!= 0){
			t2 = vec(0, t1.z, -t1.y);
		}else{
			t2 = vec(0, 1, 0);
		}

		t3 = cross(normal, t1);
		normalize(&t3);
		
		float restM[16]={t3.x, t3.y, t3.z, 0, 
				 	  t1.x, t1.y, t1.z, 0, 
					  t2.x, t2.y, t2.z, 0, 
				 	   0   , 0   , 0   ,1};

		vector rp = mult_vm(restM, ch.rest[i]);
		scalar(*rp, -1);
		make_mov4(restM, rp);
		*/


		t1 = minus(&ch.debut[i+1], &ch.debut[i]);
		normalize(&t1);
		t3 = cross(normal, t1);
		normalize(&t3);
	
		//make bone matrix m such as m*x = bone_rotation*(x - previous_bone) + previous_bone
		// = bone_rotation*x + previous_bone - bone_rotation*previous_bone
		float *rotat=&ch.mats[i*16];
		
		rotat[0]=t3.x; rotat[1] = t1.x; rotat[2]  = normal.x; rotat[3] = 0; 
		rotat[4]=t3.y; rotat[5] = t1.y; rotat[6]  = normal.y; rotat[7] = 0;
		rotat[8]=t3.z; rotat[9] = t1.z; rotat[10] = normal.z; rotat[11] = 0;
		rotat[12]=0; rotat[13] = 0; rotat[14] = 0; rotat[15] = 1;

		// rp = - bone_rotation*previous_bone
		vector rp = mult_vm(rotat, ch.rest[i]);
		// rp = previous_bone - bone_rotation*previous_bone
		rp.x =  ch.debut[i].x - rp.x;
		rp.y =  ch.debut[i].y - rp.y;
		rp.z =  ch.debut[i].z - rp.z;

		make_mov4(rotat, rp);

		//print_mat(rotat);

	}

	for (int i=1; i<=ch.n; i++){

		XSetForeground(dis,gc, 0xFF0000);
		XDrawLine(dis, win, gc, (int)ch.debut[i].x, (int)ch.debut[i].y, (int)ch.debut[i-1].x, (int)ch.debut[i-1].y);
		//XDrawLine(dis, win, gc, (int)ch.debut[i].z, (int)ch.debut[i].y, (int)ch.debut[i-1].z, (int)ch.debut[i-1].y);


		//XSetForeground(dis,gc, 0x0000FF);
		//XDrawLine(dis, win, gc, (int)bones[i].x + 100, (int)bones[i].y, (int)bones[i-1].x + 100, (int)bones[i-1].y);
	}


	/*
	XSetForeground(dis,gc, 0x880033);

	XDrawLine(dis, win, gc, 0, 200, 600, 200);		

	XDrawLine(dis, win, gc, (int)ch.target.x, (int)ch.target.y, ch.target.x + t1.x*10, ch.target.y + t1.y*10);
	XDrawLine(dis, win, gc, (int)ch.target.x, (int)ch.target.y, ch.target.x + t3.x*10, ch.target.y + t3.y*10);
	XDrawLine(dis, win, gc, (int)ch.target.x, (int)ch.target.y, ch.target.x + t2.x*10, ch.target.y + t2.y*10); 
	

	XDrawLine(dis, win, gc, (int)ch.target.z, (int)ch.target.y, ch.target.z + 50, ch.target.y);
	XDrawLine(dis, win, gc, (int)ch.target.z, (int)ch.target.y, ch.target.z, ch.target.y + 50);
	


	XSetForeground(dis,gc, 0xFF00FF);
	XDrawLine(dis, win, gc, target_p.x + 100, (int)target_p.y, target_p.x + 150, target_p.y);
	XDrawLine(dis, win, gc, target_p.x + 100, (int)target_p.y, target_p.x + 100, target_p.y + 50);
	*/
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
	
	mesh3D rect[n];

	for (int i=0; i<leg.n; i++){	
		rect[i] = make_rect(leg.debut[i].x, leg.debut[i].y, leg.debut[i].z, 10, 100, 10);		
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
						
				clock_t start = clock();

				inversKin_step1(leg);

				clock_t end = clock();

				double seconds = (double)(end - start) / CLOCKS_PER_SEC;
				printf("temps ik %f \n", seconds);


				for (int i=0; i<leg.n; i++){
					print_mat(&leg.mats[i*16]);
					mult_vm(&leg.mats[i*16], rect[i].vertics, rect[i].vertics_trans, rect[i].nb); 
					draw( rect[i].vertics_trans, rect[i].face, rect[i].nb_f);
					draw( rect[i].vertics, rect[i].face, rect[i].nb_f);

				}
				 
			}

			first = 1;					
		}
	}


	
	return 0;
}


