/* include the X library headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include "../3d.h"
#include "../bitmap.h"
#include "../obj3d.h"
#include "../anim.h"
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

typedef struct chain{
	int idx = 0;
	int n;
	vector target;
	float  phi;
}chain;


typedef struct skeleton{
	skeleton(vector *debut, int n){
		this->debut = debut;
		this->rest = new vector[n+1];
		this->lenght = new float[n+1];
		for (int i=0; i<=n; i++){
			this->rest[i]=debut[i];
			vector tp = minus(&debut[(i+1)%(n+1)], &debut[i]);
			this->lenght[i] = sqrt(dot(&tp,&tp));
		}
		this->lenght[n] = this->lenght[n-1];
		this->n = n;
		mats = new float[n*16];
	}
	float *lenght;
	float *mats;
	vector *debut;
	vector *rest;
	int n;
} skeleton;

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

void draw(mesh3D *obj, int *f, int nb_f, int bone){
	vector *v = obj->vertics_trans;
	weight *w = obj->w;

	for (int i=0; i<nb_f; i+=3){
		XSetForeground(dis,gc, 0);

		for (int k=0; k<3; k++){
			for(int j = 0; j < 3; j++){
				if (w[j*obj->nb + f[i + k]].bone == bone){
					XSetForeground(dis,gc, (int)(255*w[j*obj->nb + f[i + k]].w)<<16);
				}
			}		
			XDrawLine(dis, win, gc, v[f[i + k]].x , v[f[i + k]].y, v[f[i+(k+1)%3]].x, v[f[i+(k+1)%3]].y);
		}
	}
}

void draw(vector *v, int *f, int nb_f){

	for (int i=0; i<nb_f; i+=3){
		XDrawLine(dis, win, gc, v[f[i]].x , v[f[i]].y, v[f[i+1]].x, v[f[i+1]].y);
		XDrawLine(dis, win, gc, v[f[i+1]].x , v[f[i+1]].y, v[f[i+2]].x, v[f[i+2]].y);
		XDrawLine(dis, win, gc, v[f[i+2]].x , v[f[i+2]].y, v[f[i]].x, v[f[i]].y);
	}
}


float dispertion_l(float x){
	return 1/(1+exp(-(x)*16))+1/(1+exp((x-1)*16))-1;
}

float dispertion_r(float x){
	return 1/(1+exp((x-1)*16));
}

float dispertion_fin(float x){
	return 1/(1+exp(-x*16));
}


void print_weight(mesh3D &obj){
	for (int i=447; i< 510; i++){
		float sum = 0;
		for (int j=0; j<3; j++){
			sum+=obj.w[j*obj.nb + i].w;
			printf("w %f %d | ", obj.w[j*obj.nb + i].w, obj.w[j*obj.nb + i].bone);
		}
		printf("sum %f  \n", sum);
	}
	printf("_________________________________  \n");
}


void init_weight(mesh3D &obj, int nb_bone){
	
	obj.wf = new float[obj.nb*nb_bone];

	for (int i=0; i< obj.nb*nb_bone; i++){
			obj.wf[i] = 0;
	}
}

void compute_weight(mesh3D &obj, vector *bones, int start_bone, int nb_bones){
	vector *v = obj.vertics;
	float d_perpendicular, d_projected;	

	for (int b=start_bone; b < nb_bones-1; b++){//
		vector dir_bone = minus(&bones[b+1], &bones[b]);
		float  norm_dir_bone = normalize(&dir_bone);
		vector bx;
		vector bxXbone;

		for (int x=0; x< obj.nb; x++){
			bx = minus(&v[x], &bones[b]);
			d_projected = dot(&bx, &dir_bone)/norm_dir_bone;
			/*
			if (b == start_bone && nb_bones - start_bone > 2){
				d_projected = dispertion_r(dot(&bx, &dir_bone)/norm_dir_bone);
			}else 
			*/
			if (b == nb_bones-2){
				d_projected = dispertion_fin(dot(&bx, &dir_bone)/norm_dir_bone);
			}else{
				d_projected = dispertion_l(dot(&bx, &dir_bone)/norm_dir_bone);
			}

			bxXbone = cross(bx, dir_bone);
			//d_perpendicular = 1/(1 + sqrt(dot(&bxXbone, &bxXbone))/norm_dir_bone);
			d_perpendicular = dispertion_r(sqrt(dot(&bxXbone, &bxXbone))/norm_dir_bone*4);

			weight w;
			obj.wf[b*obj.nb + x] = d_projected*d_perpendicular;
		}
	//print_weight(obj);	
	}
}
void normalize_weight( mesh3D &obj, int nb_bone){
	for (int x=0; x< obj.nb; x++){
		float sum = 0;

		for (int b=0; b<nb_bone; b++){
			sum += obj.wf[b*obj.nb + x];
		}

		for (int b=0; b<nb_bone; b++){
			obj.wf[b*obj.nb + x]/=sum;
		}			
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


void inversKin2D(vector2D *bones, float *length, int n, vector2D &target){
	bones[n].x= target.x;
	bones[n].y= target.y;
		

	for (int converge=0; converge<5; converge++){
		bones[n].x= target.x;
		bones[n].y= target.y;
		
		for (int i=n; i>1; i--){
			vector2D bone;				
			bone.x = (bones[i-1].x + (converge==0? 10:0)) - bones[i].x;
			bone.y =  bones[i-1].y  - bones[i].y;
			normalize(&bone);
			scalar(&bone, length[i-1]);
			bones[i-1]=add(&bones[i], &bone);
		}

		for (int i=0; i<n; i++){				
			vector2D bone;
			bone.x = bones[i+1].x - bones[i].x;
			bone.y = bones[i+1].y - bones[i].y;
			normalize(&bone);
			scalar(&bone, length[i]);
			bones[i+1]=add(&bones[i], &bone);
		}			
	}	
}

//pas sur d'utiliser
matrix makeOrtho(vector *begin, vector *end, float cosP, float sinP){
	vector t1 = minus( end, begin);
	normalize(&t1);
	vector t3;
	if (t1.y==0 && t1.x== 0){
		t3 = vec(1, 0, 0);
	}else{
		t3 = vec( t1.y, -t1.x, 0);
	}
	normalize(&t3);
	normalize(&t3);
	vector t2 = cross(t3, t1);
	vector normal = vec(t2.x*cosP - t3.x*sinP, t2.y*cosP - t3.y*sinP, t2.z*cosP - t3.z*sinP); 

	float mat[16]={t3.x*cosP + t2.x*sinP, t1.x, normal.x, begin->x, 
				t3.y*cosP + t2.y*sinP, t1.y, normal.y, begin->y, 
				t3.z*cosP + t2.z*sinP, t1.z, normal.z, begin->z, 
				0   , 0   , 0   ,1};
	
	matrix *out = (matrix*)mat;

	return *out;
} 


void inversKin(chain &ch, skeleton &sk, float phi){

	vector *bones3D = &sk.debut[ch.idx];
	vector *rest = &sk.rest[ch.idx];

	vector t1 = minus( &ch.target, &bones3D[0]);
	normalize(&t1);
	vector t3;

	if (t1.y==0 && t1.x== 0){
		t3 = vec(1, 0, 0);
	}else{
		t3 = vec( t1.y, -t1.x, 0);
	}
	normalize(&t3);

	float cosP = cos(phi);
	float sinP = sin(phi);

	vector t2 = cross(t1, t3);
	vector normal = vec(t2.x*cosP - t3.x*sinP, t2.y*cosP - t3.y*sinP, t2.z*cosP - t3.z*sinP); 

	float mat[16]={t3.x*cosP + t2.x*sinP, t1.x, normal.x, bones3D[0].x, 
				t3.y*cosP + t2.y*sinP, t1.y, normal.y, bones3D[0].y, 
				t3.z*cosP + t2.z*sinP, t1.z, normal.z, bones3D[0].z, 
				0   , 0   , 0   ,1};

	//printf("normal %f %f %f \n", normal.x, normal.y, normal.z);
	//print_mat(mat);

	t1 = minus( &bones3D[ch.n], &bones3D[0]);
	normalize(&t1);
	if (t1.y==0 && t1.x== 0){
		t3 = vec(1, 0, 0);
	}else{
		t3 = vec( t1.y, -t1.x, 0);
	}
	normalize(&t3);
	t2 = cross(t1, t3);

	float mat2[16]={t3.x*cosP + t2.x*sinP, t1.x, -t3.x*sinP + t2.x*cosP, bones3D[0].x, 
				 t3.y*cosP + t2.y*sinP, t1.y, -t3.y*sinP + t2.y*cosP, bones3D[0].y, 
				 t3.z*cosP + t2.z*sinP, t1.z, -t3.z*sinP + t2.z*cosP, bones3D[0].z, 
				 0   , 0   , 0   ,1};

	vector2D target_p = inv_vm_p(mat, ch.target);
	vector2D bones[ch.n+1];

	for (int i=0; i<= ch.n; i++){
		bones[i] = inv_vm_p(mat2, bones3D[i]);	
	}

	for (int i=1; i<=ch.n; i++){
				
		XSetForeground(dis,gc, 0x808080);
		XDrawLine(dis, win, gc, (int)bones[i].x + 100, (int)bones[i].y, (int)bones[i-1].x + 100, (int)bones[i-1].y);
	}
	
	inversKin2D(bones, &sk.lenght[ch.idx] ,ch.n , target_p);

	for (int i=1; i<= ch.n; i++){	
		bones3D[i]=mult_vm(mat, bones[i]);
	}

	for (int i=0; i<ch.n; i++){

		//make rest matrix "restM": transform bone from rest to neutral position (bone at orgine pointing y axis)
		t1 = minus(&rest[i+1], &rest[i]);
		normalize(&t1);
		if (t1.y==0 && t1.x== 0){
			t3 = vec(1, 0, 0);
		}else{
			t3 = vec( t1.y, -t1.x, 0);
		}

		t2 = cross(t1, t3);
		normalize(&t2);
		
		float restM[16]={t3.x, t3.y, t3.z, 0, 
				 	  t1.x, t1.y, t1.z, 0, 
					  t2.x, t2.y, t2.z, 0, 
				 	   0   , 0   , 0   ,1};

		vector rp = mult_vm(restM, rest[i]);
		scalar(&rp, -1);
		make_mov4(restM, rp);

		t1 = minus(&bones3D[i+1], &bones3D[i]);
		normalize(&t1);
		t3 = cross(normal, t1);
		normalize(&t3);
	
		//make current matrix transform "current": transform bone from neutral position (bone at orgine pointing y axis) to current position
		float *bone_trans=&sk.mats[(i+ch.idx)*16];
		float current[16]={ t3.x, t1.x, normal.x, 0, 
						t3.y, t1.y, normal.y, 0,
						t3.z, t1.z, normal.z, 0,
						0   , 0   ,        0, 1};

		make_mov4(current, bones3D[i]);

		//make bone transform matrix "bone_trans" = current*restM, transform  bone from rest postion to current position
		mult_m(bone_trans, current, restM);
	}

	for (int i=1; i<=ch.n; i++){

		XSetForeground(dis,gc, 0xFF0000);
		XDrawLine(dis, win, gc, (int)bones3D[i].x, (int)bones3D[i].y, (int)bones3D[i-1].x, (int)bones3D[i-1].y);
	}
		
	/*
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
	float bones_legs_f[] = {300,100,0 , 300,150,0 ,300,300,0, 300,350,0, 300,500,0};
	vector *bones_legs = (vector*)bones_legs_f;

	float bones_man_f[] = {171,200,0 , 158,315,0 ,142,454,0, 232,200,0 , 243,315,0 ,260,454,0, 200,200,0,  200,-23,0};
	vector *bones_man = (vector*)bones_man_f;

	skeleton leg(bones_legs, 2);
	chain leg_c;
	leg_c.n = leg.n;
	leg_c.target = vec(100,100,100);

	skeleton man_s(bones_man, 7);
	chain leg_right_c;
	leg_right_c.idx = 0;
	leg_right_c.n = 2;
	leg_right_c.target = bones_man[2];
	leg_right_c.target.z = 100;

	chain leg_left_c;
	leg_left_c.idx = 3;
	leg_left_c.n = 2;
	leg_left_c.target = bones_man[5];

	chain spine_c;
	spine_c.idx = 6;
	spine_c.n = 1;
	spine_c.target = bones_man[7];

	transform_data entities_transform;
	entities_transform.make_entity(100,0,0 ,0,PI/6, PI/6);
	entities_transform.compute_mat();
	//print_mat(entities_transform.mats[0].data);

	int first = 0;

	mesh3D rect[leg.n];
	mesh3D gros_cylindre = make_cylinder(leg.debut[0].x, leg.debut[0].y, leg.debut[0].z, 100*leg.n, 20, 8, 40);
	mesh3D man;
	read_mesh3D(&man, "homme.obj");
	float manM[16];
	make_ide4(manM);
	vector dir=vec(2,2,0);
	make_mov4(manM, dir);
	mult_sm(manM, 100);
	print_mat(manM);
	mult_vm(manM, man.vertics,  man.nb );

	//compute_weight(gros_cylindre, leg.debut, leg.n+1);

	init_weight(man, man_s.n);
	compute_weight(man, man_s.debut, leg_left_c.idx, leg_left_c.idx+3);
	compute_weight(man, man_s.debut, leg_right_c.idx, leg_right_c.idx+3);
	compute_weight(man, man_s.debut, spine_c.idx, spine_c.idx+2);
	normalize_weight(man, man_s.n);
	//print_weight(man);


	//print_weight(man);

	float frame_list[12] = {man_s.debut[2].x-20, man_s.debut[2].y,0,0,
					   man_s.debut[2].x-20, 2*man_s.debut[0].y-man_s.debut[2].y,0, PI,
 					   man_s.debut[2].x-20, man_s.debut[2].y,0,0};

	sequence seq = make_sequence((vector4D*)frame_list, 3, 100);
	make_tangents(&seq);

	for(int i=0; i<3; i++){
		printf("seq p%f %f %f \n", seq.frames[i].p.x, seq.frames[i].p.y, seq.frames[i].p.z);
		printf(" seq t%f %f %f \n", seq.frames[i].t.x, seq.frames[i].t.y, seq.frames[i].t.z);
	}
	
	for (int i=0; i<leg.n; i++){	
		//rect[i] = make_cylinder(leg.debut[i].x, leg.debut[i].y, leg.debut[i].z, 100, 10, 8, 8);
		//rect[i] = make_rect(leg.debut[i].x, leg.debut[i].y, leg.debut[i].z, 20, 100, 20);		
	}

	int show_bone = 0;
	float phi = 0;
	int iter = 0;
	vector4D *aim = (vector4D*)&leg_right_c.target;
	
	/* look for events forever... */
	double seconds = 0;

	while(1) {		
		/* get the next event and stuff it into our event variable.
		   Note:  only events we set the mask for are detected!
		*/
		
		XNextEvent(dis, &event);

		if (event.type == KeyPress){
			if (event.xkey.keycode >= 10 && event.xkey.keycode <= 17){
				show_bone = event.xkey.keycode-10;
			}
			if (event.xkey.keycode == 86){
				phi+=0.2;
			};
			if (event.xkey.keycode == 82){
				phi-=0.2;
			};

		} 
		
		if (event.type) {
			
			/*
			leg_left_c.target.x = event.xbutton.x;
			leg_left_c.target.y = event.xbutton.y;
			*/
			leg_right_c.target.x = event.xbutton.x;
			leg_right_c.target.y = event.xbutton.y;
			
			//interpolate(aim, &seq, iter);
						

			//inversKin(leg_c, leg);

			inversKin(leg_left_c, man_s, 0);
			inversKin(spine_c, man_s, 0);
			inversKin(leg_right_c, man_s, phi);

			//printf("temps ik %f \n", seconds);
			//printf("x %d  y %d \n", event.xbutton.x, event.xbutton.y);

			//mult_vm(leg.mats, gros_cylindre.vertics, gros_cylindre.vertics_trans,  gros_cylindre.nb);
			//mult_vm(leg.mats, gros_cylindre.vertics, gros_cylindre.vertics_trans, gros_cylindre.w, gros_cylindre.nb );
			clock_t start = clock();
			mult_vm(man_s.mats, man.vertics, man.vertics_trans, man.wf, man.nb, man_s.n );
			clock_t end = clock();

			seconds += (double)(end - start) / CLOCKS_PER_SEC;
			if(iter%100==0){
				printf("temps ik %f \n", seconds/100);
				seconds = 0;
			}

			//draw( gros_cylindre.vertics_trans,  gros_cylindre.face,  gros_cylindre.nb_f );
			redraw();
			draw( &man, man.face, man.nb_f, show_bone);
			XSync( dis, 0);
			usleep(2000);
				
			XSetForeground(dis,gc, 0x808080);
			XDrawRectangle(dis, win, gc, dir.x*100  + 20, dir.y*100, 10, 10);
			XDrawRectangle(dis, win, gc, dir.x*100  - 30, dir.y*100, 10, 10);

			for (int i=0; i<leg.n; i++){
				float *m_trans = &leg.mats[i*16];
				//mult_vm(m_trans, rect[i].vertics, rect[i].vertics_trans, rect[i].nb); 
				//m_trans[3]=0;
				//m_trans[7]=0;
				//m_trans[11]=0;
				//print_mat(mult_m(m_trans, transpose(m_trans).data).data);
				//print_mat(m_trans);
				//draw( rect[i].vertics_trans, rect[i].face, rect[i].nb_f);
				//draw( rect[i].vertics, rect[i].face, rect[i].nb_f);

			}			 
					
		}
	iter++;
	}


	
	return 0;
}


