#ifndef VECTOR
#define VECTOR
#include<stdlib.h>
#include<math.h>
#define PI (3.14159265)

typedef struct vector{
	double x;
	double y;
	double z;
} vector;

typedef struct vector2D{
	double x;
	double y;
} vector2D;

typedef struct vector4D{
	double x;
	double y;
	double z;
	double phi;
} vector4D;

typedef struct weight{
	int bone;
	double w;
} weight;

typedef struct matrix{
	
	double data[16];

} matrix;



double dot(vector2D *a, vector2D *b){
	return a->x*b->x + a->y*b->y;
}

double dot(vector4D *a, vector4D *b){
	return a->x*b->x + a->y*b->y + a->z*b->z + a->phi*b->phi;
}

double dot(vector *a, vector *b){
	return a->x*b->x + a->y*b->y + a->z*b->z;
}

double acc(vector *a, vector *b){
	a->x+=b->x;
	a->y+=b->y;
	a->z+=b->z;
}

double acc(vector *a, vector *b, double alpha){
	a->x+=alpha*b->x;
	a->y+=alpha*b->y;
	a->z+=alpha*b->z;
}

double sub(vector *a, vector *b){
	a->x-=b->x;
	a->y-=b->y;
	a->z-=b->z;
}

vector2D vec2D(double x, double y){
	vector2D out;
	out.x=x;
	out.y=y;
	return out;
}


vector vec(double x, double y, double z){
	vector out;
	out.x=x;
	out.y=y;
	out.z=z;
	return out;
}

vector linear(vector *a, vector *b, double alpha, double beta){
	vector out;
	out.x = alpha*a->x+beta*b->x;
	out.y = alpha*a->y+beta*b->y;
	out.z = alpha*a->z+beta*b->z;
	return out;
}

void reflect(vector *a, vector *n){
	double p = dot(a, n);
	/*
	a->x=n->x;
	a->y=n->y;
	a->z=n->z;
	*/
	a->x-=2*p*n->x;
	a->y-=2*p*n->y;
	a->z-=2*p*n->z;

}

vector matrixMul(vector n[3], vector *b){
	vector out;
	out.x = b->x*n[0].x + b->y*n[1].x + b->z*n[2].x;
	out.y = b->x*n[0].y + b->y*n[1].y + b->z*n[2].y;
	out.z = b->x*n[0].z + b->y*n[1].z + b->z*n[2].z;
	return out;
}

vector matrixMul2D(vector2D *n1, vector2D *n2, vector2D *n3, vector *b){
	vector out;
	out.x = b->x*n1->x + b->y*n2->x + b->z*n3->x;
	out.y = b->x*n1->y + b->y*n2->y + b->z*n3->y;
	return out;
}

void scalar(vector2D *a, double alpha){
	a->x *= alpha;
	a->y *= alpha;
}

void scalar(vector *a, double alpha){

	a->x *= alpha;
	a->y *= alpha;
	a->z *= alpha;
}

void scalar(vector4D *a, double alpha){

	a->x *= alpha;
	a->y *= alpha;
	a->z *= alpha;
	a->phi*=alpha;
}

template <typename T>
double normalize(T *v){
	double norm = sqrt(dot(v,v));
	if (norm == 0) return 0;
	scalar(v, 1./norm);
	return norm;
}

vector add(vector *a, vector *b, double dt){
	vector out;
	out.x = a->x + b->x*dt;
	out.y = a->y + b->y*dt;
	out.z = a->z + b->z*dt;
	return out;
}


vector add(vector *a, vector *b){
	vector out;
	out.x = a->x + b->x;
	out.y = a->y + b->y;
	out.z = a->z + b->z;
	return out;
}

vector2D add(vector2D *a, vector2D *b){
	vector2D out;
	out.x = a->x + b->x;
	out.y = a->y + b->y;

	return out;
}

vector inv(vector *a){
	vector out;
	out.x = -a->x;
	out.y = -a->y;
	out.z = -a->z;
	return out;
}

vector minus(vector *a, vector *b){
	vector out;
	out.x = a->x - b->x;
	out.y = a->y - b->y;
	out.z = a->z - b->z;
	return out;
}

vector4D minus(vector4D *a, vector4D *b){
	vector4D out;
	out.x = a->x - b->x;
	out.y = a->y - b->y;
	out.z = a->z - b->z;
	out.phi = a->phi - b->phi;

	return out;
}

vector2D minus2D(vector2D *a, vector2D *b){
	vector2D out;
	out.x = a->x - b->x;
	out.y = a->y - b->y;
	return out;
}

vector cross(vector &a, vector &b){
	vector out;
	out.x = a.y*b.z - a.z*b.y;
	out.y = a.z*b.x - a.x*b.z;
	out.z = a.x*b.y - a.y*b.x;
	return out;
}

vector cross(vector &a, vector &b, double s){
	vector out;
	out.x = s*(a.y*b.z - a.z*b.y);
	out.y = s*(a.z*b.x - a.x*b.z);
	out.z = s*(a.x*b.y - a.y*b.x);
	return out;
}

void printVect(vector4D &v){
	printf("x:%f y:%f z:%f  phi:%f\n",v.x,v.y,v.z, v.phi);
}

void printVect(vector &v){
	printf("x:%f y:%f z:%f\n",v.x,v.y,v.z);
}

void printVect(vector2D &v){
	printf("x:%f y:%f \n",v.x,v.y);
}

void print_mat(double mat[16]){
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			printf("%2.2f |",mat[4*i + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void make_ide4(double mat[16]){

	for(int i=0; i<16; i++){
		if(i%4 != i/4){
			mat[i]=0;
		}else{
			mat[i]=1;
		}
	}	
}

vector2D inv_vm_p(double mat[16], vector &x){
	vector xx;
	vector2D y;
	xx.x = x.x - mat[3 ];
	xx.y = x.y - mat[7 ];
	xx.z = x.z - mat[11];
	double *yy=(double*)&y;

	for (int i=0; i<2; i++){
		yy[i] = mat[0 + i]*xx.x + mat[4 + i]*xx.y + mat[8 + i]*xx.z;
	}
	return y;
}

vector inv_vm(double mat[16], vector &x){
	vector xx;
	vector y;
	xx.x = x.x - mat[3 ];
	xx.y = x.y - mat[7 ];
	xx.z = x.z - mat[11];
	double *yy=(double*)&y;

	for (int i=0; i<3; i++){
		yy[i] = mat[0 + i]*xx.x + mat[4 + i]*xx.y + mat[8 + i]*xx.z;
	}
	return y;
}

void mult_sm(double mat[16], double alpha){
	for (int i=0; i<16; i++){
		mat[i]*=alpha;
	}
}

vector mult_vm(double mat[16], vector2D &x){
	vector y; 
	double *yy=(double*)&y;
	for (int i=0; i<12; i+=4){
		yy[i/4] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[3 + i];
	}
	return y;
}

vector mult_vm(double mat[16], vector &x){
	vector y; 
	double *yy=(double*)&y;
	for (int i=0; i<12; i+=4){
		yy[i/4] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[2 + i]*x.z + mat[3 + i];
	}
	return y;
}

void mult_vm(double mat[16], vector *in, int n){
	for (int k=0; k<n; k++){
		vector x=in[k]; 
		double *y=(double*)&in[k];
		for (int i=0; i<12; i+=4){
			y[i/4] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[2 + i]*x.z + mat[3 + i];
		}
	}
}

void mult_vm(double mat[16], vector *in, vector *out, int n){
	for (int k=0; k<n; k++){
		double *y=(double*)&out[k];
		for (int i=0; i<12; i+=4){
			y[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
		}
	}
}

void mult_vm(double *mat_list, vector *in, vector *out, double *W, int n, int nb_bone){
	vector tpv;
	double* tp = (double*)&tpv; 
	double w, *mat, *y;

	for (int k=0; k<n; k++){
		y=(double*)&out[k];
		mat = &mat_list[0];

		for (int i=0; i<12; i+=4){
			y[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
			y[i/4]*= W[k]; 
		}
	}
	
	for (int w_iter =1; w_iter < nb_bone; w_iter++){ 
		mat = &mat_list[w_iter*16];

		for (int k=0; k<n; k++){
			w = W[w_iter*n + k];
			if ( w > 0.001){
				for (int i=0; i<12; i+=4){
					tp[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
					tp[i/4]*= w;
				}
				acc(&out[k], &tpv);
			}
		}
	}
}

void mult_vm(double *mat_list, vector *in, vector *out, weight *W, int n){
	vector tpv;
	double* tp = (double*)&tpv; 
	////mmmmmm

	for (int k=0; k<n; k++){
		double *y=(double*)&out[k];
		double *mat = &mat_list[W[k].bone*16];

		for (int i=0; i<12; i+=4){
			y[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
			y[i/4]*= W[k].w; 
		}

		//if(y[0]==0 && y[1]==0 && y[2]==0)printf("k %d \n", k);
	}
	
	for (int w=1; w<3; w++){ 
		for (int k=0; k<n; k++){
			int idw = w*n + k;
			if (W[idw].w > 0){
				double *mat = &mat_list[W[idw].bone*16];
				for (int i=0; i<12; i+=4){
					tp[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
					tp[i/4]*= W[idw].w;
				}
				acc(&out[k], &tpv);
			}
		}
	}
}

void make_mov4(double mat[16], vector &dir){
	mat[3]=dir.x;
	mat[7]=dir.y;
	mat[11]=dir.z;
}

void make_rot4_u(double m[16], vector *u, double a){

	m[0] = cos(a) + u->x*u->x*(1-cos(a));
	m[1] = u->y*u->x*(1 - cos(a)) + u->z*sin(a);
	m[2] = u->z*u->x*(1 - cos(a)) - u->y*sin(a);
	m[3] = 0;
		
	m[4] = u->x*u->y*(1 - cos(a)) - u->z*sin(a);
	m[5] = cos(a) + u->y*u->y*(1 - cos(a));
	m[6] = u->z*u->y*(1 - cos(a)) + u->x*sin(a);
	m[7] = 0;
		
	m[8]  = u->x*u->z*(1 - cos(a)) + u->y*sin(a);
	m[9]  = u->y*u->z*(1 - cos(a)) - u->x*sin(a);
	m[10] = cos(a) + u->z*u->z*(1-cos(a));
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;
}

matrix make_rot4( double a, double b, double c){
	matrix mat;
	double *m = mat.data;
	m[0] =  cos(a)*cos(c) + sin(a)*sin(b)*sin(c);
	m[1] =  -(sin(a)*cos(c) - cos(a)*sin(b)*sin(c));
	m[2] = -cos(b)*sin(c);
	m[3] = 0;
		
	m[4]  =  sin(a)*cos(b);
	m[5]  =  -(-cos(a)*cos(b));
	m[6]  =  sin(b);
	m[7] = 0;
		
	m[8]  =  cos(a)*sin(c) - sin(a)*sin(b)*cos(c);
	m[9]  =  -(sin(a)*sin(c) + cos(a)*sin(b)*cos(c));
	m[10]  =  cos(b)*cos(c);
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;

	return mat;
}

matrix transpose(double a[16]){
	matrix c;	
	int i,j;

	for (int t=0; t<16; t++){
		i = t/4;
		j = t%4;
		c.data[4*i + j]=a[4*j + i];
	}
	return c;
}

void mult_m(double c[16], double a[16], double b[16]){	
	int i,j;

	for (int t=0; t<16; t++){
		i = t/4;
		j = t%4;
		c[4*i + j]=0;
		for (int k=0; k<4; k++){
			c[4*i + j]+=b[k*4 + j]*a[i*4 + k];
			
		}	
	}
}

matrix mult_m(double a[16], double b[16]){
	matrix c;	
	int i,j;

	for (int t=0; t<16; t++){
		i = t/4;
		j = t%4;
		c.data[4*i + j]=0;
		for (int k=0; k<4; k++){
			c.data[4*i + j]+=b[k*4 + j]*a[i*4 + k];
			
		}	
	}
	return c;
}

#endif
