#ifndef VECTOR
#define VECTOR
#include<stdlib.h>
#include<math.h>
#define PI (3.14159265)

typedef struct vector{
	union{
		struct{
			float x;
			float y;
			float z;
		};
		float el[3];
	};
} vector;

typedef struct vector2D{
	float x;
	float y;
} vector2D;

typedef struct vector4D{
	float x;
	float y;
	float z;
	float phi;
} vector4D;

typedef struct weight{
	int bone;
	float w;
} weight;

typedef struct matrix{
	
	float data[16];

} matrix;

typedef struct matrix3{
	
	float data[9];

} matrix3;

float dot(vector2D *a, vector2D *b){
	return a->x*b->x + a->y*b->y;
}

float dot(vector4D *a, vector4D *b){
	return a->x*b->x + a->y*b->y + a->z*b->z + a->phi*b->phi;
}

float dot(vector *a, vector *b){
	return a->x*b->x + a->y*b->y + a->z*b->z;
}

float quad_sol1(float a,float b,float c){
	return (-b + sqrt(b*b -4*a*c))/(2*a);
}

float quad_sol2(float a,float b,float c){
	return (-b - sqrt(b*b -4*a*c))/(2*a);
}

float acc(vector *a, vector *b){
	a->x+=b->x;
	a->y+=b->y;
	a->z+=b->z;
}

float acc(vector *a, vector *b, float alpha){
	a->x+=alpha*b->x;
	a->y+=alpha*b->y;
	a->z+=alpha*b->z;
}
float acc_one_to_many(vector *a, vector *one, float alpha, int n){
	for (int i=0; i<n; i++){
		a[i].x+= alpha*one->x;
		a[i].y+= alpha*one->y;
		a[i].z+= alpha*one->z;
	}
}

float acc(vector *a, vector *b, float alpha, int n){
	for (int i=0; i<n; i++){
		a[i].x+= alpha*b[i].x;
		a[i].y+= alpha*b[i].y;
		a[i].z+= alpha*b[i].z;
	}
}

float sub(vector *a, vector *b){
	a->x-=b->x;
	a->y-=b->y;
	a->z-=b->z;
}

float sub_overWrite_left(vector *a, vector *b){
	b->x = a->x - b->x;
	b->y = a->y - b->y;
	b->z = a->z - b->z;
}

vector2D vec2D(float x, float y){
	vector2D out;
	out.x=x;
	out.y=y;
	return out;
}


vector vec(float x, float y, float z){
	vector out;
	out.x=x;
	out.y=y;
	out.z=z;
	return out;
}

vector linear(vector *a, vector *b, float alpha, float beta){
	vector out;
	out.x = alpha*a->x+beta*b->x;
	out.y = alpha*a->y+beta*b->y;
	out.z = alpha*a->z+beta*b->z;
	return out;
}

void reflect(vector *a, vector *n){
	float p = dot(a, n);
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
	out.x = b->x*n[0].x + b->y*n[2].x + b->z*n[1].x;
	out.y = b->x*n[0].y + b->y*n[2].y + b->z*n[1].y;
	out.z = b->x*n[0].z + b->y*n[2].z + b->z*n[1].z;
	return out;
}

vector matrixMul2D(vector2D *n1, vector2D *n2, vector2D *n3, vector *b){
	vector out;
	out.x = b->x*n1->x + b->y*n2->x + b->z*n3->x;
	out.y = b->x*n1->y + b->y*n2->y + b->z*n3->y;
	return out;
}

void scalar(vector2D *a, float alpha){
	a->x *= alpha;
	a->y *= alpha;
}

void scalar(vector *a, float alpha){

	a->x *= alpha;
	a->y *= alpha;
	a->z *= alpha;
}

void scalar(vector4D *a, float alpha){

	a->x *= alpha;
	a->y *= alpha;
	a->z *= alpha;
	a->phi*=alpha;
}

template <typename T>
float normalize(T *v){
	float norm = sqrt(dot(v,v));
	if (norm == 0) return 0;
	scalar(v, 1./norm);
	return norm;
}

vector copy(vector *a){
	vector out;
	out.x = a->x;
	out.y = a->y;
	out.z = a->z;
	return out;
}

void swap( vector **a, vector **b){
	vector *tp = *a;
	*a = *b;
	*b = tp;
}

void copy(vector *pos, vector *vel, int n){
	for (int i=0; i<n; i++){
		vel[i] = pos[i];
	}
}

void set_zero(vector *vel, int n){
	for (int i=0; i<n; i++){
		vel[i] = vec(0,0,0);
	}
}


vector add(vector *a, vector *b, float dt){
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

vector cross(vector &a, vector &b, float s){
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

void printVect(vector *v,int n){
	for (int i= 0; i<n; i++){
		printVect(v[i]);
	}
}

void printVect(vector2D &v){
	printf("x:%f y:%f \n",v.x,v.y);
}

void print_mat(float mat[16]){
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			printf("%2.2f |",mat[4*i + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_mat9(float mat[9]){
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			printf("%2.2f |",mat[3*i + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void make_ide4(float mat[16]){

	for(int i=0; i<16; i++){
		if(i%4 != i/4){
			mat[i]=0;
		}else{
			mat[i]=1;
		}
	}
}

void make_scale(float mat[16], float v){

	for(int i=0; i<16; i++){
		if(i%4 != i/4){
			mat[i]=0;
		}else{
			mat[i]=v;
		}
	}
	mat[15] = 1;
}

vector2D solve_uv(vector *u, vector *v, vector *b){

	int coord0 = 0, coord1 = 1;

	float det = u->el[0] * v->el[1] - u->el[1] * v->el[0];
	float det2 = u->el[1] * v->el[2] - u->el[2] * v->el[1];

	if (det2*det2 > det*det){
		det = det2;
		coord0 = 1;
		coord1 = 2;
	}

	det2 = u->el[0] * v->el[2] - u->el[2] * v->el[0];
	if (det2*det2 > det*det){
 		det = det2;
		coord0 = 0;
		coord1 = 2;
	}

	det = 1./det;

	return vec2D( (b->el[coord0]*v->el[coord1] - b->el[coord1]*v->el[coord0])*det, (- b->el[coord0]*u->el[coord1] + b->el[coord1]*u->el[coord0])*det);
}


int plane_intersection(vector *a, vector *n, vector *o, vector *r, float *lambda_out){

	*lambda_out = 0;

	float epsilone = 0.001;

	float r_dot_n = dot(r,n);
	if (r_dot_n*r_dot_n < epsilone)return 0;

	float lambda = (dot(a,n) - dot(o,n))/r_dot_n;
	//printf("lambda %f \n", lambda);

	if (lambda > 0 + epsilone  && lambda < 1 - epsilone){
		*lambda_out = lambda;
		return 1;
	}

	return 0;
}


int triangle_intersection(vector *a, vector *b, vector *c, vector *n, vector *o, vector *r){

	float lambda = (dot(a,n) - dot(o,n))/dot(r,n);

	if (lambda >= 0 && lambda <= 1){
		vector p = vec(r->x*lambda, r->y*lambda, r->z*lambda);
		acc(&p, o);
		sub(&p, a);
		vector ab = minus(b,a);
		vector ac = minus(c,a);
		vector2D uv = solve_uv(&ab, &ac, &p);
		if (uv.x >= 0 && uv.x <= 1 && uv.y >= 0 && uv.y <= 1) return 1; 
	}

	return 0;
}

vector2D inv_vm_p(float mat[16], vector &x){
	vector xx;
	vector2D y;
	xx.x = x.x - mat[3 ];
	xx.y = x.y - mat[7 ];
	xx.z = x.z - mat[11];
	float *yy=(float*)&y;

	for (int i=0; i<2; i++){
		yy[i] = mat[0 + i]*xx.x + mat[4 + i]*xx.y + mat[8 + i]*xx.z;
	}
	return y;
}

vector inv_vm(float mat[16], vector &x){
	vector xx;
	vector y;
	xx.x = x.x - mat[3 ];
	xx.y = x.y - mat[7 ];
	xx.z = x.z - mat[11];
	float *yy=(float*)&y;

	for (int i=0; i<3; i++){
		yy[i] = mat[0 + i]*xx.x + mat[4 + i]*xx.y + mat[8 + i]*xx.z;
	}
	return y;
}

void mult_sm(float mat[16], float alpha){
	for (int i=0; i<16; i++){
		mat[i]*=alpha;
	}
}

vector mult_vm(float mat[16], vector2D &x){
	vector y; 
	float *yy=(float*)&y;
	for (int i=0; i<12; i+=4){
		yy[i/4] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[3 + i];
	}
	return y;
}


vector mult_vm9(float mat[9], vector &x){
	vector y; 
	float *yy=(float*)&y;
	for (int i=0; i<9; i+=3){
		yy[i/3] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[2 + i]*x.z;
	}
	return y;
}

vector mult_vm(float mat[16], vector &x){
	vector y; 
	float *yy=(float*)&y;
	for (int i=0; i<12; i+=4){
		yy[i/4] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[2 + i]*x.z + mat[3 + i];
	}
	return y;
}

void mult_vm(float mat[16], vector *in, int n){
	for (int k=0; k<n; k++){
		vector x=in[k]; 
		float *y=(float*)&in[k];
		for (int i=0; i<12; i+=4){
			y[i/4] = mat[0 + i]*x.x + mat[1 + i]*x.y + mat[2 + i]*x.z + mat[3 + i];
		}
	}
}

void mult_vm(float mat[16], vector *in, vector *out, int n){
	for (int k=0; k<n; k++){
		float *y=(float*)&out[k];
		for (int i=0; i<12; i+=4){
			y[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
		}
	}
}

void mult_vm(float *mat_list, vector *in, vector *out, float *W, int n, int nb_bone){
	vector tpv;
	float* tp = (float*)&tpv; 
	float w, *mat, *y;

	for (int k=0; k<n; k++){
		y=(float*)&out[k];
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

void mult_vm(float *mat_list, vector *in, vector *out, weight *W, int n){
	vector tpv;
	float* tp = (float*)&tpv; 
	////mmmmmm

	for (int k=0; k<n; k++){
		float *y=(float*)&out[k];
		float *mat = &mat_list[W[k].bone*16];

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
				float *mat = &mat_list[W[idw].bone*16];
				for (int i=0; i<12; i+=4){
					tp[i/4] = mat[0 + i]*in[k].x + mat[1 + i]*in[k].y + mat[2 + i]*in[k].z + mat[3 + i];
					tp[i/4]*= W[idw].w;
				}
				acc(&out[k], &tpv);
			}
		}
	}
}

void make_mov4(float mat[16], vector &dir){
	mat[3]=dir.x;
	mat[7]=dir.y;
	mat[11]=dir.z;
}

void make_rot4_u(float m[16], vector *u, float a){

	m[0] = ( cos(a) + u->x*u->x*(1-cos(a)) );
	m[1] = u->y*u->x*(1 - cos(a)) + u->z*sin(a);
	m[2] = u->z*u->x*(1 - cos(a)) - u->y*sin(a);
	m[3] = 0;
		
	m[4] = (  u->x*u->y*(1 - cos(a)) - u->z*sin(a) );
	m[5] =  cos(a) + u->y*u->y*(1 - cos(a));
	m[6] =  u->z*u->y*(1 - cos(a)) + u->x*sin(a);
	m[7] = 0;
		
	m[8]  = ( u->x*u->z*(1 - cos(a)) + u->y*sin(a) );
	m[9]  = u->y*u->z*(1 - cos(a)) - u->x*sin(a);
	m[10] = cos(a) + u->z*u->z*(1-cos(a));
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;
}
void scale(float m[16], float s){
	for (int i=0; i<16; i++){
		m[i]*= s;
	}	
}

matrix scale_copy(float m[16], float s){
	matrix out;
	for (int i=0; i<16; i++){
		out.data[i] = m[i]*s;
	}
	return out;
}

matrix make_rot4( float a, float b, float c){
	matrix mat;
	float *m = mat.data;
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

matrix3 transpose9(float a[9]){
	matrix3 c;	
	int i,j;

	for (int t=0; t<9; t++){
		i = t/3;
		j = t%3;
		c.data[3*i + j]=a[3*j + i];
	}
	return c;
}

matrix transpose(float a[16]){
	matrix c;	
	int i,j;

	for (int t=0; t<16; t++){
		i = t/4;
		j = t%4;
		c.data[4*i + j]=a[4*j + i];
	}
	return c;
}

matrix transpose3(float a[16]){
	matrix c;	
	int i,j;

	for (int t=0; t<12; t++){
		i = t/4;
		j = t%4;
		c.data[4*i + j]=a[4*j + i];
	}
	for (int t=12; t<15; t++){
		c.data[t]=0;
	}
	c.data[15] = 1;
	return c;
} 

matrix inv_m(float mat[16]){

	matrix out = transpose3(mat);

	for (int i=0; i<12; i += 4){
		out.data[3 + i] = -out.data[i]*mat[3] - out.data[i + 1]*mat[7]  - out.data[i + 2]*mat[11];
	}
	return out;
}


void mult_m(float c[16], float a[16], float b[16]){	
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


matrix3 mult_m3(float a[9], float b[9]){
	matrix3 c;	
	int i,j;

	for (int t=0; t<9; t++){
		i = t/3;
		j = t%3;
		c.data[3*i + j]=0;
		for (int k=0; k<3; k++){
			c.data[3*i + j]+=b[k*3 + j]*a[i*3 + k];
			
		}	
	}
	return c;
}


matrix mult_m(float a[16], float b[16]){
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

vector2D svd(vector2D *x1, vector2D *x2){
	//a factoriser peut etre, bof
	vector2D u = vec2D(dot(x1,x1), dot(x1,x2));
	vector2D v = vec2D(u.y, dot(x2,x2));
	printVect(u);
	printVect(v);

	float a = 1, b =  - u.x -v.y, c = u.x*v.y - u.y*v.x;
	vector2D eigen = vec2D(quad_sol1(a,b,c), quad_sol2(a,b,c));
	printVect(eigen);
	return eigen;
}
//pas fini bof
matrix3 gauss(float mat[9]){
	float reduce[9];
	float inverse[9] = {1,0,0, 0,1,0, 0,0,1};
	int pivot[3];

	float max = mat[0]*mat[0];
	float div = mat[0];
	pivot[0] = 0;
	for (int i = 0; i<9; i+=3){
		if( mat[i]*mat[i] > max){
			max = mat[i]*mat[i];
			pivot[0] = i/3;
		}
	}
	int p = pivot[0];
	div = mat[3*p];
	max = -1;
	pivot[1] = (p+1)%3;

	for (int i = 0; i<9; i+=3){
		if(i/3 != p){
			inverse[i    ] = inverse[i    ] - inverse[3*p    ]*mat[i]/div;
			inverse[i + 1] = inverse[i + 1] - inverse[3*p + 1]*mat[i]/div;
			inverse[i + 2] = inverse[i + 2] - inverse[3*p + 2]*mat[i]/div;

			reduce[i    ] = 0; //mat[i    ] - mat[3*p    ]*mat[i]/div;
			reduce[i + 1] = mat[i + 1] - mat[3*p + 1]*mat[i]/div;
			reduce[i + 2] = mat[i + 2] - mat[3*p + 2]*mat[i]/div;

			if (reduce[i + 1]*reduce[i + 1] > max){
				max = reduce[i + 1]*reduce[i + 1];
				pivot[1] = i/3;
			}
		}else{
			reduce[i] = mat[i ];
			reduce[i + 1] = mat[i + 1];
			reduce[i + 2] = mat[i + 2];
		}
	}
	p = pivot[1];
	div = reduce[3*p + 1];
	max = -1;
	int np = 2;

	//printf("debut\n");
	//print_mat9(reduce);
	//print_mat9(inverse);
	//printf("----\n");

	for (int i = 0; i<9; i+=3){
		if(i/3 != p && i/3 != pivot[0]){
			//printf("debug value %f, p %d i %d value %f\n ", inverse[3*p + 1], p, i/3);
			inverse[i    ] = inverse[i    ] - inverse[3*p    ]*reduce[i+1]/div;
			inverse[i + 1] = inverse[i + 1] - inverse[3*p + 1]*reduce[i+1]/div;
			inverse[i + 2] = inverse[i + 2] - inverse[3*p + 2]*reduce[i+1]/div;

			reduce[i + 2] = reduce[i + 2] - reduce[3*p + 2]*reduce[i+1]/div;
			reduce[i + 1] = 0;//reduce[i + 1] - reduce[3*p + 1]*reduce[i+1]/div;
			pivot[2] = i/3;
			np = i/3;
		}
	}
	//printf("----\n");
	//print_mat9(reduce);
	//print_mat9(inverse);

	for (int i = 0; i<9; i+=3){
		if(i/3 != np){
			inverse[i    ] = inverse[i    ] - inverse[3*np    ]*reduce[i + 2]/reduce[3*np + 2];
			inverse[i + 1] = inverse[i + 1] - inverse[3*np + 1]*reduce[i + 2]/reduce[3*np + 2];
			inverse[i + 2] = inverse[i + 2] - inverse[3*np + 2]*reduce[i + 2]/reduce[3*np + 2];

			reduce[i + 2] = 0;//reduce[i + 2] - reduce[3*np + 2]*reduce[i + 2]/reduce[3*np + 2];
		}
	}
	//printf("wtf----\n");
	//print_mat9(reduce);
	//print_mat9(inverse);
	//printf("----\n");

	inverse[3*pivot[0]    ] = inverse[3*pivot[0]    ] - inverse[3*pivot[1]    ]*reduce[3*pivot[0] + 1]/reduce[3*pivot[1] + 1];
	inverse[3*pivot[0] + 1] = inverse[3*pivot[0] + 1] - inverse[3*pivot[1] + 1]*reduce[3*pivot[0] + 1]/reduce[3*pivot[1] + 1];
	inverse[3*pivot[0] + 2] = inverse[3*pivot[0] + 2] - inverse[3*pivot[1] + 2]*reduce[3*pivot[0] + 1]/reduce[3*pivot[1] + 1];

	reduce[3*pivot[0] + 1] = 0;//reduce[3*pivot[0] + 1] - reduce[3*pivot[1] + 1]*reduce[3*pivot[0] + 1]/reduce[3*pivot[1] + 1];

	//print_mat9(reduce);
	//print_mat9(inverse);
	
	//printf("%d %d %d \n", pivot[0] ,pivot[1] ,pivot[2]);

	for (int i=0; i<3; i++){
		inverse[3*pivot[i]    ]/= reduce[3*pivot[i]+i];
		inverse[3*pivot[i] + 1]/= reduce[3*pivot[i]+i];
		inverse[3*pivot[i] + 2]/= reduce[3*pivot[i]+i];
	}
	//print_mat9(inverse);
	matrix3 out;

	for (int i=0; i<3; i++){
		out.data[3*i    ]= inverse[3*pivot[i]    ];
		out.data[3*i + 1]= inverse[3*pivot[i] + 1];
		out.data[3*i + 2]= inverse[3*pivot[i] + 2];
	}

	matrix3 ide = mult_m3(out.data, mat);
	//print_mat9(ide.data);
	return out;	
}

#endif
