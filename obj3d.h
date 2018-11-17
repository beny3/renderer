#ifndef OBJ
#define OBJ
#include <stdlib.h>
#include <stdio.h> 
#include <string.h>
#include "vector.h"
#include "bitmap.h"

typedef struct mesh3D{
	vector *vertics;
	vector *vertics_trans;  
	weight *w;
	float *wf;
	vector *normals;
	vector *normals_f;	
	vector2D *v_text;
	vector center;
	float radius;
	int nb;
	int nb_t;
	int nb_f;
	int *face;
	bitmap texture;
} mesh3D;

struct colition_hull{
	mesh3D mesh;
	int	*vect_graph_table;
	int faces_per_vertex;	
	int current_vert;
};

struct bounding_box{
	vector zero;
	float limits[3]={0,0,0};
};

struct physic{
	vector position;
	vector velocity;
	vector moment;
	vector rest_angle;
	float angular_vel;
	float angle;
	matrix mat;
};

bounding_box make_bounding_box(mesh3D *obj){
	bounding_box out;
	float min[3]={1000000,1000000,1000000};
	float max[3]={-1000000,-1000000,-1000000};

	for (int i=0; i < obj->nb; i++){
		for (int j=0; j<3; j++){
			if(obj->vertics[i].el[j] < min[j]){min[j]=obj->vertics[i].el[j];}
			if(obj->vertics[i].el[j] > max[j]){max[j]=obj->vertics[i].el[j];}
		}
	}
	out.limits[0]=max[1]-min[1];
	out.limits[1]=max[0]-min[0];
	out.limits[2]=max[2]-min[2];
	out.zero = vec(min[0], min[1], min[2]);

	return out;
}


void make_normal(mesh3D *obj){

	for (int i=0; i < obj->nb_f; i+=3){

		vector ab = minus(&obj->vertics[obj->face[i+1]], &obj->vertics[obj->face[i]]);
		vector ac = minus(&obj->vertics[obj->face[i+2]], &obj->vertics[obj->face[i]]);
		normalize(&ac);
		normalize(&ab);
		vector n = cross(ab,ac);
		
		obj->normals_f[i/3] = n;
		normalize(&obj->normals_f[i/3]);
		
		acc(&obj->normals[obj->face[i]]  , &n);
		acc(&obj->normals[obj->face[i+1]], &n);
		acc(&obj->normals[obj->face[i+2]], &n);
	
	}
	for (int i=0; i < obj->nb; i++){
		normalize(&obj->normals[i]);
		//printf("n %f \n", dot(&obj->normals[i], &obj->normals[i]));
	}
}


void euler(physic *phy, float dt){

	vector m = phy->moment;
	float theta = normalize(&m);

	matrix rot;
	make_rot4_u(rot.data, &m, phy->angle);
	matrix rot_rest = make_rot4(phy->rest_angle.x, phy->rest_angle.y, phy->rest_angle.z);
	mult_m(phy->mat.data, rot.data, rot_rest.data);
	make_mov4(phy->mat.data, phy->position);

	acc(&phy->position, &phy->velocity, dt);
	phy->angle += theta*dt;
}

void free_mesh3D(mesh3D *obj){
	delete[] obj->vertics;
	delete[] obj->vertics_trans;
	//delete[] obj->s;
	//delete[] obj->gk;
}

int read_mesh3D(mesh3D *obj, const char* name){
	FILE *file;
	char ligne[256];
	ssize_t read;
	int v=0;
	int vn=0;
	int vt=0;
	int f=0;
	double f1,f2,f3;
	
	file = fopen(name, "r");
	if (file == NULL){
		return -1;
	}
	obj->nb=0;
	obj->nb_f=0;
	obj->nb_t=0;

	while(fgets(ligne, 256, file) != NULL){
		if (ligne[0]=='v' && ligne[1]==' '){

			obj->nb++;
		}
		if (ligne[0]=='v' && ligne[1]=='t'){

			obj->nb_t++;
		}
		if (ligne[0]=='f' ){

			obj->nb_f+=3;
		}
	}
	printf("face number %d \n", obj->nb_f/3);

	obj->vertics = new vector[obj->nb];
	obj->normals = new vector[obj->nb];
	obj->vertics_trans = new vector[obj->nb];
	obj->w = new weight[obj->nb*3];

	vector *vertics_trans; 
	obj->v_text = new vector2D[obj->nb_t];
	obj->face = new int[obj->nb_f];
	fseek(file,0,SEEK_SET);
	obj->radius = 0;

	while(fgets(ligne, 256, file) != NULL){
		if (ligne[0]=='v' && ligne[1]==' '){
			sscanf(ligne, "v %f %f %f", &obj->vertics[v].x, &obj->vertics[v].y, &obj->vertics[v].z);
			float r = dot(&obj->vertics[v], &obj->vertics[v]);
			if (r > obj->radius)obj->radius = r;
			obj->normals[v] = vec(0,0,0);
			v++;	
		}
		if (ligne[0]=='v' && ligne[1]=='t'){	
			sscanf(ligne, "vt %f %f", &obj->v_text[vt].x, &obj->v_text[vt].y);
			vt++;
		}
		if (ligne[0]=='f' ){	
			sscanf(ligne, "f %d%*[^ ]%d%*[^ ]%d%*[^\n]", &obj->face[f], &obj->face[f+1], &obj->face[f+2]);
			obj->face[f]--;
			obj->face[f+1]--;
			obj->face[f+2]--;
			f+=3;
		}
	}
	fclose(file);
	//calcul normal
	obj->normals_f = new vector[f/3];
	make_normal(obj);

	return 0;

}

void compute_center(mesh3D *obj){
	vector center=vec(0,0,0);
	for (int i = 0; i< obj->nb; i++){
		acc(&center, &obj->vertics[i]);
	}
	scalar(&center, 1./obj->nb);
	obj->center = center;
}

mesh3D make_cylinder(float x, float y, float z, float sy, float  r, int nb_l, int nb_s){

	mesh3D cylinder;
	cylinder.nb = nb_l*nb_s;
	cylinder.nb_f = (nb_s-1)*nb_l*2*3;

	vector *v=new vector[cylinder.nb];
	int *f=new int[cylinder.nb_f];

	cylinder.vertics = v;
	cylinder.vertics_trans = new vector[cylinder.nb];
	cylinder.face  = f;

	int k=0;
	float stepY = sy/(nb_s-1);
	for (int i=0; i<nb_s; i++){
		for (int j=0; j<nb_l; j++){
			float theta = (2*PI*j)/nb_l;
			v[k].x = x + cos(theta)*r;
			v[k].z = z + sin(theta)*r;
			v[k].y = y + i*stepY;
			k++;
		}
	}
	k=0;
	for (int i=0; i<nb_s-1; i++){
		for (int j=0; j<nb_l; j++){
			f[k  ]=i*nb_l + j;
			f[k+1]=i*nb_l + (j+1)%nb_l;
			f[k+2]=(i+1)*nb_l + j;

			f[k+3]=(i+1)*nb_l + j;
			f[k+4]=i*nb_l + (j+1)%nb_l;
			f[k+5]=(i+1)*nb_l + (j+1)%nb_l;
			k+=6;
		}
	}
	
	for (int i=0; i<nb_s-1; i++){
	}

	compute_center(&cylinder);

	return cylinder;
	
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
	v[0] = vec(x - sx, y - sy, z-sz);
	v[1] = vec(x + sx, y - sy, z-sz);
	v[2] = vec(x + sx, y + sy, z-sz);
	v[3] = vec(x - sx, y + sy, z-sz);

	v[4] = vec(x - sx, y - sy, z + sz);
	v[5] = vec(x + sx, y - sy, z + sz);
	v[6] = vec(x + sx, y + sy, z + sz);
	v[7] = vec(x - sx, y + sy, z + sz);

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

	compute_center(&rect);

	return rect;
}


colition_hull make_graph(mesh3D obj){
	int nb = obj.nb;
	int nb_f = obj.nb_f;
	int *face = obj.face;

	///find maximum faces connected to a vertex
	int *count = new int[nb];
	memset(count, 0, nb*sizeof(int));

	for(int i=0; i< nb_f; i++){
		count[face[i]]++;
	}

	int max = 0;
	for (int i=0; i<nb; i++){
		if ( count[i] + 1 > max){
			max = count[i] + 1;
		}
	}

	delete[] count;

	int *vect_graph_table = new int[nb*max];
	for (int i=0; i < nb*max; ++i){
		vect_graph_table[i] = -1;		
	}
	
	for (int i=0; i < nb_f; i+=3){

		for (int j = 0; j<3; j++){

			int vertex_idx = face[i + j]*max;
			for (int k = 1; k<3; k++){

				int voisin = face[i + (j+k)%3];
				int is_not_there = 1;

				for (int l=0; l<max; l++){
					if (vect_graph_table[vertex_idx + l] == -1 && is_not_there){
						vect_graph_table[vertex_idx + l] = voisin;
						break;
				
					}
					if (vect_graph_table[vertex_idx + l] == voisin){
						is_not_there = 0;
						break;
					}			
				}
			}
		}
	}

 	colition_hull hull;
	hull.faces_per_vertex = max;
	hull.vect_graph_table = vect_graph_table;
	hull.mesh = obj;
	hull.current_vert = 0;
		
	return hull;
}

void print_hull(colition_hull *h){
	for (int i=0; i < h->mesh.nb; i++){
		for (int j=0; j < h->faces_per_vertex; j++){
			int *line = &h->vect_graph_table[i*h->faces_per_vertex];
			printf("%4d ", line[j] );
		}
		printf(" \n");
	}
}

#endif
