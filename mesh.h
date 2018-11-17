#ifndef MESH_H
#define MESH_H
#include <stdlib.h>
#include <stdio.h> 
#include <string.h>

typedef struct mesh{
        double *vertics;
        double *uk;
        double *gk;
		unsigned int *sem1;
		unsigned int *sem2;
        size_t nb;
        size_t nb_f;
        unsigned short *face;
        double mass;
        double X[3];
        double U[3];
	    double *s;
} mesh;


typedef struct mesh_cuda{
		size_t nb;
        double *vx;
		double *vy;
		double *vz;
		
		double *gx;
		double *gy;
		double *gz;
		
		double *delta_ux;
		double *delta_uy;
		double *delta_uz;
		
	    double *s;
} mesh_cuda;

void free_mesh(mesh *obj){
	delete[] obj->vertics;
	delete[] obj->s;
	delete[] obj->gk;
}

void write_mesh(mesh *obj){
	FILE *f = fopen("sphere.ib", "w");
	double* v = obj->vertics;
	double* s = obj->s;

	fprintf(f, "%zd \n", obj->nb);
	for (unsigned int i=0; i < obj->nb;  i+=3){
		fprintf(f, "%f %f %f\n", v[i], v[i+1], v[i+2]);
	}

	for (unsigned int i=0; i < obj->nb/3;  i+=3){
		fprintf(f, "%f %f %f\n", s[i], s[i+1], s[i+2]);
	}

	fclose(f);
}


int read_mesh(mesh *obj, const char* name){
	FILE *f;
	char ligne[256];
	
	f = fopen(name, "r");
	if (f == NULL){
		return -1;
	}

	fgets(ligne, 256, f);
	sscanf(ligne, "%zd", &obj->nb);
	obj->s = new double[obj->nb/3]();
	obj->vertics = new double[obj->nb]();
	obj->gk = new double[obj->nb]();
	memset(obj->gk, 0, obj->nb*sizeof(double));
	printf("np points %zd \n", obj->nb);

	for (size_t i=0; i<obj->nb; i+=3){
		if(fgets(ligne, 256, f) == NULL){
			return -1;
		}
		sscanf(ligne, "%lf %lf %lf", &obj->vertics[i], &obj->vertics[i + 1], &obj->vertics[i + 2]);
		//printf("%f %f %f\n", obj->vertics[i], obj->vertics[i + 1], obj->vertics[i + 2]);
	}
	
	for (size_t i=0; i<obj->nb/3; i+=3){
		if(fgets(ligne, 256, f) == NULL){
			return -1;
		}
		sscanf(ligne, "%lf %lf %lf", &obj->s[i], &obj->s[i + 1], &obj->s[i + 2]);
		//printf(ligne, "%f %f %f", obj->s[i], obj->s[i + 1], obj->s[i + 2]);
	}
	
	fclose(f);
	return 0;

}


#endif
