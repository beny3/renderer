#ifndef RENDER
#define RENDER
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bitmap.h"
#include "vector.h"
#include "obj3d.h"


void draw_tri(vector2D face[3],vector3D normals[3], vector3D normal_f, float z[3], vector2D &vt1, vector2D &vt2, vector2D &vt3, bitmap &buffer, int r, int g, int b, vector3D *light, int c){

	int maxx= face[0].x;
	int minx= face[1].x;
	int maxy= face[0].y;
	int miny= face[1].y;
	vector3D coord;
	vector3D n;
	vector2D p;
	vector2D ab = minus2D(&face[1], &face[0]);
	vector2D ac = minus2D(&face[2], &face[0]);

	int i, x, y;
	float cz;
	float tfx, tfy;
	
	for (i=0; i<3; i++){
		if (maxx < face[i].x) maxx = face[i].x+1;
		if (minx > face[i].x) minx = face[i].x;
		if (maxy < face[i].y) maxy = face[i].y+1;
		if (miny > face[i].y) miny = face[i].y;
	}
	/*		
	maxx=(maxx <= buffer.w-2 ? maxx+1: buffer.w-1);
	minx=(minx >= 1 ? minx-1:0);
	
	maxy=(maxy <= buffer.h-2 ? maxy+1: buffer.h-1);
	miny=(miny >= 1 ? miny-1:0);
	*/
	//printf("%d %d %d %d \n",minx,maxx, miny,maxy);
	
	for (x=minx; x<maxx; x++){
		for (y=miny; y<maxy; y++){
			if ( x > 0 && x < buffer.w && y > 0 && y < buffer.h){

				coord.z = ac.x*(face[0].y - y) - (face[0].x -x)*ac.y;
				coord.y = (face[0].x - x)*ab.y - ab.x*(face[0].y - y);
				cz = ab.x*ac.y - ac.x*ab.y;
				coord.z/=cz;
				coord.y/=cz;
				coord.x=1. -coord.y -coord.z;
				
				if(coord.x >= 0 && coord.y >= 0 && coord.z >= 0){
				
					cz = z[0]*coord.x+z[2]*coord.y+z[1]*coord.z;

					if (buffer.zbuffer[IDS(x,y,buffer.w)]>=cz){
						n = matrixMul(normals, &coord);
						//n = normal_f;
						normalize(&n);
						float diffuse = dot(&n, light);
						diffuse = fmax(diffuse, 0);
						diffuse += 0.2;
						diffuse = fmin(diffuse, 1);
					
						/*
						tfx = vt1.x*coord.x + vt2.x*coord.y + vt3.x*coord.z;
						tfy = vt1.y*coord.x + vt2.y*coord.y + vt3.y*coord.z;
					
						int tx = ((int)(tfx*texture.w) + texture.w)%texture.w;
						int ty = ((int)(tfx*texture.h) + texture.h)%texture.h;
						int idtex = ID3(tx, ty, 0, texture.w);
						*/
					
						int idpix = ID(x,y,0,buffer.w);
						//printf("%d \n", idpix);
						//int color =  (diffuse*255 < 205 ? diffuse*255 + 50: 255);
					
						buffer.data[idpix] = diffuse*r;//*texture.data[idtex];
						buffer.data[idpix+1] = diffuse*g;//*texture.data[idtex + 1];
						buffer.data[idpix+2] = diffuse*b;//*texture.data[idtex + 2];
						//buffer.data[idpix+3]= diffuse*255.;
						buffer.zbuffer[IDS(x,y,buffer.w)]=cz;
					}						
				}				
			}	
		}
	}	
}


void draw_tri(vector2D face[3], vector3D normals[3], float z[3], vector2D &vt1, vector2D &vt2, vector2D &vt3, bitmap &buffer, bitmap &texture, vector3D *light, int c){
	

	int maxx= face[0].x;
	int minx= face[1].x;
	int maxy= face[0].y;
	int miny= face[1].y;
	vector3D coord;
	vector3D n;
	vector2D p;
	vector2D ab = minus2D(&face[1], &face[0]);
	vector2D ac = minus2D(&face[2], &face[0]);

	int i, x, y;
	float cz;
	float tfx, tfy;
	
	for (i=0; i<3; i++){
		if (maxx < face[i].x) maxx = face[i].x+1;
		if (minx > face[i].x) minx = face[i].x;
		if (maxy < face[i].y) maxy = face[i].y+1;
		if (miny > face[i].y) miny = face[i].y;
	}
	/*		
	maxx=(maxx <= buffer.w-2 ? maxx+1: buffer.w-1);
	minx=(minx >= 1 ? minx-1:0);
	
	maxy=(maxy <= buffer.h-2 ? maxy+1: buffer.h-1);
	miny=(miny >= 1 ? miny-1:0);
	*/
	//printf("%d %d %d %d \n",minx,maxx, miny,maxy);
	
	for (x=minx; x<maxx; x++){
		for (y=miny; y<maxy; y++){
			if ( x > 0 && x < buffer.w && y > 0 && y < buffer.h){

				coord.z = ac.x*(face[0].y - y) - (face[0].x -x)*ac.y;
				coord.y = (face[0].x - x)*ab.y - ab.x*(face[0].y - y);
				cz = ab.x*ac.y - ac.x*ab.y;
				coord.z/=cz;
				coord.y/=cz;
				coord.x=1. -coord.y -coord.z;
				
				if(coord.x >= 0 && coord.y >= 0 && coord.z >= 0){
				
					cz = z[0]*coord.x+z[2]*coord.y+z[1]*coord.z;

					if (buffer.zbuffer[IDS(x,y,buffer.w)]>=cz){
						n = matrixMul(normals, &coord);	
						normalize(&n);
						float diffuse = dot(&n, light);
						diffuse = fmax(diffuse, 0);
					
						/*
						tfx = vt1.x*coord.x + vt2.x*coord.y + vt3.x*coord.z;
						tfy = vt1.y*coord.x + vt2.y*coord.y + vt3.y*coord.z;
					
						int tx = ((int)(tfx*texture.w) + texture.w)%texture.w;
						int ty = ((int)(tfx*texture.h) + texture.h)%texture.h;
						int idtex = ID3(tx, ty, 0, texture.w);
						*/
					
						int idpix = ID(x,y,0,buffer.w);
						//printf("%d \n", idpix);
						int color =  (diffuse*255 < 205 ? diffuse*255 + 50: 255);
					
						buffer.data[idpix]= color;//*texture.data[idtex];
						buffer.data[idpix+1]= color;//*texture.data[idtex + 1];
						buffer.data[idpix+2]= color;//*texture.data[idtex + 2];
						//buffer.data[idpix+3]= diffuse*255.;
						buffer.zbuffer[IDS(x,y,buffer.w)]=cz;
					}						
				}				
			}	
		}
	}	
}


void draw_dir(bitmap buffer, float mat[16]){
	double dir[]= { 
	0, -1, -1, //0
    -1,  0, -1, //1
	0,  0, -1, //2
	1,  0, -1, //3
	0,  1, -1, //4

    -1, -1,  0, //5
	0, -1,  0, //6
	1, -1,  0, //7

    -1,  0,  0, //8
	0,  0,  0, //9
	1,  0,  0, //10

    -1,  1,  0, //11
	0,  1,  0, //12
	1,  1,  0, //13

	0, -1,  1, //14
    -1,  0,  1, //15
	0,  0,  1, //16
	1,  0,  1, //17
	0,  1,  1  //18
	};
	vector2D center;
	double z = dir[2] + 500;
	center.x = (dir[0]*500)/(500 + z) +  buffer.w/2;
	center.y = (dir[1]*500)/(500 + z) +  buffer.h/2;
	

}
void draw_mesh4(mesh3D &obj, bitmap buffer, float mat[16]){
	vector2D face[3];
	vector3D normal[3];
	int sx = buffer.w/2;
	int sy = buffer.h/2;
	float z[3];
	vector3D light = vec(-1,0.5, -1);
	normalize(&light);
	//print_mat(mat);

	for(int i=0; i<obj.nb_f; i+=3){
		//printf("%d \n", i);
		for (int k=0; k<3; k++){
			int idx = obj.face[i+k];
			//printf("idx %d \n", idx);
			z[k] = obj.vertics_trans[idx].z + 500;
			face[k].x = (obj.vertics_trans[idx].x*500)/(500 + z[k]) + buffer.w/2;
			face[k].y = (obj.vertics_trans[idx].y*500)/(500 + z[k]) + buffer.h/2;

			normal[k].x = (obj.normals[idx].x*mat[0] +obj.normals[idx].y*mat[1] + obj.normals[idx].z*mat[2]);
			normal[k].y = (obj.normals[idx].x*mat[4] +obj.normals[idx].y*mat[5] + obj.normals[idx].z*mat[6]);
			normal[k].z = (obj.normals[idx].x*mat[8] +obj.normals[idx].y*mat[9] + obj.normals[idx].z*mat[10]);
			//printVect(*(normal + k));
			normalize(normal + k);
		}

		draw_tri(face,normal,z, obj.v_text[obj.face[i]], obj.v_text[obj.face[i+1]], obj.v_text[obj.face[i+2]],buffer, obj.texture, &light, i);
	}
}

void draw_mesh4(mesh3D &obj, bitmap buffer, float mat[16], char *which_shape){
	vector2D face[3];
	vector3D normal[3];
	int sx = buffer.w/2;
	int sy = buffer.h/2;
	float z[3];
	vector3D light = vec(-1,0.5, -1);
	normalize(&light);
	//print_mat(mat);

	for(int i=0; i<obj.nb_f; i+=3){
		//printf("%d \n", i);
		for (int k=0; k<3; k++){
			int idx = obj.face[i+k];
			//printf("idx %d \n", idx);
			z[k] = obj.vertics_trans[idx].z + 500;
			face[k].x = (obj.vertics_trans[idx].x*500)/(500 + z[k]) + buffer.w/2;
			face[k].y = (obj.vertics_trans[idx].y*500)/(500 + z[k]) + buffer.h/2;

			normal[k].x = (obj.normals[idx].x*mat[0] +obj.normals[idx].y*mat[1] + obj.normals[idx].z*mat[2]);
			normal[k].y = (obj.normals[idx].x*mat[4] +obj.normals[idx].y*mat[5] + obj.normals[idx].z*mat[6]);
			normal[k].z = (obj.normals[idx].x*mat[8] +obj.normals[idx].y*mat[9] + obj.normals[idx].z*mat[10]);
			//printVect(*(normal + k));
			normalize(normal + k);
		}
		/**/
		if (which_shape[i/3] || !which_shape[i/3]){
			draw_tri(face, normal, obj.normals_f[i/3],z, obj.v_text[obj.face[i]], obj.v_text[obj.face[i+1]], obj.v_text[obj.face[i+2]],buffer, 255, 0, 0, &light, i);
		}else{
			draw_tri(face, normal, obj.normals_f[i/3],z, obj.v_text[obj.face[i]], obj.v_text[obj.face[i+1]], obj.v_text[obj.face[i+2]],buffer, 50, 50, 50, &light, i);
		}
		/**/
	}

}

void draw_mesh(mesh3D &obj, bitmap buffer, float mat[9]){
	vector2D face[3];
	vector3D normal[3];
	int sx = buffer.w/2;
	int sy = buffer.h/2;
	float z[3];
	vector3D light = vec(0,0,-1);	
	
	for(int i=0; i<obj.nb_f; i+=3){
		for (int k=0; k<3; k++){
			int idx = obj.face[i+k];
			face[k].x = (obj.vertics[idx].x*mat[0] +obj.vertics[idx].y*mat[1] + obj.vertics[idx].z*mat[2])*15+sx;
			face[k].y = (obj.vertics[idx].x*mat[3] +obj.vertics[idx].y*mat[4] + obj.vertics[idx].z*mat[5])*15+sy;
			z[k] = (obj.vertics[idx].x*mat[6] +obj.vertics[idx].y*mat[7] + obj.vertics[idx].z*mat[8])*15;

			normal[k].x = (obj.normals[idx].x*mat[0] +obj.normals[idx].y*mat[1] + obj.normals[idx].z*mat[2]);
			normal[k].y = (obj.normals[idx].x*mat[3] +obj.normals[idx].y*mat[4] + obj.normals[idx].z*mat[5]);
			normal[k].z = (obj.normals[idx].x*mat[6] +obj.normals[idx].y*mat[7] + obj.normals[idx].z*mat[8]);

		}
		draw_tri(face,normal,z, obj.v_text[obj.face[i]], obj.v_text[obj.face[i+1]], obj.v_text[obj.face[i+2]],buffer, obj.texture, &light, i);
	}

}

void draw_cube(float  mat[16], bitmap buffer){
	float f = 500;
	vector3D light = vec(-1,-1,1);
	matrix invers_model = inv_m(mat);
	//print_mat( mat );
	float z, dz;
	matrix ide = mult_m(mat , invers_model.data);
	print_mat(ide.data);
	vector3D n;

	for (float y=0; y < 2; y += 1/256.){
		for (float x=0; x < 2; x += 1/256.){
			vector3D eye_o = vec(x, y, 0);
			vector3D view_line_o = vec(x/(f), y/(f), 1);
			vector3D eye  = mult_vm(invers_model.data, eye_o);
			vector3D view_line = mult_vm(invers_model.data, view_line_o);

			//printVect(view_line);
			//printVect(eye);
			z = buffer.zbuffer[IDS((int)(x*256),(int)(y*256),buffer.w)];

			float lambda = -eye.z/view_line.z;
			if (floor(eye.x + lambda*view_line.x)==0 && floor(eye.y + lambda*view_line.y)==0){
				z = eye_o.z + lambda*view_line_o.z;
				n = vec(mat[2], mat[6], mat[10]);
			}
			
			lambda = (1-eye.z)/view_line.z;
			if (floor(eye.x + lambda*view_line.x)==0 && floor(eye.y + lambda*view_line.y)==0){
				dz = eye_o.z + lambda*view_line_o.z;
				if (z > dz) {
					z = dz;
					n = vec(-mat[2], -mat[6], -mat[10]);
				}
			}

			lambda = (-eye.y)/view_line.y;
			if (floor(eye.z + lambda*view_line.z)==0 && floor(eye.x + lambda*view_line.x)==0){
				dz = eye_o.z + lambda*view_line_o.z;
				if (z > dz) {
					z = dz;
					n = vec(mat[1], mat[5], mat[9]);
				}
			}

			lambda = (1 - eye.y)/view_line.y;
			if (floor(eye.z + lambda*view_line.z)==0 && floor(eye.x + lambda*view_line.x)==0){
				dz = eye_o.z + lambda*view_line_o.z;
				if (z > dz) {
					z = dz;
					n = vec(-mat[1], -mat[5], -mat[9]);
				}
			}

			lambda = (-eye.x)/view_line.x;
			if (floor(eye.z + lambda*view_line.z)==0 && floor(eye.y + lambda*view_line.y)==0){
				dz = eye_o.z + lambda*view_line_o.z;
				if (z > dz) {
					z = dz;
					n = vec(mat[0], mat[4], mat[8]);
				}
			}

			lambda = (1 - eye.x)/view_line.x;
			if (floor(eye.z + lambda*view_line.z)==0 && floor(eye.y + lambda*view_line.y)==0){
				dz = eye_o.z + lambda*view_line_o.z;
				if (z > dz) {
					z = dz;
					n = vec(-mat[0], -mat[4], -mat[8]);
				}
			}
			
			int idpix = ID((int)(x*256),(int)(y*256),0,buffer.w);
			
			float diffuse = dot(&n, &light);
			diffuse = fmax(diffuse, 0);

			if (buffer.zbuffer[IDS((int)(x*256),(int)(y*256),buffer.w)] > z ){
				buffer.data[idpix] = (int)(255*diffuse);
				buffer.data[idpix+1] = 0;
				buffer.data[idpix+2] = (int)(255*diffuse);
			}
			
		}
	}
}

#endif
