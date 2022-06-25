#pragma once
#include "../3d.h"
#include "../obj3d.h"

struct triangle {
	vector3D a;
	vector3D u3D;
	vector3D ut3D;
	vector3D n3D;
	vector2D u;
	vector2D v;
	vector2D w;
	vector2D un;
	vector2D vn;
	vector2D wn;
};

triangle make_triangle(mesh3D *obj, int f){

	triangle tri;
	vector3D u3D = minus(obj->vertics_trans + obj->face[3*f+1], obj->vertics_trans + obj->face[3*f]);
	vector3D v3D = minus(obj->vertics_trans + obj->face[3*f+2], obj->vertics_trans + obj->face[3*f]);
	vector3D w3D = minus(&v3D, &u3D);
	tri.a = obj->vertics_trans[obj->face[3*f]];
	tri.n3D = cross(u3D, v3D);
	normalize(&tri.n3D);

	tri.u3D = u3D;
	normalize(&tri.u3D);
	tri.ut3D = cross(tri.n3D, tri.u3D);
	normalize(&tri.ut3D);

	tri.u.x = dot(&u3D, &tri.u3D);
	tri.u.y = dot(&u3D, &tri.ut3D);

	tri.v.x = dot(&v3D, &tri.u3D);
	tri.v.y = dot(&v3D, &tri.ut3D);
	
	tri.w.x = dot(&w3D, &tri.u3D);
	tri.w.y = dot(&w3D, &tri.ut3D);
	
	vector2D uva;
	vector2D medu;
	vector2D medv;
	uva.x = (tri.u.x + tri.v.x)/2.;
	uva.y = (tri.u.y + tri.v.y)/2.;
	
	medu.x = tri.u.x/2. - tri.v.x;
	medu.y = tri.u.y/2. - tri.v.y;
	
	medv.x = tri.v.x/2. - tri.u.x;
	medv.y = tri.v.y/2. - tri.u.y;	
	
	tri.un.x = -tri.u.y;
	tri.un.y =  tri.u.x;
	
	tri.vn.x = -tri.v.y;
	tri.vn.y =  tri.v.x;
	
	tri.wn.x = -tri.w.y;
	tri.wn.y = -tri.v.x;
	normalize(&tri.un);
	normalize(&tri.vn);
	normalize(&tri.wn);
	
	if (dot(&tri.wn, &uva) > 0) {
		tri.wn.x = -tri.wn.x;
		tri.wn.y = -tri.wn.y;
	}
	
	if (dot(&tri.un, &medu) > 0) {
		tri.un.x = -tri.un.x;
		tri.un.y = -tri.un.y;
	}
	
	if (dot(&tri.vn, &medv) > 0) {
		tri.vn.x = -tri.vn.x;
		tri.vn.y = -tri.vn.y;
	}
	
	return tri;
}

vector3D to3D(triangle *tri, vector2D *out){

	return vec(tri->u3D.x*out->x + tri->ut3D.x*out->y + tri->a.x,
               tri->u3D.y*out->x + tri->ut3D.y*out->y + tri->a.y,
               tri->u3D.z*out->x + tri->ut3D.z*out->y + tri->a.z);
}

vector3D proj(triangle *tri, vector3D *point){

	vector3D pointProj = minus(point, &tri->a);
	proj(&pointProj, &tri->n3D);
	return pointProj;
}
static int nb_call = 0;
//TODO broken triangle proj
vector3D triangle_in_box(triangle *tri, vector3D *pointProj, vector3D *point) {

	nb_call++;
	
	vector2D point2D;
	vector2D res;
	
	point2D.x = dot(pointProj, &tri->u3D);
	point2D.y = dot(pointProj, &tri->ut3D);
	
	//*debug = to3D(tri, &point2D);
	
	res = point2D;
	
	float mind = 400000000;
	
	//could be simplify
	vector2D pointOa;
	pointOa.x = -point2D.x;
	pointOa.y = -point2D.y;

	vector2D pointOu;
	pointOu.x = tri->u.x - point2D.x;
	pointOu.y = tri->u.y - point2D.y;	
	float f = dot(&pointOa, &tri->vn);

	if (f > 0) {

	  res.x = point2D.x + f*tri->vn.x;
	  res.y = point2D.y + f*tri->vn.y;
	  //if triangle is not degenerate v.y != 0;
	  float scale = res.y/tri->v.y;
	  
	  if (scale < 0) {
	  	res.x = 0;
	  	res.y = 0;
	  	mind = dot(&point2D, &point2D);
	  
	  } else if (scale > 1) {
	  
	  	res.x = tri->v.x;
	  	res.y = tri->v.y;
	  	mind = dist2Dsquare(&point2D, &tri->v);

	  } else {
	  	return to3D(tri, &res);
	  }
	}
	//TODO find if we are inside and return early
	f = dot(&pointOu, &tri->wn);
	
	if (f > 0) {

		vector2D tp_res;
		tp_res.x = point2D.x + f*tri->wn.x;
		tp_res.y = point2D.y + f*tri->wn.y;
		//if triangle is not degenerate w.y != 0;
		float scale = (tp_res.y - tri->u.y)/tri->w.y;

		if (scale < 0) {
		
			float d = dist2Dsquare(&point2D, &tri->u);
			if(d < mind){
				res.x = tri->u.x;
				res.y = tri->u.y;
				mind = d;
			}

		} else if (scale > 1) {
		
			float d = dist2Dsquare(&point2D, &tri->v);
			if(d < mind){
				res.x = tri->v.x;
				res.y = tri->v.y;
				mind = d;
			}

		} else {
			return to3D(tri, &tp_res);
		}
	}
	
	f = dot(&pointOa, &tri->un);
	
	if (f > 0) {

		vector2D tp_res;
		tp_res.x = point2D.x + f*tri->un.x;
		tp_res.y = point2D.y + f*tri->un.y;
		//if triangle is not degenerate u.x == 1;
		float scale = tp_res.x/tri->u.x;
		//TODO check if we tried this corner give up otherwise
		if (scale < 0) {
		
			float d = dot(&point2D, &point2D);
			if (d < mind) {
				res.x = 0;
				res.y = 0;
				mind = d;
			}

		} else if (scale > 1) {
		
			float d = dist2Dsquare(&point2D, &tri->u);
			if (d < mind) {
				res.x = tri->u.x;
				res.y = tri->u.y;
				mind = d;
			}

		} else {
			return to3D(tri, &tp_res);
		}
	}
	return to3D(tri, &res);
}
