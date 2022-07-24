#pragma once
#include "../3d.h"
#include "../obj3d.h"

inline  int ISX(const int i) {return i&1;}
inline  int ISY(const int i) {return ((i)&2)/2;}
inline  int ISZ(const int i) {return ((i)&4)/4;}

struct triangle {
	vector3D a;
	vector3D b;
	vector3D c;
	vector3D u3D;
	vector3D ut3D;
	vector3D n3D;
	vector2D u;
	vector2D v;
	vector2D w;
	vector2D un;
	vector2D vn;
	vector2D wn;
	vector2D du;
	vector2D dv;
	vector2D dw;
};

inline bool inBox(const vector3D &inter,
                  const vector3D &point, 
				  float jump) {

	float fact = 1.2;
	//if (level == octree->maxLevel) fact = 1;
	return inter.x < point.x + fact*jump && inter.x >= point.x-jump*0.1 &&
		   inter.y < point.y + fact*jump && inter.y >= point.y-jump*0.1 &&
		   inter.z < point.z + fact*jump && inter.z >= point.z-jump*0.1;

}

inline vector2D to2D(const vector3D* point,
                     const triangle* tri) {
	vector2D out;
	out.x = dot(point, &tri->u3D);
	out.y = dot(point, &tri->ut3D);
	return out;
}	

vector3D to3D(triangle *tri, vector2D *out) {
 
	return vec(tri->u3D.x*out->x + tri->ut3D.x*out->y + tri->a.x,
               tri->u3D.y*out->x + tri->ut3D.y*out->y + tri->a.y,
               tri->u3D.z*out->x + tri->ut3D.z*out->y + tri->a.z);
}

vector3D to3DO(triangle *tri, vector2D *out) {

	return vec(tri->u3D.x*out->x + tri->ut3D.x*out->y,
               tri->u3D.y*out->x + tri->ut3D.y*out->y,
               tri->u3D.z*out->x + tri->ut3D.z*out->y);
}

inline float flatner(float x) {
	
	const float eps = 0.000001;
	if( x > eps) {
		return 1;
	} else if( x < -eps) {
		return -1;
	}
	return 0;
}

triangle make_triangle(mesh3D *obj, int f){

	triangle tri;
	
	vector3D u3D = minus(obj->vertics_trans + obj->face[3*f+1],
	                     obj->vertics_trans + obj->face[3*f  ]);

	vector3D v3D = minus(obj->vertics_trans + obj->face[3*f+2],
	                     obj->vertics_trans + obj->face[3*f  ]);

	vector3D w3D = minus(&v3D, &u3D);
	tri.a        = obj->vertics_trans[obj->face[3*f]];
	tri.n3D      = cross(u3D, v3D);
	normalize(&tri.n3D);

	tri.u3D = u3D;
	normalize(&tri.u3D);
	tri.ut3D = cross(tri.n3D, tri.u3D);
	normalize(&tri.ut3D);

	tri.b = obj->vertics_trans[obj->face[3*f + 1]];
	tri.c = obj->vertics_trans[obj->face[3*f + 2]];

	tri.u = to2D(&u3D, &tri);
	tri.v = to2D(&v3D, &tri);
	tri.w = to2D(&w3D, &tri);
	
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
	//project on cardinal directions. to find the shortest path in term of L1 instead of L2
	u3D = to3DO(&tri, &tri.un); 
	u3D.x = flatner(u3D.x);
	u3D.y = flatner(u3D.y);
	u3D.z = flatner(u3D.z);
	proj(&u3D, &tri.n3D);
	tri.du = to2D(&u3D, &tri);

	v3D = to3DO(&tri, &tri.vn); 
	v3D.x = flatner(v3D.x);
	v3D.y = flatner(v3D.y);
	v3D.z = flatner(v3D.z);
	proj(&v3D, &tri.n3D);
	tri.dv = to2D(&v3D, &tri);

	w3D = to3DO(&tri, &tri.wn); 
	w3D.x = flatner(w3D.x);
	w3D.y = flatner(w3D.y);
	w3D.z = flatner(w3D.z);
	proj(&w3D, &tri.n3D);
	tri.dw = to2D(&w3D, &tri);

	return tri;
}

vector3D proj(triangle *tri, vector3D *point) {
	
	vector3D ray;
	ray.x = flatner(tri->n3D.x);
	ray.y = flatner(tri->n3D.y);
	ray.z = flatner(tri->n3D.z);
	vector3D tp = minus(&tri->a, point);
	float lambda = dot(&tp, &tri->n3D)/dot(&ray, &tri->n3D);
	tp.x = point->x + lambda*ray.x;
	tp.y = point->y + lambda*ray.y;
	tp.z = point->z + lambda*ray.z;
	return tp;
	
	//vector3D pointProj = minus(point, &tri->a);
	///proj(&pointProj, &tri->n3D);
	//return pointProj;
}
static int nb_call = 0;
//TODO broken triangle proj
vector3D triangle_in_box(triangle *tri, vector3D *pointProj) {

	nb_call++;
	
	sub(pointProj, &tri->a);
	vector2D point2D = to2D(pointProj, tri);
	vector2D res = point2D;
	
	float mind = 400000000;
	
	//could be simplify
	vector2D pointOa;
	pointOa.x = -point2D.x;
	pointOa.y = -point2D.y;

	vector2D pointOu;
	pointOu.x = tri->u.x - point2D.x;
	pointOu.y = tri->u.y - point2D.y;

	float fv = dot(&pointOa, &tri->vn);
	float fw = dot(&pointOu, &tri->wn);
	float fu = dot(&pointOa, &tri->un);
	if (fv <= 0 && fw <= 0 && fu <= 0) {
		return to3D(tri, &res);
	}

	if (fv > 0) {

			res.x = point2D.x + fv*tri->vn.x;
			res.y = point2D.y + fv*tri->vn.y;
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
	
	if (fw > 0) {

		vector2D tp_res;
		tp_res.x = point2D.x + fw*tri->wn.x;
		tp_res.y = point2D.y + fw*tri->wn.y;
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
		
	if (fu > 0) {

		vector2D tp_res;
		tp_res.x = point2D.x + fu*tri->un.x;
		tp_res.y = point2D.y + fu*tri->un.y;
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


bool point_in_triangle(triangle *tri, vector3D *pointProj) {

	sub(pointProj, &tri->a);
	vector2D point2D = to2D(pointProj, tri);
	vector2D res = point2D;
	
	//could be simplify
	vector2D pointOa;
	pointOa.x = -point2D.x;
	pointOa.y = -point2D.y;

	vector2D pointOu;
	pointOu.x = tri->u.x - point2D.x;
	pointOu.y = tri->u.y - point2D.y;

	if(dot(&pointOa, &tri->vn) > 0)return false;
	if(dot(&pointOu, &tri->wn) > 0)return false;
	if(dot(&pointOa, &tri->un) > 0)return false;

	return true;
}

bool triangle_in_box2(triangle& tri, 
                      vector3D& corner,
					  float     side, 
					  int       level) {

	//the triangle is all in
	//return true;
	//if (inBox(tri.a, corner, side))return true;
	//if (inBox(tri.b, corner, side))return true;
	//if (inBox(tri.c, corner, side))return true;
	
	//the box is all in
	vector3D middle = corner;
	middle.x += side/2;
	middle.y += side/2;
	middle.z += side/2;
	vector3D projPoint = proj(&tri, &middle);
	//acc(&projPoint, &tri.a, 1);
	if (level < 7) {
		vector3D inter = triangle_in_box(&tri, &projPoint);
		
		if (dist3Dsquare(&inter, &middle) < 3/4.*side*side)return true;
		return false;
	}

	if (inBox(projPoint, corner, side) &&
	    point_in_triangle(&tri, &projPoint))return true;
	return false;
}
