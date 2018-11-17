#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "vector.h"


inline float hermit(float p1, float p2, float r1, float r2, float t){
	float ts = t*t;
	float tc = t*ts;
	//printf("tc - ts %f \n", r2*(tc - ts));
	return p1*(2*tc - 3*ts + 1) + p2*(- 2*tc + 3*ts) + r1*(tc -2*ts + t) + r2*(tc - ts);
}

struct keyFrame{
	vector p;
	vector t;
	int temps;
};

struct keyFrame4D{
	vector4D p;
	vector4D t;
	int temps;
};


typedef struct sequence{
	int type = 1;
	union{
		keyFrame *frames;
		keyFrame4D *frames4D;
	};
	int current = 0;
	int nb = 0;
} sequence;

sequence make_sequence(vector *points, int nb, int temps){
	sequence seq;
	seq.frames = new keyFrame[nb];
	seq.nb = nb;
	int t = 0;

	for (int i=0; i<nb; i++){
		seq.frames[i].p = points[i];
		seq.frames[i].temps = t;
		t += temps;
	}
	return seq;
}

sequence make_sequence(vector4D *points, int nb, int temps){
	sequence seq;
	seq.type = 0;
	seq.frames4D = new keyFrame4D[nb];
	seq.nb = nb;
	int t = 0;

	for (int i=0; i<nb; i++){
		seq.frames4D[i].p = points[i];
		seq.frames4D[i].temps = t;
		t += temps;
	}
	return seq;
}


void make_tangents(sequence *seq){
	if (seq->type){
		keyFrame *frames = seq->frames;
		vector tangent1, tangent2, tangent3;
		float l1, l2, l3;
	
		frames[0].t = minus(&frames[1].p, &frames[0].p);
		tangent3 = minus(&frames[2].p, &frames[1].p);
		frames[1].t = minus(&frames[2].p, &frames[0].p);

		normalize(&frames[1].t);
		l1 = normalize(&frames[0].t);
		l3 = normalize(&tangent3);
		//printf("mk tangent %f %f\n",l1, l3 );

		l2 = (l1  +  l3)/2 *(0.25 + (dot(&frames[0].t, &tangent3)+1)/2*0.75);
		scalar(&frames[0].t, l1);
		scalar(&frames[1].t, l2);

		for (int i=1; i < seq->nb-2; i++){
			//printf("%d /n", i);

			tangent1 = minus(&frames[i+1].p, &frames[i].p);
			tangent3 = minus(&frames[i+2].p, &frames[i+1].p);
			frames[i + 1].t = minus(&frames[i+2].p, &frames[i].p);

			normalize(&frames[i + 1].t);
			l1 = normalize(&tangent1);
			l3 = normalize(&tangent3);

			l2 = (l1  +  l3)/2 *(0.25 + (dot(&tangent1, &tangent3)+1)/2*0.75);
			scalar(&frames[i + 1].t, l2);
		}

		frames[seq->nb-1].t = minus(&frames[seq->nb-1].p, &frames[seq->nb-2].p);
	}else{
		keyFrame4D *frames = seq->frames4D;
		vector4D tangent1, tangent2, tangent3;
		float l1, l2, l3;
	
		frames[0].t = minus(&frames[1].p, &frames[0].p);
		tangent3 = minus(&frames[2].p, &frames[1].p);
		frames[1].t = minus(&frames[2].p, &frames[0].p);

		normalize(&frames[1].t);
		l1 = normalize(&frames[0].t);
		l3 = normalize(&tangent3);
		//printf("mk tangent %f %f\n",l1, l3 );

		l2 = (l1  +  l3)/2 *(0.25 + (dot(&frames[0].t, &tangent3)+1)/2*0.75);
		scalar(&frames[0].t, l1);
		scalar(&frames[1].t, l2);

		for (int i=1; i < seq->nb-2; i++){
			//printf("%d /n", i);

			tangent1 = minus(&frames[i+1].p, &frames[i].p);
			tangent3 = minus(&frames[i+2].p, &frames[i+1].p);
			frames[i + 1].t = minus(&frames[i+2].p, &frames[i].p);

			normalize(&frames[i + 1].t);
			l1 = normalize(&tangent1);
			l3 = normalize(&tangent3);

			l2 = (l1  +  l3)/2 *(0.25 + (dot(&tangent1, &tangent3)+1)/2*0.75);
			scalar(&frames[i + 1].t, l2);
			}
		frames[seq->nb-1].t = minus(&frames[seq->nb-1].p, &frames[seq->nb-2].p);

	}
}



void interpolate(vector *v, sequence *seq, int t){
	 keyFrame *frames = seq->frames;
	 int c = seq->current;

	 t = t%(frames[seq->nb-1].temps+1);
	 t -= frames[c].temps;

	 float rt = (float)t/(frames[c+1].temps - frames[c].temps);

	 v->x = hermit(frames[c].p.x, frames[c + 1].p.x, frames[c].t.x, frames[c + 1].t.x, rt);
	 v->y = hermit(frames[c].p.y, frames[c + 1].p.y, frames[c].t.y, frames[c + 1].t.y, rt);
	 v->z = hermit(frames[c].p.z, frames[c + 1].p.z, frames[c].t.z, frames[c + 1].t.z, rt);

	 if( rt>= 1){
		seq->current= (seq->current+1)%(seq->nb-1);
		c = seq->current;
	 }
}

void interpolate(vector4D *v, sequence *seq, int t){
	 keyFrame4D *frames = seq->frames4D;
	 int c = seq->current;

	 t = t%(frames[seq->nb-1].temps+1);
	 t -= frames[c].temps;

	 float rt = (float)t/(frames[c+1].temps - frames[c].temps);

	 v->x = hermit(frames[c].p.x, frames[c + 1].p.x, frames[c].t.x, frames[c + 1].t.x, rt);
	 v->y = hermit(frames[c].p.y, frames[c + 1].p.y, frames[c].t.y, frames[c + 1].t.y, rt);
	 v->z = hermit(frames[c].p.z, frames[c + 1].p.z, frames[c].t.z, frames[c + 1].t.z, rt);
	 v->phi = hermit(frames[c].p.phi, frames[c + 1].p.phi, frames[c].t.phi, frames[c + 1].t.z, rt);

	 if( rt>= 1){
		seq->current= (seq->current+1)%(seq->nb-1);
		c = seq->current;
	 }
}




