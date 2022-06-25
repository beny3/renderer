#pragma once
#include "triangle.h"
#include <vector>

#define ISX(i) ((i)&1)
#define ISY(i) (((i)&2)/2)
#define ISZ(i) (((i)&4)/4)

const int maxVoxel = 1 << 10;
const int maxIdx = 1 << 20;

struct Node{
	vector3D corner;
	float size;
	int childs[8] = {0, 0, 0, 0, 0, 0, 0, 0};
};

struct Octree {
    int maxLevel = 10;
	std::vector<Node> nodes;
};

void octree_rec(triangle *tri, Octree *octree, int faceId, int nodeId, int level){
  
	vector3D tp = octree->nodes[nodeId].corner;
	vector3D projvec;
	float jump = octree->nodes[nodeId].size/2;  
	tp.x += jump/2;
	tp.y += jump/2;
	tp.z += jump/2;
	
	float boundingSphere = jump*jump*3;

	for(int i = 0; i < 8; ++i){

		vector3D point = tp;
		point.x += ISX(i)*jump;
		point.y += ISY(i)*jump;
		point.z += ISZ(i)*jump;
		
		projvec = proj(tri, &point);

		vector3D projvecPlusA = add(&projvec, &(tri->a));
		
		//if(dist3Dsquare(&projvecPlusA, &point) > boundingSphere)continue;
	
		vector3D inter = triangle_in_box(tri, &projvec, &point);

		//TODO jumb should be power of two?
		if( (int)(inter.x/jump) == (int)(point.x/jump) && 
			(int)(inter.y/jump) == (int)(point.y/jump) &&
			(int)(inter.z/jump) == (int)(point.z/jump)){

			if (level == octree->maxLevel){
				octree->nodes[nodeId].childs[i] = -faceId;
			} else {

				if (octree->nodes[nodeId].childs[i] == 0) {

					Node n;
					n.corner = octree->nodes[nodeId].corner;
					n.corner.x += ISX(i)*jump;
					n.corner.y += ISY(i)*jump;
					n.corner.z += ISZ(i)*jump;
					n.size = jump;
					octree->nodes.push_back(n);
					int idx = octree->nodes.size();
					octree->nodes[nodeId].childs[i] = idx;
					//TODO the index minus one appread two time maybe it should be in a function
					octree_rec(tri, octree, faceId, idx - 1, level + 1);	
				} else {
					int next_idx = octree->nodes[nodeId].childs[i] - 1;
					octree_rec(tri, octree, faceId, next_idx, level + 1);					
				}
			}	
		}
	}
}

Octree make_octree(mesh3D *obj){
    
	Octree octree;
	octree.nodes.reserve(10*maxVoxel);
	Node first;
	first.corner.x = -200;
	first.corner.y = -200;
	first.corner.z = -200;
	first.size = 400;
	octree.nodes.push_back(first);	
	triangle tri = make_triangle(obj, 0);
	octree_rec(&tri, &octree, 1, 0, 0);
	return octree;
}	

void dfsRec(std::vector<vector3D> &res, Octree *octree, Node *node){

	for(int i = 0; i < 8; ++i){
			
		if (node->childs[i] > 0) {
			 dfsRec(res, octree, &octree->nodes[node->childs[i] - 1]);
		} else if (node->childs[i] < 0)  {
			vector3D point = node->corner;
			point.x += node->size/2;
			point.y += node->size/2;
			point.z += node->size/2;
			res.push_back(node->corner);
		}
	}
}

std::vector<vector3D> dfsOctree(Octree *octree){
	std::vector<vector3D> res;
	if(octree->nodes.size() > 0){
		dfsRec(res, octree, &octree->nodes[0]);
	}
	return res;
}


