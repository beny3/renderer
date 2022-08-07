#pragma once
#include "triangle.h"
#include <vector>

const int maxVoxel = 1 << 10;
const int maxIdx = 1 << 20;

struct Node{
	vector3D corner;
	float size;
	int childs[8] = {0, 0, 0, 0, 0, 0, 0, 0};
};

struct Octree {
    int maxLevel = 8;
	std::vector<Node> nodes;
};

void octree_rec(triangle *tri, Octree *octree, int faceId, int nodeId, int level){
  
	vector3D tp = octree->nodes[nodeId].corner;
	//octree->nodes.reserve(10000000);
	vector3D projvec;
	float jump = octree->nodes[nodeId].size/2.;  
	//tp.x += jump/2.;
	//tp.y += jump/2.;
	//tp.z += jump/2.;
	
	float boundingSphere = jump*jump*3;

	for(int i = 0; i < 8; ++i){

		vector3D point = tp;

		//Loop over bottom corner of boxes
		point.x += ISX(i)*jump;
		point.y += ISY(i)*jump;
		point.z += ISZ(i)*jump;

		if (triangle_in_box2(*tri, point, jump, level)) {
			if (level == octree->maxLevel) {
				octree->nodes[nodeId].childs[i] = -faceId;
			} else {
				if (octree->nodes[nodeId].childs[i] == 0) {

					octree->nodes.emplace_back();
					int idx = octree->nodes.size();
					octree->nodes[nodeId].childs[i] = idx;

					octree->nodes[idx-1].corner = octree->nodes[nodeId].corner;
					octree->nodes[idx-1].corner.x += ISX(i)*jump;
					octree->nodes[idx-1].corner.y += ISY(i)*jump;
					octree->nodes[idx-1].corner.z += ISZ(i)*jump;
					octree->nodes[idx-1].size = jump;

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

int testSpeed(triangle *tri) {
	int total = 0;
	for(int z = 0 ; z < 200; z++) {
		for(int y = 0 ; y < 200; y++) {
			for(int x = 0 ; x < 200; x++) {
				vector3D corner = vec(x,y,z); 
				total += triangle_in_box2(*tri, corner, 2, 0);
			}
		}
	}
	return total;
}

Octree make_octree(mesh3D *obj){
    
	const float eps = 1;

	Octree octree;
	octree.nodes.reserve(10*maxVoxel);
	Node first;
	first.corner = obj->bboxa.zero;
	first.corner.x -= eps;
	first.corner.y -= eps;
	first.corner.z -= eps;
	printVect(first.corner);
    //first.corner.x = -200;
	//first.corner.y = -200;
	//first.corner.z = -213.25;
	//first.size = 450;
	first.size = obj->bboxa.edge.x + 2*eps;
	if(obj->bboxa.edge.y > first.size)first.size = obj->bboxa.edge.y + 2*eps;
	if(obj->bboxa.edge.z > first.size)first.size = obj->bboxa.edge.z + 2*eps;
	printf("size %f \n", first.size);
	octree.nodes.push_back(first);
	for(int i = 0; i < obj->nb_f/3; ++i) {
		triangle tri = make_triangle(obj, i);;
		//octree.maxLevel = testSpeed(&tri);
		octree_rec(&tri, &octree, i+1, 0, 0);		
	}
	return octree;	
}	

void dfsRec(std::vector<vector3D> &res, Octree *octree, int nodeId){

	Node *node =  &octree->nodes[nodeId];

	for(int i = 0; i < 8; ++i){
			
		if (node->childs[i] > 0) {
			 dfsRec(res, octree, node->childs[i] - 1);
		} else if (node->childs[i] < 0)  {
			//printf("%d \n", nodeId);
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
		dfsRec(res, octree, 0);
	}
	return res;
}

