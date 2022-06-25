#ifndef BITMAP
#define BITMAP
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ID(x,y,c,w) 4*((y)*w + (x)) + c
#define ID3(x,y,c,w) 3*((y)*w + (x)) + c
#define IDS(x,y,w) ((y)*w + (x))

struct bitmap{
	bitmap(){};
	bitmap(int w_, int h_):w(w_), h(w_){
		data = (unsigned char*)malloc(h_*w_*4);
		zbuffer = (float*)malloc(h_*w_*sizeof(float));
	}
	unsigned char *data = nullptr;
	float *zbuffer = nullptr;
	int w = 0;
	int h = 0;

};

void reset(bitmap *buffer){
	memset(buffer->data, 0, buffer->w*buffer->h*4);
	for (int i=0; i<buffer->w*buffer->h; i++){
		buffer->zbuffer[i]=1000;
	}
}

void readBMP(bitmap *buffer, const char *name){
	unsigned char size[2];
	

	//get size;
	FILE *file = fopen(name, "rb");
	fseek(file, 18, SEEK_SET);
	int nb_byte = fread(size,1,2,file);
	buffer->w = size[0] | (size[1] << 8);
	fseek(file, 22, SEEK_SET);
	nb_byte = fread(size,1,2,file);
	buffer->h = size[0] | (size[1] << 8);
	
	//allocat memory for pixel
	buffer->data = (unsigned char*)malloc(buffer->w*buffer->h*3);
	fseek(file, 54, SEEK_SET);
	nb_byte = fread(buffer->data,1,buffer->w*buffer->h*3,file);
}

#endif
