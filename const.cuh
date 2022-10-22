#pragma once
#pragma warning

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h> 
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
//#include <device_atomic_functions.hpp>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include <chrono>
#include <random>


#define IN_FILE_1 "MPS_physical.txt"
#define IN_FILE_2 "initial.txt"

#define Dns_Num 5
#define GST -1//�S�[�X�g���q
#define FLD 0//���̗��q
#define WLL 1//�Ǘ��q
#define OBJ 2//���Ǘ��q
#define OBJ2 3//���Ǘ��q
#define MRR 4//�~���[���q

#define Surface 0//�ǖ�
#define Edge 1//�Ǌp

#define Dns_FLD 1000.0f
#define Dns_WLL 1000.0f
#define Dns_OBJ 1000.0f
#define Dns_OBJ2 1000.0f

#define pi 3.141592f
//#define WEI(dst, re) ((re/dst) - 1.0f)
#define WEI(dst, re) ((re/dst) + (dst/re) - 2.0f)
#define WEI_grad(dst, re) ((re/dst) - (dst/re))

#define THREADS 512

#define NumMRR 10//1�̗��̗��q���琶�������~���[���q�̍ő�l(�\���ȔC�ӂ̑傫��)

typedef float real;
//typedef double real;

typedef struct {
	real x, y, z;
}treal3;

typedef struct {
	char x, y, z;
}tchar3;

typedef struct {
	real max, min;
}treal_m2;

typedef struct {
	real* x;
	real* y;
	real* z;
}areal3;

typedef struct {
	int* x;
	int* y;
	int* z;
}aint3;

typedef struct {
	real* x;
	real* y;
	real* z;
}host_vec3;


//CUDA error check
#define CHECK(call)																							\
{																														\
	const cudaError_t error = call;																		\
	if (error != cudaSuccess) {																				\
		printf_s("Error: %s:%d, ", __FILE__, __LINE__);										\
		printf_s("code:%d, reason: %s\n", error, cudaGetErrorString(error));			\
		exit(1);																										\
	}																													\
}