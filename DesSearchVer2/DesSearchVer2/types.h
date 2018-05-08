#pragma once

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#include <immintrin.h>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <bitset>
#include <algorithm>

using namespace std;

#define BLOCK_SIZE_INBIT (64)
#define ROUND_SIZE_INBIT (32)
#define E_OUTPUT_SIZE_INBIT (48)
#define SBOX_INPUT_SIZE_INBIT (6)
#define SBOX_OUTPUT_SIZE_INBIT (4)
#define SBOX_NUMBER (8)
#define SBOX_INPUTS_NUMBER (64)
#define SBOX_OUTPUTS_NUMBER (16)
#define ROUND_N (16)
#define PR_NUMBER (8)

#define WMAX_CIPHER (64)
#define WMIN_S 2
#define WMAX_S 5
#define WMIN_R 2
#define WMAX_R (40)

typedef char si8;
typedef short si16;
typedef int si32;
typedef long long si64;

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long long u64;
typedef float prType;

#define ALIGNED_(x) __declspec(align(x))
#define ALIGNED_TYPE_(t,x) t ALIGNED_(x)