#include "cipher.h"

//P置换按比特的表
int PTable[32] = { 15, 6, 19, 20, 28, 11, 27, 16,
0, 14, 22, 25, 4, 17, 30, 9,
1, 7, 23, 13, 31, 26, 2, 8,
18, 12, 29, 5, 21, 10, 3, 24
};

//扩展置换E中调用的函数，从input的index位置开始取6个比特
void pick6(u16* output, u32 input, int index){
	if (index > 26){
		int carry = 32 - index;
		*output = ((input << (index - 26)) ^ (input >> (58 - index))) & 0x3f;
	}
	else{
		*output = (input >> (26 - index)) & 0x3f;
	}
}

//扩展置换E
void Expansion(u16* output, u32 input){
	for (int i = 0; i < 8; i++){
		pick6(output + i, input, ((4 * i - 1) % 32+32)%32);
	}
}

//P置换
u32 Permutation(u32 x){
	u32 y = 0;
	u32 t = 0x80000000;
	for (int i = 0; i<32; i++, t = t >> 1){
		if (PTable[i] - i>0) y |= (x << (PTable[i] - i))&t;
		else y |= (x >> (i - PTable[i])) & t;
	}
	return y;
}

//8个4比特数组转换成32比特串
u32 SboxOutput2word(u16* input){
	u32 output=0;
	for(int i=0;i<7;i++){
		output+=*(input+i);
		output<<=4;
	}
	output+=*(input+7);
	return output;
}

//输入是S盒输出，是8个4比特数组
//输出是E的输出，是8个6比特数组
//实现线性层
void PE(u16* output, u16* input){
	u32 temp=SboxOutput2word(input);
	temp=Permutation(temp);
	Expansion(output,temp);
}
