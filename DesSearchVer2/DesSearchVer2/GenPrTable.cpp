#include "GenPrTable.h"
#include <algorithm>

int PDT[SBOX_INPUTS_NUMBER][SBOX_OUTPUTS_NUMBER];
int PDT_Number;
int PDT_0_Number[SBOX_INPUTS_NUMBER][PR_NUMBER];
int PDT_0_Offset[SBOX_INPUTS_NUMBER][PR_NUMBER][2];
int PDT_1_Number[PR_NUMBER];
int PDT_1_Offset[PR_NUMBER][2];
int PDT_1_Non0Num[PR_NUMBER];
int PDT_1_Non0Val[PR_NUMBER][SBOX_INPUTS_NUMBER];

prType PDT_MaxProb[SBOX_INPUTS_NUMBER];

ALIGNED_TYPE_(u16,8) SPE[SBOX_NUMBER][SBOX_INPUTS_NUMBER][MAX_OUTPUTDIFFS_NUMBER][SBOX_NUMBER];

#define PRINT_PDTTABLE 1
#define PRINT_PROB 1
#define PRINT_STATISTICS 1
#define PRINT_STORAGE 1
#define PRINT_PDTSTATISTICS 1

static const u8 Wti[64] = {0x0, //1
	0x1, 0x2, 0x4, 0x8, 0x10, 0x20, //6
	0x3, 0x5, 0x6, 0x9, 0xa, 0xc, 0x11, 0x12, 0x14, 0x18, 0x21, 0x22, 0x24, 0x28, 0x30, //15
	0x7, 0xb, 0xd, 0xe, 0x13, 0x15, 0x16, 0x19, 0x1a, 0x1c, 0x23, 0x25, 0x26, 0x29, 0x2a, 0x2c, 0x31, 0x32, 0x34, 0x38,//20
	0xf, 0x17, 0x1b, 0x1d, 0x1e, 0x27, 0x2b, 0x2d, 0x2e, 0x33, 0x35, 0x36, 0x39, 0x3a, 0x3c,//15
	0x1f, 0x2f, 0x37, 0x3b, 0x3d, 0x3e,//6
	0x3f//1
};
static const u8 Wto[16] = {0x0, 0x1, 0x2, 0x4, 0x8, 0x3, 0x5, 0x6, 0x9, 0xa, 0xc, 0x7, 0xb, 0xd, 0xe, 0xf};

/*static const u8 WtiForTravel[SBOX_INPUTS_NUMBER-1] = {
	0x1, 0x2, 0x4, 0x8, 0x10, 0x20, //6
	0x3, 0x5, 0x6, 0x9, 0xa, 0xc, 0x11, 0x12, 0x14, 0x18, 0x21, 0x22, 0x24, 0x28, 0x30, //15
	0x7, 0xb, 0xd, 0xe, 0x13, 0x15, 0x16, 0x19, 0x1a, 0x1c, 0x23, 0x25, 0x26, 0x29, 0x2a, 0x2c, 0x31, 0x32, 0x34, 0x38,//20
	0xf, 0x17, 0x1b, 0x1d, 0x1e, 0x27, 0x2b, 0x2d, 0x2e, 0x33, 0x35, 0x36, 0x39, 0x3a, 0x3c,//15
	0x1f, 0x2f, 0x37, 0x3b, 0x3d, 0x3e,//6
	0x3f//1
};*/
u8 WtiForTravel[SBOX_INPUTS_NUMBER-1] = {
	0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7, 0x8, 0x9 ,0xa, 0xb, 0xc, 0xd, 0xe, 0xf,
	0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19 ,0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
	0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29 ,0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
	0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39 ,0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f};

si8 WtoForTravelNumber[SBOX_INPUTS_NUMBER];
prType WtoForTravelProb[SBOX_INPUTS_NUMBER][MAX_OUTPUTDIFFS_NUMBER];

int cmp(int x,int y)                    //排序是依据a数组进行的  
{  
    return PDT_MaxProb[x]<PDT_MaxProb[y];  
}

//替换函数S
u16 Sbox[4][16] ={
14,5,7,2,11,8,1,15,0,10,9,4,6,13,12,3,
5,0,8,15,14,3,2,12,11,7,6,9,13,4,1,10,
4,9,2,14,8,7,13,0,10,12,15,1,5,11,3,6,
9,6,15,5,3,8,4,11,7,1,12,2,0,14,10,13
};
u16 Substitution(u16 input){
	u16 x=input,x1,x2;
	x1=(x&0x1)|((x>>4)&0x2);
	x2=(x>>1)&0xf;
	return Sbox[x1][x2];
}

prType Prob[8];

//生成float型概率势数数组，概率从大到小
void GenProb(){
#if (PRINT_PROB)
	FILE *fp=fopen("Prob.txt","w");
#endif
	for(int i=0;i<=7;i++){
		Prob[i]=6.0-log((prType)((8-i)*2))/log(2.0);
#if (PRINT_PROB)
		fprintf(fp,"%f\t",Prob[i]);
#endif
	}
#if (PRINT_PROB)
	fclose(fp);
#endif
}

//生成概率分布表
void GenPDT(){
	for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
		for(u16 j=0x0;j<SBOX_OUTPUTS_NUMBER;j++){
			PDT[i][j]=0x0;
		}
	}//全部初始化为0

	for(u16 id=0x0;id<SBOX_INPUTS_NUMBER;id++){
		for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
			u16 od=Substitution(i)^Substitution(i^id);
			PDT[id][od]++;
		}
	}

#if (PRINT_PDTTABLE)
	FILE *fp=fopen("PDT_Table.txt","w");
	for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
		fprintf(fp,"/*0x%02x*/{",i);
		for(u16 j=0x0;j<SBOX_OUTPUTS_NUMBER;j++){
			fprintf(fp,"%d,",PDT[i][j]);
		}
		fprintf(fp,"},\n");
	}
	fclose(fp);
#endif
}

//生成PDT统计信息，Storage()里用
void Statistics(){
	PDT_Number=0;
	memset(PDT_1_Number,0,PR_NUMBER*sizeof(int));
	memset(PDT_0_Number,0,PR_NUMBER*SBOX_INPUTS_NUMBER*sizeof(int));
	//初始化为0

	for(int i=0;i<SBOX_INPUTS_NUMBER;i++){
		for(int o=0;o<SBOX_OUTPUTS_NUMBER;o++){
			if(PDT[i][o]!=0&&PDT[i][o]!=64){
				PDT_Number++;
				PDT_1_Number[(8-PDT[i][o]/2)]++;
				PDT_0_Number[i][(8-PDT[i][o]/2)]++;
			}
		}
		PDT_0_Offset[i][0][0]=0;
		PDT_0_Offset[i][0][1]=PDT_0_Number[i][0];
		for(int p=1;p<PR_NUMBER;p++){
			PDT_0_Offset[i][p][0]=PDT_0_Offset[i][p-1][1];
			PDT_0_Offset[i][p][1]=PDT_0_Offset[i][p-1][1]+PDT_0_Number[i][p];
		}
	}
	PDT_1_Offset[0][0]=0;
	PDT_1_Offset[0][1]=PDT_1_Number[0];
	for(int p=1;p<PR_NUMBER;p++){
		PDT_1_Offset[p][0]=PDT_1_Offset[p-1][1];
		PDT_1_Offset[p][1]=PDT_1_Offset[p-1][1]+PDT_1_Number[p];
	}

#if (PRINT_STATISTICS)
	FILE *fp=fopen("Statistics_0_Table.txt","w");
	FILE *fq=fopen("Statistics_1_Table.txt","w");
	for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
		fprintf(fp,"*/0x%02x*/",i);
		for(u16 p=0x0;p<PR_NUMBER;p++){
			fprintf(fp,"%d:%d-(%d,%d)\t",p,PDT_0_Number[i][p],PDT_0_Offset[i][p][0],PDT_0_Offset[i][p][1]);
		}
		fprintf(fp,"\n");
	}
	for(u16 p=0x0;p<PR_NUMBER;p++){
		fprintf(fq,"%d:%d-(%d,%d)\t",p,PDT_1_Number[p],PDT_1_Offset[p][0],PDT_1_Offset[p][1]);
	}
	fprintf(fq,"\n");
	fclose(fp);
	fclose(fq);
#endif

}

void Storage(){
	int PDT_1_Count[PR_NUMBER];
	for(int p=0;p<PR_NUMBER;p++){
		PDT_1_Count[p]=PDT_1_Offset[p][0];
	}
	for(int i=0;i<SBOX_INPUTS_NUMBER;i++){
		int PDT_0_Count[PR_NUMBER];
		for(int p=0;p<PR_NUMBER;p++){
			PDT_0_Count[p]=PDT_0_Offset[Wti[i]][p][0];
		}
		for(int o=0;o<SBOX_OUTPUTS_NUMBER;o++){
			if(PDT[Wti[i]][Wto[o]]!=0&&PDT[Wti[i]][Wto[o]]!=64){//按照汉明重量从小到大存。真的有用吗？
				for(int si=0;si<SBOX_NUMBER;si++){
					memcpy(SPE[si][Wti[i]][PDT_0_Count[8-PDT[Wti[i]][Wto[o]]/2]],PE_Table[si][Wto[o]],SBOX_NUMBER*sizeof(u16));
				}
				PDT_1_Count[8-PDT[Wti[i]][Wto[o]]/2]++;
				PDT_0_Count[8-PDT[Wti[i]][Wto[o]]/2]++;
			}
		}
	}

#if (PRINT_STORAGE)
	FILE *fp=fopen("Storage_Table.txt","w");
	si8 si=7;
	for(u16 i=0x0;i<SBOX_INPUTS_NUMBER;i++){
		fprintf(fp,"*/0x%02x*/\t",i);
		for(u16 p=0x0;p<PR_NUMBER;p++){
			fprintf(fp,"%d:\t",p);
			for(u16 k=PDT_0_Offset[i][p][0];k<PDT_0_Offset[i][p][1];k++){
				fprintf(fp,"(");
				for(u16 l=0;l<8;l++){
					fprintf(fp,"%02x,",SPE[si][i][k][l]);
				}
				fprintf(fp,")\t");
			}
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
#endif
}

void PDTStatistics(){
	for(int p=0;p<PR_NUMBER;p++){
		PDT_1_Non0Num[p]=0;
	}
	for(int p=0;p<PR_NUMBER;p++){
		for(u16 i=0;i<SBOX_INPUTS_NUMBER;i++){
			if(PDT_0_Number[i][p]!=0){
				PDT_1_Non0Val[p][PDT_1_Non0Num[p]]=i;
				PDT_1_Non0Num[p]++;
			}
		}
	}

	PDT_MaxProb[0]=0;
	for(int p=PR_NUMBER-1;p>=0;p--){
		for(u16 i=0;i<PDT_1_Non0Num[p];i++){
			PDT_MaxProb[PDT_1_Non0Val[p][i]]=Prob[p];
		}
	}
	sort(WtiForTravel,WtiForTravel+SBOX_INPUTS_NUMBER-1,cmp);
	
	for(u16 i=1;i<SBOX_INPUTS_NUMBER;i++){
		WtoForTravelNumber[i]=PDT_0_Offset[i][PR_NUMBER-1][1];
		for(u16 p=0x0;p<PR_NUMBER;p++){
			for(u16 j=PDT_0_Offset[i][p][0];j<PDT_0_Offset[i][p][1];j++){
				WtoForTravelProb[i][j]=Prob[p];
			}
		}
	}

#if (PRINT_PDTSTATISTICS)
	FILE *fp=fopen("PDTStatistics_Table.txt","w");
	for(int p=0;p<PR_NUMBER;p++){
		fprintf(fp,"%d:%d\n",p,PDT_1_Non0Num[p]);
		for(u16 i=0;i<PDT_1_Non0Num[p];i++){
			fprintf(fp,"%02x\t",PDT_1_Non0Val[p][i]);
		}
		fprintf(fp,"\n");
	}

	fprintf(fp,"\n");
	for(u16 i=0;i<SBOX_INPUTS_NUMBER-1;i++){
		fprintf(fp,"%02x:%f\n",WtiForTravel[i],PDT_MaxProb[WtiForTravel[i]]);
	}
	fprintf(fp,"\n");
	for(u16 i=0;i<SBOX_INPUTS_NUMBER;i++){
		fprintf(fp,"%02x:%f\n",i,PDT_MaxProb[i]);
	}
	fprintf(fp,"\n");
	for(u16 i=1;i<SBOX_INPUTS_NUMBER;i++){
		fprintf(fp,"%02x:%d---",i,WtoForTravelNumber[i]);
		for(u16 j=0;j<WtoForTravelNumber[i];j++){
			fprintf(fp,"%f ",WtoForTravelProb[i][j]);
		}fprintf(fp,"\n");
	}
	fclose(fp);
#endif
}

void GenPrTable(){
	GenProb();
	GenPDT();
	Statistics();
	Storage();
	PDTStatistics();
}