#include "search.h"

si8 a[2][SBOX_NUMBER+1];//a[*][0]=-1

inline void ResetCharacter(int k,int l,int round){
	for(int i=k+1;i<l;i++){
		r_od_l[round][i]=0;
	}
}
void setAndPrint(){
	Bnc[Round-1]=r_pr[Round-1];
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_l[0][si]);
	}fprintf(fp_trails,"\t%f\n",r_pr[0]);

	/*for(int i=0;i<=SBOX_NUMBER;i++){
		fprintf(fp_trails,"%d ",a[0][i]);
	}fprintf(fp_trails,"\n");*/

	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_l[1][si]);
	}fprintf(fp_trails,"\t");
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_r[1][si]);
	}fprintf(fp_trails,"\t%f\n",r_pr[1]);

	/*for(int i=0;i<=SBOX_NUMBER;i++){
		fprintf(fp_trails,"%d ",a[1][i]);
	}fprintf(fp_trails,"\n");*/
	
	for(int r=2;r<Round-1;r++){
		for(int si=0;si<SBOX_NUMBER;si++){
			fprintf(fp_trails,"%02x ",r_od_l[r][si]);
		}fprintf(fp_trails,"\t");
		for(int si=0;si<SBOX_NUMBER;si++){
			fprintf(fp_trails,"%02x ",r_od_r[r][si]);
		}
		fprintf(fp_trails,"\t%f",r_pr[r]);
		fprintf(fp_trails,"\t%d: ",r_an[r]);
		for(int si=0;si<r_an[r];si++){
			fprintf(fp_trails,"%d ",r_ai[r][si]);
		}fprintf(fp_trails,"\n");
	}
	
	for(int si=0;si<SBOX_NUMBER;si++){
		fprintf(fp_trails,"%02x ",r_od_l[Round-1][si]);
	}
	fprintf(fp_trails,"\t%f",r_pr[Round-1]);
	fprintf(fp_trails,"\t%d: ",r_an[Round-1]);
	for(int si=0;si<r_an[Round-1];si++){
		fprintf(fp_trails,"%d ",r_ai[Round-1][si]);
	}fprintf(fp_trails,"\n");
	
	fprintf(fp_trails,"---------------\n");
}

void getInfo(int r,__m128i tmp){
	__m128i a;
	__m128i b;
	__m128i c;
	__m128i d;
	int lm;
	int la;

	//获得r_an[r]
	a=_mm_setzero_si128();
	a=_mm_cmpgt_epi16(tmp,a);//tmp非零的16比特置为0xffff
	lm=_mm_movemask_epi8(a);
	la=_mm_popcnt_u32(lm)/2;
	r_an[r]=la;

	//获得r_ai[r][]
	c=_mm_load_si128((__m128i *)(&(W16v[lm][0])));
	_mm_storeu_si128((__m128i *)(r_ai[r]),c);
}

void Round__(int r,int j,prType pr_round,__m128i tmp0){
	si8 an=r_an[r-1];
	si16 ai=r_ai[r-1][j];//第j个活跃S盒下标
	si8 an_remain=an-j-1;
	prType prob;
	si8 s;
	si8 m;
	u16 idv=r_od_l[r-1][ai];
	__m128i tmp1;
	__m128i *odp;
	odp=(__m128i *)(r_od_r[r-1]);

	for(si16 pr=0;pr<PR_NUMBER;pr++){//遍历概率
		prob=Prob[pr]+pr_round;
		if((prob+r_pr[r-2]+Bn[Round-r-1]+WMIN_S*an_remain)<(Bnc[Round-1]+1e-10)){
			s=PDT_0_Offset[idv][pr][0];
			m=PDT_0_Number[idv][pr];
			for(si16 k=s;k<s+m;k++){//遍历输出差分
				tmp1=_mm_xor_si128(tmp0,*(__m128i *)(SPE[ai][idv][k]));
				if(an_remain==0){
					r_pr[r-1]=r_pr[r-2]+prob;
					_mm_store_si128(odp,tmp1);
					Round_(r+1);
				}else{
					Round__(r,j+1,prob,tmp1);
				}
			}
		}else break;
	}
}

void Round_N_(int j,prType pr_round){
	si8 an=r_an[Round-1];
	si16 ai=r_ai[Round-1][j];//第j个活跃S盒下标
	si8 an_remain=an-j-1;
	prType prob;
	u16 idv=r_od_l[Round-1][ai];

	prob=PDT_MaxProb[idv]+pr_round;
	if((prob+r_pr[Round-2]+WMIN_S*an_remain)<(Bnc[Round-1]+1e-10)){
		if(an_remain==0){
			r_pr[Round-1]=r_pr[Round-2]+prob;
			setAndPrint();
		}else{
			Round_N_(j+1,prob);
		}
	}
}

void Round_(int r){
	
	__m128i *dp1;
	__m128i *dp2;
	__m128i *idp;
	__m128i *odp;
	dp1=(__m128i *)(r_od_l[r-3]);
	dp2=(__m128i *)(r_od_r[r-2]);
	idp=(__m128i *)(r_od_l[r-1]);
	odp=(__m128i *)(r_od_r[r-1]);

	__m128i tmp0;
	__m128i tmp1;
	tmp0=_mm_xor_si128(*dp1,*dp2);
	tmp1=_mm_setzero_si128();
	_mm_store_si128(idp,tmp0);
	getInfo(r-1,*idp);
	
	if(r_an[r-1]==0){
		r_pr[r-1]=r_pr[r-2];
		_mm_store_si128(odp,tmp1);
		if(r==Round){setAndPrint();}
		else{Round_(r+1);}
	}else{
		if((r_pr[r-2]+Bn[Round-r-1]+WMIN_S*r_an[r-1])<(Bnc[Round-1]+1e-10)){}
		else return;
	}

	if(r==Round){
		Round_N_(0,0);//1是起始活跃S盒，0是起始概率
	}else{
		Round__(r,0,0,tmp1);//r是轮数，0是起始活跃S盒，0是起始概率，tmp1是起始输出差分
	}
}

//j是当前处理第j个活跃S盒
//pr_round是当前轮在调用本次遍历之前的概率
//tmp0是当前轮在调用本次遍历之前的输出差分
void Round_2_(int j,prType pr_round,__m128i tmp0){
	prType prob;
	u16 dx;
	__m128i *odp;
	odp=(__m128i *)(r_od_r[1]);
	__m128i tmp1;
	si8 s;
	si8 m;
	u16 idv;
	si8 idv_num;
	
	for(a[1][j]=a[1][j-1]+1;a[1][j]<8;a[1][j]++){
		if(j!=1 && (r_od_l[1][a[1][j-1]]&0x3)!=0){
			if(a[1][j]!=(a[1][j-1]+1)){
				break;
			}
		}
		memset(r_od_l[1]+a[1][j-1]+1,0,sizeof(u16)*(a[1][j]-a[1][j-1]-1));
//******************************
//将第a[2][j-1]个S盒至第a[2][j]个S盒之间的S盒输入输出全部置为0。
//a[2][j-1]后两个比特不为0，那么a[2][j]只能取a[2][j-1]+1；但若j=1，不存在这个问题。
//******************************
		
//------------------------------
//a[2][j]==8
//------------------------------
		if(a[1][j]==7){
			if(j!=1||firstRoundActive==1){
//******************************
//当第一轮活跃时，第二轮才能不活跃。
//j不为1时，说明第二轮必不活跃；activeflag==1时，说明第一轮是活跃的。
//******************************

				if( 0==(r_od_l[1][6]&0x3) && 0==(r_od_l[1][0]&0x30) ){
					r_od_l[1][7]=0;
					_mm_store_si128(odp,tmp0);
					r_pr[1]=pr_round+r_pr[0];//最后一个S盒weight为0，肯定通过剪枝条件
					Round_(3);
					//setAndPrint();
				}
			}
//******************************
//a[2][j]==8且dx[2][8]==0时，dx[1][7]的后两个比特和dx[1][1]的前两个比特必须为0。
//累加概率，剪枝，因a[2][j]==8，通过剪枝则进入下一轮。
//******************************

			for(si16 pr=0;pr<PR_NUMBER;pr++){//遍历概率
				prob=Prob[pr]+pr_round;
				if((prob+r_pr[0]+Bn[Round-3])<(Bnc[Round-1]+1e-10)){
					idv_num=PDT_1_Non0Num[pr];
					for(si16 i=0;i<idv_num;i++){//遍历输入差分
						idv=PDT_1_Non0Val[pr][i];
						if((idv&0x30)==((r_od_l[1][6]&0x3)<<4) && (idv&0x3)==((r_od_l[1][0]&0x30)>>4)){
							r_od_l[1][7]=idv;
							s=PDT_0_Offset[idv][pr][0];
							m=PDT_0_Number[idv][pr];
							for(si16 k=s;k<s+m;k++){//遍历输出差分
								tmp1=_mm_xor_si128(tmp0,*(__m128i *)(SPE[7][idv][k]));
								r_pr[1]=r_pr[0]+prob;
								_mm_store_si128(odp,tmp1);
								Round_(3);
								//setAndPrint();
							}
						}
					}
				}else break;
			}
//******************************
//a[2][j]==8且dx[2][8]!=0时，遍历dx[1][8]的可能，其实这里dx[2][8]的自由度只有两个比特。
//按照概率从大到小遍历所有输入、输出差分值对。一旦发生剪枝，则结束遍历，因为之后的概率一定更小。
//累加概率，剪枝，因a[2][j]==8，通过剪枝则进入下一轮。
//******************************

//------------------------------
//a[2][j]==1
//------------------------------
		}else if(a[1][j]==0){
			for(si16 pr=0;pr<PR_NUMBER;pr++){//遍历概率
				prob=Prob[pr]+pr_round;
				if((prob+r_pr[0]+Bn[Round-3])<(Bnc[Round-1]+1e-10)){
					idv_num=PDT_1_Non0Num[pr];
					for(si16 i=0;i<idv_num;i++){//遍历输入差分
						idv=PDT_1_Non0Val[pr][i];
						r_od_l[1][0]=idv;
						s=PDT_0_Offset[idv][pr][0];
						m=PDT_0_Number[idv][pr];
						for(si16 k=s;k<s+m;k++){//遍历输出差分
							tmp1=_mm_xor_si128(tmp0,*(__m128i *)(SPE[0][idv][k]));
							//_mm_store_si128(odp,tmp1);
							Round_2_(j+1,prob,tmp1);
						}
					}
				}else break;
			}
//******************************
//a[2][j]==1时，此时也有j=1，遍历dx[2][1]的非零可能，这里dx[2][1]的自由度有六个比特。
//注：第1个S盒自由度6比特，第2~7个S盒自由度4比特，第8个S盒自由度2比特。
//按照概率从大到小遍历所有输入、输出差分值对。一旦发生剪枝，则结束遍历，因为之后的概率一定更小。
//累加概率，剪枝，因a[2][j]!=8，通过剪枝则进入下一个S盒。
//******************************

//------------------------------
//a[2][j]==2~7
//------------------------------
		}else{
			for(si16 pr=0;pr<PR_NUMBER;pr++){//遍历概率
				prob=Prob[pr]+pr_round;
				if((prob+r_pr[0]+Bn[Round-3])<(Bnc[Round-1]+1e-10)){
					idv_num=PDT_1_Non0Num[pr];
					for(si16 i=0;i<idv_num;i++){//遍历输入差分
						idv=PDT_1_Non0Val[pr][i];
						if( (idv&0x30) == ((r_od_l[1][a[1][j]-1]&0x3)<<4) ){
							r_od_l[1][a[1][j]]=idv;
							s=PDT_0_Offset[idv][pr][0];
							m=PDT_0_Number[idv][pr];
							for(si16 k=s;k<s+m;k++){//遍历输出差分
								tmp1=_mm_xor_si128(tmp0,*(__m128i *)(SPE[a[1][j]][idv][k]));
								//_mm_store_si128(odp,tmp1);
								Round_2_(j+1,prob,tmp1);
							}
						}
					}
				}else break;
			}
		}
//******************************
//a[2][j]在2至7之间时，遍历dx[2][a[2][j]]的非零可能，这里dx[2][a[2][j]]的自由度有四个比特。
//按照概率从大到小遍历所有输入、输出差分值对。一旦发生剪枝，则结束遍历，因为之后的概率一定更小。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************
	}
}

//------------------------------
//第二轮搜索
//------------------------------
void Round_2(){
	__m128i tmp;
	tmp=_mm_setzero_si128();
	Round_2_(1,0,tmp);
}

void Round_1_(int j,prType pr_round){
	prType prob;
	u16 dx;
	for(a[0][j]=a[0][j-1]+1;a[0][j]<8;a[0][j]++){

		if(j!=1 && (r_od_l[0][a[0][j-1]]&0x3)!=0){
			if(a[0][j]!=(a[0][j-1]+1)){
				break;
			}
		}
		memset(r_od_l[0]+a[0][j-1]+1,0,sizeof(u16)*(a[0][j]-a[0][j-1]-1));
		//ResetCharacter(a[0][j-1],a[0][j],0);
//******************************
//将第a[1][j-1]个S盒至第a[1][j]个S盒之间的S盒输入输出全部置为0。
//a[1][j-1]后两个比特不为0，那么a[1][j]只能取a[1][j-1]+1；但若j==1，不存在这个问题。
//注意j!=1时a[1][j]!=1。
//******************************

//------------------------------
//a[1][j]==8
//------------------------------
		if(a[0][j]==7){
			r_od_l[0][7]=0;
			if(j==1){
				firstRoundActive=0;
				r_pr[0]=0;
				Round_2();
			}else{
				if( 0==(r_od_l[0][6]&0x3) && 0==(r_od_l[0][0]&0x30) ){
					r_pr[0]=pr_round;
					Round_2();
				}
			}
//******************************
//j==1且a[1][j]==8且dx[1][8]==0时，第一轮差分为0，activeflag置为0，这时第二轮必须活跃
//a[1][j]==8且dx[1][8]==0时，dx[1][7]的后两个比特和dx[1][1]的前两个比特必须为0。
//这里j不必为1，前面j==1用于判断第一轮差分是否为0。
//累加概率，剪枝，因a[1][j]==8，通过剪枝则进入下一轮。
//******************************
			
			firstRoundActive=1;
			for(si16 i=0;i<SBOX_INPUTS_NUMBER-1;i++){
				dx=WtiForTravel[i];
 				r_od_l[0][7]=dx;
 				if( (dx&0x30)==((r_od_l[0][6]&0x3)<<4) && (dx&0x3)==((r_od_l[0][0]&0x30)>>4) ){
					prob=PDT_MaxProb[dx]+pr_round;
					if((prob+Bn[Round-2])<(Bnc[Round-1]+1e-10)){
						r_pr[0]=prob;
						Round_2();
					}else break;
				}
			}
//******************************
//a[1][j]==8且dx[1][8]!=0时，遍历dx[1][8]的可能，其实这里dx[1][8]的自由度只有两个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，因a[1][j]==8，通过剪枝则进入下一轮。
//******************************

//------------------------------
//a[1][j]==1
//------------------------------
		}else if(a[0][j]==0){
			for(si16 i=0;i<SBOX_INPUTS_NUMBER-1;i++){
				dx=WtiForTravel[i];
				r_od_l[0][0]=dx;
				prob=PDT_MaxProb[dx]+pr_round;
				if((prob+Bn[Round-2])<(Bnc[Round-1]+1e-10)){
					Round_1_(j+1,prob);
				}else break;
			}

//******************************
//a[1][j]==1时，此时也有j=1，遍历dx[1][1]的非零可能，这里dx[1][1]的自由度有六个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************

//------------------------------
//a[1][j]==2~7
//------------------------------
		}else{
			for(si16 i=0;i<SBOX_INPUTS_NUMBER-1;i++){
				dx=WtiForTravel[i];
 				r_od_l[0][a[0][j]]=dx;
				if( (dx&0x30) == ((r_od_l[0][a[0][j]-1]&0x3)<<4) ){
					prob=PDT_MaxProb[dx]+pr_round;
					if((prob+Bn[Round-2])<(Bnc[Round-1]+1e-10)){
						Round_1_(j+1,prob);
					}else break;
				}
			}
		}
//******************************
//a[1][j]在2至7之间时，遍历dx[1][a[1][j]]的非零可能，这里dx[1][a[1][j]]的自由度有四个比特。
//固定了输入差分，输出差分取概率最大的所有值。
//累加概率，剪枝，通过剪枝则进入下一个S盒。
//******************************

	}
}

void Round_1(){
	fn_trails_str=to_string((_ULonglong)Round)+"_trails.txt";
	errno_t err;
	err = fopen_s(&fp_trails, fn_trails_str.c_str(), "w" );
	a[0][0]=-1;
	a[1][0]=-1;
	//__m128i tmp0;
	//tmp0=_mm_setzero_si128();
	Round_1_(1,0);

	/*__m128i *idp;
	idp=(__m128i *)(r_od_l[0]);
	_mm_store_si128(idp,tmp0);
	r_pr[0]=0;
	firstRoundActive=0;
	Round_2();*/

	fclose(fp_trails);
}