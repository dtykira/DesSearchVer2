#pragma once
#include "types.h"

extern bool firstRoundActive;

extern ALIGNED_TYPE_(u16,8) r_od_l[ROUND_N][SBOX_NUMBER];//存储轮差分
extern ALIGNED_TYPE_(u16,8) r_od_r[ROUND_N][SBOX_NUMBER];
extern si8 r_an[ROUND_N];//存储r_od的活跃S盒个数
extern si16 r_ai[ROUND_N][SBOX_NUMBER];//存储r_od的活跃S盒的下标
extern prType r_pr[ROUND_N];

extern ALIGNED_TYPE_(u16,8) r_od_best[ROUND_N][SBOX_NUMBER];
extern si8 r_an_best[ROUND_N];
extern si8 r_ai_best[ROUND_N][SBOX_NUMBER];
extern prType r_pr_best[ROUND_N];

extern prType Bnc[ROUND_N];//Bnc[R]保存R+1轮迹当前候选weight
extern prType Bn[ROUND_N];//保存最佳weight

extern int Round;//当前搜索的轮数

extern ALIGNED_TYPE_(si16,8) W16v[256*256][8];

extern FILE *fp_trails;
extern string fn_trails_str;