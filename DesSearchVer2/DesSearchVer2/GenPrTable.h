#pragma once
#include "types.h"
#include "GenPermTable.h"

#define MAX_OUTPUTDIFFS_NUMBER (16)

extern int PDT[SBOX_INPUTS_NUMBER][SBOX_OUTPUTS_NUMBER];

extern int PDT_Number;
extern int PDT_0_Number[SBOX_INPUTS_NUMBER][PR_NUMBER];
extern int PDT_0_Offset[SBOX_INPUTS_NUMBER][PR_NUMBER][2];
extern int PDT_1_Number[PR_NUMBER];
extern int PDT_1_Offset[PR_NUMBER][2];

extern int PDT_1_Non0Num[PR_NUMBER];
extern int PDT_1_Non0Val[PR_NUMBER][SBOX_INPUTS_NUMBER];

extern prType PDT_MaxProb[SBOX_INPUTS_NUMBER];

extern ALIGNED_TYPE_(u16,8) SPE[SBOX_NUMBER][SBOX_INPUTS_NUMBER][MAX_OUTPUTDIFFS_NUMBER][SBOX_NUMBER];

extern const u8 Wti[64];
extern const u8 Wto[16];
extern u8 WtiForTravel[SBOX_INPUTS_NUMBER-1];

extern si8 WtoForTravelNumber[SBOX_INPUTS_NUMBER];
extern prType WtoForTravelProb[SBOX_INPUTS_NUMBER][MAX_OUTPUTDIFFS_NUMBER];

extern u16 *sbox;
void Substitution(u16* output,u16 input);

extern prType Prob[8];

void GenProb();
void GenPDT();
void Statistics();
void Storage();
void PDTStatistics();
void GenPrTable();