#pragma once
#include "types.h"

void Expansion(u16* output, u32 input);
u32 Permutation(u32 x);
u32 SboxOutput2word(u16* input);

void PE(u16* output, u16* input);