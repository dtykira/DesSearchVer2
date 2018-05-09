#include "types.h"
#include "GenPermTable.h"
#include "GenPrTable.h"
#include "search.h"

int main(){
	Gen_PE_Table();
	GenPrTable();
	clock_t start,end;

	Round=3;
	Bn[0]=0;
	Bn[1]=2;
	Bn[2]=4;
	Bn[3]=9.022719;
	Bn[4]=14.045439;
	Bn[5]=22.186218;
	Bn[6]=29.249580;
	Bn[7]=37.594231;
	Bn[8]=41.764153;
	Bn[9]=49.786873;
	Bn[10]=52.705261;
	Bn[11]=58.102371;
	Bn[12]=60.780441;
	Bn[13]=68.780441;

	Bnc[Round-1]=8;
	start=clock();
	Round_1();
	end = clock();
	printf("%d&%f\n",Round,(double)(end-start)/CLK_TCK);
	/*for(Round=4;Round<=14;Round++){
		Bnc[Round-1]=(int)Bn[Round-1]+1;

		start=clock();
		Round_1();
		end = clock();
		printf("%d&%f\n",Round,(double)(end-start)/CLK_TCK);
	}*/
	system("pause");
	return 0;
}