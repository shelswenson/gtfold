//Attempting to follow Ding and Lawrence as closely as possible

/**s1 -> partial_external
*  s2 -> partial_multi
*  s3 -> partial_multi2
*  u1 -> u_multi
*  f -> cond_dangle (conditional)
*  ebi -> eL
**/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "energy.h"
#include "utils.h"
#include "global.h"
#include "data.h"
#include "algorithms-partition.h"
#include "partition-dangle.h"

#define MIN_TURN 4
#define MAX_LOOP 30

//double RT = (0.00198721 * 310.15)/100.00; 

double cond_dangle(int j, int h, int l){
	if(j - 1 == l)
		return 1;
	else 
		return exp(-Ed3(h,l,l+1)/RT);
}

dangle_struct malloc_partition_arrays_d(int length){
	dangle_struct ret_struct;
	ret_struct.length = length;
	
	//Add 1 since we're not zero indexed
	ret_struct.u = mallocTwoD(length + 1, length + 1); 
	ret_struct.up = mallocTwoD(length + 1, length + 1); 
	ret_struct.upm = mallocTwoD(length + 1, length + 1); 
	ret_struct.u1 = mallocTwoD(length + 1, length + 1); 
	ret_struct.s1 = mallocTwoD(length + 1, length + 1); 
	ret_struct.s2 = mallocTwoD(length + 1, length + 1); 
	ret_struct.s3 = mallocTwoD(length + 1, length + 1); 

	int i,j;

	for(i = 0; i < length + 1; i++){
		for(j = 0; j < length + 1; j++){
			ret_struct.u[i][j] = 1;
			ret_struct.up[i][j] = 0;
			ret_struct.upm[i][j] = 0;
			ret_struct.u1[i][j] = 0;
			ret_struct.s1[i][j] = 0;
			ret_struct.s2[i][j] = 0;
			ret_struct.s3[i][j] = 0;
		}
	}
	return ret_struct;
}

void free_partition_arrays_d(struct partition_d part)
{
	freeTwoD(part.u, part.length, part.length);
	freeTwoD(part.up, part.length, part.length);
	freeTwoD(part.upm, part.length, part.length);
	freeTwoD(part.u1, part.length, part.length);
	freeTwoD(part.s1, part.length, part.length);
	freeTwoD(part.s2, part.length, part.length);
	freeTwoD(part.s3, part.length, part.length);
}

void fill_partition_arrays_d(dangle_struct part_struct){
	//So that we don't have the /100 in every exponent
	//Seeing as we missed some last time

	double ** partial_external = part_struct.s1;
	double ** partial_multi = part_struct.s2;
	double ** partial_multi2 = part_struct.s3;
	double ** u_multi = part_struct.u1;
	double ** u = part_struct.u;
	double ** up = part_struct.up;
	double ** upm = part_struct.upm;
	int len = part_struct.length;

	int seg_length;
	volatile int i,j,l;
	double second_half = 1; //used to calculate conditional terms

	for(seg_length = MIN_TURN; seg_length <= len; seg_length++){
		printf("Length %d\n", seg_length);
		//Insert parallelism here.
		for(i = 1; i < len - seg_length; i++){
			j = i + seg_length - 1;

			if(canPair(RNA[i],RNA[j])){
				for(l = i+2; l < j; l++){

					double temp_term = up[i+1][l] * 
								exp(-(Ea + 2 * Ec + auPenalty(i + 1,l)/RT));
					if(l + 2 < j)
					{
						temp_term *= (exp(-Ed3(i + 1,l,l+1)/RT) * 
									u_multi[l + 2][j - 1] + 
									u_multi[l + 1][j - 1] -
									u_multi[l + 2][j - 1]);
					}
					upm[i][j] += temp_term; 

					if(l != i+2){ //goes from i + 2 < l < j
						temp_term = up[i + 2][l] * 
								exp(-(Ea + 2 * Ec + Eb + Ed3(j,i,i + 1) + 
								auPenalty(i + 2, l))/RT);

						if(l + 2 < j - 1){
							temp_term *= (exp(-Ed3(i + 2,l,l + 1)/RT)*
									u_multi[l + 2][j - 1] + 
									u_multi[l + 1][j - 1]  - 
									u_multi[l + 2][j - 1]);
						}
					if(temp_term < -0.01){
						printf("Second temp_term is negative\n");
						exit(-1);
					}

						upm[i][j] += temp_term; 

						if(l != j - 1){
							//Changed h in the paper to l it makes no difference
							//other than unifying loops
							temp_term = exp(-Ed3(j,i,i + 1)/RT)*
									exp(-(Ea + 2 * Ec + (l - i - 1) * Eb)/RT) * 
									partial_multi[l][j];

					if(temp_term < -0.01){
						printf("Third temp_term is negative\n");
						exit(-1);
					}
							upm[i][j] += temp_term;  

						}
					}
				}
				if(upm[i][j]<0)
				printf("What the hell?! why is upm[%d][%d] negative? %f\n", i,j,upm[i][j]);
				


				up[i][j] += exp(-eH(i,j)/RT) + 
					exp(-eS(i,j)/RT) * up[i + 1][j - 1] + upm[i][j];

				for(l = j - 1; l > i + 1; l--){
					if((j - l) - 2> MAX_LOOP){
						break;
					}
					int h;
					for(h = i + 1 ; h < l; h++){
						if(!(i == h - 1 && j == l - 1)){ //If this is true we have a stack
							up[i][j] += exp(-eL(i,j,h,l)/RT);
						}
						if((j - l) + (h - i) - 2 > MAX_LOOP){
							break;
						}
					}
				}
				if(up[i][j] < 0){
					printf("What the hell?! why is up[%d][%d] negative? %f\n", i,j,up[i][j]);
				}
			}//End checkpair conditional

			for(l = i + 1; l <= j ; l++){
				double temp_term = up[i][l] * exp(-(Ec + auPenalty(i,l))/RT); 
				if(l + 2 < j)
				{
					temp_term *= (cond_dangle(j + 1,i,l) * 
							exp(-(j - l) * Eb/RT) + 
							exp(-Ed3(i,l,l + 1)/RT) * u_multi[l + 2][j] + 
							u_multi[l + 1][j] -
							u_multi[l + 2][j]);
				}

				u_multi[i][j] += temp_term;

				if(temp_term < 0){
					printf("First temp_term in u_multi[%d][%d] is negative: %f\n",
							i,j,temp_term);
					exit(-1);
				}

				if(l != i+1){ //this term runs from i + 2 to j	
					temp_term = up[i + 1][l]*
							exp(-(Ec + Eb + auPenalty(i + 1,l))/RT);
					if(l + 2 < j)
					{
						temp_term *= (cond_dangle(j + 1, i + 1, l) * 
								exp(-(j - l)*Eb/RT) + 
								exp(-Ed3(i + 1,l,l + 1)/RT) * u_multi[l + 2][j] + 
								u_multi[l + 1][j] - 
								u_multi[l + 2][j]);
					}
					u_multi[i][j] += temp_term;

					if(temp_term < 0){
						printf("Second temp_term in u_multi[%d][%d] is negative: %f\n",
						i,j,temp_term);
						exit(-1);
					}

					if(l != j){ //The third summation only runs from i+2 to j-1
						u_multi[i][j] += exp(-(Ec + (l - i) * Eb)/RT) * 
										partial_multi2[l][j];

					}
				}
			}
			if(u_multi[i][j] < 0){
				printf("u_multi[%d][%d] is %f\n", i, j, u_multi[i][j]);
				exit(-1);
			}


			int h = i; //To stay consistant with the notation in the paper.

			for(l = h+1; l < j; l++){
				double temp_ext, temp_mul, temp_mul2;
				temp_ext = up[h][l] * 
							exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT);

				temp_mul = up[h][l] * 
							exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT);

				temp_mul2 = up[h][l] * 
							exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT);

				if(l + 2 < j){
					temp_ext *= 
							(exp(-Ed3(h,l,l + 1)/RT)*u[l + 2][j] + 
								u[l + 1][j] -
								u[l + 2][j]);

					temp_mul *= 
							(exp(-Ed3(h,l,l + 1)/RT)*u_multi[l + 2][j - 1] + 
								u_multi[l + 1][j - 1] -
								u_multi[l + 2][j - 1]);

					temp_mul2 *=
							(cond_dangle(j + 1, h, l)*exp(-(j - l) * Eb / RT) + 
							exp(-Ed3(h,l,l + 1)/RT)*u_multi[l + 2][j] + 
								u_multi[l + 1][j] - 
								u_multi[l + 2][j]);
				}
				
				partial_external[h][j] += temp_ext;

				partial_multi[h][j] += temp_mul;

				partial_multi2[h][j] += temp_mul2;
							
			}
			/** partial_multi2 goes up to j **/
			//Since l = j here there's no need for a second_half
			partial_multi2[h][j] += up[h][l] * 
						exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT);

			//Finally we do the u matrix
			u[i][j] = 1 + up[i][j]*exp(-auPenalty(i,j) /RT);

			for(l = i + 1; l < j; l++){
				//Replaced l with h in the next line to unify loops
				u[i][j] += up[l][j] * 
						exp(-(Ed5(l, j , l - 1) + auPenalty(h,j))/RT);

				if(l + 2 < j){
					second_half = (exp(-Ed3(i,l,l + 1)/RT) * 
						u[l + 2][j] +
						u[l + 1][j] - 
						u[l + 2][j]);
				}
				u[i][j] += exp(-(auPenalty(i,l))/RT)*second_half;

				if(l != j - 1) {
					u[i][j] += partial_external[l][j];
				}
			}
		}//end of minor (paralellizeable) for loop
	}//End of major for loop
	printf("Total partition function: %f \n", u[1][len]);
}
