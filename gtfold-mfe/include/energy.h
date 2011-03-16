#ifndef _ENERGY_TABLES_H_
#define _ENERGY_TABLES_H_

extern int *V; 
extern int *W; 
extern int **VBI; 
extern int **VM; 
extern int **WM; 
extern int *indx; 


#define V(i,j) V[indx[i]+j]
//#define VM(i,j) VM[indx[i]+j]
//#define WM(i,j) WM[indx[i]+j]
//#define VBI(i,j) VBI[indx[i]+j]


#define auPen(i, j) ((( (i)==BASE_U || (j)==BASE_U ) && ( (i)==BASE_A || (i)==BASE_G || (j)==BASE_A || (j)==BASE_G )) ? auend : 0)

#ifdef __cplusplus
extern "C" {
#endif
int Ed3(int i, int j, int k);
int Ed5(int i, int j, int k);
int auPenalty(int i, int j);

int Ec();
int Eb();
int Ea(); 

int eS(int i, int j);
int eH(int i, int j);
int eL(int i, int j, int ip, int jp);

void create_tables(int len);
void init_tables(int len);
void free_tables(int len);
#ifdef __cplusplus
}
#endif

#endif
