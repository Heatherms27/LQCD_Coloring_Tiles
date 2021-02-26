#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Helper functions -------------------------------------------

#define MAX(x, y) ((x) > (y)) ? (x) : (y)
#define MIN(x, y) ((x) < (y)) ? (x) : (y)

void index2coord(int ix, short d[], short c[], short n);
int coord2index(short c[], short d[], short n);
short posMod(short a, short n);

void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

void randshuffle(int *arr, int N) {
    int i;
    for(i = N-1; i > 0; i--) {
        int j = rand() % (i+1);
        swap(&arr[i], &arr[j]);
    }
}

void randperm(int *perm, int N) {
   int i;
   for(i = 0; i < N; ++i) perm[i] = i;

   // Random permutation the order
   randshuffle(perm, N);
}

void RBRand(int *perm, short pts[]) {
   int N = pts[0]*pts[1]*pts[2]*pts[3];
   short c[4];
   short dims = 4;

   int* odds = (int*)calloc(N, sizeof(int));

   int numOdds = 0;
   int numEvens = 0;
   int sumc, i;
   for (i = 0; i < N; ++i) {
      index2coord(i, pts, c, dims);
      sumc = c[0]+c[1]+c[2]+c[3];
      if(sumc % 2) 
         perm[numEvens++] = i;
      else
         odds[numOdds++] = i;
   }
   randshuffle(perm, numEvens);
   for(i = numEvens; i < N; i++)
      perm[i] = odds[i-numEvens];
   randshuffle(perm + numEvens, numOdds);
   
   free(odds);
//p = [evens(randperm(length(evens)));odds(randperm(length(odds)))];
}

void RB(int *perm, short pts[]) {
   int N = pts[0]*pts[1]*pts[2]*pts[3];
   short c[4];
   short dims = 4;

   int* odds = (int*)calloc(N, sizeof(int));

   int numOdds = 0;
   int numEvens = 0;
   int sumc, i;
   for (i = 0; i < N; ++i) {
      index2coord(i, pts, c, dims);
      sumc = c[0]+c[1]+c[2]+c[3];
      if(sumc % 2) 
         perm[numEvens++] = i;
      else
         odds[numOdds++] = i;
   }
   for(i = numEvens; i < N; i++)
      perm[i] = odds[i-numEvens];
   
   free(odds);
//p = [evens(randperm(length(evens)));odds(randperm(length(odds)))];
}

// Find maximum value in array. Assume N > 0
int findMax(short C[], int N){
	short max = C[0];
	int i;
	for(i = 1; i < N; ++i)
		if(C[i] > max)
			max = C[i];
	return max;
}

// a % n returns a negative number when a < 0, need it to return a postive number
short posMod(short a, short n){
	return (n + (a % n)) % n;
}

// ------------------------------------------------------------

// Convert lattice coordinate into array index
int coord2index(short c[], short d[], short n){

   int ix = c[0];
   int temp = d[0];
   short i;

   for(i = 1; i < n; ++i){
      ix += temp*c[i];
      temp = temp * d[i];
   }
   return ix;

}

// Convert array index into lattice coordinate
void index2coord(int ix, short d[], short c[], short n){
	short i;
	for(i = 0; i < n; ++i){
		c[i] = ix % d[i];
		ix = (int)(ix/d[i]);
	}
}

// Find the L2 norm of a coordinate point
double L2Norm(short *x){
   
   int i;
   double norm = 0.0;
   for(i = 0; i < 4; ++i){
      norm += (double)x[i] * (double)x[i];
   }

   return sqrt(norm);

}

short colorNode(int neighbors[], int nCount, short colors[], short colorArray[]){
   
   int i, j;
   for(i = 0; i < nCount; ++i)
      colorArray[colors[neighbors[i]]] = 1;

   for(i = 1; i < 10000; ++i)
      if(colorArray[i] == 0)
         break;
    
   for(j = 0; j < nCount; ++j)
      colorArray[colors[neighbors[j]]] = 0;     

   return i;

}

short addToNeighbors(int index, int neighbors[], int nCount){

   int i;
   for(i = 0; i < nCount; ++i)
      if(index == neighbors[i] || index == 0)
         return 0;

   return 1;

}

int createStencil(short stencil[][4], short p, short numk, short ks[], int maxN, short d[]){

   int nCount = 0;
   int index;
   short i, k, j1, j2, j3, j4, s1, s2, s3;
   int *neighbors = (int*)malloc(maxN * sizeof(int)); 
   
   for(i = 0; i < numk; ++i){
      k = ks[i];
      for(j1 = -p; j1 <= p; ++j1){
         for(j2 = -p; j2 <= p; ++j2){
            for(j3 = -p-k; j3 <= p-k; ++j3){
               for(j4 = -p; j4 <= p; ++j4){
                  short coord[4] = {j1, j2, j3, j4};
                  if(L2Norm(coord) <= p){
                     unsigned short coord2[4] = {posMod(j1, d[0]), posMod(j2, d[1]), posMod(j3, d[2]), posMod(j4, d[3])};
                     index = coord2index(coord2, d, 4);
                     if(addToNeighbors(index, neighbors, nCount)){
                        neighbors[nCount] = index;
                        stencil[nCount][0] = coord2[0];
                        stencil[nCount][1] = coord2[1];
                        stencil[nCount][2] = coord2[2];
                        stencil[nCount][3] = coord2[3];
                        ++nCount;
                     }  // End if           
                  }  // End if    
               }  // End j4
            }  // End j3
            for(j3 = -p+k; j3 <= p+k; ++j3){
               for(j4 = -p; j4 <= p; ++j4){
                  short coord[4] = {j1, j2, j3, j4};
                  if(L2Norm(coord) <= p){
                     unsigned short coord2[4] = {posMod(j1, d[0]), posMod(j2, d[1]), posMod(j3, d[2]), posMod(j4, d[3])};
                     index = coord2index(coord2, d, 4);
                     if(addToNeighbors(index, neighbors, nCount)){
                        neighbors[nCount] = index;
                        stencil[nCount][0] = coord2[0];
                        stencil[nCount][1] = coord2[1];
                        stencil[nCount][2] = coord2[2];
                        stencil[nCount][3] = coord2[3];
                        ++nCount;
                     }  // End if           
                  }  // End if    
               }  // End j4
            } // End j3
         }  // End j2
      }  //End j1

   } // End i
   free(neighbors);
   return nCount;
}

// MAX NEIGHBORHOOD SIZE FORMULAS ----------------------------------------------------------

// Formula: 1/3 * (2n^4 + 4n^3 + 10n^2 + 8n + 3) for each neighborhood
// Because we are in 4D, eight neighborhoods in total, plus the original point.
int maxNeighbors4D(short p){
	return pow(2*p+1, 4);	
}

/******************************************************************************************************
                                             4D Lattices
******************************************************************************************************/

void NaturalOrdering(short numk, short k[], short p, short d[], FILE *file, short *colors, short *numColors){

   fprintf(stderr, "Beginning Natural Ordering\n");
   // Find number of nodes in Lattice
   int N = d[0] * d[1] * d[2] * d[3];

   // allocate array for permutation *********************
   int x;
   int* perm = (int*)calloc(N, sizeof(int));
   for(x = 0; x < N; x++)
      perm[x] = x;            // NATURAL ORDERING
   // ****************************************************  


   // Creating the Stencil ******************************* 
   // No coloring occurs here incase node order changes 
   fprintf(stderr, "Creating Stencil.\n");
   int maxN = maxNeighbors4D(p) * 2 * numk;
   int nCount;
   short stencil[maxN][4];
   nCount = createStencil(stencil, p, numk, k, maxN, d);
   fprintf(stderr, "Stencil created with %d elements. Beginning coloring.\n", nCount);
   // *************************************************** 

   int *Neighbors = (int*)malloc(nCount * sizeof(int));
   short *colorArray = (short*)calloc(10000, sizeof(short));
   short *colors2 = (short*)calloc(N, sizeof(short));
   int ix, j;
   short coord[4], coord2[4];
   for(x = 0; x < N; ++x){
      ix = perm[x];
      index2coord(ix, d, coord, 4);
      for(j = 0; j < nCount; ++j){

         coord2[0] = (coord[0] + stencil[j][0]) % d[0];
         coord2[1] = (coord[1] + stencil[j][1]) % d[1];
         coord2[2] = (coord[2] + stencil[j][2]) % d[2];
         coord2[3] = (coord[3] + stencil[j][3]) % d[3];
         Neighbors[j] = coord2index(coord2, d, 4);

      } // For each stencil item ended.
      
      colors2[ix] = colorNode(Neighbors, nCount, colors2, colorArray);
   
   } // For each node ended.

   int nc = findMax(colors2, N);
   fprintf(stderr, "Number of colors: %d\n", nc);

   if(nc < (*numColors)){
      (*numColors) = nc;
      for(x = 0; x < N; ++x)
         colors[x] = colors2[x];
   }

   fprintf(file, "Natural Ordering. Number of colors: %d\n", findMax(colors2, N));

   free(perm);
   free(Neighbors);
   free(colorArray);
   free(colors2);   

}

void RandomOrdering(short numk, short k[], short p, short d[], FILE *file, short *colors, short *numColors){

   fprintf(stderr, "Beginning Random Ordering\n");
   // Find number of nodes in Lattice
   int N = d[0] * d[1] * d[2] * d[3];
   
   short RandBest = 10000;
   short i;
   short *colors2 = (short*)calloc(N, sizeof(short));
   for(i = 0; i < 1000; ++i){

      // allocate array for permutation *********************
      int x;
      int* perm = (int*)calloc(N, sizeof(int));
      randperm(perm, N);   
      // ****************************************************  

      // Creating the Stencil ******************************* 
      // No coloring occurs here incase node order changes 
      fprintf(stderr, "Creating Stencil.\n");
      int maxN = maxNeighbors4D(p) * 2 * numk;
      int nCount;
      short stencil[maxN][4];
      nCount = createStencil(stencil, p, numk, k, maxN, d);
      fprintf(stderr, "Stencil created with %d elements. Beginning coloring.\n", nCount);
      // *************************************************** 

      int *Neighbors = (int*)malloc(nCount * sizeof(int));
      short *colorArray = (short*)calloc(10000, sizeof(short));
      int ix, j;
      short coord[4], coord2[4];
      for(x = 0; x < N; ++x){
         ix = perm[x];
         index2coord(ix, d, coord, 4);
         for(j = 0; j < nCount; ++j){

            coord2[0] = (coord[0] + stencil[j][0]) % d[0];
            coord2[1] = (coord[1] + stencil[j][1]) % d[1];
            coord2[2] = (coord[2] + stencil[j][2]) % d[2];
            coord2[3] = (coord[3] + stencil[j][3]) % d[3];
            Neighbors[j] = coord2index(coord2, d, 4);

         } // For each stencil item ended.
      
         colors2[ix] = colorNode(Neighbors, nCount, colors2, colorArray);
   
      } // For each node ended.

      int nc = findMax(colors2, N);

      if(nc < (*numColors)){
         (*numColors) = nc;
         for(x = 0; x < N; ++x)
            colors[x] = colors2[x];
      }

      if(nc < RandBest)
         RandBest = nc;

      free(Neighbors);
      free(colorArray);
      free(perm);
   }

   fprintf(file, "Random Ordering. Best Number of colors: %d\n", RandBest);

   free(colors2);

}

void RBOrdering(short numk, short k[], short p, short d[], FILE *file, short *colors, short *numColors){

   fprintf(stderr, "Beginning Red-Black Ordering\n");
   // Find number of nodes in Lattice
   int N = d[0] * d[1] * d[2] * d[3];

   // allocate array for permutation *********************
   int x;
   int* perm = (int*)calloc(N, sizeof(int));
   RB(perm, d);
   // ****************************************************  

   // Creating the Stencil ******************************* 
   // No coloring occurs here incase node order changes 
   fprintf(stderr, "Creating Stencil.\n");
   int maxN = maxNeighbors4D(p) * 2 * numk;
   int nCount;
   short stencil[maxN][4];
   nCount = createStencil(stencil, p, numk, k, maxN, d);
   fprintf(stderr, "Stencil created with %d elements. Beginning coloring.\n", nCount);
   // *************************************************** 

   int *Neighbors = (int*)malloc(nCount * sizeof(int));
   short *colorArray = (short*)calloc(10000, sizeof(short));
   short *colors2 = (short*)calloc(N, sizeof(short));
   int ix, j;
   short coord[4], coord2[4];
   for(x = 0; x < N; ++x){
      ix = perm[x];
      index2coord(ix, d, coord, 4);
      for(j = 0; j < nCount; ++j){

         coord2[0] = (coord[0] + stencil[j][0]) % d[0];
         coord2[1] = (coord[1] + stencil[j][1]) % d[1];
         coord2[2] = (coord[2] + stencil[j][2]) % d[2];
         coord2[3] = (coord[3] + stencil[j][3]) % d[3];
         Neighbors[j] = coord2index(coord2, d, 4);

      } // For each stencil item ended.
      colors2[ix] = colorNode(Neighbors, nCount, colors2, colorArray);
   
   } // For each node ended.

   int nc = findMax(colors2, N);
   fprintf(stderr, "Number of colors: %d\n", nc);

   if(nc < (*numColors)){
      (*numColors) = nc;
      for(x = 0; x < N; ++x)
         colors[x] = colors2[x];
   }

   fprintf(file, "Red-Black Ordering. Number of colors: %d\n", findMax(colors2, N));

   free(perm);
   free(Neighbors);
   free(colorArray);
   free(colors2);
   
}

void RBRandomOrdering(short numk, short k[], short p, short d[], FILE *file, short *colors, short *numColors){

   fprintf(stderr, "Beginning Random Red-Black Ordering\n");
   // Find number of nodes in Lattice
   int N = d[0] * d[1] * d[2] * d[3];
   short RBRandBest = 10000;   
   short i;
   short *colors2 = (short*)calloc(N, sizeof(short));
   for(i = 0; i < 1000; ++i){

      // allocate array for permutation *********************
      int x;
      int* perm = (int*)calloc(N, sizeof(int));
      RBRand(perm, d);   
      // ****************************************************  

      // Creating the Stencil ******************************* 
      // No coloring occurs here incase node order changes 
      fprintf(stderr, "Creating Stencil.\n");
      int maxN = maxNeighbors4D(p) * 2 * numk;
      int nCount;
      short stencil[maxN][4];
      nCount = createStencil(stencil, p, numk, k, maxN, d);
      fprintf(stderr, "Stencil created with %d elements. Beginning coloring.\n", nCount);
      // *************************************************** 

      int *Neighbors = (int*)malloc(nCount * sizeof(int));
      short *colorArray = (short*)calloc(10000, sizeof(short));
      int ix, j;
      short coord[4], coord2[4];
      for(x = 0; x < N; ++x){
         ix = perm[x];
         index2coord(ix, d, coord, 4);
         for(j = 0; j < nCount; ++j){

            coord2[0] = (coord[0] + stencil[j][0]) % d[0];
            coord2[1] = (coord[1] + stencil[j][1]) % d[1];
            coord2[2] = (coord[2] + stencil[j][2]) % d[2];
            coord2[3] = (coord[3] + stencil[j][3]) % d[3];
            Neighbors[j] = coord2index(coord2, d, 4);

         } // For each stencil item ended.
      
         colors2[ix] = colorNode(Neighbors, nCount, colors2, colorArray);
   
      } // For each node ended.

      int nc = findMax(colors2, N);

      if(nc < (*numColors)){
         (*numColors) = nc;
         for(x = 0; x < N; ++x)
            colors[x] = colors2[x];
      }
      if(nc < RBRandBest)
         RBRandBest = nc;

      free(perm);
      free(Neighbors);
      free(colorArray);

   }

   fprintf(file, "RB Random Ordering. Best Number of colors: %d\n", RBRandBest);
   free(colors2);
}

int main(int argc, char* argv[]){

   if(argc < 8){
      printf("Invalid number of arguments.\nTo run: '<executable> <Size of Dims> <Number of Displacements> <List of Displacements> <Coloring Distance>'\n Ex: './LatticeColor 16 16 16 16 1 2 5'\n"); 
      return 1;
   }

   short d[4] = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])};
   short numk = atoi(argv[5]);
   short k[numk];   
   short i;

   for(i = 6; i < 6 + numk; ++i)
      k[i-6] = atoi(argv[i]);

   short p = atoi(argv[6 + numk]);

   clock_t start, end;
   start = clock();
   srand(time(NULL));

   int N = d[0] * d[1] * d[2] *d[3];
   short *colors = (short*)calloc(N, sizeof(short));
   FILE *file = fopen("L2Tiles.txt", "a");
   fprintf(file, "d = [%d %d %d %d], distance = %d, k = [", d[0], d[1], d[2], d[3], p);
   for(i = 0; i < numk-1; ++i)
      fprintf(file, "%d, ", k[i]);
   fprintf(file, "%d]\n", k[numk-1]);
   
   short numColors = 10000;

   NaturalOrdering(numk, k, p, d, file, colors, &numColors);
   RBOrdering(numk, k, p, d, file, colors, &numColors);
   //RandomOrdering(numk, k, p, d, file, colors, &numColors);
   //RBRandomOrdering(numk, k, p, d, file, colors, &numColors);

   fprintf(file, "Best number of colors: %d\n\n", findMax(colors, N));
   fclose(file);

   char filename[100];
   sprintf(filename, "L2tile%d_%d_%d_%dk", d[0], d[1], d[2], d[3]);
   char temp[20];
   for(i = 0; i < numk-1; ++i){
      sprintf(temp, "%d_", k[i]);
      strcat(filename, temp);
   }
   sprintf(temp, "%dp%dc%d.txt", k[numk-1], p, findMax(colors, N));
   strcat(filename, temp);

   fprintf(stderr, "%s\n", filename);
   FILE *file2 = fopen(filename, "w");

   int x;
   for(x = 0; x < N; ++x)
      fprintf(file2, "%d\n", colors[x]);
   fclose(file2);

   end = clock();
   fprintf(stderr, "Coloring saved to %s\n", filename);
   fprintf(stderr, "Total run time: %4f\n", ((double)(end-start)/CLOCKS_PER_SEC));


   return 0;
}
