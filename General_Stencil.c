/******************************************************************
 * Heather Switzer
 * hmswitzer@email.wm.edu
 *
 * To compile:
 *    gcc General_Stencil.c -o GS -O3 -lm
 * To run:
 *    ./GS <options>
 * To get a full list of options, run './GS --help' 
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Helper functions -------------------------------------------

#define MAX(x, y) ((x) > (y)) ? (x) : (y)
#define MIN(x, y) ((x) < (y)) ? (x) : (y)

// Convert lattice coordinate into array index
int coord2index(short *c, short *d, short n)
{

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
void index2coord(int ix, short d[], short c[], short n)
{
	short i;
	for(i = 0; i < n; ++i){
		c[i] = ix % d[i];
		ix = (int)(ix/d[i]);
	}
}

void Find_Neighbor_Locations_L1(short *record, short *temp, short nk, short *k, short dk, short p, short dim_num, short d, short *dims)
{
   
   int i;
   if(dim_num == d)
   {
      for(i = 0; i < nk; i++)
      {
         temp[dk] = (temp[dk] + k[i]) % dims[dk];
         record[coord2index(temp, dims, d)] = 1;
         temp[dk] = ((temp[dk] - 2*k[i]) + dims[dk]) % dims[dk];
         record[coord2index(temp, dims, d)] = 1;
         temp[dk] = (temp[dk] + k[i]) % dims[dk];
      }
      return;
   }

   for(i = -p; i <= p; i++)
   {
      temp[dim_num] = (i+dims[dim_num]) % dims[dim_num];
      Find_Neighbor_Locations_L1(record, temp, nk, k, dk, p-abs(i), dim_num+1, d, dims);
   }

}

void Find_Neighbor_Locations_L2(short *record, short *temp, short nk, short *k, short dk, short p, short dim_num, short d, short *dims, int sum)
{
   
   int i;
   if(dim_num == d)
   {
      if(sqrt(sum) <= p)
         for(i = 0; i < nk; i++)
         {
            temp[dk] = (temp[dk] + k[i]) % dims[dk];
            record[coord2index(temp, dims, d)] = 1;
            temp[dk] = ((temp[dk] - 2*k[i]) + dims[dk]) % dims[dk];
            record[coord2index(temp, dims, d)] = 1;
            temp[dk] = (temp[dk] + k[i]) % dims[dk];
         }
      return;
   }

   for(i = -p; i <= p; i++)
   {
      temp[dim_num] = (i+dims[dim_num]) % dims[dim_num];
      Find_Neighbor_Locations_L2(record, temp, nk, k, dk, p, dim_num+1, d, dims, sum+pow(i, 2));
   }

}

int Create_Stencil(short **Stencil, short *record, short nk, short *k, short dk, short p, short d, short *dims, int num_nodes, short norm)
{
   int i, j, count;
   short coord[d];

   short *temp = (short *)malloc(d * sizeof(short));
   if(norm == 1)
      Find_Neighbor_Locations_L1(record, temp, nk, k, dk, p, 0, d, dims);
   else
      Find_Neighbor_Locations_L2(record, temp, nk, k, dk, p, 0, d, dims, 0);

   record[0] = -1;
   count = 0;
   for(i = 1; i < num_nodes; ++i)
   {
      if(record[i] == 1)
      {
         index2coord(i, dims, coord, d);
         for(j = 0; j < d; j++)
            Stencil[count][j] = coord[j];
         count++;
      }
      record[i] = -1;
   }

   free(temp);

   return count;

}

/* Find the number of neighbors each node would have in an infinite lattice. Used to allocate memory */
void Max_NeighborsL1(short p, short dim_num, short d, int *max_Neighbors)
{

   if(dim_num == d)
   {
      (*max_Neighbors)++;
      return;
   }

   int i;
   for(i = -p; i <= p; i++)
      Max_NeighborsL1(p-abs(i), dim_num+1, d, max_Neighbors);

}

/* Find the number of neighbors each node would have in an infinite lattice. Used to allocate memory */
void Max_NeighborsL2(short p, short dim_num, short d, int sum, int *max_Neighbors)
{

   if(dim_num == d)
   {
      if((short)sqrt(sum) <= p)
         (*max_Neighbors)++;
      return;
   }

   int i;
   for(i = -p; i <= p; i++)
      Max_NeighborsL2(p, dim_num+1, d, sum+pow(i, 2), max_Neighbors);

}

/* Swap function for random ordering */
void swap(int *a, int *b)
{
   int temp = *a;
   *a = *b;
   *b = temp;
   return;
}

/* Randomly shuffle an array */
void shuffle(int *perm, int N)
{
   int i, j;

   for(i = N-1; i > 0; i--)
      swap(&perm[i], &perm[rand()%(i+1)]);

   return;
}

void ColorLattice(short *colors, short **Stencil, int count, short d, short *dims, int num_nodes, short ordering)
{

   int *perm               = (int*)malloc(num_nodes * sizeof(int));     /* Permutation array for coloring nodes */

   bool *avail_colors      = (bool*)malloc(count+1 * sizeof(bool));
   short *coord            = (short*)malloc(d * sizeof(short));         /* For index2coord */
   short *new_coord        = (short*)malloc(d * sizeof(short));         /* Use for applying stencil to current node location */

   int *neighbor_colors    = (int*)malloc(count * sizeof(int));         /* To make it easier to reset avail_colors for each iteration. More storage, but less computations */
   short num_unique_cols   = 0;  

   int i, j, l, ix;            /* Loop variables */
   short color;

   /* Ordering: 1 = Natural, 2 = Random, 3 = Red-Black, 4 = Red-Black+Random */
   int numOdds    = 0;
   int numEvens   = 0;
   int midpoint   = num_nodes/2;
   switch(ordering)
   {
      case 1:     /* Natural Ordering */

         for(i = 0; i < num_nodes; i++)
            perm[i] = i;
         break;

      case 2:     /* Random Ordering */

         for(i = 0; i < num_nodes; i++)
            perm[i] = i;
         shuffle(perm, num_nodes);
         break;

      case 3:     /* Red-Black Ordering */

         for(i = 0; i < num_nodes; i++)
         {
            l = 0;
            index2coord(i, dims, coord, d);
            for(j = 0; j < d; j++)
               l += coord[j];
            if(l % 2 == 0)
               perm[numEvens++] = i;
            else
               perm[midpoint+(numOdds++)] = i;
         }
         break;

      case 4:     /* Red-Black + Random Ordering */
         for(i = 0; i < num_nodes; i++)
         {
            l = 0;
            index2coord(i, dims, coord, d);
            for(j = 0; j < d; j++)
               l += coord[j];
            if(l % 2 == 0)
               perm[numEvens++] = i;
            else
               perm[midpoint+(numOdds++)] = i;
         }   
         shuffle(perm, numEvens);
         shuffle(perm+numEvens, numOdds);

   } /* End Permutation Creation */

   for(i = 0; i <= count; i++)
      avail_colors[i] = true;

   for(ix = 0; ix < num_nodes; ix++)
   {

      i = perm[ix];

      /* Reset boolean array storing which colors are already in use */
      for(j = 0; j < num_unique_cols; j++)
         avail_colors[neighbor_colors[j]] = true;
      num_unique_cols = 0;

      /* Apply Stencil to node i and record color of neighbors */
      index2coord(i, dims, coord, d);
      for(j = 0; j < count; j++)
      {
         for(l = 0; l < d; l++)
            new_coord[l] = (coord[l] + Stencil[j][l]) % dims[l];
         if((color = colors[coord2index(new_coord, dims, d)]) > -1 && avail_colors[color])
         {
            neighbor_colors[num_unique_cols++] = color;
            avail_colors[color] = false;
         }
      }      
      
      /* Find first available color not in use */
      for(j = 0; j <= count; j++)
         if(avail_colors[j])
         {
            colors[i] = j;
            break;
         }

   }

   free(avail_colors);
   free(neighbor_colors);
   free(coord);
   free(new_coord);

}

int main(int argc, char* argv[])
{

   int i, j;

   /****************************************************************
    * SETTING UP FUNCTION PARAMETERS
    ****************************************************************/

   short d                 = -1;
   short nk                = -1;
   short dk                = -1;
   short p                 = -1;
   short ordering          =  1;
   short norm              =  1;
   short *dims;
   short *k;
   bool print_to_file      = false;
   char filename[100];
   
   for(i = 1; i < argc; i++)
   {
      if(strcmp(argv[i], "--help") == 0)
      {
         printf("Lattice Coloring with Displacements Usage:\n");
         printf("   -d <num_dim>             : Number of dimensions this lattice has. d must be greater than 0, and come before '-n' option (Default 4)\n");
         printf("   -n <nx> <ny>    ...      : List of dimension sizes equal to the number of dimensions (Default 16^4)\n");
         printf("   -nk <num_dis>            : Number of displacements (Default 0)\n");
         printf("   -k <k1> <k2>    ...      : List of displacements being applied equal to the number of displacements given (Default 0)\n");
         printf("   -k_dim <kd>              : Which dimension the displacements will be placed in (default MIN(3, <num_dim>))\n");
         printf("   -p <distance>            : L1 Coloring Distance (Default 1)\n");
         printf("   -ordering <num>          : What order the visitation of nodes will occur (Default 1)\n");
         printf("      Options:\n");
         printf("         1: Natural Ordering\n");
         printf("         2: Random Ordering\n");
         printf("         3: Red-Black Ordering\n");
         printf("         4: Red-Black+Random Ordering\n");
         printf("   -norm <1 or 2>           : Specify whether the coloring distance is Manhattan or Euclidian (Default 1: Manhattan)\n");
         printf("   -outfile <filname>       : Output coloring to filename (Default prints to screen)");
         printf("\nDefault parameters: -d 4 -n 10 10 10 10 -nk 1 -k 0 -k_dim 3 -p 1 -ordering 1\n");
         return 0;
      }
      else if(strcmp(argv[i], "-d") == 0)
      {
         if(argc < i+2 || (d = atoi(argv[i+1])) < 1)
         {
            printf("Must specify a number greater than 0 after '-d' specifying the number of dimensions.\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
         dims = (short*)malloc(d * sizeof(short));
         for(j = 0; j < d; ++j)
            dims[j] = 16;
      }
      else if(strcmp(argv[i], "-n") == 0)
      {
         if(d == -1)
         {
            printf("Please specify the number of dimensions first using '-d'\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
         for(j = 1; j <= d; j++)
            dims[j-1] = (short)atoi(argv[i+j]);
      }
      else if(strcmp(argv[i], "-nk") == 0)
      {
         if(argc < i+2 || (nk = atoi(argv[i+1])) < 0)
         {
            printf("Must specify a number greater than -1 after '-nk' specifying the number of displacements.\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
         if(nk > 0)
            k = (short*)malloc(nk * sizeof(short));
         else
         {
            k = (short*)malloc(sizeof(short));
            k[1] = 0;
         }
      }
      else if(strcmp(argv[i], "-k") == 0)
      {
         if(nk == -1)
         {
            printf("Please specify the number of displacements first using '-nk'\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
         for(j = 1; j <= nk; j++)
            k[j-1] = (short)atoi(argv[i+j]);
      }
      else if(strcmp(argv[i], "-p") == 0)
      {
         if(argc < i+2 || (p = atoi(argv[i+1])) < 1)
         {
            printf("Must specify a number greater than 0 after '-p' specifying the coloring distance.\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
      }
      else if(strcmp(argv[i], "-k_dim") == 0)
      {
         if(argc < i+2 || (dk = atoi(argv[i+1])-1) < 0)
         {
            printf("Must specify a number greater than 0 after '-k_dim' specifying which dimension the displacements should be placed in (1-based indexing).\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
      }
      else if(strcmp(argv[i], "-ordering") == 0)
      {
         if(argc < i+2 || (ordering = atoi(argv[i+1])) < 1 || ordering > 4)
         {
            printf("Must specify a number between 1 and 4 after '-ordering' specifying the node visitation order.\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
      }
      else if(strcmp(argv[i], "-outfile") == 0)
      {
         if(argc < i+2)
         {
            printf("Must specify filename after '-outfile' specifying the coloring file to print to.\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
         print_to_file = true;
         strcpy(filename, argv[i+1]);
      }
      else if(strcmp(argv[i], "-norm") == 0)
      {
         if(argc < i+2 || ((norm = atoi(argv[i+1])) != 1 && norm != 2))
         {
            printf("Must specify norm type after '-norm' for finding the distance between nodes. \nMust specify either '1' for Manhattan distance or '2' for Euclidian.\n");
            printf("Please use '--help' for more information.\n");
            return 1;
         }
      }
   }

   if(d == -1)
   {
      d = 4;
      dims = (short*)malloc(4 * sizeof(short));
      for(i = 0; i < 4; i++)
         dims[i] = 16;
   }
   if(nk == -1)
   {
      nk = 1;
      k = (short*)malloc(sizeof(short));
      k[0] = 0;
   }
   if(p == -1)
      p = 1;
   if(dk == -1)
      dk = MIN(d-1, 2);

   /* Display function parameters */

   fprintf(stderr, "|*****************************|\n");
   fprintf(stderr, "|      Coloring Settings      |\n");
   fprintf(stderr, "|*****************************|\n");
   fprintf(stderr, " Number of Dimensions: %d\n", d);
   fprintf(stderr, " Dimensions: ");
   for(i = 0; i < d; i++)
      fprintf(stderr, "%d ", dims[i]);
   fprintf(stderr, "\n Coloring Distance: %d\n", p);
   fprintf(stderr, " Number of Displacements: %d\n", nk);
   fprintf(stderr, " Displaced Dimension: %d\n", dk+1);
   if(nk > 0)
   {
      fprintf(stderr, " Displacements: ");
      for(i = 0; i < nk; i++)
         fprintf(stderr, "%d ", k[i]);
   }
   else
      fprintf(stderr, "0");
   switch(ordering)
   {
      case 1:
         fprintf(stderr, "\n Ordering: Natural\n");
         break;
      case 2:
         fprintf(stderr, "\n Ordering: Random\n");
         break;
      case 3:
         fprintf(stderr, "\n Ordering: Red-Black\n");
         break;
      case 4:
         fprintf(stderr, "\n Ordering: RB + Random\n");
   }
   if(print_to_file)
      fprintf(stderr, " Output: %s\n", filename);
   else
      fprintf(stderr, " Output: Terminal\n");
   fprintf(stderr, " Norm Type: %d\n", norm);

   fprintf(stderr, "|*****************************|\n");
   
   fprintf(stderr, "\nBeginning Coloring.\n");
   clock_t start = clock();

   /****************************************************************
   * CREATING STENCIL
   ****************************************************************/
   /* Find number of nodes in lattice */
   int num_nodes = 1;
   for(i = 0; i < d; i++)
      num_nodes *= dims[i];

   int max_Neighbors = 0;
   if(norm == 1)
      Max_NeighborsL1(p, 0, d, &max_Neighbors);
   else
      Max_NeighborsL2(p, 0, d, 0, &max_Neighbors);
   max_Neighbors = MIN(max_Neighbors*2*nk, num_nodes);

   short **Stencil;
   Stencil = (short**)malloc(max_Neighbors * sizeof(short *));
   for(i = 0; i < max_Neighbors; i++)
      Stencil[i] = (short*)malloc(d * sizeof(short));

   short *colors = (short*)calloc(num_nodes, sizeof(short));
   int Num_Neighbors = Create_Stencil(Stencil, colors, nk, k, dk, p, d, dims, num_nodes, norm);

   ColorLattice(colors, Stencil, Num_Neighbors, d, dims, num_nodes, ordering);

   /* Print coloring to screen or to a specified file */
   if(print_to_file)
   {
      FILE *f = fopen(filename, "w");
      for(i = 0; i < num_nodes; i++)
         fprintf(f, "%d\n", colors[i]);
      fclose(f);
   }
   else
      for(i = 0; i < num_nodes; i++)
         printf("%d\n", colors[i]);

   short max = 0;
   for(i = 0; i < num_nodes; i++)
      if(colors[i] > max)
         max = colors[i];

   /* Free local variables */
   free(k);
   free(dims);
   for(i = 0; i < max_Neighbors; i++)
      free(Stencil[i]);
   free(Stencil);
   free(colors);

   clock_t end = clock();   

   fprintf(stderr, "Number of colors used: %d\nRun time = %3f seconds.\nExiting.\n", max, (double)(end-start)/CLOCKS_PER_SEC);

   return 0;

}
