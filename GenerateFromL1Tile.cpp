#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include <dirent.h>


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

int main(int argc, char* argv[]){
/* Accepts a desired 4d size (pts[4]) which are all multiples of  
 * and a filename for the coloring to be used to tile it.
 * Finally, type in the 4d size of the tile at <filename>
 * Ex:   
 * ./Tiling 32 32 32 64 TileFileName 16 16 16 16 >> NewColoringFile
 * To take the 16^4 coloring from 'TileFileName", and make it into a 32^3 x 64 tile
 * with coloring located in 'NewColoringFile'
 */

   if(argc < 10)
   {
      fprintf(stderr, "Too few arguments. Aborting. \n");
      return 1;
   }

   short pts[4] = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])};
   char filename[100];
   strcpy(filename, argv[5]);
   printf("Tile coloring from: %s\n",filename);

   FILE* file = NULL;
   if((file = fopen(filename, "r")) == NULL){
      fprintf(stderr, "Unable to open tiling file. Aborting.\n");
      return 1;
   }

   int color;
   size_t len = 0;
   short coord[4];

   short tile[4] = {atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9])};

   short colors[tile[0]][tile[1]][tile[2]][tile[3]];

   int index = 0;
   while(fscanf(file, "%d", &color) != EOF){
      index2coord(index, tile, coord, 4);
      colors[coord[0]][coord[1]][coord[2]][coord[3]] = color;
      ++index;
   }
   fclose(file);

 
   short d1, d2, d3, d4;
   for(d4 = 0; d4 < pts[3]; ++d4)
      for(d3 = 0; d3 < pts[2]; ++d3)
         for(d2 = 0; d2 < pts[1]; ++d2)
            for(d1 = 0; d1 < pts[0]; ++d1){
               printf("%d\n", colors[d1 % tile[0]][d2 % tile[1]][d3 % tile[2]][d4 % tile[3]]);
            }

return 0;
}
