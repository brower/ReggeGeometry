/**********

r1 = (0,1,1), r2 = (1,0,1), r3 = (1,1,0)

Xindex(nx,ny,nz) = n1 * r1 + n2*r2 + n3 *r3
n1 + n2 + n3 <= L

integer n's gerated even lattice
n's all half integer are odd lattice


Nlayer[L_] := Sum[l + 1, {l, 0, L}] = 1/2 (1 + L) (2 + L)
NumVert[L_] := Sum[1/2  (1 + l)  (2 + l), {l, 0, L}]  = 1/6 (1 + L) (2 + L) (3 + L)

**********/

#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include <set>
#include "S3.h"
using namespace std;


int main( int argc, char *argv[])
{
  param p;
  p.L = 4;
  
  if(argc ==1)
    printf("In program  %s no input default %d \n",argv[0],p.L);
  else
    p.L = atoi(argv[1]);
  p.N = p.L*p.L*p.L;
   
  printf("L = %d  N = %d  \n", p.L,p.N);

  p.EvenNum = ((p.L)*(p.L+1)*(p.L +2))/6;

  cout << "  p.EvenNum "<<  p.EvenNum << endl;
  
  int tetraVertex[p.EvenNum][3];
  Eigen::Vector3d r1 = {0,0,0}; 
  Eigen::Vector3d r2 = {1,1,0}; 
  Eigen::Vector3d r3 = {1,0,1}; 
  Eigen::Vector3d r4 = {0,1,1}; 
  
  printf("test\n");
  std::vector<int> test = {1,2,3,4};
  std::vector<int> object = {3,1,2,4}; 
  std::vector<int> res = MatchIndices(test, object); //gives indices to match object[i] to target
  //res = (1,2,0,3)


  // int nx,ny, nz;

  int n = 0;
  
  for(int lev =0; lev < p.L; lev++)
      {
	for(int n1 = 0; n1<= lev; n1++)
	  for(int n2 = 0; n2<= lev; n2++)
	    for(int n3 = 0; n3<= lev; n3++)
	      {	
		if (n1 + n2 + n3 == lev)
		  { // nx = n2 + n3; ny = n1 + n3;  nz = n1 + n2;
		    printf(" (%d,%d,%d) ",  n2 + n3, n1 + n3, n1 + n2);
		    // This convention for printing out the lattice coordinates: 
        //tetraVertex[n][0] =  n2 + n3; tetraVertex[n][1] =  n1 + n3; tetraVertex[n][2] =  n1 + n2;
		    tetraVertex[n][0] =  n1; tetraVertex[n][1] =  n2; tetraVertex[n][2] =  n3;
        n++; 
		  }
	      }
	printf("\n \n");
      }

  printEven(tetraVertex, p);
  std::vector<int> orbitlabels = FindOrbitLabels(tetraVertex,p);
  for (int i = 0; i<p.EvenNum; i++)
  {std::cout<<orbitlabels[i]<<" "<< tetraVertex[i][0]<<" "<<tetraVertex[i][1]<<" "<<tetraVertex[i][2]<<"\n"; }
  printf("\n"); 
  std::vector<std::vector<int>> orbit_reps = FindOrbitRepresentative(orbitlabels, tetraVertex, p); 
  for (size_t i = 0; i < orbit_reps.size(); ++i) 
  {  // Loop over rows
    std::cout << i << " ";
    for (size_t j = 0; j < orbit_reps[i].size(); ++j) 
    {  // Loop over columns
        std::cout << orbit_reps[i][j] << " ";
    }
    std::cout << std::endl;  // New line after each row    
  }
  for (int i=0; i<p.EvenNum; i++){
    std::cout << i << " "<< orbitlabels[i] << "\n";
  }

  std::cout << "Printing Tetrahedrons...\n"; 
  std::vector<T> tetra= FindTetrahedronRepresentative(orbit_reps, tetraVertex,p); 
  for (size_t i = 0; i < tetra.size(); ++i) {  // Loop over rows
        std::cout << i << " ";
        std::cout << " wt: "<< tetra[i].wt<<"\n";
    for (size_t j = 0; j < tetra[i].Vert.size(); ++j) {  // Loop over columns
        std::cout << tetra[i].Vert[j] << " ";
    }
    std::cout << std::endl;
    for (size_t j = 0; j < tetra[i].ClassLabels.size(); ++j) {  // Loop over columns
        std::cout << tetra[i].ClassLabels[j] << " ";
    }
    std::cout << std::endl;  // New line after each row
}  
  return 0;
}


