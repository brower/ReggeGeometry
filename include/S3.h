#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <set>
#include <Eigen/Core> // Core Eigen functionality
#include <Eigen/Sparse>
#include <Eigen/Dense>


template <typename T>
using Tensor2D = std::vector<std::vector<T>>;

struct param{
  int L;
  int N;      
  int EvenNum;
};      


struct T{
  std::vector<int> Vert; 
  std::vector<int> ClassLabels; 
  int wt; 
};

double TetrahedronVol(Eigen::Vector4d& r1, Eigen::Vector4d& r2, Eigen::Vector4d& r3, Eigen::Vector4d& r4)
{
  Eigen::Vector4d e1 = r4-r1;
  Eigen::Vector4d e2 = r3-r1;
  Eigen::Vector4d e3 = r2-r1;  
  return std::abs(e1.dot(e2.cross(e3)))/6.0; 
}



//Compute Nabla_x V for x a vertex of the tetrahedron (Nabla_x V  = 1/3 A hat{N} where A is the area of the triangle
//facing the vertex, hat{N} is the inwards normal)
Eigen::Vector4d VertexDerivative(Eigen::Vector4d& r1, Eigen::Vector4d& r2, Eigen::Vector4d& r3, Eigen::Vector4d& r4)
{
  Eigen::Vector4d e1 = r4-r2;  
  Eigen::Vector4d e2 = r3-r2; //e1 e2 defines vectors on the opposing triangle
  Eigen::Vector4d e3 = r2-r1; //e3 is used to ensure the direction of Nabla_xV is inwards normal
  Eigen::Vector4d Area = e1.cross(e2)/2.0;  //Need to generalize the way to find normal direction
  if (Area.dot(e3)>=0){return -Area/3.0;}
  else {return Area/3.0;}
}

Tensor2D<Eigen::Vector4d> FindDerivative(std::vector<Eigen::Vector4d>& r, Tensor2D<int>& orbitReps)
{
  Tensor2D<Eigen::Vector4d> results; 
  for (size_t i=0; i< orbitReps.size(); i++)
  {
    Eigen::Vector4d nd={0.0, 0.0, 0.0, 0.0};
    std::vector<int> vec = orbitReps[i]; 
    auto it = std::find_if(vec.begin(), vec.end(), [](int i) { return i != 0; });
    int distance = std::distance(vec.begin(), it); //Number of steps from the begining to first non-zero
    
    if (distance!=0)
    {//Not interior points     
      std::set<int> unique_labels(it, vec.end());
      int n_distinct = unique_labels.size(); //Number of distinct points for the non-zero part

      if (distance==3)
      {results.push_back({nd, nd, nd, nd});} //We are at the vertex of the tetrahedron, not a dof
      
      else if (distance==2)
      {//We are at an edge
          if (n_distinct==1){results.push_back({nd, nd, nd});}
          else{results.push_back({nd, nd, r[2]-r[3]});}
      }

      else if (distance==1)//We are at a face
      {
        if (n_distinct==1){results.push_back({nd, nd, nd});} //Center of Triangle face

        else if (n_distinct==2)
          {
            if (vec[distance]==vec[distance+1])
            {//xi^2 = xi^3 constrant//There is only one d.o.f.
              Eigen::Vector4d v1 = r[1]+r[2]-2*r[3]; 
              results.push_back({nd, v1, nd});
            } 
            else
              {//xi^3 = xi^4 constrant
              Eigen::Vector4d v1 = r[1]-2*r[2]-2*r[3]; 
              results.push_back({nd, v1, nd});
              }
          }
        else
          {
            Eigen::Vector4d v1 = r[1]-r[3]; 
            Eigen::Vector4d v2 = r[2]-r[3]; 
            results.push_back({nd, v1, v2});           
          }
      }
    }

    else
    {//We are at an interior point 
      std::set<int> unique_labels(vec.begin(), vec.end());
      int n_distinct = unique_labels.size(); //Number of distinct points for the non-zero part
      if (n_distinct==1){ results.push_back({nd, nd, nd}); } //Skip at center of tetrahedron
      else if (n_distinct==2)
      {
        if (vec[0]==vec[1])//patterns (a,a,a,b) (a,a,b,b)
        {
          if (vec[0]==vec[2]) //pattern (a,a,a,b)
          {
            Eigen::Vector4d v0 = r[0]+r[1]+r[2]-3*r[3]; 
            results.push_back({v0, nd, nd});     
          }
          else //pattern (a,a,b,b)
          {
            Eigen::Vector4d v0 = r[0]+r[1]-r[2]-r[3]; 
            results.push_back({v0, nd, nd});    
          }
        }
      else //patterns (a,b,b,b) (a,b,b,a)
      {
        if (vec[1]==vec[3])//(a,b,b,b) 
        {
          Eigen::Vector4d v0 = r[0]-(r[1]+r[2]+r[3])/3; 
          results.push_back({v0, nd, nd});     
        }
        else //(a,b,b,a)
        {
          Eigen::Vector4d v0 = r[0]-r[1]-r[2]+r[3]; 
          results.push_back({v0, nd, nd});    
        }   
      }  
      }
      else if (n_distinct ==3)
      {
        if (vec[0]==vec[1]) //(a,a,b,c) and 
        {
          Eigen::Vector4d v0 = r[0]+r[1]-2*r[3]; 
          Eigen::Vector4d v1 = r[2]-r[3]; 
          results.push_back({v0, nd, v1});    
        }
        else if (vec[1]==vec[2])//(a,b,b,c)
        {
          Eigen::Vector4d v0 = r[0]-r[3]; 
          Eigen::Vector4d v1 = r[1]+r[2]-2*r[0]; 
          results.push_back({v0, v1, nd});    
        }
        else //(a,b,c,c)
        {
          Eigen::Vector4d v0 = r[0]-(r[2]+r[3])*0.5; 
          Eigen::Vector4d v1 = r[1]-(r[2]+r[3])*0.5; 
          results.push_back({v0, v1, nd});    
        }
      }
      else
      {
        Eigen::Vector4d v0 = r[0]-r[3]; 
        Eigen::Vector4d v1 = r[1]-r[3]; 
        Eigen::Vector4d v2 = r[2]-r[3]; 
        results.push_back({v0, v1, v2});       
      }
    }
  }
  return results; 
}





Tensor2D<int> FindOrbitRepresentative(std::vector<int>& OrbitLabels, int tetraVertex[][3], param p){
  std::set<int> unique_labels(OrbitLabels.begin(), OrbitLabels.end());
  int n_orbit = unique_labels.size(); 
  std::vector<std::vector<int>> result; 
  for (int o=0; o<n_orbit; o++){
    auto it = std::find(OrbitLabels.begin(), OrbitLabels.end(), o);
    int index = std::distance(OrbitLabels.begin(), it);
    std::vector<int> pos = {tetraVertex[index][0],tetraVertex[index][1],tetraVertex[index][2],
                     p.L-1-tetraVertex[index][0]-tetraVertex[index][1]-tetraVertex[index][2]};
    result.push_back(pos); 
  }
  return result;
}

//gives indices to reorder target like objective
std::vector<int> MatchIndices(std::vector<int> &objective, std::vector<int> &target)
{
    std::vector<int> indices;
    std::unordered_map<int, std::vector<int>> index_map;
    for (int i = 0; i < target.size(); ++i) {
        index_map[target[i]].push_back(i);
    }
    std::unordered_map<char, int> used_count;
    for (const auto& elem : target) 
    {
        used_count[elem] = 0;
    }
    for (int element: objective)
    {
        int idx = used_count[element];
        indices.push_back(index_map[element][idx]);
        used_count[element]++;
    }
  return indices; 
    
}


int FindRep(Tensor2D<int>& orbit_reps, int n1, int n2, int n3, param p){
  std::vector<int> v = {n1, n2, n3, p.L-1-n1-n2-n3}; 
  std::sort(v.begin(), v.end()); 
  auto it = std::find(orbit_reps.begin(), orbit_reps.end(), v);
  int index = std::distance(orbit_reps.begin(), it);
  return index;
}

int FindInd(int tetraVertex[][3], int n1, int n2, int n3, param p){
  int indx = 0; 
  for (int i =0; i<p.EvenNum; i++){
      if (n1==tetraVertex[i][0]&&n2==tetraVertex[i][1]&&n3==tetraVertex[i][2]){
        indx = i; 
        break;
      }
  }
  return indx;
}

//Check if we get correct total number of tetrahedrons!!
//std::vector<std::vector<int>> 
std::vector<T> FindTetrahedronRepresentative(Tensor2D<int>& orbit_reps, int tetraVertex[][3] ,param p){
  std::vector<std::vector<int>> tmp; 
  std::vector<T> result; 
  for(int lev =0; lev < p.L-1; lev++)
    {
	for(int n1 = 0; n1<= lev; n1++)
	  for(int n2 = 0; n2<= lev; n2++)
	    for(int n3 = 0; n3<= lev; n3++)
	      {	
		if (n1 + n2 + n3 == lev)
		  { 
        int v1 = FindRep(orbit_reps, n1, n2, n3, p); 
        int v2 = FindRep(orbit_reps, n1, n2, n3+1, p);
        int v3 = FindRep(orbit_reps, n1, n2+1, n3, p);
        int v4 = FindRep(orbit_reps, n1+1, n2, n3, p);
        std::vector<int> tetrahedron = {v1,v2,v3,v4}; 
        std::sort(tetrahedron.begin(), tetrahedron.end()); 
        auto it = std::find(tmp.begin(), tmp.end(), tetrahedron);
        if (it==tmp.end()){
          int ind1 = FindInd(tetraVertex, n1, n2, n3, p);    
          int ind2 = FindInd(tetraVertex, n1, n2, n3+1,p); 
          int ind3 = FindInd(tetraVertex, n1, n2+1, n3, p); 
          int ind4 = FindInd(tetraVertex, n1+1, n2, n3,p); 
          // std::vector<int> indList = {ind1, ind2, ind3, ind4}; 
          // std::sort(indList.begin()+1, indList.end());
          T tetra = {{ind1, ind2, ind3, ind4}, {v1,v2,v3, v4}, 1};
          tmp.push_back(tetrahedron);
          result.push_back(tetra); }
        else{
          int index = std::distance(tmp.begin(), it);
          result[index].wt+=1;
        }
		  }
	      }

      }
  return result; 
}

std::vector<int> FindOrbitLabels(int tetraVertex[][3], param p){

  std::vector<int> counter(p.EvenNum); //to keep track of points already labeled
  std::vector<int> classlabels(p.EvenNum); //a list that returns the label for all points 
  int dof_counter = 0; 
  for(int n = 0; n<p.EvenNum; n++){
      if (counter[n]!=0){continue;}
      else{
        classlabels[n]=dof_counter;
        counter[n] = 1; 
        int ncounted = 1; 
        std::vector<int> xi= {tetraVertex[n][0],tetraVertex[n][1],tetraVertex[n][2],
                    p.L-1-tetraVertex[n][0]-tetraVertex[n][1]-tetraVertex[n][2]};
        std::sort(xi.begin(), xi.end()); 
      // std::cout << xi[0] << ' ' << xi[1] << ' ' << xi[2]<< ' ' << xi[3] << '\n';
       // std::cout << "begin Check....\n";
        for (int n1=0; n1<p.EvenNum; n1++){
            if (ncounted==24){break;}
            if (counter[n1]!=0){continue;}
            else{
            std::vector<int>  xi2 = {tetraVertex[n1][0],tetraVertex[n1][1],tetraVertex[n1][2],
                    p.L-1-tetraVertex[n1][0]-tetraVertex[n1][1]-tetraVertex[n1][2]}; 
            std::sort(xi2.begin(),xi2.end());

            //std::cout << xi2[0] << ' ' << xi2[1] << ' ' << xi2[2]<< ' ' << xi2[3] << '\n';

            if (xi==xi2){
              classlabels[n1]=dof_counter;
              counter[n1] = 1; 
              ncounted += 1; 
              //std::cout << "Success....\n";
            }
            }
        }
        //std::cout << "end Check....\n";
        dof_counter+=1; 
      }
  }

  return classlabels;}


void printEven(int tetraVertex[][3], param p)
{
  FILE* fptr = NULL;  // C style
  char out_name[64];
  //  sprintf(out_name,"data/CGstate_%d_%d.dat",N,frame); // filename
  sprintf(out_name,"data/EvenVert_%d.dat",p.L); // filename
  fptr = fopen(out_name,"w");
  if(fptr == NULL)
    {
      printf("Error!");   
      exit(1);             
    }
  
  for(int n = 0; n<p.EvenNum; n++)
	{
	  fprintf(fptr, " %d ", n);
	  for(int k = 0; k<3; k++)
	    {
	      fprintf(fptr,"  %d ", tetraVertex[n][k]);
	    }
	  fprintf(fptr,"\n");
	}
  fclose(fptr);
}
