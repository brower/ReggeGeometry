/**********************************
THESE ARE THE  BASIC TRIANGLE IDENTITIES:

Triangle(a,b,c) edge lenths squared  x = a^2, y = b^2 , z = c^2
I THINK WE SHOULD GO BACK TO a,b,c NOT the Squares! Particularly since
lstar  internal can be negative.

4 A R = \sqrt{x y z}   4 A R = a b c
4 A[x,y,z] = Sqrt[(2 x y + 2 y z + 2 z x - x x - y y - z z)]
R[x,y,z] = Sqrt[x y z]/Sqrt[(2 x y + 2 y z + 2 z x - x x - y y - z z)]
D[R[x,y,z],x] = R^3[x,y,z] ( x^2 - (y-z)^2)/(2 x^2 y z)


***************************************/

#pragma once
#include <iostream>
#include <fstream>
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
#include <unsupported/Eigen/IterativeSolvers>


#define TWOPI  6.283185307179586

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
template <typename T>
using Tensor3D = std::vector<std::vector<std::vector<T>>>;

template <typename T>
using Tensor2D = std::vector<std::vector<T>>;

struct Triangle{
  std::vector<int> Vert; 
  std::vector<int> ClassLabels; 
  int wt; 
};

struct Lattice{
  Tensor2D<int> Vertices; //size = L where L number of points in lattice
  std::vector<Eigen::Vector4d> rvec;   //size = L
  std::vector<Eigen::Vector4d> xvec;   //size = L
  std::vector<Eigen::Vector4d> xivec;  //size = L
  std::vector<int> VertexLabels;       //size = L given by FindOrbitLabels
  std::vector<int> OrbitPos;           //size = Number of Orbits
  Tensor2D<int> Basis;                 //size = (Number of Orbits, 2)
  int ndof;                            //Number of Orbits
  int q; //q =  600 cell
};

double returnRMS(Eigen::VectorXd& vec){
  double mean = vec.array().mean();
  double sq_mean = vec.array().square().mean();
  return sq_mean - mean*mean;
}

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


double TetrahedronVol(Eigen::Vector4d& r1, Eigen::Vector4d& r2, Eigen::Vector4d& r3, Eigen::Vector4d& r4)
{ //Using definition that relies on the Cayley Menger determinant
  double d12 = (r1-r2).norm(); 
  double d13 = (r1-r3).norm(); 
  double d14 = (r1-r4).norm(); 
  double d23 = (r2-r3).norm(); 
  double d24 = (r2-r4).norm(); 
  double d34 = (r3-r4).norm(); 
  MatrixXd cm;
  cm << 0,1,1,1,1,
        1,0,d12*d12, d13*d13,d14*d14,
        1,d12*d12, 0, d23*d23,d24*d24,
        1,d13*d13, d23*d23, 0,d34*d34,
        1,d14*d14, d24*d24, d34*d34,0;
  return sqrt(cm.determinant()/288.0); 
}

//Take derivative of Area w.r.t r1, where r2, r3 are the remaining coordinates of 
Eigen::Vector3d VertexDerivative(Eigen::Vector3d& r1, Eigen::Vector3d& r2, Eigen::Vector3d& r3)
{
  // Eigen::Vector3d e1 = r3-r2; //e1 defines vector opposite to vertex
  // Eigen::Vector3d e2 = r2-r1; //e2 is used to ensure the direction of Nabla_r1 A is inwards normal
  // Eigen::Vector3d normalDirection = (e2.cross(e1)).cross(e1);
  // return 0.5*e1.norm()*normalDirection/(normalDirection.norm()); //Old Code I think it does not work
  Eigen::Vector3d e1 = r1-r3; 
  Eigen::Vector3d e2 = r1-r2;
  Eigen::Vector3d e3 = r3-r2;

  double l1 = e1.dot(e1); 
  double l2 = e2.dot(e2); 
  double l3 = e3.dot(e3); 
  double l1s = Compute_ls(l1, l2,l3, l1);
  double l2s = Compute_ls(l1, l2,l3, l2);
  return e1*l1s/e1.norm()+e2*l2s/e2.norm();
}


Tensor2D<int> TriangleVertexMapping(Tensor2D<int> &Vertices, int L){
    int rowcounter = 0; 
    int arrPosCounter = 0; 
    Tensor2D<int> arr; 
    for (int layer = L; layer>0; layer--){
        
        int numOftriangle = 1+(layer-1)*2;
        // printf("Layer: %d, Num of Triangle: %d\n", layer, numOftriangle);
        int nx = 0; 
        int ny = rowcounter; 
        int poscounter = 0; 
        // int oddcounter = 0; 
        // int evencounter = 0; 
    for (int point = 0; point<numOftriangle; point++){
        if (point%2==0){
            // int x1 = nx+evencounter; 
            // int y1 = ny+evencounter; 
            auto i1 = std::find(Vertices.begin(), Vertices.end(), std::vector<int>{nx, ny, L-nx-ny});
            auto i2 = std::find(Vertices.begin(), Vertices.end(), std::vector<int>{nx+1, ny, L-nx-ny-1});
            auto i3 = std::find(Vertices.begin(), Vertices.end(), std::vector<int>{nx, ny+1, L-nx-ny-1});
            int index1 = std::distance(Vertices.begin(), i1);
            int index2 = std::distance(Vertices.begin(), i2);
            int index3 = std::distance(Vertices.begin(), i3);

            arr.push_back({index1, index2, index3});
            ny+=1; 
        }
        else{
            auto i1 = std::find(Vertices.begin(), Vertices.end(), std::vector<int>{nx, ny, L-nx-ny});
            auto i2 = std::find(Vertices.begin(), Vertices.end(), std::vector<int>{nx+1, ny, L-nx-ny-1});
            auto i3 = std::find(Vertices.begin(), Vertices.end(), std::vector<int>{nx+1, ny-1, L-nx-ny});
            int index1 = std::distance(Vertices.begin(), i1);
            int index2 = std::distance(Vertices.begin(), i2);
            int index3 = std::distance(Vertices.begin(), i3);
            arr.push_back({index1, index2, index3});
            ny-=1; 
            nx+=1;
        }
        arrPosCounter +=1; 
    }
        rowcounter+=1; 
    }
  return arr; 
}





std::vector<int> FindOrbitLabels(Tensor2D<int> &Vertices){

  std::vector<int> counter(Vertices.size()); //to keep track of points already labeled
  std::vector<int> classlabels(Vertices.size()); //a list that returns the label for all points 
  int dof_counter = 0; 
  for(int n = 0; n< Vertices.size(); n++){
      if (counter[n]!=0){continue;}
      else{
        classlabels[n]=dof_counter;
        counter[n] = 1; 
        int ncounted = 1; 
        std::vector<int> xi= {Vertices[n][0],Vertices[n][1],Vertices[n][2]};
        std::sort(xi.begin(), xi.end()); 
      // std::cout << xi[0] << ' ' << xi[1] << ' ' << xi[2]<< ' ' << xi[3] << '\n';
       // std::cout << "begin Check....\n";
        for (int n1=0; n1< Vertices.size(); n1++){
            if (ncounted==24){break;}
            if (counter[n1]!=0){continue;}
            else{
            std::vector<int>  xi2 = {Vertices[n1][0],Vertices[n1][1],Vertices[n1][2]}; 
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




void InitializeOrbitInformation(Lattice &lattice)
{
  //Setting up the orbit labels for each point in the lattice
  lattice.VertexLabels = FindOrbitLabels(lattice.Vertices); 

  //Finding the position in Vertices of each orbit representative
  //Orbit rep is conventionally choosen to be the sorted barycentric coordinates
  std::vector<int> OrbitPos; 
  std::set<int> unique_labels(lattice.VertexLabels.begin(), lattice.VertexLabels.end());
  int n_orbit = unique_labels.size(); 
  Tensor2D<int> result; 
  for (int o=0; o<n_orbit; o++){
    auto it = std::find(lattice.VertexLabels.begin(), lattice.VertexLabels.end(), o);
    int index = std::distance(lattice.VertexLabels.begin(), it);
    std::vector<int> pos = {lattice.Vertices[index][0],lattice.Vertices[index][1],lattice.Vertices[index][2]};
    std::sort(pos.begin(), pos.end()); 
    auto it2 = std::find(lattice.Vertices.begin(), lattice.Vertices.end(), pos);
    OrbitPos.push_back(std::distance(lattice.Vertices.begin(), it2));
  }
  lattice.OrbitPos = OrbitPos; 

  // Finding the basis

  Tensor2D<int> results; //Basis 
  int dof = 0; 
  for (int i=0; i<n_orbit; i++)
  {
    std::vector<int> vec = lattice.Vertices[lattice.OrbitPos[i]];
    auto it = std::find_if(vec.begin(), vec.end(), [](int i) { return i != 0; });
    int distance = std::distance(vec.begin(), it);
    if (distance!=0)
    {
      if (distance==2){results.push_back({-1, 0});}
      else
      {  
      std::set<int> unique_labels(it, vec.end());
      int n_distinct = unique_labels.size(); 
      if (n_distinct==1){results.push_back({-1, 0});}
      else
      {
        results.push_back({dof, 1});
        dof+=1; 
      }

      }
    }
    else
    {
      std::set<int> unique_labels(vec.begin(), vec.end());
      int n_distinct = unique_labels.size(); 
      if (n_distinct==1){results.push_back({-1, 0});}
      else if (n_distinct==2){results.push_back({dof, 1});dof+=1;}
      else{results.push_back({dof,2});dof+=2;}
    }
  }
  lattice.ndof = dof; 
  lattice.Basis = results; 
}




Lattice GenerateLattice(int L, std::vector<Eigen::Vector3d>& r, int q)
{
  Lattice lattice; 
  lattice.q = q; 
  Tensor2D<int> Vertices; 
  std::vector<Eigen::Vector3d> rvec;   
  std::vector<Eigen::Vector3d> xvec;  
  std::vector<Eigen::Vector3d> xivec; 

   for(int ny = 0; ny<= L; ny++) //Store Vertex information
    {
      for(int nx = 0; nx <= L; nx++)
      {   
        for (int nz = 0; nz<=L; nz++)
        {
          if(nx+ny+nz==L)
          {
            Vertices.push_back({nx, ny, nz}); 
          }
        }
      }
    }
  lattice.Vertices = Vertices; 
  double R = r[0].norm(); 
  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    Eigen::Vector3d vertex(Vertices[i][0], Vertices[i][1],Vertices[i][2]);
    Eigen::Vector3d xpos  = Vertices[i][0]*r[0]+ Vertices[i][1]*r[1]+ Vertices[i][2]*r[2];
    vertex/=L;
    xpos/=L;
    Eigen::Vector3d rpos  = R*xpos/xpos.norm(); 
    xivec.push_back(vertex); 
    xvec.push_back(xpos);
    rvec.push_back(rpos);
  }
  lattice.xivec = xivec; 
  lattice.xvec  = xvec; 
  lattice.rvec  = rvec; 
  InitializeOrbitInformation(lattice); 
  return lattice; 
}



//Not helpful since I realized that even "identical triangles" can have distinct area depending on placement of the vertices.
std::vector<Triangle> FindTriangleRepresentative(Lattice &lattice, int L){
  std::vector<Triangle> result; 
  Tensor2D<int> tmp; 
  Tensor2D<int> TriangleList = TriangleVertexMapping(lattice.Vertices, L); 
  std::vector<int> PointLabels = lattice.VertexLabels; 
  for (std::vector<int> TT: TriangleList)
  {
    int v1 = PointLabels[TT[0]];
    int v2 = PointLabels[TT[1]];
    int v3 = PointLabels[TT[2]];
    std::vector<int> elem = {v1,v2,v3};  //Triangle written in terms of point representatives  
    std::vector<int> original = elem; 
    std::sort(elem.begin(), elem.end()); 
    std::vector<int> indices = MatchIndices(elem, original); 

    auto it = std::find(tmp.begin(), tmp.end(), elem);
    if (it==tmp.end()){
      Triangle A = {{TT[indices[0]], TT[indices[1]], TT[indices[2]]}, elem, 1};
      tmp.push_back(elem);
      result.push_back(A);  
    }
    else{
      int index = std::distance(tmp.begin(), it);
      result[index].wt+=1;
    }
  }
  return result; 
}


Tensor2D<Eigen::Vector3d> FindDerivative(std::vector<Eigen::Vector3d>& r, Lattice &lattice)
{
  Tensor2D<Eigen::Vector3d> results; 
  std::vector<int> OrbitLabels = lattice.VertexLabels;

  for (size_t i=0; i< lattice.Vertices.size(); i++)
  {
    Eigen::Vector3d nd = {0.0, 0.0, 0.0}; 
    int orbit_pos = lattice.OrbitPos[lattice.VertexLabels[i]]; 
    std::vector<int> vec = lattice.Vertices[orbit_pos]; 
    auto it = std::find_if(vec.begin(), vec.end(), [](int i) { return i != 0; });
    int distance = std::distance(vec.begin(), it); //Number of steps from the begining to first non-zero

    if (distance!=0)
    {
      std::set<int> unique_labels(it, vec.end());
      int n_distinct = unique_labels.size(); //Number of distinct points for the non-zero part

      if (distance == 2){results.push_back({nd, nd});}//We are at a vertex of the triangle
      else if (distance==1)//We are at an edge of the triangle
      { 
        if (n_distinct==1){results.push_back({nd, nd});}//We are at the midpoint of the edge
        else
        {
        std::vector<int> Indices = MatchIndices(vec, lattice.Vertices[i]); //vec[j] = Vertices[i][Indices[j]]
        results.push_back({r[Indices[1]]-r[Indices[2]], nd});
        }
      }  
    }
    else//Interior of triangle
    {
      std::vector<int> Indices = MatchIndices(vec, lattice.Vertices[i]);
      std::set<int> unique_labels(vec.begin(), vec.end());
      int n_distinct = unique_labels.size(); 
      if (n_distinct==1){ results.push_back({nd, nd}); } //Skip at center of tetrahedron
      else if (n_distinct==2) // (a,a,b) or (a, b, b) patterns
      {
        if (vec[0]==vec[1]) // (a,a,b)
        {
          results.push_back({r[Indices[0]]+r[Indices[1]]-2*r[Indices[2]], nd});
        }
        else //(a, b, b)
        {
          results.push_back({-2*r[Indices[0]]+r[Indices[1]]+r[Indices[2]], nd});
        }
      }
      else//Interior point with two d.o.f
      {
        results.push_back({ r[Indices[1]]-r[Indices[0]], r[Indices[2]]-r[Indices[0]]});
      }

    }

  }
  return results; 
}


SpMat AreaOperator(Lattice &lattice, Tensor2D<Eigen::Vector3d> &DerivativeList,Tensor2D<int> &TList, double R)
{
  std::vector<T> tripletList; 
  Tensor2D<int> basis = lattice.Basis; 
  int nT = 0; 
  for (std::vector<int> TT: TList)//iterate over triangles, to be weighed by multiplicities
  {
    //printf("Triange is: %d %d %d\n", TT[0],TT[1],TT[2]);
    for (int i = 0; i<3; i++)
    //Iterate over all the vertices in the triangle
    {
      int ind = lattice.VertexLabels[TT[i]]; 
      int indV = TT[i]; 
      //printf("Considering Verx %d with label %d and basis %d\n", indV, ind, basis[ind][0]); 
      if (basis[ind][0]==-1){
        continue; }
      else
      {
        std::vector<Eigen::Vector3d> r = {lattice.rvec[TT[0]], lattice.rvec[TT[1]], lattice.rvec[TT[2]]};
        std::rotate(r.begin(), r.begin()+i, r.end());
        Eigen::Vector3d PointGradient = VertexDerivative(r[0], r[1], r[2]); 

        //std::cout<< "PointGradient is: "<<PointGradient<<"\n";
        Eigen::Vector3d XiGradient = DerivativeList[indV][0]/lattice.xvec[indV].norm();
        XiGradient-=lattice.xvec[indV]/pow(lattice.xvec[indV].norm(), 3)*(DerivativeList[indV][0].dot(lattice.xvec[indV]));
        //printf("target derivative is %.12f\n", DerivativeList[indV][0].dot(PointGradient)); 
        tripletList.push_back(Eigen::Triplet<double>(basis[ind][0], nT, R*XiGradient.dot(PointGradient))); 
        if(basis[ind][1]==2)
        {
        XiGradient = DerivativeList[indV][1]/lattice.xvec[indV].norm();
        XiGradient-= lattice.xvec[indV]/pow(lattice.xvec[indV].norm(), 3)*(DerivativeList[indV][1].dot(lattice.xvec[indV]));
        tripletList.push_back(Eigen::Triplet<double>(basis[ind][0]+1, nT, R*XiGradient.dot(PointGradient))); 
      
        }
      }
    }
    nT+=1; 
  }

  SpMat sparseMat(lattice.ndof, TList.size());
  sparseMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return sparseMat; 
}

Eigen::VectorXd returnCurrentArea(Lattice &lattice, Tensor2D<int> &TList)
{
  std::vector<Eigen::Vector3d> rvec = lattice.rvec; 
  Eigen::VectorXd result(TList.size()); 
  int n = 0; 
  for (std::vector<int> TT: TList)
  { 
    Eigen::Vector3d r1 = rvec[TT[0]];
    Eigen::Vector3d r2 = rvec[TT[1]];
    Eigen::Vector3d r3 = rvec[TT[2]];
    result[n] = TriangleArea(r1, r2, r3);
    n++; 
  }
  return result; 
}

Eigen::VectorXd returnCurrentDualArea(Lattice &lattice, Tensor2D<int> &TList, int L)
{
  std::vector<double> result(lattice.Vertices.size(), 0.0); 
  int q = lattice.q; 
  
  for (int i=0;i<TList.size();i++)
  {
    std::vector<Eigen::Vector3d> r = {lattice.rvec[TList[i][0]],lattice.rvec[TList[i][1]], lattice.rvec[TList[i][2]]};
    std::vector<Eigen::Vector3d> edges = {r[0]-r[1], r[0]-r[2], r[1]-r[2]};

    double l1 = edges[0].dot(edges[0]); 
    double l2 = edges[1].dot(edges[1]); 
    double l3 = edges[2].dot(edges[2]); 
    double l1s = Compute_ls(l1, l2, l3, l1);
    double l2s = Compute_ls(l1, l2, l3, l2);
    double l3s = Compute_ls(l1, l2, l3, l3);
    result[TList[i][0]] += (sqrt(l1)*l1s+sqrt(l2)*l2s)/2;
    result[TList[i][1]] += (sqrt(l1)*l1s+sqrt(l3)*l3s)/2;
    result[TList[i][2]] += (sqrt(l2)*l2s+sqrt(l3)*l3s)/2;

  }
  for (int i=0; i<lattice.Vertices.size(); i++)
  {
    int nx = lattice.Vertices[i][0]; 
    int ny = lattice.Vertices[i][1]; 
    int nz = lattice.Vertices[i][2];
    if (nx==L||ny==L||nz==L){result[i]*=double(q);}
    else if (nx==0||ny==0||nz==0){result[i]*=2; }
  }
  Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>(result.data(), result.size());
  return vec;
}

Eigen::VectorXd returnCurrentDeficit(Lattice &lattice, Tensor2D<int> &TList, int L)
{
  std::vector<double> result(lattice.Vertices.size(), 0.0); 
  int q= lattice.q;
  for (int i=0; i<TList.size(); i++)
  {
    std::vector<Eigen::Vector3d> r = {lattice.rvec[TList[i][0]],lattice.rvec[TList[i][1]], lattice.rvec[TList[i][2]]};
    std::vector<Eigen::Vector3d> edges = {r[0]-r[1], r[0]-r[2], r[1]-r[2]};

    double l1 = edges[0].norm(); 
    double l2 = edges[1].norm(); 
    double l3 = edges[2].norm();
    double p = l1+l2+l3; 
    double ratio; 

    ratio = (p-2*l1)*(p-2*l2)/(p*(p-2*l3));
    result[TList[i][0]] -= 2.0*atan(sqrt(ratio));
    ratio = (p-2*l1)*(p-2*l3)/(p*(p-2*l2));
    result[TList[i][1]] -=2.0*atan(sqrt(ratio));
    ratio = (p-2*l2)*(p-2*l3)/(p*(p-2*l1));
    result[TList[i][2]] -=2.0*atan(sqrt(ratio));
  }
  for (int i=0; i<lattice.Vertices.size(); i++)
  {
    int nx = lattice.Vertices[i][0]; 
    int ny = lattice.Vertices[i][1]; 
    int nz = lattice.Vertices[i][2];
    if (nx==L||ny==L||nz==L){result[i]*=double(q);}
    else if (nx==0||ny==0||nz==0){result[i]*=2; }
    result[i]+=TWOPI;
  }

  Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>(result.data(), result.size());
  return vec;
}


double EffectiveLatticeSpacing(Lattice &lattice, Tensor2D<int> &TList, int L)
{
  std::vector<Eigen::Vector3d> rvec = lattice.rvec; 
  for (int i=0; i<rvec.size(); i++)
  {
    rvec[i]/=rvec[i].norm();
  }

  std::vector<double> result(lattice.Vertices.size(), 0.0); 
  int q= lattice.q;
  for (int i=0; i<TList.size(); i++)
  {
    std::vector<Eigen::Vector3d> r = {rvec[TList[i][0]],rvec[TList[i][1]], rvec[TList[i][2]]};
    std::vector<Eigen::Vector3d> edges = {r[0]-r[1], r[0]-r[2], r[1]-r[2]};

    double l1 = edges[0].norm(); 
    double l2 = edges[1].norm(); 
    double l3 = edges[2].norm();
    double p = l1+l2+l3; 
    double ratio; 

    ratio = (p-2*l1)*(p-2*l2)/(p*(p-2*l3));
    result[TList[i][0]] -= 2.0*atan(sqrt(ratio));
    ratio = (p-2*l1)*(p-2*l3)/(p*(p-2*l2));
    result[TList[i][1]] -=2.0*atan(sqrt(ratio));
    ratio = (p-2*l2)*(p-2*l3)/(p*(p-2*l1));
    result[TList[i][2]] -=2.0*atan(sqrt(ratio));
  }
  for (int i=0; i<lattice.Vertices.size(); i++)
  {
    int nx = lattice.Vertices[i][0]; 
    int ny = lattice.Vertices[i][1]; 
    int nz = lattice.Vertices[i][2];
    if (nx==L||ny==L||nz==L){result[i]*=double(q);}
    else if (nx==0||ny==0||nz==0){result[i]*=2; }
    result[i]+=TWOPI;
  }

  Eigen::VectorXd D = Eigen::Map<Eigen::VectorXd>(result.data(), result.size());

  double dmean = 0.0;  
  //Multiplicity of points on the traingle is different depending on the polyhedra
  std::vector<double> factors(3); 
  int nVert; 
  if (lattice.q==5)
  {factors[0]=4.0; factors[1]=10.0;factors[2]=20.0; nVert=2+10*L*L;}
  else if (lattice.q ==4)
  {factors[0]=2.0;  factors[1]=4.0; factors[2]=8.0; nVert=2+4*L*L;}
  else
  {factors[0]=4/3;  factors[1]=2.0; factors[2]=4.0; nVert=2+2*L*L;}

  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    if (lattice.Vertices[i][0]==L|| lattice.Vertices[i][1]==L||lattice.Vertices[i][2]==L)
    {dmean += factors[0]*D[i];}
    else if  (lattice.Vertices[i][0]==0|| lattice.Vertices[i][1]==0||lattice.Vertices[i][2]==0)
    {dmean += factors[1]*D[i];}
    else
    {dmean += factors[2]*D[i];}
  }
  dmean/=nVert;  
  return dmean; 
}

Eigen::VectorXd setUpGMRESsolver(SpMat& Operator, Eigen::VectorXd& val,  int iter_num){ 
    Eigen::VectorXd sol(Operator.rows());
    SpMat A = -1*Operator; 
    SpMat AtA = A*(A.transpose());
    Eigen::VectorXd target = A*val;
    Eigen::GMRES<SpMat> solver; 
    solver.compute(AtA);
    solver.setMaxIterations(iter_num);
    sol = solver.solve(target);
    return sol; 
}

void UpdateLattice(Lattice &lattice, std::vector<Eigen::Vector3d> &r, Eigen::VectorXd& update)
{
  //Update the primary sites
  Tensor2D<int> Basis = lattice.Basis; 
  for (int i =0; i<lattice.Basis.size(); i++)
  {
    if (Basis[i][0]==-1){continue; }
    else
    {
      int n = lattice.OrbitPos[i];
      std::vector<int> orbit = lattice.Vertices[n];
      if (Basis[i][1]==2){
        Eigen::Vector3d updatevec = {-update[Basis[i][0]]-update[Basis[i][0]+1], update[Basis[i][0]],update[Basis[i][0]+1]}; 
        lattice.xivec[n]+= updatevec; }
      else
      {
        if(orbit[0]==0)//We are at an edge
        {
        //std::cout<<"Before Update "<<lattice.xivec[n]<<"\n"; 
        Eigen::Vector3d updatevec = {0,update[Basis[i][0]],-update[Basis[i][0]]}; 
        lattice.xivec[n]+= updatevec;
        //std::cout<<"After Update "<<lattice.xivec[n]<<"\n"; 
        }
        else //We are at an equal xi_i xi_j line in the interior
        {
          if (orbit[0]==orbit[1])
          {
          Eigen::Vector3d updatevec = {update[Basis[i][0]],update[Basis[i][0]],-2*update[Basis[i][0]]}; 
          lattice.xivec[n]+= updatevec;
          }
          else
          {
          Eigen::Vector3d updatevec = {-2*update[Basis[i][0]],update[Basis[i][0]],update[Basis[i][0]]}; 
          lattice.xivec[n]+= updatevec;             
          }
        }
      }
    }
  }
#if 1
  //Update all the sites using the primary sites
  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    int rep = lattice.OrbitPos[lattice.VertexLabels[i]]; //Position of the Corresponding Orbit Rep 
    if (rep==i){continue;} //We are at a primary site
    else
    {
      std::vector<int> indices = MatchIndices(lattice.Vertices[i],lattice.Vertices[rep]);
      lattice.xivec[i] = {lattice.xivec[rep][indices[0]],lattice.xivec[rep][indices[1]],lattice.xivec[rep][indices[2]]}; 
    }  
  }
  //Update xvec and rvec
  double R = r[0].norm(); 
  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    lattice.xvec[i] = r[0]*lattice.xivec[i][0]+r[1]*lattice.xivec[i][1]+r[2]*lattice.xivec[i][2]; 
    lattice.rvec[i] = R*lattice.xvec[i]/lattice.xvec[i].norm(); 
  }
#endif
}

Eigen::VectorXd returnGradient(SpMat& Operator, Eigen::VectorXd &Area)
{
  return 2*(Operator*Area).rowwise().sum(); 
}


Eigen::VectorXd UpdateArea(Lattice &lattice, std::vector<Eigen::Vector3d> &r, Eigen::VectorXd& update, Tensor2D<int> &TList)
{
  std::vector<Eigen::Vector3d> xivec = lattice.xivec; 
  std::vector<Eigen::Vector3d> xvec = lattice.xvec; 
  std::vector<Eigen::Vector3d> rvec = lattice.rvec; 
  Tensor2D<int> Basis = lattice.Basis; 
  for (int i =0; i<lattice.Basis.size(); i++)
  {
    if (Basis[i][0]==-1){continue; }
    else
    {
      int n = lattice.OrbitPos[i];
      std::vector<int> orbit = lattice.Vertices[n];
      if (Basis[i][1]==2){
        Eigen::Vector3d updatevec = {-update[Basis[i][0]]-update[Basis[i][0]+1], update[Basis[i][0]],update[Basis[i][0]+1]}; 
        xivec[n]+= updatevec; }
      else
      {
        if(orbit[0]==0)//We are at an edge
        {
        //std::cout<<"Before Update "<<lattice.xivec[n]<<"\n"; 
        Eigen::Vector3d updatevec = {0,update[Basis[i][0]],-update[Basis[i][0]]}; 
        xivec[n]+= updatevec;
        //std::cout<<"After Update "<<lattice.xivec[n]<<"\n"; 
        }
        else //We are at an equal xi_i xi_j line in the interior
        {
          if (orbit[0]==orbit[1])
          {
          Eigen::Vector3d updatevec = {update[Basis[i][0]],update[Basis[i][0]],-2*update[Basis[i][0]]}; 
          xivec[n]+= updatevec;
          }
          else
          {
          Eigen::Vector3d updatevec = {-2*update[Basis[i][0]],update[Basis[i][0]],update[Basis[i][0]]}; 
          xivec[n]+= updatevec;             
          }
        }
      }
    }
  }
  //Update all the sites using the primary sites
  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    int rep = lattice.OrbitPos[lattice.VertexLabels[i]]; //Position of the Corresponding Orbit Rep 
    if (rep==i){continue;} //We are at a primary site
    else
    {
      std::vector<int> indices = MatchIndices(lattice.Vertices[i],lattice.Vertices[rep]);
      xivec[i] = {xivec[rep][indices[0]],xivec[rep][indices[1]],xivec[rep][indices[2]]}; 
    }  
  }
  //Update xvec and rvec
  double R = r[0].norm(); 
  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    xvec[i] = r[0]*xivec[i][0]+r[1]*xivec[i][1]+r[2]*xivec[i][2]; 
    rvec[i] = R*xvec[i]/xvec[i].norm(); 
  }

  Eigen::VectorXd result(TList.size()); 
  int n = 0; 
  for (std::vector<int> TT: TList)
  { 
    Eigen::Vector3d r1 = rvec[TT[0]];
    Eigen::Vector3d r2 = rvec[TT[1]];
    Eigen::Vector3d r3 = rvec[TT[2]];
    result[n] = TriangleArea(r1, r2, r3);
    n++; 
  }
  return result; 
}


double BacktrackingLineSearch(
 Lattice &lattice,
 Tensor2D<int> &TList,
 Eigen::VectorXd & Gradient, std::vector<Eigen::Vector3d> &r, Eigen::VectorXd& update, double epsilon_max)
{
  Eigen::VectorXd avec = returnCurrentArea(lattice, TList);
  double AG_condition= -update.dot(Gradient)/2.0;
  double ratio = 0.8; //This ratio technically can be finetuned too
  double action = avec.dot(avec); 
  int condition = 0;   
  int counter = 0; 
  while (condition!=1){
      double epsilon = epsilon_max*pow(ratio, counter); 
      Eigen::VectorXd testupdate = update*epsilon; 
      Eigen::VectorXd newArea = UpdateArea(lattice, r, testupdate, TList);
      if (action-newArea.dot(newArea)>=AG_condition*epsilon){           
        epsilon_max = epsilon;
        condition=1;
        }
        else{
          counter+=1; 
        }
    } 
  return epsilon_max;
}


void PrintGeometry(Lattice &lattice, Tensor2D<int> &TList, int L)
{
  Eigen::VectorXd P(TList.size()); 
  Eigen::VectorXd Cr(TList.size()); 
  for (int i=0; i<TList.size(); i++)
  {
    Eigen::Vector3d r1 = lattice.rvec[TList[i][0]];
    Eigen::Vector3d r2 = lattice.rvec[TList[i][1]];
    Eigen::Vector3d r3 = lattice.rvec[TList[i][2]];
    P[i] = Perimeter(r1,r2,r3); 
    Cr[i] = Circumradius(r1,r2,r3); 
  }
  Eigen::VectorXd A = returnCurrentArea(lattice, TList); 
  Eigen::VectorXd D = returnCurrentDualArea(lattice, TList, L);
  Eigen::VectorXd DA = returnCurrentDeficit(lattice, TList, L);

  //Calculate the average quantities across the whole sphere
  double dmean = 0.0;  
  double dmean_sq = 0.0; 
  double defmean = 0.0; 
  double defmean_sq = 0.0; 

  //Multiplicity of points on the traingle is different depending on the polyhedra
  std::vector<double> factors(3); 
  int nVert; 
  if (lattice.q==5)
  {factors[0]=4.0; factors[1]=10.0;factors[2]=20.0; nVert=2+10*L*L;}
  else if (lattice.q ==4)
  {factors[0]=2.0;  factors[1]=4.0; factors[2]=8.0; nVert=2+4*L*L;}
  else
  {factors[0]=4/3;  factors[1]=2.0; factors[2]=4.0; nVert=2+2*L*L;}

  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    if (lattice.Vertices[i][0]==L|| lattice.Vertices[i][1]==L||lattice.Vertices[i][2]==L)
    {
      dmean      += factors[0]*D[i]; 
      dmean_sq   += factors[0]*D[i]*D[i]; 
      defmean    += factors[0]*DA[i]; 
      defmean_sq += factors[0]*DA[i]*DA[i]; 
    }
    else if  (lattice.Vertices[i][0]==0|| lattice.Vertices[i][1]==0||lattice.Vertices[i][2]==0)
    {
      dmean      += factors[1]*D[i]; 
      dmean_sq   += factors[1]*D[i]*D[i]; 
      defmean    += factors[1]*DA[i]; 
      defmean_sq += factors[1]*DA[i]*DA[i]; 
    }
    else
    {
      dmean      += factors[2]*D[i]; 
      dmean_sq   += factors[2]*D[i]*D[i]; 
      defmean    += factors[2]*DA[i]; 
      defmean_sq += factors[2]*DA[i]*DA[i]; 
    }
  }
  
  dmean/=nVert;  
  dmean_sq/=nVert; 
  defmean/=nVert; 
  defmean_sq/=nVert; 
  double dRMS = dmean_sq-dmean*dmean; 
  double defRMS = defmean_sq-defmean*defmean; 



  printf("%.16f %.16f %.16f %.16f %.16f ",
  A.array().mean(), P.array().mean(), Cr.array().mean(), dmean, defmean);
  printf("%.16f %.16f %.16f %.16f %.16f ",
  A.array().square().mean(), P.array().square().mean(), Cr.array().square().mean(), 
  dmean_sq, defmean_sq);
   printf("%.16f %.16f %.16f %.16f %.16f\n", returnRMS(A), returnRMS(P), returnRMS(Cr),
   dRMS, defRMS); 
}

void PrintRenormalizedGeometry(Lattice &lattice, Tensor2D<int> &TList, int L)
{
  Eigen::VectorXd P(TList.size()); 
  Eigen::VectorXd Cr(TList.size()); 
  for (int i=0; i<TList.size(); i++)
  {
    Eigen::Vector3d r1 = lattice.rvec[TList[i][0]];
    Eigen::Vector3d r2 = lattice.rvec[TList[i][1]];
    Eigen::Vector3d r3 = lattice.rvec[TList[i][2]];
    P[i] = Perimeter(r1,r2,r3); 
    Cr[i] = Circumradius(r1,r2,r3); 
  }
  Eigen::VectorXd A = returnCurrentArea(lattice, TList); 
  Eigen::VectorXd D = returnCurrentDualArea(lattice, TList, L);
  Eigen::VectorXd DA = returnCurrentDeficit(lattice, TList, L);

  //Calculate the average quantities across the whole sphere
  double dmean = 0.0;  
  double dmean_sq = 0.0; 
  double defmean = 0.0; 
  double defmean_sq = 0.0; 

  //Multiplicity of points on the traingle is different depending on the polyhedra
  std::vector<double> factors(3); 
  int nVert; 
  if (lattice.q==5)
  {factors[0]=4.0; factors[1]=10.0;factors[2]=20.0; nVert=2+10*L*L;}
  else if (lattice.q ==4)
  {factors[0]=2.0;  factors[1]=4.0; factors[2]=8.0; nVert=2+4*L*L;}
  else
  {factors[0]=4/3;  factors[1]=2.0; factors[2]=4.0; nVert=2+2*L*L;}

  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    if (lattice.Vertices[i][0]==L|| lattice.Vertices[i][1]==L||lattice.Vertices[i][2]==L)
    {
      dmean      += factors[0]*D[i]; 
      dmean_sq   += factors[0]*D[i]*D[i]; 
      defmean    += factors[0]*DA[i]; 
      defmean_sq += factors[0]*DA[i]*DA[i]; 
    }
    else if  (lattice.Vertices[i][0]==0|| lattice.Vertices[i][1]==0||lattice.Vertices[i][2]==0)
    {
      dmean      += factors[1]*D[i]; 
      dmean_sq   += factors[1]*D[i]*D[i]; 
      defmean    += factors[1]*DA[i]; 
      defmean_sq += factors[1]*DA[i]*DA[i]; 
    }
    else
    {
      dmean      += factors[2]*D[i]; 
      dmean_sq   += factors[2]*D[i]*D[i]; 
      defmean    += factors[2]*DA[i]; 
      defmean_sq += factors[2]*DA[i]*DA[i]; 
    }
  }
  
  dmean/=nVert;  
  dmean_sq/=nVert; 
  defmean/=nVert; 
  defmean_sq/=nVert; 
  double dRMS = dmean_sq-dmean*dmean; 
  double defRMS = defmean_sq-defmean*defmean; 



  printf("%.16f %.16f %.16f %.16f %.16f ",
  A.array().mean(), P.array().mean(), Cr.array().mean(), dmean, defmean);
  printf("%.16f %.16f %.16f %.16f %.16f ",
  A.array().square().mean(), P.array().square().mean(), Cr.array().square().mean(), 
  dmean_sq, defmean_sq);
   printf("%.16f %.16f %.16f %.16f %.16f\n", returnRMS(A)/pow(A.array().mean(),2), returnRMS(P)/pow(P.array().mean(),2), returnRMS(Cr)/pow(Cr.array().mean(),2),
   dRMS/pow(dmean,2), defRMS/pow(defmean,2)); 
}

void SaveBarycentric(Lattice &lattice, FILE* file)
{
  for (int i=0; i<lattice.OrbitPos.size();i++)
  {
  std::vector<double> xi = {lattice.xivec[lattice.OrbitPos[i]][0],lattice.xivec[lattice.OrbitPos[i]][1],lattice.xivec[lattice.OrbitPos[i]][2]};
   std::sort(xi.begin(), xi.end(), std::greater<double>());
   fprintf(file, "%d %.16f %.16f\n", i, xi[0], xi[1]);
  }
}
void SaveTriangle(Tensor2D<int> &TList, FILE* file)
{
  for (int i=0; i<TList.size();i++)
  {
   fprintf(file, "%d %d %d\n",TList[i][0], TList[i][1], TList[i][2]);
  }
}

void SavePosition(Lattice &lattice, FILE* file)
{
  for (int i=0; i<lattice.Vertices.size();i++)
  {
   fprintf(file, "%d %.16f %.16f %.16f\n", i, lattice.rvec[i][0], lattice.rvec[i][1], lattice.rvec[i][2]);
  }
}

void ReadBarycentric(Lattice &lattice, const std::string& filePath, double R, std::vector<Eigen::Vector3d> &r)
{
  std::ifstream file(filePath);
  std::string line;
  int row = 0;
  for (int classn=0; classn<lattice.OrbitPos.size(); classn++)
  {
    if (!getline(file, line))
    {
      std::cerr << "Error reading line for class number " << classn << std::endl;
      return;
    }
    std::istringstream iss(line);
    int orbitID; 
    iss>>orbitID; 
    for (int i=0; i<lattice.OrbitPos.size(); i++)
    {   
        if (i==orbitID){
        double xi1, xi2, xi3; 
        iss>>xi1; 
        iss>>xi2;
        xi3 = 1-xi1-xi2; 
        //printf("%.12f %.12f %.12f\n", xi1, xi2, xi3); 
        std::vector<double> xi = {xi1, xi2, xi3}; 
        std::sort(xi.begin(), xi.end()); 
        int o = lattice.OrbitPos[i];
        Eigen::Vector3d xivec = {xi[0], xi[1], xi[2]};
        lattice.xivec[o] =  xivec; }
    }
  }
  for (int i=0; i<lattice.Vertices.size(); i++)
  {
    int o = lattice.OrbitPos[lattice.VertexLabels[i]];
    if (o==i){continue;}
    else
    {
      std::vector<int> indices = MatchIndices(lattice.Vertices[i],lattice.Vertices[o]);
      Eigen::Vector3d xivec = {lattice.xivec[o][indices[0]], lattice.xivec[o][indices[1]], lattice.xivec[o][indices[2]]};
      lattice.xivec[i] = xivec;
    }
  }
  #if 1
  for (int i =0; i<lattice.Vertices.size(); i++)
  {
    lattice.xvec[i] = r[0]*lattice.xivec[i][0]+r[1]*lattice.xivec[i][1]+r[2]*lattice.xivec[i][2]; 
    lattice.rvec[i] = R*lattice.xvec[i]/lattice.xvec[i].norm(); 
  }
  #endif 
}