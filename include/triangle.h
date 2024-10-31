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
  std::vector<Eigen::Vector3d> rvec;   //size = L
  std::vector<Eigen::Vector3d> xvec;   //size = L
  std::vector<Eigen::Vector3d> xivec;  //size = L
  std::vector<int> VertexLabels;       //size = L given by FindOrbitLabels
  std::vector<int> OrbitPos;           //size = Number of Orbits
  Tensor2D<int> Basis;                 //size = (Number of Orbits, 2)
  int ndof;                            //Number of Orbits
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

//double CHECK!!
inline double Compute_ls(double a, double b, double c, double edge_sq){
  double l_sq_sum = a+b+c;
  double A = sqrt(2*a*b + 2*b*c + 2*a*c - a*a - b*b - c*c)/4;
  return sqrt(edge_sq)*(l_sq_sum-2.0*edge_sq)/(8*A);
}
double TriangleArea(Eigen::Vector3d& r1, Eigen::Vector3d& r2, Eigen::Vector3d& r3)
{
  Eigen::Vector3d e1 = r3-r1;
  Eigen::Vector3d e2 = r2-r1;  
  return std::abs((e1.cross(e2)).norm())/2.0; 
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




Lattice GenerateLattice(int L, std::vector<Eigen::Vector3d>& r)
{
  Lattice lattice; 
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
    Eigen::Vector3d rpos  = xpos; 
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
        tripletList.push_back(Eigen::Triplet<double>(basis[ind][0], nT, DerivativeList[indV][0].dot(PointGradient))); 
        if(basis[ind][1]==2)
        {
        XiGradient = DerivativeList[indV][1]/lattice.xvec[indV].norm();
        XiGradient-= lattice.xvec[indV]/pow(lattice.xvec[indV].norm(), 3)*(DerivativeList[indV][1].dot(lattice.xvec[indV]));
        tripletList.push_back(Eigen::Triplet<double>(basis[ind][0]+1, nT, DerivativeList[indV][1].dot(PointGradient))); 
      
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


Eigen::VectorXd setUpGMRESsolver(SpMat& Operator, Eigen::VectorXd& Area,  int iter_num){ 
    Eigen::VectorXd sol(Operator.rows());
    SpMat A = -1*Operator; 
    SpMat AtA = A*(A.transpose());
    Eigen::VectorXd target = A*Area;
    std::cout<<"target: "<<target<<"\n"; 
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
    lattice.rvec[i] = lattice.xvec[i]; 
  }
#endif
}







//Old Code. maybe not needed?

// SpMat AreaOperator(Lattice &lattice, Tensor2D<Eigen::Vector3d> &DerivativeList,std::vector<Triangle> &TList, double R)
// {
//   std::vector<T> tripletList; 

//   std::vector<Triangle> FindTriangleRepresentative(Lattice &lattice, int L)
//   Tensor2D<int> basis = lattice.Basis; 
//   int nT = 0; 
//   for (Triangle TT: TList)//iterate over triangles, to be weighed by multiplicities
//   {
//     //printf("Triange is: %d %d %d\n", TT.Vert[0],TT.Vert[1],TT.Vert[2]);
//     for (int i = 0; i<3; i++)
//     //Iterate over all the vertices in the triangle
//     {

//       int ind = lattice.VertexLabels[i]; 
//       int indV = TT.Vert[i]; 
//       if (basis[ind][0]==-1){continue; }
//       else
//       {
//         //printf("Considering Verx %d....\n", indV); 
//         std::vector<Eigen::Vector3d> r = {lattice.rvec[TT.Vert[0]], lattice.rvec[TT.Vert[1]], lattice.rvec[TT.Vert[2]]};
//         std::rotate(r.begin(), r.begin()+i, r.end());
//         Eigen::Vector3d PointGradient = VertexDerivative(r[0], r[1], r[2]); 

//         //std::cout<< "PointGradient is: "<<PointGradient<<"\n";
//         Eigen::Vector3d XiGradient = DerivativeList[indV][0]/lattice.xvec[indV].norm();
//         XiGradient-=lattice.xvec[indV]/pow(lattice.xvec[indV].norm(), 3)*(DerivativeList[indV][0].dot(lattice.xvec[indV]));
//         //printf("target derivative is %.12f\n", DerivativeList[indV][0].dot(PointGradient)); 
//         tripletList.push_back(Eigen::Triplet<double>(basis[ind][0], nT, sqrt(TT.wt)*DerivativeList[indV][0].dot(PointGradient))); 
//         if(basis[ind][1]==2)
//         {
//         XiGradient = DerivativeList[indV][1]/lattice.xvec[indV].norm();
//         XiGradient-= lattice.xvec[indV]/pow(lattice.xvec[indV].norm(), 3)*(DerivativeList[indV][1].dot(lattice.xvec[indV]));
//         tripletList.push_back(Eigen::Triplet<double>(basis[ind][0]+1, nT, sqrt(TT.wt)*DerivativeList[indV][1].dot(PointGradient))); 
      
//         }
//       }
//     }
//     nT+=1; 
//   }

//   SpMat sparseMat(lattice.ndof, TList.size());
//   sparseMat.setFromTriplets(tripletList.begin(), tripletList.end());
//   return sparseMat; 
// }