#include "triangle.h"


using namespace std;


#define L 4
#define TWOPI  6.283185307179586
#define Root2 1.4142135623730951
#define Two 2
#define Three 3
#define Debug 0
#define Zero 0.0
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

//Implement size as a command line parameter
int main(int argc, char *argv[])
{
    int N_iter=100; 
  
    //double z = sqrt((7 + 3 * sqrt(5))/8); 
    double z = 0; 
    Eigen::Vector3d r1 =  { sqrt(3)/2.0 , -1/2.0,  z };
    Eigen::Vector3d r2 =  {           0,   1.0 ,   z };
    Eigen::Vector3d r3 =  { -sqrt(3)/2.0, -1/2.0,  z };

    // r1*=(L/r1.norm());
    // r2*=(L/r2.norm());
    // r3*=(L/r3.norm());

    double R= r1.norm(); 
    std::vector<Eigen::Vector3d> r={r1,r2,r3}; 

    Lattice lattice = GenerateLattice(L, r);
    std::vector<Triangle> TList = FindTriangleRepresentative(lattice, L); 
    Tensor2D<Eigen::Vector3d> D = FindDerivative(r, lattice); 
    Tensor2D<int> TTList = TriangleVertexMapping(lattice.Vertices, L);
    
    Eigen::VectorXd displace(2);
    displace<< 0.2, 0.1;

    UpdateLattice(lattice, r, displace); 
    for (int n=0; n< N_iter; n++)
    {
    double epsilon =1; 
    
    SpMat Op = AreaOperator(lattice, D, TTList, R); 
    Eigen::VectorXd Area = returnCurrentArea(lattice, TTList); 
    double rms = returnRMS(Area); 
    Eigen::VectorXd sol = setUpGMRESsolver(Op, Area, 1);
    printf("iter = %d, rms = %.12f,grad  = %.12f, action = %.12f\n", n, rms, sol.norm(), Area.norm()); 

    sol*=epsilon;  
    UpdateLattice(lattice, r, sol); 
  
    }
   
  
return 0;
}

