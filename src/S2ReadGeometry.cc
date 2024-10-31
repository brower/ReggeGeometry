#include "S2.h"
using namespace std;

int main(int argc, char *argv[])
{
    int L = 4; 
    int q = 5; 
    std::string data_dir = "";
    if (argc==1)
    { std::cout << "Using default values: L = " << L << ", q = " << q  << std::endl;}
    else{
    if (argc<=3&&argc>1)
    {
    L= atoi(argv[1]);
    q= atoi(argv[2]);
    }
    else
    {
    L= atoi(argv[1]);
    q= atoi(argv[2]);
    data_dir = argv[3];   
    }}

    double z = z = sqrt((7 + 3 * sqrt(5))/8);
    if (q==4){z = 1/sqrt(2); }
    if (q==3){z = 1/(2*sqrt(2)); }
    Eigen::Vector3d r1 =  { sqrt(3)/2.0 , -1/2.0,  z };
    Eigen::Vector3d r2 =  {           0,   1.0 ,   z };
    Eigen::Vector3d r3 =  { -sqrt(3)/2.0, -1/2.0,  z };

    r1*=(L/r1.norm());
    r2*=(L/r2.norm());
    r3*=(L/r3.norm());

    double R= r1.norm(); 
    std::vector<Eigen::Vector3d> r={r1,r2,r3}; 

    Lattice lattice = GenerateLattice(L, r, q);
    Tensor2D<int> TTList = TriangleVertexMapping(lattice.Vertices, L);
    
    if (data_dir.empty())
    {
    //printf("Initialize lattice-----\n"); 
    double alat = EffectiveLatticeSpacing(lattice, TTList, L); 
    printf("%d %.16f ", L, alat); 
    PrintGeometry(lattice, TTList, L);
    }

    else
    {
    //printf("Downloading Orbit-----\n"); 
    ReadBarycentric(lattice, data_dir, double(L), r); 
    double alat = EffectiveLatticeSpacing(lattice, TTList, L); 
    printf("%d %.16f ", L, alat); 
    PrintRenormalizedGeometry(lattice, TTList, L);
    
    size_t pos = data_dir.find_last_of('/');
    std::string Directory;
    if (pos != std::string::npos) {
        // Resize the string to remove the characters after the last backslash
        Directory = data_dir.substr(0, pos);
    }
    
    FILE* testFull = NULL;  
    char test[128];
    sprintf(test, "%s/rvec_%d.dat", Directory.c_str(), L); // filename
    testFull = fopen(test,"w");
    SavePosition(lattice, testFull);
    fclose(testFull);
    
    sprintf(test, "%s/Triangle_%d.dat", Directory.c_str(), L); // filename
    testFull = fopen(test,"w");
    SaveTriangle(TTList, testFull);
    fclose(testFull);
    }

return 0;
}

