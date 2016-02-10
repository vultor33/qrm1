#include "MatrixRotation.h"

#include<math.h>
#include <vector>

#include "EularAngle.h"
#include "Molecule.h"

using namespace std;

vector< vector<double> > MatrixRotation::CalcRotatingMatrix(const Molecule& mol, int A, int B)
{
	vector< vector<double> > rotatingMatrix;
	int size=4;
	rotatingMatrix.resize(size);
	for(int i=0; i<size;i++)
	{
		rotatingMatrix[i].resize(size);
	}

	double x = mol.MolXYZ[B].x - mol.MolXYZ[A].x;
	double y = mol.MolXYZ[B].y - mol.MolXYZ[A].y;
    double z = mol.MolXYZ[B].z - mol.MolXYZ[A].z;

	EularAngle eularAngle(x, y, z);
	double alpha = eularAngle.GetAlpha();
	double beta  = eularAngle.GetBeta();

	// meu s e s beleza
	// x = 1
	// y = 2
	// z = 3

   // rotating matrix for s-function
   rotatingMatrix[0][0] = 1.0e0;
   rotatingMatrix[0][1] = 0.0e0;
   rotatingMatrix[0][2] = 0.0e0;
   rotatingMatrix[0][3] = 0.0e0;
   rotatingMatrix[1][0] = 0.0e0;
   rotatingMatrix[2][0] = 0.0e0;
   rotatingMatrix[3][0] = 0.0e0;

   // rotating matrix for p-function
   // dMatrix is (53) with gamma=0 in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
   rotatingMatrix[1][1] = cos(alpha)*cos(beta);
   rotatingMatrix[1][2] = -1.0e0*sin(alpha);
   rotatingMatrix[1][3] = cos(alpha)*sin(beta);

   rotatingMatrix[2][1] = sin(alpha)*cos(beta);
   rotatingMatrix[2][2] = cos(alpha);
   rotatingMatrix[2][3] = sin(alpha)*sin(beta);

   rotatingMatrix[3][1] = -1.0e0*sin(beta);
   rotatingMatrix[3][2] = 0.0e0;
   rotatingMatrix[3][3] = cos(beta);

   return rotatingMatrix;
/*
   // rotating matrix for d-function
   // dMatrix is (37) in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
   double dMatrix[OrbitalType_end][OrbitalType_end];
   dMatrix[dzz][dzz] = 0.5*(3.0*(cos(beta)*cos(beta)) - 1.0);
   dMatrix[dxxyy][dxxyy] = cos(0.5*beta)*cos(0.5*beta)*cos(0.5*beta)*cos(0.5*beta);
   dMatrix[dzx][dzx] = (2.0*cos(beta)-1.0)*cos(0.5*beta)*cos(0.5*beta);
   dMatrix[dxxyy][dzx] = -2.0*sin(0.5*beta)*cos(0.5*beta)*cos(0.5*beta)*cos(0.5*beta);
   dMatrix[dxxyy][dzz] = sqrt(6.0)*(sin(0.5*beta)*sin(0.5*beta))*((cos(0.5*beta))*(cos(0.5*beta)));
   dMatrix[dxxyy][dyz] = -2.0*sin(0.5*beta)*sin(0.5*beta)*sin(0.5*beta)*cos(0.5*beta);
   dMatrix[dxxyy][dxy] = sin(0.5*beta)*sin(0.5*beta)*sin(0.5*beta)*sin(0.5*beta);
   dMatrix[dzx][dzz] = -sqrt(6.0)*cos(beta)*cos(0.5*beta)*sin(0.5*beta);
   dMatrix[dzx][dyz] = (2.0*cos(beta)+1.0)*(sin(0.5*beta)*sin(0.5*beta));

   rotatingMatrix[dxy][dxy] = cos(2.0*alpha)*            (dMatrix[dxxyy][dxxyy] - dMatrix[dxxyy][dxy]);
   rotatingMatrix[dxy][dyz] = cos(2.0*alpha)*            (-1.0*dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);
   rotatingMatrix[dxy][dzz] = sqrt(2.0)*sin(2.0*alpha)*  dMatrix[dxxyy][dzz];
   rotatingMatrix[dxy][dzx] = sin(2.0*alpha)*            (-1.0*dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
   rotatingMatrix[dxy][dxxyy] = sin(2.0*alpha)*          (dMatrix[dxxyy][dxxyy] + dMatrix[dxxyy][dxy]);

   rotatingMatrix[dyz][dxy] = cos(alpha)*                (dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
   rotatingMatrix[dyz][dyz] = cos(alpha)*                (dMatrix[dzx][dzx] + dMatrix[dzx][dyz]);
   rotatingMatrix[dyz][dzz] = -1.0*sqrt(2.0)*sin(alpha)* dMatrix[dzx][dzz];
   rotatingMatrix[dyz][dzx] = sin(alpha)*                (dMatrix[dzx][dzx] - dMatrix[dzx][dyz]);
   rotatingMatrix[dyz][dxxyy] = sin(alpha)*              (dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);

   rotatingMatrix[dzz][dxy] = 0.0;
   rotatingMatrix[dzz][dyz] = 0.0;
   rotatingMatrix[dzz][dzz] = dMatrix[dzz][dzz];
   rotatingMatrix[dzz][dzx] = sqrt(2.0)*dMatrix[dzx][dzz];
   rotatingMatrix[dzz][dxxyy] = sqrt(2.0)*dMatrix[dxxyy][dzz];

   rotatingMatrix[dzx][dxy] = -1.0*sin(alpha)*           (dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
   rotatingMatrix[dzx][dyz] = -1.0*sin(alpha)*           (dMatrix[dzx][dzx] + dMatrix[dzx][dyz]);
   rotatingMatrix[dzx][dzz] = -1.0*sqrt(2.0)*cos(alpha)* dMatrix[dzx][dzz];
   rotatingMatrix[dzx][dzx] = cos(alpha)*                (dMatrix[dzx][dzx] - dMatrix[dzx][dyz]);
   rotatingMatrix[dzx][dxxyy] = cos(alpha)*              (dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);

   rotatingMatrix[dxxyy][dxy] = -1.0*sin(2.0*alpha)*     (dMatrix[dxxyy][dxxyy] - dMatrix[dxxyy][dxy]);
   rotatingMatrix[dxxyy][dyz] = -1.0*sin(2.0*alpha)*     (-1.0*dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);
   rotatingMatrix[dxxyy][dzz] = sqrt(2.0)*cos(2.0*alpha)*dMatrix[dxxyy][dzz];
   rotatingMatrix[dxxyy][dzx] = cos(2.0*alpha)*          (-1.0*dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
   rotatingMatrix[dxxyy][dxxyy] = cos(2.0*alpha)*        (dMatrix[dxxyy][dxxyy] + dMatrix[dxxyy][dxy]);
*/
   
}




vector< vector<double> > MatrixRotation::CalcRotatingMatrix_over(const Molecule& mol, int A, int B)
{
	vector< vector<double> > rotatingMatrix;
	int size=4;
	rotatingMatrix.resize(size);
	for(int i=0; i<size;i++)
	{
		rotatingMatrix[i].resize(size);
	}

	double x = mol.MolXYZ[B].x - mol.MolXYZ[A].x;
	double y = mol.MolXYZ[B].y - mol.MolXYZ[A].y;
    double z = mol.MolXYZ[B].z - mol.MolXYZ[A].z;

	EularAngle eularAngle(x, y, z);
	double alpha = eularAngle.GetAlpha();
	double beta  = eularAngle.GetBeta();

	// meu s e s beleza
	// x = 1
	// y = 2
	// z = 3

   // rotating matrix for s-function
   rotatingMatrix[0][0] = 1.0e0;
   rotatingMatrix[0][1] = 0.0e0;
   rotatingMatrix[0][2] = 0.0e0;
   rotatingMatrix[0][3] = 0.0e0;
   rotatingMatrix[1][0] = 0.0e0;
   rotatingMatrix[2][0] = 0.0e0;
   rotatingMatrix[3][0] = 0.0e0;

   // rotating matrix for p-function -- overlap definition
   // dMatrix is (53) with gamma=0 in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
   rotatingMatrix[3][3] = cos(alpha)*cos(beta);
   rotatingMatrix[3][1] = -1.0e0*sin(alpha);
   rotatingMatrix[3][2] = cos(alpha)*sin(beta);

   rotatingMatrix[1][3] = sin(alpha)*cos(beta);
   rotatingMatrix[1][1] = cos(alpha);
   rotatingMatrix[1][2] = sin(alpha)*sin(beta);

   rotatingMatrix[2][3] = -1.0e0*sin(beta);
   rotatingMatrix[2][1] = 0.0e0;
   rotatingMatrix[2][2] = cos(beta);


   return rotatingMatrix;
   
}




/*
Os overlaps diatomicos são construidos assim:
A - B
     s    px    py    pz
s    
px
py
pz

Esssa é a rotação deles

R*Overlap*R

*/



























  /*
  2 CENTER 2 CORE ROTATION
   // rotate (slow algorithm)
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               matrix[mu][nu][lambda][sigma] = 0.0;
               for(int i=0; i<dxy; i++){
                  for(int j=0; j<dxy; j++){
                     for(int k=0; k<dxy; k++){
                        for(int l=0; l<dxy; l++){
                           matrix[mu][nu][lambda][sigma] += oldMatrix[i][j][k][l] 
                                                            *rotatingMatrix[mu][i] 
                                                            *rotatingMatrix[nu][j] 
                                                            *rotatingMatrix[lambda][k] 
                                                            *rotatingMatrix[sigma][l];
                        }
                     }
                  }
               }
            }
         }
      }
   }
   */