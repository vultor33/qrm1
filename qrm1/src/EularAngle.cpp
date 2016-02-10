#include"EularAngle.h"

#include <cmath>

using namespace std;

EularAngle::EularAngle(){
   // e the [BFB_1997] for defenitions of alpha, beta, gamma;
   alpha = 0.0; // (= "phi" in P25 in J. A. Pople book)
   beta  = 0.0; // (= "theta" in P25 in J. A. Pople book)
   gamma = 0.0;
}

EularAngle::EularAngle(double x, double y, double z){
   double r = 0.0;

   // calc. beta
   r = sqrt( pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) );
   beta = acos(z/r);

   // calc. alpha
   if(x==0.0 && y==0.0){
      alpha = 0.0;
   }
   else{ 
      r = sqrt( pow(x, 2.0) + pow(y, 2.0) );
      alpha = atan2(y/r, x/r);
   }

   // set gamma
   gamma = 0.0;
   
}

EularAngle::EularAngle(double* angles){
   alpha = angles[0];
   beta  = angles[1];
   gamma = angles[2];
}

