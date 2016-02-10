#ifndef EULARANGLE_H
#define EULARANGLE_H

class EularAngle{
public:
   EularAngle();
   EularAngle(double x, double y, double z);
   explicit EularAngle(double* angles);
   double GetAlpha() const{return this->alpha;}
   double GetBeta()  const{return this->beta;}
   double GetGamma() const{return this->gamma;}
   void SetAlpha(double alpha){this->alpha = alpha;}
   void SetBeta (double beta) {this->beta  = beta;}
   void SetGamma(double gamma){this->gamma = gamma;}
private:
   double alpha;
   double beta;
   double gamma;
};

#endif
