#ifndef READMLPDD2011UNI_H 
#define READMLPDD2011UNI_H 1

// Include files

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#include "IClassifierReader.h"

/** @class ReadMLPDD2011Uni ReadMLPDD2011Uni.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-12-11
 */
class ReadMLPDD2011Uni : public IClassifierReader {
public: 
  /// Standard constructor
  ReadMLPDD2011Uni( ); 
  ReadMLPDD2011Uni( std::vector<std::string>& theInputVars ) ;
  
  virtual ~ReadMLPDD2011Uni( ); ///< Destructor

  double GetMvaValue( const std::vector<double>& inputValues ) const;
  
protected:

private:
  
  // method-specific destructor
   void Clear();

   // input variable transformation

   int nvar;

   double  cumulativeDist[20][3][1829];
   double  X[20][3][1829];
   double xMin[20][3];
   double xMax[20][3];
   int    nbins[20][3];
   void InitTransform_1();
   void Transform_1( std::vector<double> & iv, int sigOrBgd ) const;
   void InitTransform();
   void Transform( std::vector<double> & iv, int sigOrBgd ) const;

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   const bool fIsNormalised;
   bool IsNormalised() const { return fIsNormalised; }
   double fVmin[20];
   double fVmax[20];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[20];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)

   double ActivationFnc(double x) const;
   double OutputActivationFnc(double x) const;

   int fLayers;
   int fLayerSize[3];
   double fWeightMatrix0to1[26][21];   // weight matrix from layer 0 to 1
   double fWeightMatrix1to2[1][26];   // weight matrix from layer 1 to 2

   double * fWeights[3];
};
#endif // READMLPDD2011UNI_H
