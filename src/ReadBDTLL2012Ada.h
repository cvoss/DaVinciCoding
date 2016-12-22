#ifndef READBDTLL2012ADA_H 
#define READBDTLL2012ADA_H 1

#define NN_LL new BDTLL2012AdaNode

// Include files
#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#include "IClassifierReader.h"

/** @class ReadBDTLL2012Ada ReadBDTLL2012Ada.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-11-25
 */

class BDTLL2012AdaNode {
   
public:
   
   // constructor of an essentially "empty" node floating in space
   BDTLL2012AdaNode ( BDTLL2012AdaNode* left,BDTLL2012AdaNode* right,
                          int selector, double cutValue, bool cutType, 
                          int nodeType, double purity, double response ) :
   fLeft         ( left         ),
   fRight        ( right        ),
   fSelector     ( selector     ),
   fCutValue     ( cutValue     ),
   fCutType      ( cutType      ),
   fNodeType     ( nodeType     ),
   fPurity       ( purity       ),
   fResponse     ( response     ){
   }

   virtual ~BDTLL2012AdaNode();

   // test event if it decends the tree at this node to the right
   virtual bool GoesRight( const std::vector<double>& inputValues ) const;
   BDTLL2012AdaNode* GetRight( void )  {return fRight; };

   // test event if it decends the tree at this node to the left 
   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;
   BDTLL2012AdaNode* GetLeft( void ) { return fLeft; };   

   // return  S/(S+B) (purity) at this node (from  training)

   double GetPurity( void ) const { return fPurity; } 
   // return the node type
   int    GetNodeType( void ) const { return fNodeType; }
   double GetResponse(void) const {return fResponse;}

private:

   BDTLL2012AdaNode*   fLeft;     // pointer to the left daughter node
   BDTLL2012AdaNode*   fRight;    // pointer to the right daughter node
   int                     fSelector; // index of variable used in node selection (decision tree)   
   double                  fCutValue; // cut value appplied on this node to discriminate bkg against sig
   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise
   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal 
   double                  fPurity;   // Purity of node from training
   double                  fResponse; // Regression response value of node
};

//_______________________________________________________________________
class ReadBDTLL2012Ada : public IClassifierReader {
public: 
  /// Standard constructor
  ReadBDTLL2012Ada( ); 
  ReadBDTLL2012Ada( std::vector<std::string>& theInputVars ) ;
  
  virtual ~ReadBDTLL2012Ada( ); ///< Destructor

  double GetMvaValue( const std::vector<double>& inputValues ) const;

protected:

private:

  // method-specific destructor
  void Clear();
  
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
  std::vector<BDTLL2012AdaNode*> fForest;       // i.e. root nodes of decision trees
  std::vector<double>            fBoostWeights; // the weights applied in the individual boosts
};
#endif // READBDTLL2012ADA_H
