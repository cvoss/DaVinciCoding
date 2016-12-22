// Include files 



// local
#include "ReadMLPDD2012.h"

//-----------------------------------------------------------------------------
// Implementation file for class : ReadMLPDD2012
//
// 2014-11-25 : Christian Voss
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ReadMLPDD2012::ReadMLPDD2012( std::vector<std::string>& theInputVars ) 
  : IClassifierReader()
  , fClassName( "ReadMLPDD2012" )
  , fNvars( 20 )
  , fIsNormalised( false )
{
  // the training input variables
  const char* inputVars[] = { 
    "log_pi_IP",
    "log_p_IP",
    "log_Lpi_IP",
    "log_Lp_IP",
    "log_BIP",
    "log_LIP",
    "log_BVtxchi2",
    "FitProb",
    "log_B_DauSumchi2",
    "log_B_FLBchi2",
    "log_Lambda_FDchi2",
    "B_PT",
    "B_Eta",
    "log_B_ctau", 
    "log_L_ctau", 
    "BWerner",
    "LWerner", 
    "ppiDist", 
    "LAngle", 
    "Pointing" 
  };
  
  // sanity checks
  if (theInputVars.size() <= 0) {
    std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
    fStatusIsClean = false;
  }
  
  if (theInputVars.size() != fNvars) {
    std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
              << theInputVars.size() << " != " << fNvars << std::endl;
    fStatusIsClean = false;
  }
  
  // validate input variables
  for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
    if (theInputVars[ivar] != inputVars[ivar]) {
      std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
      fStatusIsClean = false;
    }
  }
  
  // initialize min and max vectors (for normalisation)
  fVmin[0] = -1;
  fVmax[0] = 1;
  fVmin[1] = -1;
  fVmax[1] = 1;
  fVmin[2] = -1;
  fVmax[2] = 0.99999988079071;
  fVmin[3] = -1;
  fVmax[3] = 1;
  fVmin[4] = -1;
  fVmax[4] = 1;
  fVmin[5] = -1;
  fVmax[5] = 1;
  fVmin[6] = -1;
  fVmax[6] = 1;
  fVmin[7] = -1;
  fVmax[7] = 1;
  fVmin[8] = -1;
  fVmax[8] = 0.99999988079071;
  fVmin[9] = -1;
  fVmax[9] = 1;
  fVmin[10] = -1;
  fVmax[10] = 1;
  fVmin[11] = -1;
  fVmax[11] = 1;
  fVmin[12] = -1;
  fVmax[12] = 1;
  fVmin[13] = -1;
  fVmax[13] = 1;
  fVmin[14] = -1;
  fVmax[14] = 1;
  fVmin[15] = -1;
  fVmax[15] = 0.99999988079071;
  fVmin[16] = -1;
  fVmax[16] = 1;
  fVmin[17] = -1;
  fVmax[17] = 1;
  fVmin[18] = -1;
  fVmax[18] = 1;
  fVmin[19] = -1;
  fVmax[19] = 0.99999988079071;
  
  // initialize input variable types
  fType[0] = 'F';
  fType[1] = 'F';
  fType[2] = 'F';
  fType[3] = 'F';
  fType[4] = 'F';
  fType[5] = 'F';
  fType[6] = 'F';
  fType[7] = 'F';
  fType[8] = 'F';
  fType[9] = 'F';
  fType[10] = 'F';
  fType[11] = 'F';
  fType[12] = 'F';
  fType[13] = 'F';
  fType[14] = 'F';
  fType[15] = 'F';
  fType[16] = 'F';
  fType[17] = 'F';
  fType[18] = 'F';
  fType[19] = 'F';
  
  // initialize constants
  Initialize();
  
  // initialize transformation
  InitTransform();
}
//=============================================================================
// Destructor
//=============================================================================
ReadMLPDD2012::~ReadMLPDD2012()
{
  Clear(); // method-specific 
}
//=============================================================================
double ReadMLPDD2012::GetMvaValue__( const std::vector<double>& inputValues ) const
{
  if (inputValues.size() != (unsigned int)fLayerSize[0]-1) {
    std::cout << "Input vector needs to be of size " << fLayerSize[0]-1 << std::endl;
    return 0;
  }
  
  for (int l=0; l<fLayers; l++)
    for (int i=0; i<fLayerSize[l]; i++) fWeights[l][i]=0;
  
  for (int l=0; l<fLayers-1; l++)
    fWeights[l][fLayerSize[l]-1]=1;
  
  for (int i=0; i<fLayerSize[0]-1; i++)
    fWeights[0][i]=inputValues[i];
  
  // layer 0 to 1
  for (int o=0; o<fLayerSize[1]-1; o++) {
    for (int i=0; i<fLayerSize[0]; i++) {
      double inputVal = fWeightMatrix0to1[o][i] * fWeights[0][i];
      fWeights[1][o] += inputVal;
    }
    fWeights[1][o] = ActivationFnc(fWeights[1][o]);
  }
  // layer 1 to 2
  for (int o=0; o<fLayerSize[2]; o++) {
    for (int i=0; i<fLayerSize[1]; i++) {
      double inputVal = fWeightMatrix1to2[o][i] * fWeights[1][i];
      fWeights[2][o] += inputVal;
    }
    fWeights[2][o] = OutputActivationFnc(fWeights[2][o]);
  }
  
  return fWeights[2][0];
}
//=============================================================================
double ReadMLPDD2012::ActivationFnc(double x) const {
  // hyperbolic tan
  return tanh(x);
}
//=============================================================================
double ReadMLPDD2012::OutputActivationFnc(double x) const {
  // identity
  return x;
}
//=============================================================================
// Clean up
inline void ReadMLPDD2012::Clear() 
{
  // clean up the arrays
  for (int lIdx = 0; lIdx < 3; lIdx++) {
    delete[] fWeights[lIdx];
  }
}
//=============================================================================
double ReadMLPDD2012::GetMvaValue( const std::vector<double>& inputValues ) const
{
  // classifier response value
  double retval = 0;
  
  // classifier response, sanity check first
  if (!IsStatusClean()) {
    std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
              << " because status is dirty" << std::endl;
    retval = 0;
  }
  else {
    if (IsNormalised()) {
      // normalise variables
      std::vector<double> iV;
      iV.reserve(inputValues.size());
      int ivar = 0;
      for (std::vector<double>::const_iterator varIt = inputValues.begin();
           varIt != inputValues.end(); varIt++, ivar++) {
        iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));
      }
      Transform( iV, -1 );
      retval = GetMvaValue__( iV );
    }
    else {
      std::vector<double> iV;
      int ivar = 0;
      for (std::vector<double>::const_iterator varIt = inputValues.begin();
           varIt != inputValues.end(); varIt++, ivar++) {
        iV.push_back(*varIt);
      }
      Transform( iV, -1 );
      retval = GetMvaValue__( iV );
    }
  }
  
  return retval;
}
//_______________________________________________________________________
void ReadMLPDD2012::InitTransform_1()
{
  // Normalization transformation, initialisation
   fMin_1[0][0] = -1.46886312962;
   fMax_1[0][0] = 1.04128730297;
   fMin_1[1][0] = -1.54731798172;
   fMax_1[1][0] = 0.760807096958;
   fMin_1[2][0] = -1.54731798172;
   fMax_1[2][0] = 1.04128730297;
   fMin_1[0][1] = -1.41660022736;
   fMax_1[0][1] = 0.696034431458;
   fMin_1[1][1] = -1.55906772614;
   fMax_1[1][1] = 0.598295807838;
   fMin_1[2][1] = -1.55906772614;
   fMax_1[2][1] = 0.696034431458;
   fMin_1[0][2] = 0.233619987965;
   fMax_1[0][2] = 2.0077559948;
   fMin_1[1][2] = 0.0131395924836;
   fMax_1[1][2] = 1.92900812626;
   fMin_1[2][2] = 0.0131395924836;
   fMax_1[2][2] = 2.0077559948;
   fMin_1[0][3] = -0.428629815578;
   fMax_1[0][3] = 1.35506546497;
   fMin_1[1][3] = -0.481942236423;
   fMax_1[1][3] = 1.23818778992;
   fMin_1[2][3] = -0.481942236423;
   fMax_1[2][3] = 1.35506546497;
   fMin_1[0][4] = -3.83516454697;
   fMax_1[0][4] = -0.883335828781;
   fMin_1[1][4] = -4.45256614685;
   fMax_1[1][4] = -0.735439062119;
   fMin_1[2][4] = -4.45256614685;
   fMax_1[2][4] = -0.735439062119;
   fMin_1[0][5] = -2.40407109261;
   fMax_1[0][5] = 0.706591188908;
   fMin_1[1][5] = -3.2453725338;
   fMax_1[1][5] = 0.674849152565;
   fMin_1[2][5] = -3.2453725338;
   fMax_1[2][5] = 0.706591188908;
   fMin_1[0][6] = 0.0551982149482;
   fMax_1[0][6] = 4.9922785759;
   fMin_1[1][6] = 0.0982858240604;
   fMax_1[1][6] = 4.96775722504;
   fMin_1[2][6] = 0.0551982149482;
   fMax_1[2][6] = 4.9922785759;
   fMin_1[0][7] = 1.14490721899e-05;
   fMax_1[0][7] = 0.999765396118;
   fMin_1[1][7] = 1.23313402582e-05;
   fMax_1[1][7] = 0.998426258564;
   fMin_1[2][7] = 1.14490721899e-05;
   fMax_1[2][7] = 0.999765396118;
   fMin_1[0][8] = 1.29851078987;
   fMax_1[0][8] = 5.05440044403;
   fMin_1[1][8] = 1.51460409164;
   fMax_1[1][8] = 4.71614408493;
   fMin_1[2][8] = 1.29851078987;
   fMax_1[2][8] = 5.05440044403;
   fMin_1[0][9] = 0.83063018322;
   fMax_1[0][9] = 2.81794953346;
   fMin_1[1][9] = 0.790921747684;
   fMax_1[1][9] = 2.43504452705;
   fMin_1[2][9] = 0.790921747684;
   fMax_1[2][9] = 2.81794953346;
   fMin_1[0][10] = 0.762708127499;
   fMax_1[0][10] = 3.33264160156;
   fMin_1[1][10] = 0.752305984497;
   fMax_1[1][10] = 3.16218924522;
   fMin_1[2][10] = 0.752305984497;
   fMax_1[2][10] = 3.33264160156;
   fMin_1[0][11] = 1512.6138916;
   fMax_1[0][11] = 48625.7695312;
   fMin_1[1][11] = 1492.55407715;
   fMax_1[1][11] = 42515.5351562;
   fMin_1[2][11] = 1492.55407715;
   fMax_1[2][11] = 48625.7695312;
   fMin_1[0][12] = 2.39824223518;
   fMax_1[0][12] = 5.55322599411;
   fMin_1[1][12] = 2.31465554237;
   fMax_1[1][12] = 5.45913934708;
   fMin_1[2][12] = 2.31465554237;
   fMax_1[2][12] = 5.55322599411;
   fMin_1[0][13] = -1.22993326187;
   fMax_1[0][13] = 0.496834427118;
   fMin_1[1][13] = -1.23079657555;
   fMax_1[1][13] = 0.49398842454;
   fMin_1[2][13] = -1.23079657555;
   fMax_1[2][13] = 0.496834427118;
   fMin_1[0][14] = 0.811039686203;
   fMax_1[0][14] = 2.48783326149;
   fMin_1[1][14] = 0.851882338524;
   fMax_1[1][14] = 2.37544369698;
   fMin_1[2][14] = 0.811039686203;
   fMax_1[2][14] = 2.48783326149;
   fMin_1[0][15] = 4.11361133956e-06;
   fMax_1[0][15] = 0.0299081020057;
   fMin_1[1][15] = 2.77817107417e-05;
   fMax_1[1][15] = 0.0299989879131;
   fMin_1[2][15] = 4.11361133956e-06;
   fMax_1[2][15] = 0.0299989879131;
   fMin_1[0][16] = -3.84514093399;
   fMax_1[0][16] = -0.35182890296;
   fMin_1[1][16] = -3.85285067558;
   fMax_1[1][16] = -0.340868204832;
   fMin_1[2][16] = -3.85285067558;
   fMax_1[2][16] = -0.340868204832;
   fMin_1[0][17] = 3.96544828618e-06;
   fMax_1[0][17] = 0.129571795464;
   fMin_1[1][17] = 8.18750777398e-06;
   fMax_1[1][17] = 0.129972323775;
   fMin_1[2][17] = 3.96544828618e-06;
   fMax_1[2][17] = 0.129972323775;
   fMin_1[0][18] = 0.0188157968223;
   fMax_1[0][18] = 3.07901358604;
   fMin_1[1][18] = 0.0470530167222;
   fMax_1[1][18] = 3.13562917709;
   fMin_1[2][18] = 0.0188157968223;
   fMax_1[2][18] = 3.13562917709;
   fMin_1[0][19] = 4.84630254505e-05;
   fMax_1[0][19] = 0.391300797462;
   fMin_1[1][19] = 0.000690708809998;
   fMax_1[1][19] = 0.474036633968;
   fMin_1[2][19] = 4.84630254505e-05;
   fMax_1[2][19] = 0.474036633968;
}
//_______________________________________________________________________
void ReadMLPDD2012::Transform_1( std::vector<double>& iv, int cls) const
{
   // Normalization transformation
   if (cls < 0 || cls > 2) {
   if (2 > 1 ) cls = 2;
      else cls = 2;
   }
   const int nVar = 20;

   // get indices of used variables

   // define the indices of the variables which are transformed by this transformation
   static std::vector<int> indicesGet;
   static std::vector<int> indicesPut;

   if ( indicesGet.empty() ) { 
      indicesGet.reserve(fNvars);
      indicesGet.push_back( 0);
      indicesGet.push_back( 1);
      indicesGet.push_back( 2);
      indicesGet.push_back( 3);
      indicesGet.push_back( 4);
      indicesGet.push_back( 5);
      indicesGet.push_back( 6);
      indicesGet.push_back( 7);
      indicesGet.push_back( 8);
      indicesGet.push_back( 9);
      indicesGet.push_back( 10);
      indicesGet.push_back( 11);
      indicesGet.push_back( 12);
      indicesGet.push_back( 13);
      indicesGet.push_back( 14);
      indicesGet.push_back( 15);
      indicesGet.push_back( 16);
      indicesGet.push_back( 17);
      indicesGet.push_back( 18);
      indicesGet.push_back( 19);
   } 
   if ( indicesPut.empty() ) { 
      indicesPut.reserve(fNvars);
      indicesPut.push_back( 0);
      indicesPut.push_back( 1);
      indicesPut.push_back( 2);
      indicesPut.push_back( 3);
      indicesPut.push_back( 4);
      indicesPut.push_back( 5);
      indicesPut.push_back( 6);
      indicesPut.push_back( 7);
      indicesPut.push_back( 8);
      indicesPut.push_back( 9);
      indicesPut.push_back( 10);
      indicesPut.push_back( 11);
      indicesPut.push_back( 12);
      indicesPut.push_back( 13);
      indicesPut.push_back( 14);
      indicesPut.push_back( 15);
      indicesPut.push_back( 16);
      indicesPut.push_back( 17);
      indicesPut.push_back( 18);
      indicesPut.push_back( 19);
   } 

   static std::vector<double> dv;
   dv.resize(nVar);
   for (int ivar=0; ivar<nVar; ivar++) dv[ivar] = iv[indicesGet.at(ivar)];
   for (int ivar=0;ivar<20;ivar++) {
      double offset = fMin_1[cls][ivar];
      double scale  = 1.0/(fMax_1[cls][ivar]-fMin_1[cls][ivar]);
      iv[indicesPut.at(ivar)] = (dv[ivar]-offset)*scale * 2 - 1;
   }
}

//_______________________________________________________________________
void ReadMLPDD2012::InitTransform()
{
   InitTransform_1();
}

//_______________________________________________________________________
void ReadMLPDD2012::Transform( std::vector<double>& iv, int sigOrBgd ) const
{
   Transform_1( iv, sigOrBgd );
}
//_______________________________________________________________________
void ReadMLPDD2012::Initialize()
{
   // build network structure
   fLayers = 3;
   fLayerSize[0] = 21; fWeights[0] = new double[21]; 
   fLayerSize[1] = 26; fWeights[1] = new double[26]; 
   fLayerSize[2] = 1; fWeights[2] = new double[1]; 
   // weight matrix from layer 0 to 1
   fWeightMatrix0to1[0][0] = -0.547621023118515;
   fWeightMatrix0to1[1][0] = 1.71837360380379;
   fWeightMatrix0to1[2][0] = 0.657502394227848;
   fWeightMatrix0to1[3][0] = 1.53500883537181;
   fWeightMatrix0to1[4][0] = -1.59340200047872;
   fWeightMatrix0to1[5][0] = -1.2006175088291;
   fWeightMatrix0to1[6][0] = -0.526694744453009;
   fWeightMatrix0to1[7][0] = 1.76522830792282;
   fWeightMatrix0to1[8][0] = -1.30120114672178;
   fWeightMatrix0to1[9][0] = -1.20337722348957;
   fWeightMatrix0to1[10][0] = -1.30969354413382;
   fWeightMatrix0to1[11][0] = -0.240596462142571;
   fWeightMatrix0to1[12][0] = -0.956000412309808;
   fWeightMatrix0to1[13][0] = -0.163748359848593;
   fWeightMatrix0to1[14][0] = 0.538947813150695;
   fWeightMatrix0to1[15][0] = 1.14556709387676;
   fWeightMatrix0to1[16][0] = -0.413641711454501;
   fWeightMatrix0to1[17][0] = 1.90737652450903;
   fWeightMatrix0to1[18][0] = 0.648139264522359;
   fWeightMatrix0to1[19][0] = 1.74021970028867;
   fWeightMatrix0to1[20][0] = -1.43716814341997;
   fWeightMatrix0to1[21][0] = -0.883679760233042;
   fWeightMatrix0to1[22][0] = 0.950291577875211;
   fWeightMatrix0to1[23][0] = 0.117455813784386;
   fWeightMatrix0to1[24][0] = -1.14095256610572;
   fWeightMatrix0to1[0][1] = -0.22193891322223;
   fWeightMatrix0to1[1][1] = 1.50692889915256;
   fWeightMatrix0to1[2][1] = -0.745387073877011;
   fWeightMatrix0to1[3][1] = -1.79900004184736;
   fWeightMatrix0to1[4][1] = -1.33888467155573;
   fWeightMatrix0to1[5][1] = 0.58009838617049;
   fWeightMatrix0to1[6][1] = 1.62144094583041;
   fWeightMatrix0to1[7][1] = 0.0742387348354247;
   fWeightMatrix0to1[8][1] = -0.307960938247377;
   fWeightMatrix0to1[9][1] = -0.132703106054574;
   fWeightMatrix0to1[10][1] = -0.725146970069841;
   fWeightMatrix0to1[11][1] = -1.42288263025115;
   fWeightMatrix0to1[12][1] = 1.7717246726747;
   fWeightMatrix0to1[13][1] = -0.610582753776211;
   fWeightMatrix0to1[14][1] = -0.338088718611676;
   fWeightMatrix0to1[15][1] = 1.69583411111088;
   fWeightMatrix0to1[16][1] = 0.840828779347867;
   fWeightMatrix0to1[17][1] = 1.85776123066739;
   fWeightMatrix0to1[18][1] = 1.21867956150574;
   fWeightMatrix0to1[19][1] = -0.375357610658757;
   fWeightMatrix0to1[20][1] = -1.23164414169139;
   fWeightMatrix0to1[21][1] = 0.575388869058029;
   fWeightMatrix0to1[22][1] = -0.230498167865036;
   fWeightMatrix0to1[23][1] = 1.22995240258836;
   fWeightMatrix0to1[24][1] = 1.409691600091;
   fWeightMatrix0to1[0][2] = 1.1625518007926;
   fWeightMatrix0to1[1][2] = 1.42335391948674;
   fWeightMatrix0to1[2][2] = -1.55231292138203;
   fWeightMatrix0to1[3][2] = 1.38083503895744;
   fWeightMatrix0to1[4][2] = -0.165435043025712;
   fWeightMatrix0to1[5][2] = -0.986591725550627;
   fWeightMatrix0to1[6][2] = -1.5079654801373;
   fWeightMatrix0to1[7][2] = -1.38781335141556;
   fWeightMatrix0to1[8][2] = 1.64995128057307;
   fWeightMatrix0to1[9][2] = 0.88916720883909;
   fWeightMatrix0to1[10][2] = -1.74940056230217;
   fWeightMatrix0to1[11][2] = 0.25717404029907;
   fWeightMatrix0to1[12][2] = -0.488260312735659;
   fWeightMatrix0to1[13][2] = 0.682369294712246;
   fWeightMatrix0to1[14][2] = 0.524529706782914;
   fWeightMatrix0to1[15][2] = -0.416578344347059;
   fWeightMatrix0to1[16][2] = 0.226067815947851;
   fWeightMatrix0to1[17][2] = -1.36951729416033;
   fWeightMatrix0to1[18][2] = 0.689461261294795;
   fWeightMatrix0to1[19][2] = -0.472881384423014;
   fWeightMatrix0to1[20][2] = -0.269990884555645;
   fWeightMatrix0to1[21][2] = -1.43826890402407;
   fWeightMatrix0to1[22][2] = 0.841222358997791;
   fWeightMatrix0to1[23][2] = 0.585876434277006;
   fWeightMatrix0to1[24][2] = 1.21379331357021;
   fWeightMatrix0to1[0][3] = -0.694549505939446;
   fWeightMatrix0to1[1][3] = -1.90077309542615;
   fWeightMatrix0to1[2][3] = -1.93976204377455;
   fWeightMatrix0to1[3][3] = 1.04470202159538;
   fWeightMatrix0to1[4][3] = 0.0602973955931599;
   fWeightMatrix0to1[5][3] = 1.89701613114047;
   fWeightMatrix0to1[6][3] = 0.40566365169723;
   fWeightMatrix0to1[7][3] = 0.445981617092708;
   fWeightMatrix0to1[8][3] = -0.61385191866728;
   fWeightMatrix0to1[9][3] = -0.777198981102614;
   fWeightMatrix0to1[10][3] = -0.456324303261828;
   fWeightMatrix0to1[11][3] = 1.15769576858258;
   fWeightMatrix0to1[12][3] = -1.13988620831092;
   fWeightMatrix0to1[13][3] = -1.50362740455717;
   fWeightMatrix0to1[14][3] = -0.658920158646953;
   fWeightMatrix0to1[15][3] = 0.0192098779883943;
   fWeightMatrix0to1[16][3] = -0.407999139738286;
   fWeightMatrix0to1[17][3] = 1.45814656414314;
   fWeightMatrix0to1[18][3] = -1.28483911639713;
   fWeightMatrix0to1[19][3] = -0.152047443774279;
   fWeightMatrix0to1[20][3] = 0.322117466957915;
   fWeightMatrix0to1[21][3] = -0.895034251987954;
   fWeightMatrix0to1[22][3] = -1.40236860977579;
   fWeightMatrix0to1[23][3] = -1.77112103146611;
   fWeightMatrix0to1[24][3] = 0.0502541155617683;
   fWeightMatrix0to1[0][4] = -1.74509434922991;
   fWeightMatrix0to1[1][4] = -1.65639241379124;
   fWeightMatrix0to1[2][4] = 0.587313776656462;
   fWeightMatrix0to1[3][4] = 1.79675880633034;
   fWeightMatrix0to1[4][4] = -1.07266512034219;
   fWeightMatrix0to1[5][4] = 1.5245850383255;
   fWeightMatrix0to1[6][4] = -1.06733175639899;
   fWeightMatrix0to1[7][4] = -0.406774439046625;
   fWeightMatrix0to1[8][4] = 0.222771261456085;
   fWeightMatrix0to1[9][4] = 0.0161467863602446;
   fWeightMatrix0to1[10][4] = -1.4523175105185;
   fWeightMatrix0to1[11][4] = 0.329805684660783;
   fWeightMatrix0to1[12][4] = 0.181823352602571;
   fWeightMatrix0to1[13][4] = 0.742045950384608;
   fWeightMatrix0to1[14][4] = -0.113868088071574;
   fWeightMatrix0to1[15][4] = 1.67909879934666;
   fWeightMatrix0to1[16][4] = 0.19176195757348;
   fWeightMatrix0to1[17][4] = -0.93359093072189;
   fWeightMatrix0to1[18][4] = 0.551594811588662;
   fWeightMatrix0to1[19][4] = -1.19718981964744;
   fWeightMatrix0to1[20][4] = -3.11147905282195;
   fWeightMatrix0to1[21][4] = 0.0611693438276768;
   fWeightMatrix0to1[22][4] = -0.320573800103744;
   fWeightMatrix0to1[23][4] = 1.50229933531426;
   fWeightMatrix0to1[24][4] = 0.603279807773776;
   fWeightMatrix0to1[0][5] = -0.28889404533001;
   fWeightMatrix0to1[1][5] = -0.317777848934413;
   fWeightMatrix0to1[2][5] = 0.116735602975199;
   fWeightMatrix0to1[3][5] = -1.61985775518731;
   fWeightMatrix0to1[4][5] = -0.284028739977464;
   fWeightMatrix0to1[5][5] = 0.418332122936449;
   fWeightMatrix0to1[6][5] = -1.93623741537961;
   fWeightMatrix0to1[7][5] = 0.607386644672658;
   fWeightMatrix0to1[8][5] = -1.46592579294099;
   fWeightMatrix0to1[9][5] = -0.351665714685886;
   fWeightMatrix0to1[10][5] = 0.859175837693608;
   fWeightMatrix0to1[11][5] = 1.6805882150392;
   fWeightMatrix0to1[12][5] = -2.01240256841516;
   fWeightMatrix0to1[13][5] = 0.328701309440203;
   fWeightMatrix0to1[14][5] = 1.33643426051556;
   fWeightMatrix0to1[15][5] = 1.29294798392738;
   fWeightMatrix0to1[16][5] = 0.304905536161487;
   fWeightMatrix0to1[17][5] = -1.63789858583986;
   fWeightMatrix0to1[18][5] = -0.286088594897248;
   fWeightMatrix0to1[19][5] = -1.27844332940104;
   fWeightMatrix0to1[20][5] = 1.02506832256963;
   fWeightMatrix0to1[21][5] = 1.25939438330297;
   fWeightMatrix0to1[22][5] = -0.498928797717137;
   fWeightMatrix0to1[23][5] = -0.886032184042065;
   fWeightMatrix0to1[24][5] = -0.0892370254301625;
   fWeightMatrix0to1[0][6] = -0.888492939094713;
   fWeightMatrix0to1[1][6] = 0.442630188484819;
   fWeightMatrix0to1[2][6] = 1.52489723211494;
   fWeightMatrix0to1[3][6] = -0.746196003850723;
   fWeightMatrix0to1[4][6] = -0.601603504314664;
   fWeightMatrix0to1[5][6] = -0.0594841272163887;
   fWeightMatrix0to1[6][6] = 1.40498644821677;
   fWeightMatrix0to1[7][6] = 0.071377964048619;
   fWeightMatrix0to1[8][6] = 0.745671934964292;
   fWeightMatrix0to1[9][6] = -1.26312815991601;
   fWeightMatrix0to1[10][6] = 1.62305674620482;
   fWeightMatrix0to1[11][6] = 0.90077024805302;
   fWeightMatrix0to1[12][6] = 1.03769407900979;
   fWeightMatrix0to1[13][6] = 0.0176601509711533;
   fWeightMatrix0to1[14][6] = 0.495585625380941;
   fWeightMatrix0to1[15][6] = 0.221714337321144;
   fWeightMatrix0to1[16][6] = -0.03149910327636;
   fWeightMatrix0to1[17][6] = 1.19415342978529;
   fWeightMatrix0to1[18][6] = -1.24264874426282;
   fWeightMatrix0to1[19][6] = -0.293041057114824;
   fWeightMatrix0to1[20][6] = 0.73320620627574;
   fWeightMatrix0to1[21][6] = 0.945292520802145;
   fWeightMatrix0to1[22][6] = -0.0567046701725098;
   fWeightMatrix0to1[23][6] = 1.86241612866469;
   fWeightMatrix0to1[24][6] = 2.13714916035519;
   fWeightMatrix0to1[0][7] = -2.19779995694342;
   fWeightMatrix0to1[1][7] = 0.667684878829621;
   fWeightMatrix0to1[2][7] = 0.355981338134968;
   fWeightMatrix0to1[3][7] = 0.210005849733254;
   fWeightMatrix0to1[4][7] = -0.369219243966199;
   fWeightMatrix0to1[5][7] = -0.621733736817813;
   fWeightMatrix0to1[6][7] = 1.31365066470436;
   fWeightMatrix0to1[7][7] = 1.20708374485288;
   fWeightMatrix0to1[8][7] = -1.2192279176121;
   fWeightMatrix0to1[9][7] = -0.424521788797842;
   fWeightMatrix0to1[10][7] = -0.393624145018861;
   fWeightMatrix0to1[11][7] = 0.239810260041538;
   fWeightMatrix0to1[12][7] = 1.43874700300601;
   fWeightMatrix0to1[13][7] = -0.179803526486604;
   fWeightMatrix0to1[14][7] = 0.0441547397471651;
   fWeightMatrix0to1[15][7] = -0.785423561146083;
   fWeightMatrix0to1[16][7] = 0.491117983942854;
   fWeightMatrix0to1[17][7] = 1.89031905574359;
   fWeightMatrix0to1[18][7] = 1.70781305893537;
   fWeightMatrix0to1[19][7] = 0.182797677871198;
   fWeightMatrix0to1[20][7] = 0.131903078558961;
   fWeightMatrix0to1[21][7] = -1.90483344123134;
   fWeightMatrix0to1[22][7] = 0.782098341028125;
   fWeightMatrix0to1[23][7] = 1.07075847895096;
   fWeightMatrix0to1[24][7] = -1.27201742353062;
   fWeightMatrix0to1[0][8] = -0.724196688538807;
   fWeightMatrix0to1[1][8] = 0.107952440306348;
   fWeightMatrix0to1[2][8] = 0.233565549671535;
   fWeightMatrix0to1[3][8] = -2.15049842706771;
   fWeightMatrix0to1[4][8] = 2.12770251124501;
   fWeightMatrix0to1[5][8] = 0.993494564583748;
   fWeightMatrix0to1[6][8] = -0.186462487175718;
   fWeightMatrix0to1[7][8] = 0.00177558952621251;
   fWeightMatrix0to1[8][8] = 1.43634433291277;
   fWeightMatrix0to1[9][8] = 0.905356831038628;
   fWeightMatrix0to1[10][8] = 0.424908051923684;
   fWeightMatrix0to1[11][8] = -0.0779352279828249;
   fWeightMatrix0to1[12][8] = -1.68623149402059;
   fWeightMatrix0to1[13][8] = -0.573671539930479;
   fWeightMatrix0to1[14][8] = 0.856157596470528;
   fWeightMatrix0to1[15][8] = 1.14491478679373;
   fWeightMatrix0to1[16][8] = 0.648672129983564;
   fWeightMatrix0to1[17][8] = 0.31196035583702;
   fWeightMatrix0to1[18][8] = 2.01053950537287;
   fWeightMatrix0to1[19][8] = -1.05865910686065;
   fWeightMatrix0to1[20][8] = -0.462607633178128;
   fWeightMatrix0to1[21][8] = -0.455232851145832;
   fWeightMatrix0to1[22][8] = -1.52006205180573;
   fWeightMatrix0to1[23][8] = 0.548671828907541;
   fWeightMatrix0to1[24][8] = 1.80062381178379;
   fWeightMatrix0to1[0][9] = -1.23512237117731;
   fWeightMatrix0to1[1][9] = 0.539290040526383;
   fWeightMatrix0to1[2][9] = -0.103154667193285;
   fWeightMatrix0to1[3][9] = -2.02685500644276;
   fWeightMatrix0to1[4][9] = -2.17527292703121;
   fWeightMatrix0to1[5][9] = 1.11545039354359;
   fWeightMatrix0to1[6][9] = 0.461247480409866;
   fWeightMatrix0to1[7][9] = 0.997583236277286;
   fWeightMatrix0to1[8][9] = -1.40552925389116;
   fWeightMatrix0to1[9][9] = 1.38656762459319;
   fWeightMatrix0to1[10][9] = -2.0754439135979;
   fWeightMatrix0to1[11][9] = 1.0519498716039;
   fWeightMatrix0to1[12][9] = 0.448830534102526;
   fWeightMatrix0to1[13][9] = -1.37259039686897;
   fWeightMatrix0to1[14][9] = -1.87774054774654;
   fWeightMatrix0to1[15][9] = -1.77297123629501;
   fWeightMatrix0to1[16][9] = -0.424532043208213;
   fWeightMatrix0to1[17][9] = -1.70173009786989;
   fWeightMatrix0to1[18][9] = -1.44024336892222;
   fWeightMatrix0to1[19][9] = -1.01567670385323;
   fWeightMatrix0to1[20][9] = 0.847896762949791;
   fWeightMatrix0to1[21][9] = -1.33759480453972;
   fWeightMatrix0to1[22][9] = -1.58140481101846;
   fWeightMatrix0to1[23][9] = 1.31072282009118;
   fWeightMatrix0to1[24][9] = -0.701753193858298;
   fWeightMatrix0to1[0][10] = 0.565040046475798;
   fWeightMatrix0to1[1][10] = -0.972574939343861;
   fWeightMatrix0to1[2][10] = 0.461835807688123;
   fWeightMatrix0to1[3][10] = 1.7168196147201;
   fWeightMatrix0to1[4][10] = 0.231917999745079;
   fWeightMatrix0to1[5][10] = -0.855300964312458;
   fWeightMatrix0to1[6][10] = -1.22363408550187;
   fWeightMatrix0to1[7][10] = -0.262274785537653;
   fWeightMatrix0to1[8][10] = -0.906342212405445;
   fWeightMatrix0to1[9][10] = -0.622660824298265;
   fWeightMatrix0to1[10][10] = -0.301628417832071;
   fWeightMatrix0to1[11][10] = -0.182815391247393;
   fWeightMatrix0to1[12][10] = 1.89354662412271;
   fWeightMatrix0to1[13][10] = -0.420972158723918;
   fWeightMatrix0to1[14][10] = -0.766593636134799;
   fWeightMatrix0to1[15][10] = 0.167641661541172;
   fWeightMatrix0to1[16][10] = -1.35523245458955;
   fWeightMatrix0to1[17][10] = 1.28063056212566;
   fWeightMatrix0to1[18][10] = 1.14300083645964;
   fWeightMatrix0to1[19][10] = -0.902219947957295;
   fWeightMatrix0to1[20][10] = -0.559691858205672;
   fWeightMatrix0to1[21][10] = 0.343466007395962;
   fWeightMatrix0to1[22][10] = 0.998969352000417;
   fWeightMatrix0to1[23][10] = 0.047036971037785;
   fWeightMatrix0to1[24][10] = -0.388872105864265;
   fWeightMatrix0to1[0][11] = 0.252922467524328;
   fWeightMatrix0to1[1][11] = 1.29586242661655;
   fWeightMatrix0to1[2][11] = 0.912415021353475;
   fWeightMatrix0to1[3][11] = 0.513154271884854;
   fWeightMatrix0to1[4][11] = -2.8541444312072;
   fWeightMatrix0to1[5][11] = 0.270910399582985;
   fWeightMatrix0to1[6][11] = 0.259137252869647;
   fWeightMatrix0to1[7][11] = -1.45329531111277;
   fWeightMatrix0to1[8][11] = -0.943607988886961;
   fWeightMatrix0to1[9][11] = -2.11605801229057;
   fWeightMatrix0to1[10][11] = -1.68091403760999;
   fWeightMatrix0to1[11][11] = -1.3460674254422;
   fWeightMatrix0to1[12][11] = 2.25142050736153;
   fWeightMatrix0to1[13][11] = -2.26139078439658;
   fWeightMatrix0to1[14][11] = 1.90547438115173;
   fWeightMatrix0to1[15][11] = -1.51972460036703;
   fWeightMatrix0to1[16][11] = 1.0426521169034;
   fWeightMatrix0to1[17][11] = -0.696684743495346;
   fWeightMatrix0to1[18][11] = 0.62204940562633;
   fWeightMatrix0to1[19][11] = 0.986196560118107;
   fWeightMatrix0to1[20][11] = -1.70826530296939;
   fWeightMatrix0to1[21][11] = 0.08052124160443;
   fWeightMatrix0to1[22][11] = 0.51697709741637;
   fWeightMatrix0to1[23][11] = -1.15810313933478;
   fWeightMatrix0to1[24][11] = -1.15560389570385;
   fWeightMatrix0to1[0][12] = -1.0955553618341;
   fWeightMatrix0to1[1][12] = -0.715145817499623;
   fWeightMatrix0to1[2][12] = 1.6315375426903;
   fWeightMatrix0to1[3][12] = -0.626419643409583;
   fWeightMatrix0to1[4][12] = -0.875818377635211;
   fWeightMatrix0to1[5][12] = -0.835061938925823;
   fWeightMatrix0to1[6][12] = -1.12205640036511;
   fWeightMatrix0to1[7][12] = -2.20358958683089;
   fWeightMatrix0to1[8][12] = -1.03974465138247;
   fWeightMatrix0to1[9][12] = -1.74520157099147;
   fWeightMatrix0to1[10][12] = 0.623502883026691;
   fWeightMatrix0to1[11][12] = -1.25364604694864;
   fWeightMatrix0to1[12][12] = -1.05789000918101;
   fWeightMatrix0to1[13][12] = 1.67205906241318;
   fWeightMatrix0to1[14][12] = 0.859093219442392;
   fWeightMatrix0to1[15][12] = 0.00916636421601022;
   fWeightMatrix0to1[16][12] = 2.24126889766074;
   fWeightMatrix0to1[17][12] = -1.57535999085099;
   fWeightMatrix0to1[18][12] = 0.977932407358949;
   fWeightMatrix0to1[19][12] = 1.05691774058469;
   fWeightMatrix0to1[20][12] = -0.372997610669923;
   fWeightMatrix0to1[21][12] = -0.400303906485029;
   fWeightMatrix0to1[22][12] = -0.0675512414327478;
   fWeightMatrix0to1[23][12] = -0.593244891963829;
   fWeightMatrix0to1[24][12] = 0.502260606756272;
   fWeightMatrix0to1[0][13] = 1.76137719613233;
   fWeightMatrix0to1[1][13] = 1.15519278275431;
   fWeightMatrix0to1[2][13] = 1.47468063943002;
   fWeightMatrix0to1[3][13] = -1.40708397111567;
   fWeightMatrix0to1[4][13] = 1.89446484135468;
   fWeightMatrix0to1[5][13] = -1.77282747590463;
   fWeightMatrix0to1[6][13] = 0.720439736193274;
   fWeightMatrix0to1[7][13] = -1.87518546851057;
   fWeightMatrix0to1[8][13] = -1.71739041207869;
   fWeightMatrix0to1[9][13] = -0.521109038032418;
   fWeightMatrix0to1[10][13] = 1.12587820109469;
   fWeightMatrix0to1[11][13] = 0.650842355642048;
   fWeightMatrix0to1[12][13] = 1.2526196400305;
   fWeightMatrix0to1[13][13] = 0.503622131151346;
   fWeightMatrix0to1[14][13] = 0.245676960308865;
   fWeightMatrix0to1[15][13] = -0.199805459736988;
   fWeightMatrix0to1[16][13] = -0.0444468263311241;
   fWeightMatrix0to1[17][13] = 1.96132872168993;
   fWeightMatrix0to1[18][13] = -1.58193083791436;
   fWeightMatrix0to1[19][13] = 1.09873105904035;
   fWeightMatrix0to1[20][13] = -1.68658918702239;
   fWeightMatrix0to1[21][13] = -0.672930034331926;
   fWeightMatrix0to1[22][13] = 0.043929017444397;
   fWeightMatrix0to1[23][13] = -0.181547313845415;
   fWeightMatrix0to1[24][13] = -1.19367472609898;
   fWeightMatrix0to1[0][14] = 0.57616930879542;
   fWeightMatrix0to1[1][14] = 1.58301771821053;
   fWeightMatrix0to1[2][14] = 0.433142475680594;
   fWeightMatrix0to1[3][14] = -0.288516643574946;
   fWeightMatrix0to1[4][14] = -0.602828158834229;
   fWeightMatrix0to1[5][14] = -1.66069450700488;
   fWeightMatrix0to1[6][14] = -1.90650653127566;
   fWeightMatrix0to1[7][14] = 0.833871816320766;
   fWeightMatrix0to1[8][14] = -0.389639618901728;
   fWeightMatrix0to1[9][14] = -1.48257860041873;
   fWeightMatrix0to1[10][14] = 0.138064530784293;
   fWeightMatrix0to1[11][14] = -1.26331153517135;
   fWeightMatrix0to1[12][14] = -1.05861235358888;
   fWeightMatrix0to1[13][14] = -2.41545925822389;
   fWeightMatrix0to1[14][14] = 0.933030247233472;
   fWeightMatrix0to1[15][14] = -0.111816886693405;
   fWeightMatrix0to1[16][14] = -2.13282836613175;
   fWeightMatrix0to1[17][14] = -0.588726819349272;
   fWeightMatrix0to1[18][14] = -0.601754594526975;
   fWeightMatrix0to1[19][14] = -0.680772724998647;
   fWeightMatrix0to1[20][14] = 0.0374085513827674;
   fWeightMatrix0to1[21][14] = -1.62677568884729;
   fWeightMatrix0to1[22][14] = -0.879980940034411;
   fWeightMatrix0to1[23][14] = 1.3971438057492;
   fWeightMatrix0to1[24][14] = 0.575140396632942;
   fWeightMatrix0to1[0][15] = 0.76985624731496;
   fWeightMatrix0to1[1][15] = -0.118085446306233;
   fWeightMatrix0to1[2][15] = 1.09138332258836;
   fWeightMatrix0to1[3][15] = 1.32092047454222;
   fWeightMatrix0to1[4][15] = -0.58204260117175;
   fWeightMatrix0to1[5][15] = -1.06856473989557;
   fWeightMatrix0to1[6][15] = -0.573400419225186;
   fWeightMatrix0to1[7][15] = -2.07669018392238;
   fWeightMatrix0to1[8][15] = 0.818635209954417;
   fWeightMatrix0to1[9][15] = 0.687429845588071;
   fWeightMatrix0to1[10][15] = 0.556505014367859;
   fWeightMatrix0to1[11][15] = 1.0148686165649;
   fWeightMatrix0to1[12][15] = 1.25347987503805;
   fWeightMatrix0to1[13][15] = 1.02945247121446;
   fWeightMatrix0to1[14][15] = -0.33488500706783;
   fWeightMatrix0to1[15][15] = 0.658263704690837;
   fWeightMatrix0to1[16][15] = 0.741384296492645;
   fWeightMatrix0to1[17][15] = -1.79179692495038;
   fWeightMatrix0to1[18][15] = -0.0672068857649077;
   fWeightMatrix0to1[19][15] = -1.1625888452283;
   fWeightMatrix0to1[20][15] = 1.67290111714554;
   fWeightMatrix0to1[21][15] = 0.128162704632534;
   fWeightMatrix0to1[22][15] = 1.4164465611626;
   fWeightMatrix0to1[23][15] = 1.73524086208627;
   fWeightMatrix0to1[24][15] = -0.20393402864811;
   fWeightMatrix0to1[0][16] = 1.88009149387646;
   fWeightMatrix0to1[1][16] = 1.67687253393851;
   fWeightMatrix0to1[2][16] = 0.343437076029305;
   fWeightMatrix0to1[3][16] = -0.999161081469966;
   fWeightMatrix0to1[4][16] = 0.601375508583672;
   fWeightMatrix0to1[5][16] = -1.97097168353231;
   fWeightMatrix0to1[6][16] = 0.406782740522347;
   fWeightMatrix0to1[7][16] = 1.53220873570151;
   fWeightMatrix0to1[8][16] = -0.565860369220689;
   fWeightMatrix0to1[9][16] = -0.555917785611085;
   fWeightMatrix0to1[10][16] = -0.114528090752429;
   fWeightMatrix0to1[11][16] = -1.2424876953009;
   fWeightMatrix0to1[12][16] = 0.419752256272072;
   fWeightMatrix0to1[13][16] = -1.06229741183362;
   fWeightMatrix0to1[14][16] = 0.255781423732214;
   fWeightMatrix0to1[15][16] = 0.28953180187208;
   fWeightMatrix0to1[16][16] = 1.03584530153604;
   fWeightMatrix0to1[17][16] = 1.15000534607886;
   fWeightMatrix0to1[18][16] = 1.35064826699129;
   fWeightMatrix0to1[19][16] = -0.646702451833889;
   fWeightMatrix0to1[20][16] = -0.628761829650464;
   fWeightMatrix0to1[21][16] = 0.368829882201678;
   fWeightMatrix0to1[22][16] = 1.67560370876423;
   fWeightMatrix0to1[23][16] = -1.46161591463925;
   fWeightMatrix0to1[24][16] = -1.34865269604276;
   fWeightMatrix0to1[0][17] = 0.0854612737172713;
   fWeightMatrix0to1[1][17] = -1.45281490635253;
   fWeightMatrix0to1[2][17] = -1.68230118377725;
   fWeightMatrix0to1[3][17] = -1.30293419033639;
   fWeightMatrix0to1[4][17] = 0.272397228517648;
   fWeightMatrix0to1[5][17] = -0.0343232817616589;
   fWeightMatrix0to1[6][17] = -1.65916849028504;
   fWeightMatrix0to1[7][17] = -0.898128869239394;
   fWeightMatrix0to1[8][17] = -0.211145947858921;
   fWeightMatrix0to1[9][17] = 1.01738350037384;
   fWeightMatrix0to1[10][17] = 1.31336414336704;
   fWeightMatrix0to1[11][17] = 1.26601312085536;
   fWeightMatrix0to1[12][17] = -0.197997948730542;
   fWeightMatrix0to1[13][17] = -0.85418884072155;
   fWeightMatrix0to1[14][17] = 0.218482738660716;
   fWeightMatrix0to1[15][17] = -0.265748158480688;
   fWeightMatrix0to1[16][17] = -0.661672719206526;
   fWeightMatrix0to1[17][17] = -0.68816233366662;
   fWeightMatrix0to1[18][17] = 0.0273907116558413;
   fWeightMatrix0to1[19][17] = -1.55624190608493;
   fWeightMatrix0to1[20][17] = 0.691084300849322;
   fWeightMatrix0to1[21][17] = 1.18581613963527;
   fWeightMatrix0to1[22][17] = 0.0154144629474532;
   fWeightMatrix0to1[23][17] = 0.235024973181913;
   fWeightMatrix0to1[24][17] = -0.552348944122127;
   fWeightMatrix0to1[0][18] = 1.8175244652084;
   fWeightMatrix0to1[1][18] = -0.56403517568444;
   fWeightMatrix0to1[2][18] = -0.0652113846345195;
   fWeightMatrix0to1[3][18] = 1.5666639138416;
   fWeightMatrix0to1[4][18] = 0.458972501217469;
   fWeightMatrix0to1[5][18] = -0.580088773980882;
   fWeightMatrix0to1[6][18] = 1.34209637183923;
   fWeightMatrix0to1[7][18] = -2.60006226078823;
   fWeightMatrix0to1[8][18] = 1.53226013101455;
   fWeightMatrix0to1[9][18] = 1.06768796277788;
   fWeightMatrix0to1[10][18] = -0.580811627235158;
   fWeightMatrix0to1[11][18] = -0.201477034225059;
   fWeightMatrix0to1[12][18] = 0.402499936342413;
   fWeightMatrix0to1[13][18] = -1.06962055503249;
   fWeightMatrix0to1[14][18] = 2.29978059772625;
   fWeightMatrix0to1[15][18] = 0.401057582667547;
   fWeightMatrix0to1[16][18] = -1.11932794534545;
   fWeightMatrix0to1[17][18] = -1.55257401947988;
   fWeightMatrix0to1[18][18] = -0.609750571072594;
   fWeightMatrix0to1[19][18] = 2.10741181840975;
   fWeightMatrix0to1[20][18] = -0.753164761325668;
   fWeightMatrix0to1[21][18] = -0.473304102414098;
   fWeightMatrix0to1[22][18] = -1.70797238943855;
   fWeightMatrix0to1[23][18] = 2.04921066663664;
   fWeightMatrix0to1[24][18] = 0.329876723066813;
   fWeightMatrix0to1[0][19] = -0.420203683933854;
   fWeightMatrix0to1[1][19] = -0.466657352424775;
   fWeightMatrix0to1[2][19] = -1.30617777852092;
   fWeightMatrix0to1[3][19] = 0.885599254208418;
   fWeightMatrix0to1[4][19] = 3.79992641037773;
   fWeightMatrix0to1[5][19] = -0.984026265834853;
   fWeightMatrix0to1[6][19] = -0.886916827020692;
   fWeightMatrix0to1[7][19] = 1.41855536762031;
   fWeightMatrix0to1[8][19] = 1.40228707215129;
   fWeightMatrix0to1[9][19] = 0.819675759315498;
   fWeightMatrix0to1[10][19] = 0.783430776209761;
   fWeightMatrix0to1[11][19] = -1.31884876587377;
   fWeightMatrix0to1[12][19] = -0.684039876440252;
   fWeightMatrix0to1[13][19] = 0.98580075837508;
   fWeightMatrix0to1[14][19] = 1.02994944351666;
   fWeightMatrix0to1[15][19] = 0.288790393127533;
   fWeightMatrix0to1[16][19] = -1.73179265235188;
   fWeightMatrix0to1[17][19] = -1.50736126014748;
   fWeightMatrix0to1[18][19] = -1.9936855228445;
   fWeightMatrix0to1[19][19] = -0.611955633645413;
   fWeightMatrix0to1[20][19] = -1.03743824836981;
   fWeightMatrix0to1[21][19] = 0.516059760069911;
   fWeightMatrix0to1[22][19] = -0.209250142229375;
   fWeightMatrix0to1[23][19] = 2.09443517470248;
   fWeightMatrix0to1[24][19] = 0.638516144169301;
   fWeightMatrix0to1[0][20] = -2.61985828497914;
   fWeightMatrix0to1[1][20] = 1.72570059736651;
   fWeightMatrix0to1[2][20] = 2.06555135539767;
   fWeightMatrix0to1[3][20] = 1.70389513476596;
   fWeightMatrix0to1[4][20] = -0.876244979305502;
   fWeightMatrix0to1[5][20] = -0.996466925082449;
   fWeightMatrix0to1[6][20] = 1.29647048550008;
   fWeightMatrix0to1[7][20] = -1.52278569174567;
   fWeightMatrix0to1[8][20] = 0.59109966219761;
   fWeightMatrix0to1[9][20] = -1.44171226186803;
   fWeightMatrix0to1[10][20] = 1.32100168350215;
   fWeightMatrix0to1[11][20] = -1.38061078641469;
   fWeightMatrix0to1[12][20] = -0.435465683351922;
   fWeightMatrix0to1[13][20] = 0.0945795384217279;
   fWeightMatrix0to1[14][20] = -1.57017054656852;
   fWeightMatrix0to1[15][20] = -0.671077890349129;
   fWeightMatrix0to1[16][20] = -0.822072480316464;
   fWeightMatrix0to1[17][20] = 0.40185659969794;
   fWeightMatrix0to1[18][20] = 0.293483239507991;
   fWeightMatrix0to1[19][20] = 0.62412424756977;
   fWeightMatrix0to1[20][20] = 0.902717009275546;
   fWeightMatrix0to1[21][20] = -1.41329662590025;
   fWeightMatrix0to1[22][20] = 1.44724109281529;
   fWeightMatrix0to1[23][20] = -2.72430893853328;
   fWeightMatrix0to1[24][20] = 0.535503932684931;
   // weight matrix from layer 1 to 2
   fWeightMatrix1to2[0][0] = -0.052653091206938;
   fWeightMatrix1to2[0][1] = 0.00259285250896161;
   fWeightMatrix1to2[0][2] = -0.0224187704750374;
   fWeightMatrix1to2[0][3] = -0.0069126794247282;
   fWeightMatrix1to2[0][4] = -0.416060587643179;
   fWeightMatrix1to2[0][5] = -0.039485700654896;
   fWeightMatrix1to2[0][6] = -0.0313640068264573;
   fWeightMatrix1to2[0][7] = -0.057457658940992;
   fWeightMatrix1to2[0][8] = -0.0352390687694822;
   fWeightMatrix1to2[0][9] = -0.0660417659025874;
   fWeightMatrix1to2[0][10] = -0.0772409582493431;
   fWeightMatrix1to2[0][11] = 0.0434953325344349;
   fWeightMatrix1to2[0][12] = 0.16837166712782;
   fWeightMatrix1to2[0][13] = -0.081905641504242;
   fWeightMatrix1to2[0][14] = -0.228444900468299;
   fWeightMatrix1to2[0][15] = 0.158658708515225;
   fWeightMatrix1to2[0][16] = 0.0835367033748966;
   fWeightMatrix1to2[0][17] = 0.00330530233711109;
   fWeightMatrix1to2[0][18] = 0.0515674148040513;
   fWeightMatrix1to2[0][19] = 0.067781547922357;
   fWeightMatrix1to2[0][20] = 0.182851813251914;
   fWeightMatrix1to2[0][21] = 0.019151184571152;
   fWeightMatrix1to2[0][22] = 0.0292825578580646;
   fWeightMatrix1to2[0][23] = 0.461532345000694;
   fWeightMatrix1to2[0][24] = -0.00843330184552723;
   fWeightMatrix1to2[0][25] = 0.833008127309007;
}
