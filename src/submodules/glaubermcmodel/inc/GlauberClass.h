/* Header File for GlauberClass
   This contains the class declaration and all the declarations
   of its data and function members

   Author: Chris Flores
   Date: July 10, 2014

*/

//Compiler Directives - If this file has already been compiled don't
//don't bother compiling it again. (This is good practice and is most
//important if you are working with a large codebase.)
#ifndef GLAUBERCLASS_H
#define GLAUBERCLASS_H

//We need to tell the compiler what external classes we are going to use
//This is called a forward declaration.
class TVector3;
class TRandom3;
class TF1;
class TF2;

enum nucleonDistributions {kUniformSphere, kWoodsSaxon};

//GlauberClass Declaration
class GlauberClass : public TObject{

  //DATA MEMBERS
  //These are typically private as they are here, but they don't have to be.
 private:

  Double_t crossSection;      //The nucleon-nucleon inelastic cross section
  Double_t impactParameter;   //Distance between the centers of the two nuclei along the x-axis
  Double_t nucleonHardRadius; //2*nucleonHardRadius is the minimum distance required between
                              //nucleons in the same nucleus (default = 0 E.G. nucleons are allowed full overlap)
  
  Int_t nBinaryCollisions;    //The number of nucleon-nucleon binary collisions (computed)
  Int_t nParticipants;        //The number of nucleons with at least one binary collision (computed)
  
  //User option to pass in parameters for WoodsSaxon - Does not need to be normalized
  //     useful for systematic error studies.
  Double_t *userWSPars_NucA;                 //!
  Double_t *userWSPars_NucB;                 //!

  //Hard->Soft Scaling Factor
  Double_t xHardSoft;
  
  //Lets list the variables for each of the two nuclei seperately
  //NucleusA
  Int_t nNucleons_NucA;
  std::vector<TVector3> nucleonPosition_NucA; //!

  //NucleusB
  Int_t nNucleons_NucB;
  std::vector<TVector3> nucleonPosition_NucB; //!

  //Random Number Generator
  TRandom3 *rand; //!

  //Distributions
  TF1 *negativeBinomialDist; //!
  TF1 *impactParameterDist;  //!
  TF1 *uniformSphereDist;    //!
  TF2 *woodsSaxonDist_NucA;  //!
  TF2 *woodsSaxonDist_NucB;  //!
  TF2 *woodsSaxon2DPDF_NucA; //!
  TF2 *woodsSaxon2DPDF_NucB; //!
  
  //Woods-Saxon Utilities
  void RotateNucleons(Int_t nucleusIndex);
  void CreateWoodsSaxonFunctions(Int_t nucleusIndex);
  Bool_t IsOverlapping(Int_t n, Int_t nucleusIndex);
  Double_t WoodsSaxonFunc(Double_t *x, Double_t *par);
  Double_t WoodsSaxon2D_PDF_Func(Double_t *x, Double_t *par);
  Double_t WoodsSaxonNormalizedFunc(Double_t *x, Double_t *par);
  Double_t WoodsSaxonNonNormalizedFunc(Double_t *x, Double_t *par);
  Double_t WoodsSaxonVolumeFunc(Double_t *x, Double_t *par);
  Double_t ComputeWoodsSaxonNormalizationParameter(Int_t nNucleons);
  Double_t WoodsSaxonGivenRFunc(Double_t *x, Double_t *par);
  Double_t WoodsSaxonGivenThetaFunc(Double_t *x, Double_t *par);
  
  //For the Particle Production Model
  Double_t ChooseFunction(Double_t n, Double_t k); //!
  Double_t NegativeBinomialFunction(Double_t *x, Double_t *par); //!
  
  //Particle Production
  Double_t nParticles; //!

  //MEMBER FUNCTIONS
 public:

  //The constructor(s) and destructor - ABSOLUTELY NECESSARY
  GlauberClass(); //Default Constructor
  GlauberClass(Int_t nNucleonsA, Int_t nNucleonsB, Double_t nnCrossSection, Double_t nHardRadius=0);
  ~GlauberClass();

  
  //Get Member Functions
  Int_t GetNNucleonsNucA(){return nNucleons_NucA;}
  Int_t GetNNucleonsNucB(){return nNucleons_NucB;}
  Double_t GetInelasticCrossSection(){return crossSection;}
  Double_t GetImpactParameter(){return impactParameter;}
  Double_t GetHardSoftScale(){return xHardSoft;}
  std::vector<double> GetWoodsSaxonParameters(Int_t nNucleons,Bool_t norm=true);
  Int_t GetNParticipants(){return nParticipants;}
  Int_t GetNBinaryCollisions(){return nBinaryCollisions;}
  TF2 *GetWoodsSaxonNormalized(Int_t nucleusIndex);
  TF2 *GetWoodsSaxon2DPDF(Int_t nucleusIndex);
  TF1 *GetWoodsSaxon1DGivenR(Double_t r, Int_t nucleusIndex);
  TF1 *GetWoodsSaxon1DGivenTheta(Double_t theta, Int_t nucleusIndex);
  TF1 *GetNegativeBinomialDistribution(){return negativeBinomialDist;}
  Double_t GetNuclearRadius(Int_t nNucleons);
  
  //Set Member Functions
  void SetNNucleonsNucA(Int_t val){nNucleons_NucA=val;}
  void SetNNucleonsNucB(Int_t val){nNucleons_NucB=val;}
  void SetHardSoftScale(Double_t val){xHardSoft = val;}
  void SetImpactParameter(Double_t val);
  void SetImpactParameter(Double_t min, Double_t max);
  void SetWoodsSaxonParameters(Double_t *wsPars_NucA, Double_t *wsPars_NucB = NULL);
  void SetNegativeBinomialParameters(Double_t Npp, Double_t K);
  void SetTotalParticlesProduced(Double_t val){nParticles=val;}

  //Other Useful Memeber Functions
  void PrintClassInfo();
  void GenerateNucleonPosition(Double_t &x, Double_t &y, Double_t &z,
			       Int_t nucleusIndex, Int_t methodIndex); //Single Nucleon
  void GenerateNucleonPositions(Int_t methodIndex); //Fill Nucleus
  void CountNpartNcoll();
  Double_t GetEfficiency(Double_t idealMultiplicity);
  Double_t ProduceParticles(Bool_t useTrackingEfficiency=false);
  Double_t ProduceParticles(Double_t x, Bool_t useTrackingEfficiency=false);
  Double_t ProduceParticles(Int_t nPart, Int_t nColl, Double_t x=-1, Bool_t useTrackingEfficiency=false);
  void Reset();
  
  ClassDef(GlauberClass,1); //This is a ROOT specific directive


}; //A very easy to forget semicolon! ;)


#endif
