/* Implementation File for GlauberClass
   This contains all of the function definitions for the member
   functions of Glauber Class.

   Author: Chris Flores
   Date: July 10, 2014

*/

//We need to include our class declaration and any declarations of 
//classes we are going to use in this implementation file.
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TRotation.h"

#include "GlauberClass.h"

using namespace std;

ClassImp(GlauberClass); //Another ROOT Necessary Directive

// Default Constructor ___________________________________________________
GlauberClass::GlauberClass(){

  //NULL Constructor - Not all that useful

  crossSection = 0;
  impactParameter = 0;
  nBinaryCollisions = 0;
  nParticipants = 0;
  nucleonHardRadius = 0;

  userWSPars_NucA = NULL;
  userWSPars_NucB = NULL;

  xHardSoft = 0.12; //Default from Physical Review C, 70(2):021902 2004
  
  nNucleons_NucA = 0;
  nNucleons_NucB = 0;

  nParticles = 0;

  rand = new TRandom3(0);

  negativeBinomialDist = new TF1("negBinomial",this,&GlauberClass::NegativeBinomialFunction,0,100,2);
  impactParameterDist  = new TF1("impactParameterDist","2.0*3.14*x",0,100);

  woodsSaxonDist_NucA  = NULL;
  woodsSaxonDist_NucB  = NULL;
  
  woodsSaxon2DPDF_NucA = NULL;
  woodsSaxon2DPDF_NucB = NULL;

}

// Optional Constructor __________________________________________________
GlauberClass::GlauberClass(Int_t nNucleonsA, Int_t nNucleonsB, Double_t nnCrossSection, Double_t nHardRadius){

  crossSection = nnCrossSection;
  nucleonHardRadius = nHardRadius;
  impactParameter = 0;
  nBinaryCollisions = 0;
  nParticipants = 0;

  xHardSoft = 0.12; //Default from Physical Review C, 70(2):021902 2004

  userWSPars_NucA = NULL;
  userWSPars_NucB = NULL;

  nNucleons_NucA = nNucleonsA;
  nucleonPosition_NucA.resize(nNucleons_NucA);

  nNucleons_NucB = nNucleonsB;
  nucleonPosition_NucB.resize(nNucleons_NucB);

  nParticles = 0;

  rand = new TRandom3(0);

  //Distributions for Particle Production Model
  negativeBinomialDist = new TF1("negBinomial",this,&GlauberClass::NegativeBinomialFunction,0,100,2);

  //Impact Parameter
  impactParameterDist = new TF1("impactParameterDist","2.0*3.14*x",0,100);

  //Nucleon Distributions
  woodsSaxonDist_NucA  = NULL;
  woodsSaxonDist_NucB  = NULL;
  woodsSaxon2DPDF_NucA = NULL;
  woodsSaxon2DPDF_NucB = NULL;

  //Create the Woods Saxon Functions for each Nucleus
  CreateWoodsSaxonFunctions(0);
  CreateWoodsSaxonFunctions(1);
  
}

// Destructor ___________________________________________________________
GlauberClass::~GlauberClass(){

  //If we call "new" in our class then we are responsible for deleting
  //that memory allocation.  We won't do that in our class, but here
  //is the right place do do it with al call like:

  delete rand;
  delete negativeBinomialDist;
  delete impactParameterDist;
  delete woodsSaxonDist_NucA;
  delete woodsSaxon2DPDF_NucA;
  delete woodsSaxonDist_NucB;
  delete woodsSaxon2DPDF_NucB;

}

//_________________________________________________________________________
void GlauberClass::CreateWoodsSaxonFunctions(Int_t nucleusIndex){

  //Create the Woods-Saxon Distribution of nucleus nucleusIndex,
  //nucleusIndex = 0 for nucleusA
  //nucleusIndex = 1 for nucleusB

  //Check the nucleusIndex
  if (nucleusIndex != 0 && nucleusIndex != 1){
    cout <<"ERROR: GlauberClass::CreateWoodsSaxonFunctions - nucleusIndex is not 0 or 1" <<endl;
    exit (EXIT_FAILURE);
  }

  //Nucleus A
  if (nucleusIndex == 0){
    woodsSaxonDist_NucA  = new TF2("woodsSaxonNormalized_NucA",this,&GlauberClass::WoodsSaxonNormalizedFunc,
				   0.0,TMath::Pi(),0,18,8,"GlauberClass","WoodsSaxonNormalizedFunc");
    woodsSaxon2DPDF_NucA = new TF2("woodsSaxon2DPDF_NucA",this,&GlauberClass::WoodsSaxon2D_PDF_Func,
				   0.0,TMath::Pi(),0,18,8,"GlauberClass","WoodsSaxon2D_PDF_Func");
    
    woodsSaxonDist_NucA->SetParameters(&GetWoodsSaxonParameters(nNucleons_NucA)[0]);
    woodsSaxon2DPDF_NucA->SetParameters(&GetWoodsSaxonParameters(nNucleons_NucA)[0]);

    woodsSaxonDist_NucA->SetNpx(10000);
  }
  
  //Nucleus B
  else {
    woodsSaxonDist_NucB  = new TF2("woodsSaxonNormalized_NucB",this,&GlauberClass::WoodsSaxonNormalizedFunc,
				   0.0,TMath::Pi(),0,18,8,"GlauberClass","WoodsSaxonNormalizedFunc");
    woodsSaxon2DPDF_NucB = new TF2("woodsSaxon2DPDF_NucB",this,&GlauberClass::WoodsSaxon2D_PDF_Func,
				   0.0,TMath::Pi(),0,18,8,"GlauberClass","WoodsSaxon2D_PDF_Func");

    woodsSaxonDist_NucB->SetParameters(&GetWoodsSaxonParameters(nNucleons_NucB)[0]);
    woodsSaxon2DPDF_NucB->SetParameters(&GetWoodsSaxonParameters(nNucleons_NucB)[0]);

    woodsSaxonDist_NucB->SetNpx(10000);
  }
  
}

// Print Class Info _______________________________________________________
void GlauberClass::PrintClassInfo(){

  cout <<"nNucleons NucA: " <<nNucleons_NucA <<endl
       <<"nNucleons NucB: " <<nNucleons_NucB <<endl
       <<"size of Position Vecs(A,B): " <<nucleonPosition_NucA.size() <<" " 
       <<nucleonPosition_NucB.size() <<endl
       <<"Impact Parameter: " <<impactParameter <<endl
       <<"Npart: " <<nParticipants <<endl
       <<"Ncoll: "<<nBinaryCollisions <<endl;
  
}

// Set Impact Parameter _________________________________________________
void GlauberClass::SetImpactParameter(Double_t val){

  impactParameter = val;
  
}
//_______________________________________________________________________
void GlauberClass::SetImpactParameter(Double_t min, Double_t max){

  if (min>max){
    cout <<"WARNING - GlauberClass::SetImpactParameter() - min range is GREATER than max setting impact parameter to zero\n";
    impactParameter = 0;
  }
  
  impactParameter = impactParameterDist->GetRandom(min,max);
  
}

//________________________________________________________________________
TF2 *GlauberClass::GetWoodsSaxonNormalized(Int_t nucleusIndex){

  if (nucleusIndex != 0 && nucleusIndex != 1){
    cout <<"WARNING - GlauberClass::GetWoodsSaxonNormalized() - nucleusIndex must be 0 or 1! Returning NULL!\n";
    return NULL;
  }
  
  if (nucleusIndex == 0){
    if (woodsSaxonDist_NucA == NULL)
      cout <<"WARNING - GlauberClass::GetWoodsSaxonNormalized() - WoodsSaxon is not defined for nucleus A. Did you use the correct constructor? Returning NULL\n";
    
    return woodsSaxonDist_NucA;
  }
  
  if (woodsSaxonDist_NucB == NULL)
    cout <<"WARNING - GlauberClass::GetWoodsSaxonNormalized() - WoodsSaxon is not defined for nucleus B. Did you use the correct constructor? Returning NULL\n";
  
  return woodsSaxonDist_NucB;
  
}

//________________________________________________________________________
TF2 *GlauberClass::GetWoodsSaxon2DPDF(Int_t nucleusIndex){

  if (nucleusIndex != 0 && nucleusIndex != 1){
    cout <<"WARNING - GlauberClass::GetWoodsSaxon2DPDF() - nucleusIndex must be 0 or 1! Returning NULL!\n";
    return NULL;
  }

  if (nucleusIndex == 0){
    if (woodsSaxon2DPDF_NucA == NULL)
      cout <<"WARNING - GlauberClass::GetWoodsSaxon2DPDF() - WoodsSaxon2DPDF is not defined for nucleus A. Did you use the correct constructor? Returning NULL\n";
    
    return woodsSaxon2DPDF_NucA;
  }
  
  if (woodsSaxon2DPDF_NucB == NULL)
    cout <<"WARNING - GlauberClass::GetWoodsSaxon2DPDF() - WoodsSaxon2DPDF is not defined for nucleus B. Did you use the correct constructor? Returning NULL\n";
  
  return woodsSaxon2DPDF_NucB;
  
}

//_________________________________________________________________________
void GlauberClass::GenerateNucleonPosition(Double_t &x, Double_t &y, Double_t &z,
					   Int_t nucleusIndex, Int_t methodIndex){

  //Generates the location of a single nucleon in a nucleus
  Double_t r(0),theta(0), phi(0);

  //---Uniform Sphere Option
  if (methodIndex == kUniformSphere){
    
    //Compute the Maximum Nuclear Radius based on its A
    Double_t rMax(0);
    if (nucleusIndex == 0)
      rMax = GetNuclearRadius(nNucleons_NucA);
    else
      rMax = GetNuclearRadius(nNucleons_NucB);
    
    //Get a random radial value
    r = rMax * pow(rand->Uniform(0,1),1.0/3.0);
    
    //Sample the surface of the sphere with radius r randomly and return
    //the cartesian coordinates of the nucleon
    rand->Sphere(x,y,z,r);
    return;
  }
  
  //---Woods-Saxon Option
  
  //Generate an azimuthal angle uniformly in phi=[0,2pi)
  phi = rand->Uniform(0,2.0*TMath::Pi());
  
  //The radial and polar angle are generated by sampling the 2D Woods-Saxon PDF
  //where theta=[0,pi) and r=[0,18]
  //Remember the coordiates of the Woods-Saxon are (theta,r) for ROOT reasons
  if (nucleusIndex == 0)
    woodsSaxon2DPDF_NucA->GetRandom2(theta,r);
  else
    woodsSaxon2DPDF_NucB->GetRandom2(theta,r);
  
  //Convert to Cartesian
  x = r * TMath::Sin(theta) * TMath::Cos(phi);
  y = r * TMath::Sin(theta) * TMath::Sin(phi);
  z = r * TMath::Cos(theta);
  
  return;
  
}

//_________________________________________________________________________
Bool_t GlauberClass::IsOverlapping(Int_t n, Int_t nucleusIndex){

  //Compare the distance between the center of the nTh nucleon to the centers
  //of all the nucleons in nucleonArr, If the distance between any centers
  //is less than 2*nucleonHardRadius then return true, otherwise return false

  //If nucleonHardRadius == 0 then overlapping does not matter, return false
  if (nucleonHardRadius == 0)
    return false;
  
  for (Int_t i=0; i<n; i++){

    Double_t distance(0);
    
    if (nucleusIndex == 0){
      distance = sqrt(pow(nucleonPosition_NucA[n].X()-nucleonPosition_NucA[i].X(),2) +
		      pow(nucleonPosition_NucA[n].Y()-nucleonPosition_NucA[i].Y(),2) +
		      pow(nucleonPosition_NucA[n].Z()-nucleonPosition_NucA[i].Z(),2));
    }

    else {
      distance = sqrt(pow(nucleonPosition_NucB[n].X()-nucleonPosition_NucB[i].X(),2) +
		      pow(nucleonPosition_NucB[n].Y()-nucleonPosition_NucB[i].Y(),2) +
		      pow(nucleonPosition_NucB[n].Z()-nucleonPosition_NucB[i].Z(),2));
    }
    
    if (distance <= 2.0*nucleonHardRadius)
      return true;

  }

  return false;

}

// Fill Nucleus with Nucleons _____________________________________________
void GlauberClass::GenerateNucleonPositions(Int_t methodIndex){

  //NOTE: This Uses either a Uniform Sphere Distribution methodIndex = 0
  //           or a Woods-Saxon Distribution methodIndex = 1.

  //NOTE: Continue to generate nucleon location
  //      until the new nucleon location does not overlap with current ones
  
  Double_t  x(0), y(0), z(0);  

  //Loop over the Nucleus A and Fill it with Nucleons
  for (Int_t iNucleon=0; iNucleon<nNucleons_NucA; iNucleon++){

    do {
      
      GenerateNucleonPosition(x,y,z,0,methodIndex);
      nucleonPosition_NucA[iNucleon].SetXYZ(x,y,z);

    } while(IsOverlapping(iNucleon,0));
    
  }

  //Loop over the Nucleus B and Fill it with Nucleons
  for (Int_t iNucleon=0; iNucleon<nNucleons_NucB; iNucleon++){

    do {
      
    GenerateNucleonPosition(x,y,z,1,methodIndex);
    nucleonPosition_NucB[iNucleon].SetXYZ(x,y,z);
    } while(IsOverlapping(iNucleon,1));
    
  }

  //Rotate the Nuclei unless using kUniformSphere method
  if (methodIndex != kUniformSphere){
    RotateNucleons(0);
    RotateNucleons(1);
  }
  
  //Translate Nucleus B Along the x-Direction according to its impact parameter
  for (Int_t iNucleon=0; iNucleon<nNucleons_NucB; iNucleon++){

    nucleonPosition_NucB[iNucleon].SetXYZ(nucleonPosition_NucB[iNucleon].X()+impactParameter,
					  nucleonPosition_NucB[iNucleon].Y(),
					  nucleonPosition_NucB[iNucleon].Z());
        
  }

}

//Rotate Nucleons_______________________________________________________
void GlauberClass::RotateNucleons(Int_t nucleusIndex){

  if (nucleusIndex != 0 && nucleusIndex != 1)
    return;

  //Rotating all the nucleons is computationally expensive, so only do it if
  //there are non-zero deformation Parameters
  Double_t B40(0), B20(0);
  if (nucleusIndex == 0){
    B20 = woodsSaxon2DPDF_NucA->GetParameter(4);
    B40 = woodsSaxon2DPDF_NucA->GetParameter(5);
  }
  else {
    B20 = woodsSaxon2DPDF_NucB->GetParameter(4);
    B40 = woodsSaxon2DPDF_NucB->GetParameter(5);
  }

  if (B20 == 0 && B40 == 0)
    return; //No need to rotate the nucleus is spherically symmetric

  //Generate the rotation angles - these are defined in the x-convention as
  //described in Goldstein Mechanics discussion of Euler angles
  Double_t theta = TMath::ACos(2*rand->Uniform(0,1)-1);
  Double_t phi   = 2.0*TMath::Pi()*rand->Uniform(0,1);
  Double_t psi   = 2.0*TMath::Pi()*rand->Uniform(0,1);

  //The Rotation Matrix
  TRotation r;
  r.SetXEulerAngles(phi,theta,psi);

  //Apply the rotation to all the nucleons
  if (nucleusIndex == 0){
    for (Int_t iNucleon=0; iNucleon<nNucleons_NucA; iNucleon++)
      nucleonPosition_NucA[iNucleon].Transform(r);
  }
  else {
    for (Int_t iNucleon=0; iNucleon<nNucleons_NucB; iNucleon++)
      nucleonPosition_NucB[iNucleon].Transform(r);
  }
  
}

//Get the Nuclear Radius________________________________________________
Double_t GlauberClass::GetNuclearRadius(Int_t nNucleons){

  return 1.25 * pow(nNucleons,1./3.);
  
}

//_____________________________________________________________________________
void GlauberClass::SetWoodsSaxonParameters(Double_t *wsPars_NucA, Double_t *wsPars_NucB){

  //Set the parameters that should be used for the Woods-Saxon Distribution.
  //Once called by the user the parameters here supercede the parameters stored in the code

  //NOTE: The use does not need to worry about obtaining the normalizion parameter p0. The
  //      Wood-Saxon distribution defined by the remaining parameters will be normalized
  //      automatically by ComputeWoodsSaxonNormalizationParameter()
  
  //NOTE: By default wsPars_NucB = wsPars_NucA for ease of use when doing symmetric collision,
  //      but this can be changed by simply passing in an array for wsPars_NucB.

  userWSPars_NucA = wsPars_NucA;

  if (wsPars_NucB != NULL)
    userWSPars_NucB = wsPars_NucB;
  else
    userWSPars_NucB = userWSPars_NucA;

}

// ____________________________________________________________________________
std::vector<double> GlauberClass::GetWoodsSaxonParameters(Int_t nNucleons, Bool_t norm){

  //Sets the Woods-Saxon Parameters based on the type of Nucleus
  //This could probably be done more elegently

  //If adding a new nucleus please use the include the source and list them
  //in Accending order of nucleon number for ease of reading!

  //Parameter Definitions
  //p0 - The user should set this to 1 in all arrays below, it is automatically
  //     normalized by GlauberClass::ComputeWoodsSaxonNormalizationParameter().
  //     It is constructed so that the volume integral of the WS
  //     in three space is equal to the number of nucleons for the nucleus.
  //w  - Controls the central supression of nucleons to account for coulomb
  //     repulsion of protons
  //c  - Half Density Radius. If it is not given in the source, it can be computed
  //     by doing c = 1.07*(nNucleons)^{1/3} per Ramona's Text Book
  //a  - Diffuseness/Skin-Depth. Somtimes model dependent, but if it is necessary
  //     to guess a good starting point is usually .5 fm.
  //B20 - Quadrupole deformation parameter - for spherical set it to 0.0
  //B40 - Hexadecapole deformation parameter - for spherical set it to 0.0
  //gamma2 and gamma4 - entirely model dependent
  //      gamma2 = gamma4 = 0 for Woods-Saxon and modified-c Fermi Distribution
  //      gamma2 = gamma4 = 1 for deformed Fermi Distribution
  //      gamma2 = gamma4 = (c/2)ln(3) for hard-core distribution
  //      for more see: http://dx.doi.org/10.1103/PhysRevC.8.1466 ong page 1475 

  //If you are particularly annal or ambitious or bored have a look  at
  //doi:10.1006/adnd.1995.1002 for more than you ever wanted to know about
  //nuclear shapes and deformation parameters
  
  //Parameters for Various Elements must be of the form:
  //              {p0,     w,     c(fm),  a(fm),   B20,    B40,    gamma2,  gamma4}

  //Unless otherwise noted all values of w, c, and a come from:
  //H. DE Vries et al. Nuclear Charge-Density Distribution Parameters,
  //  Atomeic Data and Nuclear Data Tables, Vol. 36, No. 3 May 1987
  //  doi:10.1016/0092-640X(87)90013-1

  //Only the most abundant natural isotope is listed

  //N 14   :
  Double_t n[]  = {1.0,   0.0,    2.570,  0.5052,  0.0,    0.0,    0.0,     0.0};
  
  //0  16  : 
  Double_t o[]  = {1.0,   -0.051, 2.608,  0.513,   0.0,    0.0,    0.0,     0.0};

  //F 19   :
  Double_t f[]  = {1.0,   0.0,    2.59,   0.564,   0.0,    0.0,    0.0,     0.0};

  //Ne 20  :
  Double_t ne[] = {1.0,   -0.168, 2.791,  0.698,   0.0,    0.0,    0.0,     0.0};

  //Mg 24  :
  Double_t mg[] = {1.0,   -0.163, 3.108,  0.607,   0.0,    0.0,    0.0,     0.0};

  //Al 27  : 
  Double_t al[] = {1.0,    0.0,   3.07,   0.519,   0.0,    0.0,    0.0,     0.0};
  
  //Si 28  : 
  Double_t si[] = {1.0,   -0.233, 3.340,  0.580,   0.0,    0.0,    0.0,     0.0};

  //P 31   :
  Double_t p[]  = {1.0,   -0.173, 3.369,  0.582,   0.0,    0.0,    0.0,     0.0};
  
  //Cl 35  :
  Double_t cl[] = {1.0,   -0.10,  3.476,  0.599,   0.0,    0.0,    0.0,     0.0};
  
  //K 39   :
  Double_t k[]  = {1.0,   -0.201, 3.743,  0.585,   0.0,    0.0,    0.0,     0.0};

  //Ca 40  : 
  Double_t ca[] = {1.0,   -0.161, 3.766,  0.586,   0.0,    0.0,    0.0,     0.0};

  //Ti 48  :
  Double_t ti[] = {1.0,   0.0,    3.843,  0.588,   0.0,    0.0,    0.0,     0.0};
  
  //V 51   :
  Double_t v[]  = {1.0,   0.0,    3.94,   0.505,   0.0,    0.0,    0.0,     0.0};

  //Cr 52  :
  Double_t cr[] = {1.0,   0.0,    4.01,   0.497,   0.0,    0.0,    0.0,     0.0};

  //Mn 55  :
  Double_t mn[] = {1.0,   0.0,    3.89,   0.567,   0.0,    0.0,    0.0,     0.0};
  
  //Fe 56  :
  Double_t fe[] = {1.0,   0.0,    4.106,  0.519,   0.0,    0.0,    0.0,     0.0};
  
  //Co 59  :
  Double_t co[] = {1.0,   0.0,    4.08,   0.569,   0.0,    0.0,    0.0,     0.0};

  //Ni 58  :
  Double_t ni[] = {1.0,   -0.131, 4.3092, 0.5169,  0.0,    0.0,    0.0,     0.0};
  
  //Cu 63  : R. Vogt - Ultra Relativistic Heavy-Ion Collisions
  Double_t cu[] = {1.0,   0.0,    4.214,  0.586,   0.0,    0.0,    0.0,     0.0};

  //Zn 64  :
  Double_t zn[] = {1.0,   0.0,    4.285,  0.584,   0.0,    0.0,    0.0,     0.0};

  //Ge 70  :
  Double_t ge[] = {1.0,   0.0,    4.44,   0.585,   0.0,    0.0,    0.0,     0.0};
  
  //Sr 88  :
  Double_t sr[] = {1.0,   0.0,    4.83,   0.496,   0.0,    0.0,    0.0,     0.0};
  
  //Y 89   :
  Double_t y[]  = {1.0,   0.0,    4.76,   0.571,   0.0,    0.0,    0.0,     0.0};

  //Nb 93  :
  Double_t nb[] = {1.0,   0.0,    4.87,   0.573,   0.0,    0.0,    0.0,     0.0};

  //Cd 114 :
  Double_t cd[] = {1.0,   0.0,    5.40,   0.537,   0.0,    0.0,    0.0,     0.0};

  //In 113 and 115 :
  Double_t in[] = {1.0,   0.0,    5.357,  0.563,   0.0,    0.0,    0.0,     0.0};

  //Sn 120 :
  Double_t sn[] = {1.0,   0.0,    5.315,  0.576,   0.0,    0.0,    0.0,     0.0};

  //Sb 121 and 123 :
  Double_t sb[] = {1.0,   0.0,    5.357,  0.563,   0.0,    0.0,    0.0,     0.0};

  //La 139 :
  Double_t sa[] = {1.0,   0.0,    5.71,   0.535,   0.0,    0.0,    0.0,     0.0};
  
  //Nd 144 :
  Double_t nd[] = {1.0,   0.0,    5.626,  0.6178,  0.0,    0.0,    0.0,     0.0};
  
  //Sm 152 :
  Double_t sm[] = {1.0,   0.0,    5.804,  0.581,   0.0,    0.0,    0.0,     0.0};
  
  //Ho 165 :
  Double_t ho[] = {1.0,   0.0,    6.18,   0.57,    0.0,    0.0,    0.0,     0.0};
  
  //Er 166 :
  Double_t er[] = {1.0,   0.19,   5.98,   0.446,   0.0,    0.0,    0.0,     0.0};

  //Ta 181 :
  Double_t ta[] = {1.0,   0.0,    6.38,   0.64,    0.0,    0.0,    0.0,     0.0};

  //W 184  :
  Double_t w[]  = {1.0,   0.0,    6.51,   0.535,   0.0,    0.0,    0.0,     0.0};
  
  //Au 197 : 
  Double_t au[] = {1.0,    0.0,   6.38,   0.535,   0.0,    0.0,    0.0,     0.0};

  //Pb 208 : R. Vogt - Ultra Relativistic Heavy-Ion Collisions
  Double_t pb[] = {1.0,    0.0,   6.624,  0.549,   0.0,    0.0,    0.0,     0.0};

  //Bi 209 :
  Double_t bi[] = {1.0,   0.0,    0.39,   6.315,   2.881,  0.0,    0.0,     0.0};

  //Th 232 : http://dx.doi.org/10.1103/PhysRevC.16.1686
  Double_t th[] = {1.0,   0.5,    6.76,   0.5,     0.216,  0.111,  0.0,     0.0};
  
  //U 238  : http://dx.doi.org/10.1103/PhysRevC.16.1686
  Double_t u[]  = {1.0,   0.5,    6.82,   0.5,     0.254,  0.052,  0.0,     0.0};

  //Pu 244 : http://dx.doi.org/10.1103/PhysRevC.16.1686
  Double_t pu[] = {1.0,   0.5,    6.87,   0.5,     0.272,  -0.010, 0.0,     0.0};


  Double_t *pars = NULL;

  //If the user has passed in an array for the parameters then use it,
  //otherwise determine which array of parameters to use based on nucleon number
  if (userWSPars_NucA && nNucleons == nNucleons_NucA){
    pars = userWSPars_NucA;
  }
  else if(userWSPars_NucB && nNucleons == nNucleons_NucB){
    pars = userWSPars_NucB;
  }
  else {
    switch (nNucleons) {
    case 14:  pars = n;  break;
    case 16:  pars = o;  break;
    case 19:  pars = f;  break;
    case 20:  pars = ne; break;
    case 24:  pars = mg; break;
    case 27:  pars = al; break;
    case 28:  pars = si; break;
    case 31:  pars = p;  break;
    case 35:  pars = cl; break;
    case 39:  pars = k;  break;
    case 40:  pars = ca; break;
    case 48:  pars = ti; break;
    case 51:  pars = v;  break;
    case 52:  pars = cr; break;
    case 55:  pars = mn; break;
    case 56:  pars = fe; break;
    case 59:  pars = co; break;
    case 58:  pars = ni; break;
    case 63:  pars = cu; break;
    case 64:  pars = zn; break;
    case 70:  pars = ge; break;
    case 88:  pars = sr; break;
    case 89:  pars = y;  break;
    case 93:  pars = nb; break;
    case 114: pars = cd; break;
    case 115: pars = in; break;
    case 120: pars = sn; break;
    case 121: pars = sb; break;
    case 139: pars = sa; break;
    case 144: pars = nd; break;
    case 152: pars = sm; break;
    case 165: pars = ho; break;
    case 166: pars = er; break;
    case 181: pars = ta; break;
    case 184: pars = w;  break;      
    case 197: pars = au; break;      
    case 208: pars = pb; break;
    case 209: pars = bi; break;
    case 132: pars = th; break;
    case 238: pars = u;  break;
    case 244: pars = pu; break;
    default: {
      cout <<"ERROR - GlauberClass::GetWoodsSaxonParameters() - Unable to set Woods Saxon Parameters!\n";
      exit (EXIT_FAILURE);
    }
    }//end Switch
    
  }//end else

  //Create a vector to return
  std::vector<double> paramVec;
  for (int i=0; i<8; i++)
    paramVec.push_back(pars[i]);
  
  //If the normalization parameter is not requested then return the array as is
  if (norm == false){
    return paramVec;
  }
  
  //Set the Normalization Parameter
  paramVec.at(0) = pars[0] = ComputeWoodsSaxonNormalizationParameter(nNucleons);
  
  return paramVec;
}

// Cound Participants and Binary Collisions _______________________________
void GlauberClass::CountNpartNcoll(){

  //Compute nnCrossSection Effective Radius
  Double_t rEff = sqrt(crossSection/(10.0*3.14));

  Double_t nnDistance =0;
  std::vector<int> nucAParticipants (nNucleons_NucA,0);
  std::vector<int> nucBParticipants (nNucleons_NucB,0);

  //Loop Over all the Nucleons in each of the nuclei
  for (Int_t iNucleonA=0; iNucleonA<nNucleons_NucA; iNucleonA++){
    for (Int_t iNucleonB=0; iNucleonB<nNucleons_NucB; iNucleonB++){


      //Compute the Distance Between iNucleonA and iNucleonB, if the distance
      //is less than rEff the nucleons collide
      nnDistance = sqrt(pow(nucleonPosition_NucB[iNucleonB].X()-
			                    nucleonPosition_NucA[iNucleonA].X(),2)+
			            pow(nucleonPosition_NucB[iNucleonB].Y()-
			                    nucleonPosition_NucA[iNucleonA].Y(),2));

      if (nnDistance < rEff){
	nBinaryCollisions++;
	nucAParticipants[iNucleonA] = 1;
	nucBParticipants[iNucleonB] = 1;
	
      }

    }//End Loop over NucA
  }//End Loop over NucB

  nParticipants = std::accumulate(nucAParticipants.begin(),nucAParticipants.end(),0);
  nParticipants = std::accumulate(nucBParticipants.begin(),nucBParticipants.end(),nParticipants);

}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxonFunc(Double_t *x, Double_t *par){

  //Three parameter woods-saxon as developed in:
  //Y.N. Kim, Mesic Atoms and Nuclear Structure, Chapter 3 and
  //W. T. Milner, C. E. Bemis, Jr., and F. K. McGowan, Phys. Rev. C 16, 1686(R) 

  //Radial and Angular Position
  Double_t r      = x[1];
  Double_t theta  = x[0];
  
  //Parameters
  Double_t p0     = par[0]; //Normaliztion Factor 
  Double_t w      = par[1]; //Controls central density
  Double_t c      = par[2]; //Half Density Radius
  Double_t a      = par[3]; //Diffuseness/Skin-Depth (z in Ramona's book)
  Double_t B20    = par[4]; //Quadrupole Deformation Parameter l=2, m=0
  Double_t B40    = par[5]; //Hexadecaople Deformation Parameter l=4, m=0
  Double_t gamma2 = par[6]; //Model dependent Quadrupole shape parameter
  Double_t gamma4 = par[7]; //Model dependent Hexadecapole shape parameter
 
  //Argument of Exponent
  Double_t Y20(0), Y40(0);
  Double_t X(0);
  
  //Spherical Harmonics - only computed if needed
  if (B20 != 0)
    Y20 = 0.25 * sqrt(5./(TMath::Pi())) * (3.0*pow(TMath::Cos(theta),2) - 1.0);
  if (B40 != 0)
    Y40 = (3./16.) * sqrt(1./TMath::Pi()) * (35*pow(TMath::Cos(theta),4) - 30*pow(TMath::Cos(theta),2) + 3.0);
      
  //Define X
  X = (r-c) * (1.0 + B20*Y20 + B40*Y40);
  X = X / (a*(1.0 + B20*gamma2*Y20 + B40*gamma4*Y40));

  Double_t numerator(0);
  if (r < c)
    numerator = 1.0 + w*pow(r/c,2);
  else
    numerator = 1.0 + w;
  
  Double_t denomenator(0);
  denomenator = 1.0 + TMath::Exp(X);

  return p0 * (numerator / denomenator);
  
}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxonNonNormalizedFunc(Double_t *x, Double_t *par){

  //Wrapper for the Woods-Saxon with p0=1;
  par[0] = 1.0;

  return WoodsSaxonFunc(x,par);
  
}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxonNormalizedFunc(Double_t *x, Double_t *par){

  //Wrapper for the Woods-Saxon with p0=correct normalization;

  return WoodsSaxonFunc(x,par);
  
}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxon2D_PDF_Func(Double_t *x, Double_t *par){

  //Function which can be sampled to obtain r and theta components.
  //NOTE the coordinate pairs will be of the form (theta, r),
  //where r=[0,18] and theta=[0,Pi]

  Double_t theta = x[0];
  Double_t r = x[1];

  return TMath::Sin(theta)*r*r*WoodsSaxonFunc(x,par);
}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxonVolumeFunc(Double_t *x, Double_t *par){

  //Function that can be integrated to compute the r and theta
  //components of the volume
  
  Double_t theta = x[0];
  Double_t r = x[1];
  
  return TMath::Sin(theta)*r*r*WoodsSaxonNonNormalizedFunc(x,par);

}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxonGivenRFunc(Double_t *x, Double_t *par){

  //Evaluate the woods saxon distribution of the nucleus according to
  //nucleusIndex=par[1] as a function of theta=x[0] at fixed r=par[1];
  
  Double_t theta = x[0];
  Double_t r = par[0];
  Int_t nucleusIndex = (Int_t)par[1];

  if (nucleusIndex == 0)
    return woodsSaxonDist_NucA->Eval(theta,r);
  
  return woodsSaxonDist_NucB->Eval(theta,r);
  
}

//______________________________________________________________________
TF1 *GlauberClass::GetWoodsSaxon1DGivenR(Double_t r, Int_t nucleusIndex){

  //Return a pointer to a new TF1 object which is the Woods-Saxon Distribution
  //of the nucleus according to nucleusIndex as a function of theta, evaluated
  //at r

  if (nucleusIndex != 0 && nucleusIndex != 1){
    cout <<"WARNING - GlauberClass::GetWoodsSaxon1DGivenR() - nucleus index must be 0 or 1! Returning NULL!";
    return NULL;
  }
  
  TF1 *f = new TF1(Form("WoodsSaxonGivenR_%g",r),this,&GlauberClass::WoodsSaxonGivenRFunc,0,TMath::Pi(),2);
  f->SetParameter(0,r);
  f->SetParameter(1,nucleusIndex);

  return f;

}

//______________________________________________________________________
Double_t GlauberClass::WoodsSaxonGivenThetaFunc(Double_t *x, Double_t *par){

  //Evaluate the Woods-Saxon Distirbution of the nucleus according to
  //nucleusIndex=par[1] as a function of r(x[0]) for a fixed value
  //of theta(par[0]);

  Double_t r = x[0];
  Double_t theta = par[0];
  Int_t nucleusIndex = (Int_t)par[1];

  if (nucleusIndex == 0)
    return woodsSaxonDist_NucA->Eval(theta,r);
  
  return woodsSaxonDist_NucB->Eval(theta,r);
  
}

//______________________________________________________________________
TF1 *GlauberClass::GetWoodsSaxon1DGivenTheta(Double_t theta, Int_t nucleusIndex){

  //Return a pointer to a new TF1 object which is the Woods Saxon distribution
  //of the nucleus according to nucleusIndex as a function of r, evaluated
  //at theta

  if (theta < 0 || theta >= TMath::Pi()){
    cout <<"WARNING - GlauberClass::GetWoodsSaxon1DGivenTheta() - theta bust bin in the range [0,pi). Setting theta to 0\n";
    theta = 0.0;
  }

  if (nucleusIndex != 0 && nucleusIndex != 1){
    cout <<"WARNING - GlauberClass::GetWoodsSaxon1DGivenTheta() - nucleus index must be 0 or 1! Returning NULL!";
    return NULL;
  }
  
  TF1 *f = new TF1(Form("WoodsSaxonGivenTheta_%d_%g",nucleusIndex,theta),this,
		   &GlauberClass::WoodsSaxonGivenThetaFunc,0,18,2);
  f->SetParameter(0,theta);
  f->SetParameter(1,(Double_t)nucleusIndex);
  return f;

}

//______________________________________________________________________
Double_t GlauberClass::ComputeWoodsSaxonNormalizationParameter(Int_t nNucleons){

  //Computes the p0 parameter such that
  // p0 = nNucleons/VolIntegral

  TF2 temp("woodsSaxonNonNormalized",this,&GlauberClass::WoodsSaxonVolumeFunc,
	   0,TMath::Pi(),0,18,8,"GlauberClass","WoodsSaxonVolumeFunc");
  temp.SetParameters(&GetWoodsSaxonParameters(nNucleons,false)[0]);

  Double_t thetaMin, rMin, thetaMax, rMax;
  temp.GetRange(thetaMin,rMin,thetaMax,rMax);

  Double_t integral = temp.Integral(thetaMin,thetaMax,rMin,rMax);
  integral = integral*2.0*TMath::Pi();

  Double_t p0 = ((double)nNucleons) / integral;

  return p0;
  
}

//______________________________________________________________________
Double_t GlauberClass::ChooseFunction(Double_t n, Double_t k){

  //Using the Gamma Function Definition from Wikipedia

  return TMath::Exp(TMath::LnGamma(n+1)-TMath::LnGamma(k+1)-TMath::LnGamma(n-k+1));
}

//______________________________________________________________________
Double_t GlauberClass::NegativeBinomialFunction(Double_t *x, Double_t *par){

  Double_t Nch   = x[0];   //Number of Charged Particles by pp collision
  Double_t npp   = par[0]; //Average number of produced particles in pp collision
  Double_t k     = par[1]; //Related to the variance of the number of particles in pp collision

  return ChooseFunction(Nch+k-1,k-1) * pow((npp/k)/(1.+(npp/k)),Nch) * pow(1.+(npp/k),-k);
}

//______________________________________________________________________
void GlauberClass::SetNegativeBinomialParameters(Double_t Npp, Double_t K){

  //Set the Parmeters of the Negative Binomial Distribution
  negativeBinomialDist->SetParameters(Npp,K);
}

//______________________________________________________________________
Double_t GlauberClass::GetEfficiency(Double_t idealMultiplicity){

  //Return a multiplicity depentent value of detector efficiency.
  //This efficiency changes as a function of multiplicity to model
  //the worsening tracking efficiency of the STAR TPC with increasing
  //occupancy. The value goes from 90% for low multiplicity events
  //which is consistent with the observed efficiency in pp collisions,
  //to roughly 70% which is the observed efficiency in central AuAu
  //events at 200 GeV.

  //Since this value only makes sense for idealMultiplicities up to
  //560 return 70% for idealMultiplicities that are greater
  if (idealMultiplicity > 560)
    return .70;

  //Return the Efficiency
  return 0.9 * (1.0 - (idealMultiplicity/2520));
  
}

//______________________________________________________________________
Double_t GlauberClass::ProduceParticles(Bool_t useTrackingEfficiency){

  //Simulate the observed multiplicity for this event. To do this we
  //sample the NB distribution m times where m controls the hard->soft
  //scaling and is defined as:
  Double_t m = xHardSoft*nBinaryCollisions + (1.0-xHardSoft)*nParticipants;
  //Where the value xHardSoft = .12 (default) is a good description for BES energies
  //per Physical Review C, 70(2):021902 2004.  

  //Generate the ideal multiplicity where ideal means that there is
  //no detector efficiency drop off due to occupancy effects
  Double_t totalParticlesProduced(0);
  for (Int_t iBinColl=0; iBinColl<TMath::CeilNint(m); iBinColl++){
    totalParticlesProduced +=negativeBinomialDist->GetRandom(0,35);
  }

  //Use tracking efficiency if requested
  if (useTrackingEfficiency){
    
    Double_t idealMultiplicity = totalParticlesProduced;
    Double_t efficiency = GetEfficiency(idealMultiplicity);

    //For each particle produced choose a random number, if
    //the number is less or equal to the efficiency count it,
    //otherwise reduce the totalParticles Produced by one.
    for (Int_t iParticle=0; iParticle < TMath::Nint(idealMultiplicity); iParticle++){
      if (rand->Uniform() > efficiency)
	totalParticlesProduced--;
    }

  }
  
  SetTotalParticlesProduced(totalParticlesProduced);

  return totalParticlesProduced;
  
}

//__________________________________________________________________________
Double_t GlauberClass::ProduceParticles(Double_t x, Bool_t useTrackingEfficiency){

  //Produce Particles with a user defined hard->Soft scale
  if (x > 0 && x < 1.0){
    xHardSoft = x;
  }
  else if (x != -1){
    cout <<"WARNING: GlauberClass::ProduceParticles() - The Hard->Soft Scaling factor must be in [0,1]. NO PARTICLES PRODUCED!" <<endl;
    return 0;
  }

  
  return ProduceParticles(useTrackingEfficiency);
  
}

//________________________________________________________________________
Double_t GlauberClass::ProduceParticles(Int_t nPart, Int_t nColl, Double_t x, Bool_t useTrackingEfficiency){

  //This function can be used do the particle production for any nPart and nColl, and x combination.
  //Its useful if you already know nPart and nColl and just want to do a quick particle production for
  //either testing or some parallel operations.

  //Set the class's nBinaryCollisions and nParticipants
  nBinaryCollisions = nColl;
  nParticipants     = nPart;

  return ProduceParticles(x,useTrackingEfficiency);
  
}

// Reset the Necessary Variables _________________________________________
void GlauberClass::Reset(){

  impactParameter = 0;
  nBinaryCollisions = 0;
  nParticipants = 0;

}


