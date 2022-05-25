
#ifndef kinematics_hpp
#define kinematics_hpp

//C and C++ libraries.
#include <cmath>
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"

struct Particle {
    double M,InvMass,KE,E,Ecm;
    double X,Y,Z;
    double Theta,Phi;
    double Px,Py,Pz,P;
    TLorentzVector LV;
    TVector3 Pos;
};

class Kinematics{
protected:
  // double mass[6];
  double FragEst=0.0;
  double BeamSigma=0.0, MeanBeamEnergy=0.0;
  TRandom3 *Rndm = new TRandom3();
  
public:
  Kinematics(){};
  ~Kinematics(){};

  Particle Beam,Target,Ejectile,Fragment,Decay_Light,Decay_Heavy;

  void SetReaction(const std::string& beam,const std::string& target,const std::string& ejectile,const std::string& fragment,const std::string& product,const std::string& heavy);
 
  double QValue = 0.0;
  double DecayQValueEst = 0.0;
  double DecayQValue = 0.0;
  double ExcEnergy = 0.0;
  double InvMassExcEnergy = 0.0;
  double RecoilQValue = 0.0;
  double ThetaCM = 0.0;
  double Sp = 0.0;
  double gamma,term1,term2,term3;
  double Ecm,Pcm;

  void SetMeanBeamEnergy(const double&);
  void SetBeamSigma(const double&);
  void SetExcitationEnergy(const double&); 
  double SeparationEnergy();

  double r,s,Pcm_in,Pcm_out,rapidity;
  bool Calc2Body(const double&);
  bool Calc3Body();

};


#endif
