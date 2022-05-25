
//==================================================================================
// Physics Calculations 
// 0 + 1 -> 2 + 3* -> 2 + (4+5)
// Beam + Target -> Ejectile + Fragment* -> Ejectile + (Decay_Light + Decay_Heavy)
// All Everything should be in MeV via natural units
//
// Author: Eli Temanson
// Date: May 25, 2022
//==================================================================================

#ifndef kinematics_cpp
#define kinematics_cpp

#include "kinematics.hpp"
#include <json/json.h>

using namespace std;

void Kinematics::SetMeanBeamEnergy(const double& e)
{ 
  MeanBeamEnergy = e;
}
void Kinematics::SetBeamSigma(const double& esig)
{ 
  BeamSigma=esig;
}
void Kinematics::SetExcitationEnergy(const double& e)
{ 
  ExcEnergy = e; 
}
double Kinematics::SeparationEnergy()
{ 
  Sp = Decay_Heavy.M + Decay_Light.M - Fragment.M; 
  return Sp;
}

//Set Primary Reaction and automatically add one decay channel
void Kinematics::SetReaction(const string &beam,const string &target,const string &ejectile,const string &fragment,const string &light,const string &heavy)
{  
  std::cout<<"Setting Reaction"<<std::endl;
  
  std::ifstream file_masses("files/nuclear_masses.json");
  if (!file_masses.is_open())
  {
    std::cerr<<"Unable to load configuration in: files/nuclear_masses.json, check that it exists"<<std::endl;
    exit(0);
  }
  
  Json::Reader reader; // Reader 
  Json::Value nuc_data; // The value of value can be any object
  
  if (!reader.parse(file_masses, nuc_data))
  {
    std::cerr << "Error parsing: files/nuclear_masses.json \n" << reader.getFormattedErrorMessages() << std::endl;
    exit(0);
  } else {
      std::cout << "Parsing successful!" << std::endl;
  }

  // double conv = 931.49410242*1e-6; // 1 amu to MeV/c2 
  Beam.M = nuc_data[beam]["isotope_mass"].asDouble();
  Target.M = nuc_data[target]["isotope_mass"].asDouble();
  Ejectile.M = nuc_data[ejectile]["isotope_mass"].asDouble();
  Fragment.InvMass = nuc_data[fragment]["isotope_mass"].asDouble();
  Decay_Heavy.M = nuc_data[heavy]["isotope_mass"].asDouble();

  if (light=="p")
  { 
    std::cerr << "Error!!! Use 1H not p for protons" << std::endl; 
    std::cout << "Replacing it with 1H and proceeding..." << std::endl;
    Decay_Light.M = nuc_data["1H"]["isotope_mass"].asDouble();
  } else {
    Decay_Light.M = nuc_data[light]["isotope_mass"].asDouble();
  }
  
  Fragment.M = Fragment.InvMass;
}

bool Kinematics::Calc2Body(const double& th)
{ 
  ThetaCM = th;
  QValue = Beam.M + Target.M - (Fragment.M + Ejectile.M + ExcEnergy);
  Beam.KE = Rndm->Gaus(MeanBeamEnergy, BeamSigma); 

  term1 = Beam.KE*Beam.M*Ejectile.M;
  term2 = Fragment.M*(Ejectile.M + Fragment.M)*QValue;
  term3 = Fragment.M*(Fragment.M + Ejectile.M - Beam.M)*Beam.KE;
  gamma = sqrt(term1/(term2+term3));

  // Lab Frame Ejectile Theta, Iliadis(C.38)
  Ejectile.Theta = acos((gamma + cos(ThetaCM))/sqrt(1 + gamma*gamma + 2*gamma*cos(ThetaCM)));
  r = sqrt(Beam.KE*Beam.M*Ejectile.M)*cos(Ejectile.Theta)/(Ejectile.M+Fragment.M); 
  s = (Beam.KE*(Fragment.M-Beam.M)+Fragment.M*QValue)/(Ejectile.M+Fragment.M);

  if ((r+sqrt(r*r+s)) > 0) // Iliadis(C.5-6)
  {
    Ejectile.KE = (r+sqrt(r*r+s))*(r+sqrt(r*r+s));
  }else{
    Ejectile.KE = (r-sqrt(r*r+s))*(r-sqrt(r*r+s));
  }

  // Ejectile.P = sqrt(Ejectile.KE*(Ejectile.KE + 2.0*Ejectile.M));
  Ejectile.E = Ejectile.KE + Ejectile.M;
  Ejectile.P = sqrt(Ejectile.E*Ejectile.E - Ejectile.M*Ejectile.M);
  Ejectile.Phi = Rndm->Uniform(0.0,2*TMath::Pi());
  Ejectile.Px=Ejectile.P*cos(Ejectile.Phi)*sin(Ejectile.Theta);
  Ejectile.Py=Ejectile.P*sin(Ejectile.Phi)*sin(Ejectile.Theta);
  Ejectile.Pz=Ejectile.P*cos(Ejectile.Theta);
  Ejectile.LV.SetPxPyPzE(Ejectile.Px, Ejectile.Py, Ejectile.Pz, Ejectile.E);

  // Assign Beam assuming it's on axis
  Beam.E = Beam.KE + Beam.M;
  Beam.P = sqrt(Beam.E*Beam.E - Beam.M*Beam.M);
  Beam.LV.SetPxPyPzE(0,0,Beam.P,Beam.E);
  
  // Assign Stationary target
  Target.LV.SetPxPyPzE(0,0,0,Target.M);

  // There are a few ways to now calculate the kinematics of the fragment 
  Fragment.Theta = atan(sin(Ejectile.Theta)/(sqrt(Beam.M*Beam.KE/(Ejectile.M*Ejectile.KE))-cos(Ejectile.Theta)));
  Fragment.KE=(Ejectile.M*Ejectile.KE/Fragment.M)*pow(sin(Ejectile.Theta)/sin(Fragment.Theta),2);
  Fragment.E = Fragment.KE + Fragment.M;
  Fragment.LV = Beam.LV + Target.LV - Ejectile.LV;
  Fragment.Theta = Fragment.LV.Theta();
  Fragment.Phi = Fragment.LV.Phi();

  // Fragment.E = Fragment.LV.E();
  // Fragment.KE = Fragment.E - Fragment.M - ExcEnergy;
  // Fragment.Phi = Rndm->Uniform(0.0,2*TMath::Pi());
  // Fragment.P = sqrt(Fragment.KE*(Fragment.KE + 2.0*Fragment.M));
  // Fragment.Px=Fragment.P*cos(Fragment.Phi)*sin(Fragment.Theta);
  // Fragment.Py=Fragment.P*sin(Fragment.Phi)*sin(Fragment.Theta);
  // Fragment.Pz=Fragment.P*cos(Fragment.Theta);
  // Fragment.LV.SetPxPyPzE(Fragment.Px, Fragment.Py, Fragment.Pz, Fragment.E);

  return true;
}


// Now calculate 3 body break up of the fragment in CM frame
bool Kinematics::Calc3Body()
{
  TVector3 boost = Fragment.LV.BoostVector();
  Fragment.LV.Boost(-boost); //Boost the Fragment to it's center of mass
  
  Ecm = (Fragment.LV.E()*Fragment.LV.E() + 
        Decay_Light.M*Decay_Light.M - 
        Decay_Heavy.M*Decay_Heavy.M) /(2.0*Fragment.LV.E());

  Pcm = sqrt(Ecm*Ecm - Decay_Light.M*Decay_Light.M);

  Decay_Light.Theta = acos(Rndm->Uniform(-1,1)); //in CM frame
  Decay_Light.Phi = Rndm->Uniform(0,2*TMath::Pi());

  Decay_Light.LV.SetPx(Pcm*cos(Decay_Light.Phi)*sin(Decay_Light.Theta));
  Decay_Light.LV.SetPy(Pcm*sin(Decay_Light.Phi)*sin(Decay_Light.Theta));
  Decay_Light.LV.SetPz(Pcm*cos(Decay_Light.Theta));
  Decay_Light.LV.SetE(Ecm);

  Decay_Heavy.LV = Fragment.LV - Decay_Light.LV;

  Fragment.LV.Boost(boost); // Boost back to the lab frame
  Decay_Light.LV.Boost(boost); // Boost back to the lab frame
  Decay_Heavy.LV.Boost(boost); // Boost back to the lab frame

  Decay_Light.E = Decay_Light.LV.E();
  Decay_Light.KE = Decay_Light.E - Decay_Light.M;
  Decay_Light.Theta = Decay_Light.LV.Theta();
  Decay_Light.Phi = Decay_Light.LV.Phi();

  Decay_Heavy.E = Decay_Heavy.LV.E();
  Decay_Heavy.KE = Decay_Heavy.E - Decay_Heavy.M;
  Decay_Heavy.Theta = Decay_Heavy.LV.Theta();
  Decay_Heavy.Phi = Decay_Heavy.LV.Phi();

  Decay_Light.Pos.SetMagThetaPhi(1.0,Decay_Light.Theta,Decay_Light.Phi);

  return true;
}






#endif
