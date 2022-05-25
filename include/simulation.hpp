//---------------------------------------------------------
// Declare anything that is used specifically in the analysis
// such as referencing the raw data, I reference global ones
// in the RN namespace in RNCore.cpp/hpp
//---------------------------------------------------------

#ifndef simulation_hpp
#define simulation_hpp

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TCutG.h"
#include "TF1.h"
#include "TSpline.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "global.hpp"
#include "kinematics.hpp"


class Simulation{

  struct SimDet{
    double X,Y,Z,Rin,Rout,zrot;
    TVector3 Pos;
  };

private:
  TRandom3 *Rndm = new TRandom3();
  double rxnType=0;
  std::string state,outfile_name, DetConfig_path;
  std::string rxnString[6]; 

public:
  Simulation(TFile*,THashTable*);
  Simulation(const char*);
  ~Simulation();
  
  Kinematics kin;
  SimDet SiA,SiB;

  THashTable *table1 = new THashTable();
  TFile *out_file = nullptr;

  std::vector<double> ad,angle;
  int nEvents;
  double ThetaCM;  
  TVector3 BeamSpot;

  void CalcEfficiency(int,double,double,double);
  bool CheckDetector();
  bool InSiDet(const SimDet, const TVector3, const TVector3);
  void Run(const std::string &);
  void Reset();
  void ProcessFill();
  void Loop();
  void Close();


  TGraph *angdis = new TGraph();
  double GetAngDis(double x){
    return angdis->Eval(x);
  }

  void GetInvKinTheta(){
    do{ //from dwba, randomize neutron events going from regular to inverse kin (180-theta)
      ThetaCM = acos(Rndm->Uniform(-1,1)); // create a uniform distribution in radians 
    }while (Rndm->Uniform(0,1) > GetAngDis(180.0-ThetaCM*TMath::RadToDeg())); // Weight that distribution by the Normalized Ang Dis. 
  }

  void SetAngDis(const std::string &filename){
    std::ifstream input(filename.c_str()); 
    if(!input.is_open()){
      std::cout<<"Could not Load Angular Distribution from: "<< filename << "\n";
      exit(0);
    }else{
      std::cout<<"Loaded Angular Distribution from: "<< filename << "\n";
    }
    std::vector<double> vtemp_angle, vtemp_ad;
    double temp_angle,temp_ad;
    double max_ad=0.0;
    while(input >> temp_angle >> temp_ad){
      if(temp_ad > max_ad){ max_ad = temp_ad; }
      vtemp_angle.push_back(temp_angle);
      vtemp_ad.push_back(temp_ad);
    }
    ad.clear();
    angle.clear();
    for(unsigned int i=0;i<vtemp_ad.size();i++){
      ad.push_back(vtemp_ad[i]/max_ad);
      angle.push_back(vtemp_angle[i]);

      angdis->SetPoint(i+1,angle[i],ad[i]);
    }
    input.close();
  }

}; // End of Simulation Class



// /*2D histogram fill wrapper*/
// inline void MyFill(THashTable *rootObj,string name,int binsx,double minx,double maxx,double valuex,int binsy,double miny,double maxy,double valuey){
//   TH2F *histo = (TH2F*) rootObj->FindObject(name.c_str());
//   if(histo != NULL) {
//     histo->Fill(valuex, valuey);
//   } else {
//     TH2F *h = new TH2F(name.c_str(), name.c_str(), binsx, minx, maxx, binsy, miny, maxy);
//     h->Fill(valuex, valuey);
//     rootObj->Add(h);
//   }
// }
// /*1D histogram fill wrapper*/
// inline void MyFill(THashTable *rootObj,string name,int binsx,double minx,double maxx,double valuex){
//   TH1F *histo = (TH1F*) rootObj->FindObject(name.c_str());
//   if(histo != NULL) {
//     histo->Fill(valuex);
//   } else {
//     TH1F *h = new TH1F(name.c_str(), name.c_str(), binsx, minx, maxx);
//     h->Fill(valuex);
//     rootObj->Add(h);
//   }
// }
// inline bool InsideGate(double val,double low,double high){
//   if(val > low && val < high){ return true; }
//   else{ return false; }
// }

#endif
