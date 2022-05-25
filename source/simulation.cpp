
//========================================================================
// Simuliation
// 
// Includes:
//    * 2- and 3-body kinematic calculations
//    * detector check for efficiency calculations
//      * silicon telescope and 0-degree ion chamber currently included
//
// Units: 
//    * distance (milimeter)
//    * energy (MeV)
// 
// Author: Eli Temanson
// Date: May 25, 2022
//=================================================================

#ifndef simulation_cpp
#define simulation_cpp

#include "simulation.hpp"
#include <json/json.h>
#include <iomanip> // for std::setprecision

using namespace std;

int main (int argc,char** argv)
{
  if (argc<2) 
  {
    std::cerr<<"Incorrect number of arguments!"<<std::endl;
    std::cout<<"If you need help, add -h or --help for assistance"<<std::endl;
    return 1;
  } else if (strcmp(argv[1], "-h")==0 || strcmp(argv[1], "--help")==0) {
    std::cout<<"Welcome to the inverse kinematics rxn simulation!"<<std::endl;
    std::cout<<"To run this code you need to provide a *.json file as an argument"<<std::endl;
    std::cout<<"Example: ./build/sim input/test.json"<<std::endl;
    std::cout<<"See the README.md for more details"<<std::endl;
    return 1;
  }

  Simulation run01(argv[1]); 
  run01.Run(""); // Put a string inside here to give the histrogram names a different ending
  run01.CalcEfficiency(10000,2.0,5.0,0.01);
  return 0;
}

//=================================================================

Simulation::Simulation(const char* filename)
{
  std::ifstream input(filename);
  if (!input.is_open())
  {
    std::cerr<<"Unable to load configuration in "<<filename<<", check that it exists"<<std::endl;
    exit(0);
  }
  
  Json::Reader reader; // Reader 
  Json::Value sim_input; // The value of value can be any object
  
  if (!reader.parse(input, sim_input))
  {
    std::cerr << "Error parsing: " << filename << "\n" << reader.getFormattedErrorMessages() << std::endl;
    exit(0);
  } else {
      std::cout << "Parsing successful!" << std::endl;
  }

  out_file = TFile::Open(sim_input["out_filename"].asString().c_str(),"RECREATE");
  
  SetAngDis(sim_input["input_dwba"].asString());

  rxnType = sim_input["reaction_type"].asInt();
  kin.SetReaction(sim_input["beam"].asString(),
                  sim_input["target"].asString(),
                  sim_input["ejectile"].asString(),
                  sim_input["fragment"].asString(),
                  sim_input["decay_light"].asString(),
                  sim_input["decay_heavy"].asString());

  nEvents = sim_input["samples"].asInt();

  kin.SetMeanBeamEnergy(sim_input["beam_energy"].asDouble());
  kin.SetBeamSigma(sim_input["beam_energy_sigma"].asDouble());
  kin.SetExcitationEnergy(sim_input["excitation_energy"].asDouble());

  // Close input file
  input.close();

  SiA.Pos.SetXYZ(0.2,-2.0,69.0);
  SiA.Rin=11.53; //11.53
  SiA.Rout=31.0; //35
  SiB.Pos.SetXYZ(0.2,-2.0,82.0);
  SiB.Rin=11.53;
  SiB.Rout=35.0;
}

Simulation::~Simulation()
{
  out_file->cd();
  table1->Write("",TObject::kOverwrite);
  std::cout<<"\nOutput File: "<<out_file->GetName()<<std::endl;
  out_file->Write(out_file->GetName(), TObject::kOverwrite);
  out_file->Close(); 
}

void Simulation::Run(const std::string &s)
{
  state = s;
  Loop();
}

void Simulation::Loop()
{ 
  std::cout << "Number of Events to loop: " << nEvents << "\n";
  for(Long64_t event=0; event < nEvents; event++)
  {  
    // kin.Calc2Body(acos(Rndm->Uniform(-1,1)));
    GetInvKinTheta(); //Sets ThetaCM
    kin.Calc2Body(ThetaCM);

    kin.Calc3Body();

    if (CheckDetector())
    {
      ProcessFill();
    }

    if (event%(nEvents/100) == 0)
    {
      std::cout<<"Percent Completed: "<<(int)100*event/nEvents<<"\r"<<std::flush;
    }
    event+=1;
  }
}

bool Simulation::CheckDetector()
{
  double beam_spot_r = 2.0*TMath::Sqrt(Rndm->Uniform(0,1.0));
  double beam_spot_theta = 2*TMath::Pi()*Rndm->Uniform(0,1.0);

  BeamSpot.SetXYZ(beam_spot_r*TMath::Cos(beam_spot_theta), 
                  beam_spot_r*TMath::Sin(beam_spot_theta), 0);

  if(!InSiDet(SiA,kin.Decay_Light.Pos, BeamSpot)) return false;
  if(!InSiDet(SiB,kin.Decay_Light.Pos, BeamSpot)) return false;
  if(InsideGate(kin.Decay_Light.Phi*TMath::RadToDeg(),66.0,88.0)) return false;
  if(InsideGate(kin.Decay_Light.Phi*TMath::RadToDeg(),-108.0,-89.0)) return false;
  if(!InsideGate(kin.Decay_Light.KE,3.0,11.0)) return false;
  // if(kin.Decay_Heavy.KE < 78.0) return false;

  return true;
}

bool Simulation::InSiDet(const SimDet Silicon, const TVector3 pv, const TVector3 bs)
{
  TVector3 normVect = Silicon.Pos.Unit();
  // normVect.RotateZ(90.0*TMath::DegToRad());

  // first see if it intersects the plane of the detector
  double udotnormv = (pv.Unit()).Dot(normVect);
  if (udotnormv <= 0.) return false;
  // Find distance from target to interesection point of detector plane 
  // use line-plane intersection formula
  double vdist = normVect.Dot(Silicon.Pos)/udotnormv;
  if (vdist <= 0.) return false;
  // vector from target to interesection point of detector plane with
  // magnitude equal to distance
  TVector3 v_to_det(pv);
  v_to_det.SetMag(vdist);
  // create vector from detector origin to interaction point on plane of detector
  TVector3 ch_vect = v_to_det + bs - Silicon.Pos;
  // radial component from center of detector
  double ch_vect_mag = ch_vect.Mag();
  // see if it falls in the detector window
  if (!InsideGate(ch_vect_mag,Silicon.Rin,Silicon.Rout)) return false;

  return true;
}

//=================================================================
// Efficiency Calculation
void Simulation::CalcEfficiency(int n,double elow, double emax, double de)
{ 
  kin.ExcEnergy = elow;
  double passedEvents=0.0;

  std::ofstream eff_data;
  eff_data.open("data/efficieny_data.txt");

  // TGraph* Sim_Eff_Graph = new TGraph();
  TGraphErrors* Sim_Eff_Graph = new TGraphErrors();

  int j=0; 
  while(kin.ExcEnergy < emax)
  {
    passedEvents = 0.0;
    for(Long64_t i=0; i<n; i++)
    {
      GetInvKinTheta(); //Sets ThetaCM
      kin.Calc2Body(ThetaCM);
      // kin.Calc2Body(acos(Rndm->Uniform(-1,1)));
      kin.Calc3Body();
      if(CheckDetector()){
        passedEvents=passedEvents+1.0;
      }
    } // End of for loop

    Sim_Eff_Graph->SetPoint(j, kin.ExcEnergy, passedEvents/n);
    Sim_Eff_Graph->SetPointError(j,0,sqrt(1/passedEvents + 1/n)/n);    
    eff_data<<kin.ExcEnergy<<" "<<passedEvents/n<<" "<<sqrt(1/passedEvents + 1/n)/n<<std::endl; 
    
    kin.ExcEnergy = kin.ExcEnergy + de;
    j = j + 1;

    std::cout<<"Percent Completed: "<<(int)100*(kin.ExcEnergy-elow)/(emax-elow)<<"\r"<<std::flush;
  } // End of while loop

  Sim_Eff_Graph->SetMaximum(1.0);
  Sim_Eff_Graph->SetMinimum(0.0);
  Sim_Eff_Graph->Write("Sim_Eff_Graph");
  // Sim_Eff_Graph->SetMarkerColor(4);
  // Sim_Eff_Graph->SetMarkerStyle(21);
  // Sim_Eff_Graph->SetFillColor(4);
  // Sim_Eff_Graph->SetFillStyle(3005);
  // Sim_Eff_Graph->SetDrawOption("a4"); //ACP

  eff_data.close();
}



void Simulation::ProcessFill()
{
  MyFill(table1,"beam_spot"+state,200,-10,10,BeamSpot.X(),200,-10,10,BeamSpot.Y());

  MyFill(table1,"qval"+state,5000,-4,1,kin.QValue);
  MyFill(table1,"theta_cm"+state,360,0,180,kin.ThetaCM*TMath::RadToDeg());
  MyFill(table1,"ejec_theta_lab"+state,360,0,180,kin.Ejectile.Theta*TMath::RadToDeg());
  MyFill(table1,"frag_theta_lab"+state,2000,0,10,kin.Fragment.Theta*TMath::RadToDeg());
  MyFill(table1,"frag_ke"+state,80,70,110,kin.Fragment.KE);
  MyFill(table1,"ejec_ke"+state,400,0,20,kin.Ejectile.KE);
  MyFill(table1,"frag_ejec_ke_corr"+state,400,60,120,kin.Fragment.KE,400,0,20,kin.Ejectile.KE);

  MyFill(table1,"ejec_ke_angle"+state,720,0,180,kin.Ejectile.Theta*TMath::RadToDeg(),2000,0,20,kin.Ejectile.KE);
  MyFill(table1,"frag_ke_angle"+state,2000,0,10,kin.Fragment.Theta*TMath::RadToDeg(),2000,0,20,kin.Fragment.KE/26.0);
  MyFill(table1,"frag_ke_ejec-angle"+state,2000,70,120,kin.Fragment.KE,720,0,180,kin.Ejectile.Theta*TMath::RadToDeg());
  MyFill(table1,"theta_cm_frag_ke"+state,720,0,180,kin.ThetaCM*TMath::RadToDeg(),400,70,110,kin.Fragment.KE);
  MyFill(table1,"ejec_ke_angle"+state,720,0,180,kin.Ejectile.Theta*TMath::RadToDeg(),2000,0,20,kin.Ejectile.KE);

  MyFill(table1,"decay_light_theta"+state,720,0,180,kin.Decay_Light.Theta*TMath::RadToDeg());
  MyFill(table1,"decay_light_phi"+state,720,-200,200,kin.Decay_Light.Phi*TMath::RadToDeg());
  MyFill(table1,"decay_light_kin"+state,720,0,180,kin.Decay_Light.Theta*TMath::RadToDeg(),800,0,20,kin.Decay_Light.KE);

  MyFill(table1,"decay_heavy_kin"+state,100,0,10,kin.Decay_Heavy.Theta*TMath::RadToDeg(),1200,0,120,kin.Decay_Heavy.KE);
  MyFill(table1,"decay_sum_ke"+state,1200,0,120,kin.Decay_Light.KE+kin.Decay_Heavy.KE);
  MyFill(table1,"decay_sum_ke_n_theta"+state,720,0,180,kin.Ejectile.Theta*TMath::RadToDeg(),1200,0,120,kin.Decay_Light.KE+kin.Decay_Heavy.KE);
  MyFill(table1,"decay_sum_ke_theta_cm"+state,720,0,180,kin.ThetaCM*TMath::RadToDeg(),1200,0,120,kin.Decay_Light.KE+kin.Decay_Heavy.KE);

}
















#endif
