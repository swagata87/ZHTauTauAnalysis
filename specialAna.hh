#ifndef specialAna_hh
#define specialAna_hh

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

///clean up the header!!!
#include "Pxl/Pxl/interface/pxl/core.hh"
#include "Pxl/Pxl/interface/pxl/hep.hh"
#include "Tools/PXL/Sort.hh"
//#include "Tools/Tools.hh"
#include "Tools/MConfig.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TLorentzVector.h"
#include <TFile.h>

#include "Main/Systematics.hh"
#include "Main/EventSelector.hh"


//----------------------------------------------------------------------
using namespace std;

//class Systematics;

class specialAna : public pxl::AnalysisProcess  {
public:
  specialAna( const Tools::MConfig &config,
	      Systematics *syst_shifter,
	      EventSelector *selector);
  virtual ~specialAna();
  
  virtual void endJob(const Serializable*);
  virtual void analyseEvent( const pxl::Event* event );
  bool Check_Muo_ID_Strict(pxl::Particle* muon, bool do_pt_cut = true, bool do_eta_cut = true);
  bool Check_Muo_ID_Mild(pxl::Particle* muon, bool do_pt_cut = true, bool do_eta_cut = true);
  bool Check_Tau_ID(pxl::Particle* tau, bool do_pt_cut = true, bool do_eta_cut = true);
  bool TriggerSelector(const pxl::Event* event);

  TFile* file1;
  //    TH1D* MC_Z_m_Gen; // Zu einer Instanz von specialAna gehoert ein Histogramm
  
  ofstream eventDisplayFile;
  std::stringstream eventsAfterCuts;
  std::stringstream eventsAfterCutsEvents;
  
  
  ///void initEvent( const pxl::Event* event );
  void endEvent( const pxl::Event* event );

  EventSelector *m_eventSelector;
  Systematics *m_syst_shifter;
  
  pxl::EventView *m_RecEvtView;
  pxl::EventView *m_GenEvtView;
  pxl::EventView *m_TrigEvtView;
  pxl::EventView *m_FilterEvtView;

  bool runOnData;
  bool useSyst;
  std::string m_dataPeriod;
  string const m_JetAlgo, m_BJets_algo, m_METType, m_TauType;
  const std::vector< std::string > m_trigger_string;
  std::unordered_set< std::string > triggers;


  std::vector< pxl::Particle* > * EleList;
  std::vector< pxl::Particle* > * MuonList;
  std::vector< pxl::Particle* > * TauList;	

};

#endif
