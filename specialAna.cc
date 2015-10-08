#include "specialAna.hh"
#include "HistClass.hh"
#include "Tools/Tools.hh"
#include <csignal>
#include "TVector3.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "boost/format.hpp"
#pragma GCC diagnostic pop

specialAna::specialAna( const Tools::MConfig &cfg,
						Systematics *syst_shifter,
						EventSelector *selector
						) :
   runOnData(       cfg.GetItem< bool >( "General.RunOnData" ) ),
   //   m_TauType(       cfg.GetItem< string >( "Tau.Type.Rec" ) ),
   m_trigger_string(Tools::splitString< std::string >(cfg.GetItem< std::string >("ZH.trigger_list"))),
   m_eventSelector( selector )
{
    string safeFileName = "ZHTauTau.root";
    file1 = new TFile(safeFileName.c_str(), "RECREATE");
    file1->cd();
	
    //histo = new TH1D("zmass","zmass",100,0,100);
    HistClass::CreateHisto("nTau",   8, -0.5, 7.5, "nTau  (Reco)");
    HistClass::CreateHisto("nMuon",  8, -0.5, 7.5, "nMuon (Reco)");
    HistClass::CreateHisto("nEle",   8, -0.5, 7.5, "nEle  (Reco)");
    HistClass::CreateHisto("nLep",   8, -0.5, 7.5, "nLep  (Reco)");

    HistClass::CreateHisto("Zmass_mu",  400, 0, 200, "M_{Z} [GeV]");
    HistClass::CreateHisto("Zmass_ele", 400, 0, 200, "M_{Z} [GeV]");
    HistClass::CreateHisto("Zmass_lep", 400, 0, 200, "M_{Z} [GeV]");

    HistClass::CreateHisto("Hmass_tau",  100, 0, 300,  "M_{H} [GeV]");
    HistClass::CreateHisto("ZPrimemass", 100, 0, 4000, "M_{X} [GeV]");

    // number of events, saved in a histogram
    HistClass::CreateHistoUnchangedName("h_counters", 10, 0, 11, "N_{events}");
    
    for(auto syst : syst_shifter->m_activeSystematics){
        if(syst->m_particleType=="met"){
            syst->m_particleType=m_METType;
        }
    }
}

specialAna::~specialAna() {
}

void specialAna::analyseEvent( const pxl::Event* event ) {
  //~ initEvent( event );
  std::cout << "\nEVT" << std::endl;

  pxl::EventView *RecEvtView = event->getObjectOwner().findObject< pxl::EventView >( "Rec" );
  bool isRec = (RecEvtView->getUserRecord("Type").asString() == "Rec");
  
  //m_eventSelector;
  std::map< std::string, std::vector< pxl::Particle* > > particleMap = m_eventSelector->getParticleLists( RecEvtView, isRec );
  std::vector< pxl::Particle* > muoList = particleMap["Muon"];
  std::vector< pxl::Particle* > eleList = particleMap["Ele"];
  std::vector< pxl::Particle* > tauList = particleMap["Tau"];
  
  //std::cout << "muoList.size() " << muoList.size() << std::endl;
  HistClass::Fill("nTau",  tauList.size(), 1.);
  HistClass::Fill("nMuon", muoList.size(), 1.);
  HistClass::Fill("nEle",  eleList.size(), 1.);
  HistClass::Fill("nLep",  (eleList.size()+muoList.size()), 1.);

  //   std::vector< pxl::Particle* > zCandidateList;
  pxl::Particle*  zCandidate_Mu = 0;
  pxl::Particle*  zCandidate_Ele = 0;
  pxl::Particle*  zCandidate = 0;
  pxl::Particle*  HiggsCandidate_Tau = 0;
  pxl::Particle*  ZPrime_Candidate = 0;


  double zmass = 91.2;
  double bestdiff = 10000000.0;
  double Hmass = 125.0;
 
  if (TriggerSelector(event)) {
    
    /////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if (muoList.size() > 1){
      //~ for( auto& muo : muoList ){
      for( vector< pxl::Particle* >::const_iterator muon1 = muoList.begin(); muon1 != muoList.end(); muon1++ ){
	pxl::Particle* mu1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *muon1 );
	for( vector< pxl::Particle* >::const_iterator muon2 = muon1+1; muon2 != muoList.end(); ++muon2 ){
	  pxl::Particle* mu2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *muon2 );
	  if ( ((Check_Muo_ID_Strict(mu1)) &&  (Check_Muo_ID_Mild(mu2))) || ((Check_Muo_ID_Strict(mu2)) &&  (Check_Muo_ID_Mild(mu1)))  ) {
	    //std::cout << "mu1->getCharge() = " << mu1->getCharge() << std::endl;
	    if ( (mu1->getCharge()+mu2->getCharge())==0 ) {
	      pxl::Particle* zcand_mu = RecEvtView->getObjectOwner().create< pxl::Particle >( *muon1 );
	      zcand_mu->addP4( (*muon2)->getVector() );
	      // zCandidate_MuList.push_back( zcand );
	      if( fabs( zcand_mu->getMass() - zmass) < bestdiff) {
		zCandidate_Mu = zcand_mu;
		bestdiff = fabs( zcand_mu->getMass() - zmass);
	      }
	    }
	  }		
	}
      }
      if (zCandidate_Mu != 0) {
	//std::cout << "**Found z mass Muon** " << zCandidate_Mu->getMass() << std::endl;
	HistClass::Fill("Zmass_mu", zCandidate_Mu->getMass(), 1.);
      }  
    }
    
    // Electron //
    if (eleList.size() > 1){
      for( vector< pxl::Particle* >::const_iterator ele1 = eleList.begin(); ele1 != eleList.end(); ele1++ ){
	pxl::Particle* e1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *ele1 );
	if ( ((e1->getUserRecord("ISOfailed").asBool()) == true) ||  ((e1->getUserRecord("IDpassed").asBool()) == true)   ) {
	  for( vector< pxl::Particle* >::const_iterator ele2 = ele1+1; ele2 != eleList.end(); ++ele2 ){
	    pxl::Particle* e2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *ele2 );
	    if ( ((e2->getUserRecord("ISOfailed").asBool()) == true) || ((e2->getUserRecord("IDpassed").asBool()) == true) ) {
	      if ( (e1->getCharge()+e2->getCharge())==0 ){
		pxl::Particle* zcand_ele = RecEvtView->getObjectOwner().create< pxl::Particle >( *ele1 );
		zcand_ele->addP4( (*ele2)->getVector() );
		if( fabs( zcand_ele->getMass() - zmass) < bestdiff) {
		  zCandidate_Ele = zcand_ele;
		  bestdiff = fabs( zcand_ele->getMass() - zmass);
		}
	      }
	    }		
	  }
	}
      }
      if (zCandidate_Ele != 0) {
	// std::cout << "**Found z mass ele** " << zCandidate_Ele->getMass() << std::endl;
	HistClass::Fill("Zmass_ele", zCandidate_Ele->getMass(), 1.);
      }
    }
    
    //  std::cout << "zCandidate_Mu=" << zCandidate_Mu << "   zCandidate_Ele=" << zCandidate_Ele << std::endl; 
    if (zCandidate_Mu && zCandidate_Ele) {
      //  std::cout << "z cand found from both mu and ele collection. Have to choose one." << std::endl;
      if  ( (fabs(zmass-zCandidate_Mu->getMass())) < (fabs(zmass-zCandidate_Ele->getMass())) )  zCandidate=zCandidate_Mu;
      else zCandidate=zCandidate_Ele;
    }
    
    else if (zCandidate_Mu && !zCandidate_Ele) {
      zCandidate=zCandidate_Mu;
    }
    else if (!zCandidate_Mu && zCandidate_Ele) {
      zCandidate=zCandidate_Ele;
    }
    
    if (zCandidate != 0){
      HistClass::Fill("Zmass_lep", zCandidate->getMass(), 1.);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    if (tauList.size() > 1){
      for( vector< pxl::Particle* >::const_iterator tau1 = tauList.begin(); tau1 != tauList.end(); tau1++ ){
	pxl::Particle* t1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau1 );
	if ( Check_Tau_ID(t1) ) {
	  for( vector< pxl::Particle* >::const_iterator tau2 = tau1+1; tau2 != tauList.end(); ++tau2 ){
	    pxl::Particle* t2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau2 );
	    if ( Check_Tau_ID(t2) ) {
	      pxl::Particle* Hcand_tau = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau1 );
	      Hcand_tau->addP4( (*tau2)->getVector() );
	      if( fabs( Hcand_tau->getMass() - Hmass) < bestdiff) {
		HiggsCandidate_Tau = Hcand_tau;
		bestdiff = fabs( Hcand_tau->getMass() - Hmass);
	      }
	    }
	  }		
	}
      }
      if (HiggsCandidate_Tau != 0) {
	// std::cout << "**Found H mass Tau** " << HiggsCandidate_Tau->getMass() << std::endl;
	HistClass::Fill("Hmass_tau", HiggsCandidate_Tau->getMass(), 1.);
      }  
    }
    
    if ( zCandidate  &&  HiggsCandidate_Tau ) {
      std::cout << "SIGNAL EVENT FOUND !!!!!" << std::endl ;
      ZPrime_Candidate = zCandidate;
      ZPrime_Candidate->addP4(HiggsCandidate_Tau->getVector() );
      std::cout << "Zprime mass=" << ZPrime_Candidate->getMass() << std::endl;
      HistClass::Fill("ZPrimemass", ZPrime_Candidate->getMass(), 1.);
      
    }
  }
}


void specialAna::endJob( const Serializable* ) { 
  //std::cout << "Inside endJob" << std::endl;		
  
  file1->mkdir("RECO");
  file1->cd("RECO/");
  HistClass::Write("nTau");
  HistClass::Write("nMuon");
  HistClass::Write("nEle");
  HistClass::Write("nLep");

  HistClass::Write("Zmass_mu");
  HistClass::Write("Zmass_ele");
  HistClass::Write("Zmass_lep");

  HistClass::Write("Hmass_tau");
  HistClass::Write("ZPrimemass");

  file1->Close();
  
  delete file1;
  
}


bool specialAna::Check_Muo_ID_Strict(pxl::Particle* muon, bool do_pt_cut, bool do_eta_cut) {
  bool muon_ID = muon->getUserRecord("isTightMuon").asBool() ? true : false;
  //bool muon_ISO = muon -> getUserRecord("IsoR3SumPt").asFloat() / muon -> getPt() < 0.1 ? true : false;
  bool muon_eta = TMath::Abs(muon -> getEta()) < 2.1 ? true : false;
  bool muon_pt = muon -> getPt() > 25. ? true : false;
  if (not do_pt_cut) {
    muon_pt = true;
  }
  if (not do_eta_cut) {
    muon_eta = true;
  }
  if (muon_ID && muon_eta && muon_pt) return true;
  return false;
}

bool specialAna::Check_Muo_ID_Mild(pxl::Particle* muon, bool do_pt_cut, bool do_eta_cut) {
  bool muon_ID = false;
  if ( (muon->getUserRecord("isTrackerMuon").asBool()) || (muon->getUserRecord("isGlobalMuon").asBool()) ) muon_ID=true;
  //bool muon_ISO = muon -> getUserRecord("IsoR3SumPt").asFloat() / muon -> getPt() < 0.1 ? true : false;
  bool muon_eta = TMath::Abs(muon -> getEta()) < 2.1 ? true : false;
  bool muon_pt = muon -> getPt() > 25. ? true : false;
  if (not do_pt_cut) {
    muon_pt = true;
  }
  if (not do_eta_cut) {
    muon_eta = true;
  }
  if (muon_ID && muon_eta && muon_pt) return true;
  return false;
}

bool specialAna::Check_Tau_ID(pxl::Particle* tau, bool do_pt_cut, bool do_eta_cut) {
  bool tau_ID = tau->getUserRecord("decayModeFinding").asFloat() >= 1 ? true : false;
  ////bool tau_ISO = tau->getUserRecord("byLooseIsolationMVA3oldDMwLT").asFloat() >= 1 ? true : false;
  bool tau_ELE = tau->getUserRecord("againstElectronLooseMVA5").asFloat() >= 1 ? true : false;
  bool tau_MUO = tau->getUserRecord("againstMuonTight3").asFloat() >= 1 ? true : false;
  bool tau_eta = TMath::Abs(tau -> getEta()) < 3.0 ? true : false;
  bool tau_pt = tau -> getPt() > 25. ? true : false;
  
  if (tau_ID && tau_ELE && tau_MUO && tau_eta && tau_pt) return true;
  return false;
}


bool specialAna::TriggerSelector(const pxl::Event* event) {
  bool triggered = false;
  for (auto const it : m_trigger_string) {
    // triggered = m_TrigEvtView->hasUserRecord(it) ? true : false;
    for (auto us : m_TrigEvtView->getUserRecords()) {
      if (std::string::npos != us.first.find(it)) {
	triggered = true;
	triggers.insert(us.first);
      }
    }
  }
  return (triggered);
}
