#ifndef HelperFunction_H
#define HelperFunction_H
// -*- C++ -*-
//
// Package:     Subsystem/Package
// Class  :     HelperFunction
// 
/**\class HelperFunction HelperFunction.h "HelperFunction.h"

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author:  Tongguang Cheng
//         Created:  Mon, 21 Dec 2015 12:47:33 GMT
//

// system include files

// user include files

// forward declarations

#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <iomanip>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TFile.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CommonTools/CandUtils/interface/CenterOfMassBooster.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include "CommonTools/CandUtils/interface/cloneDecayTree.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include <cmath>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"
#include <TMatrixD.h>

// fit result covariance matrix
#include <TMatrixDSym.h>

namespace reco { class Candidate; class Muon; class GsfElectron; class Track; class PFCandidate; }
namespace edm { class EventSetup; }

#include <vector>
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include <TMatrixDSym.h>
#include <boost/shared_ptr.hpp>

using namespace std;


class HelperFunction
{

   public:
      HelperFunction();
      virtual ~HelperFunction();

      void setdebug(int d){debug_= d;};

      //ForZ
      double pterr(reco::Candidate *c, bool isData);

      //double pterr(pat::Electron electron, bool isData);
      //double pterr(pat::Muon muon, bool isData);
      double pterr(TLorentzVector fsrPhoton);

      double pterr(reco::GsfElectron* electron, bool isData);
      double pterr(reco::Muon* muon, bool isData);

      double masserror(std::vector<TLorentzVector> p4s, std::vector<double> pTErrs);

      double masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix);

      //double masserror(std::vector<TLorentzVector> p4s, )


      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------

   private:

      HelperFunction(const HelperFunction&); // stop default

      const HelperFunction& operator=(const HelperFunction&); // stop default

      int debug_;

      boost::shared_ptr<TFile>     fmu;
      boost::shared_ptr<TFile>     fel;
      boost::shared_ptr<TH2F>      muon_corr_data;
      boost::shared_ptr<TH2F>      muon_corr_mc;
      boost::shared_ptr<TH2F>      electron_corr_data;
      boost::shared_ptr<TH2F>      electron_corr_mc;

      // ---------- member data --------------------------------

};


#endif
