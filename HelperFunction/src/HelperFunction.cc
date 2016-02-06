// -*- C++ -*-
//
// Package:     Subsystem/Package
// Class  :     HelperFunction
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Tongguang Cheng
//         Created:  Mon, 21 Dec 2015 12:47:33 GMT
//

#ifndef HelperFunction_CC
#define HelperFunction_CC
// system include files

// user include files
#include "KinZfitter/HelperFunction/interface/HelperFunction.h"

// fileinPath
#include "FWCore/ParameterSet/interface/FileInPath.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HelperFunction::HelperFunction()
{

        //declarations
        debug_ = 0;

        TString fmu_s = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/ebeOverallCorrections.Legacy2013.v0.root" ).fullPath());
        TString fel_s = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/ebeOverallCorrections.Legacy2013.v0.root" ).fullPath());

        fmu = boost::shared_ptr<TFile>( new TFile(fmu_s));
        fel = boost::shared_ptr<TFile>( new TFile(fel_s));
                
        muon_corr_data = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fmu->Get( "mu_reco53x" )->Clone() )) );
        muon_corr_mc = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fmu->Get( "mu_mc53x" )->Clone() )) );
                
        electron_corr_data = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fel->Get( "el_reco53x" )->Clone() )) );
        electron_corr_mc = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fel->Get( "el_mc53x" )->Clone() )) );


}

// HelperFunction::HelperFunction(const HelperFunction& rhs)
// {
//    // do actual copying here;
// }

HelperFunction::~HelperFunction()
{
}

//
// assignment operators
//
// const HelperFunction& HelperFunction::operator=(const HelperFunction& rhs)
// {
//   //An exception safe implementation is
//   HelperFunction temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

double HelperFunction:: masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix){

        int ndim = 3*p4s.size();
        if(debug_) cout<<""<<endl;

        TMatrixD jacobian(1,ndim);

        double e = 0; double mass = 0;
        double px = 0; double py = 0; double pz = 0;
        for (unsigned int ip = 0; ip < p4s.size(); ip++) {
         
            e = e + p4s[ip].E();
            px = px + p4s[ip].Px();
            py = py + p4s[ip].Py();
            pz = pz + p4s[ip].Pz();
        }

        mass = TMath::Sqrt(e*e-px*px-py*py-pz*pz);

        for (unsigned int i = 0, o = 0; i < p4s.size(); i++, o += 3) {

                double pxi = p4s[i].Px();
                double pyi = p4s[i].Py();
                double pzi = p4s[i].Pz();
                double ei = p4s[i].E();

                jacobian(0, o+0) = (e*(pxi/ei) - px)/mass;
                jacobian(0, o+1) = (e*(pyi/ei) - py)/mass;
                jacobian(0, o+2) = (e*(pzi/ei) - pz)/mass;
        }

        TMatrixDSym massCov = covMatrix.Similarity(jacobian);

        double dm2 = massCov(0,0);
        return (dm2 > 0 ? std::sqrt(dm2) : 0.0);

}


double HelperFunction::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){
        // if(Lep.size()!= pterr.size()!=4) {std::cout<<" Lepsize="<<Lep.size()<<", "<<pterr.size()<<std::endl;}
        TLorentzVector compositeParticle ;
        for(unsigned int i=0; i<Lep.size(); i++){
                compositeParticle+=Lep[i];
        }
        double mass  =  compositeParticle.M();

        double masserr = 0;

        for(unsigned int i=0; i<Lep.size(); i++){
                TLorentzVector variedLep; // = Lep[i];

                variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
                TLorentzVector compositeParticleVariation ;
                for(unsigned int j=0; j<Lep.size(); j++){
                        if(i!=j)compositeParticleVariation+=Lep[j];
                        else compositeParticleVariation+=variedLep;
                }

                masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
        }

        return sqrt(masserr);
}


double HelperFunction::pterr( reco::Candidate *c, bool isData){

  reco::GsfElectron *gsf; reco::Muon *mu;
  reco::PFCandidate *pf;

  double pterrLep = 0.0;

  if ((gsf = dynamic_cast<reco::GsfElectron *> (&(*c)) ) != 0)
  {
    pterrLep=pterr(gsf, isData);
  }
  else if ((mu = dynamic_cast<reco::Muon *> (&(*c)) ) != 0)
  {
    pterrLep=pterr(mu, isData);
  }
  else if ((pf = dynamic_cast<reco::PFCandidate *> (&(*c)) ) != 0)
  { 
    pterrLep=pterr(c, isData);
  }


  return pterrLep;

}

double HelperFunction::pterr( reco::Muon* mu, bool isData){

        double pterr = mu->muonBestTrack()->ptError();

        return pterr;
}

/*
double HelperFunction::pterr( pat::Muon muon, bool isData){

        TH2F* mu_corr;
        if(isData) mu_corr = dynamic_cast<TH2F*> (muon_corr_data->Clone());
        else mu_corr = dynamic_cast<TH2F*> (muon_corr_mc->Clone());
        
        TAxis* x_mupTaxis = mu_corr->GetXaxis(); TAxis* y_muetaaxis = mu_corr->GetYaxis();
        double maxPt = x_mupTaxis->GetXmax(); double minPt = x_mupTaxis->GetXmin();
                
        const pat::Muon *mu = &muon;
        double scaleFactormu = 1.0;
        int xbin = x_mupTaxis->FindBin(mu->pt()); int ybin = y_muetaaxis->FindBin(fabs(mu->eta()));
        if(mu->pt()>minPt && mu->pt()<maxPt){  scaleFactormu = mu_corr->GetBinContent(xbin,ybin);  }
        double pterr = mu->muonBestTrack()->ptError() * scaleFactormu;

        return pterr;
}
*/

double HelperFunction::pterr( reco::GsfElectron * elec, bool isData ){

        if(debug_) cout<<"reco:gsfelectron pt err"<<endl; 

        double perr = elec->p();

        if (elec->ecalDriven()){
           perr = elec->p4Error(reco::GsfElectron::P4_COMBINATION);         
        }
        else{

                 double ecalEnergy = elec->correctedEcalEnergy() ;

                 if(debug_)cout<<"ecalEnergy "<<ecalEnergy<<endl;
                 double err2 = 0.0;
                 if (elec->isEB()) {
                        err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                        err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                        err2 += 1.00e-02*1.00e-02;
                 } else if (elec->isEE()) {
                        err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                        err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                        err2 += 1.94e-03*1.94e-03;
                 }
                 perr = ecalEnergy * sqrt(err2);

        }
         
        double pterr = perr*elec->pt()/elec->p();

        return pterr;

}

/*
double HelperFunction::pterr( pat::Electron electron, bool isData ){

        TH2F* el_corr;
        if(isData) el_corr = dynamic_cast<TH2F*>(electron_corr_data->Clone());
        else el_corr = dynamic_cast<TH2F*>(electron_corr_mc->Clone());
        TAxis* x_elpTaxis = el_corr->GetXaxis(); TAxis* y_eletaaxis = el_corr->GetYaxis();
        double maxPt = x_elpTaxis->GetXmax(); double minPt = x_elpTaxis->GetXmin();

        const pat::Electron *elec = &(electron);
        double perr = 0.;
        if (elec->ecalDriven()) {

              perr = elec->p4Error(reco::GsfElectron::P4_COMBINATION);

        }
        else {
                 // Parametrization from Claude Charlot, 
                 // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CJLST/ZZAnalysis/AnalysisStep/src/ZZMassErrors.cc?revision=1.2&view=markup
#if CMSSW_VERSION<500
                        double ecalEnergy = elec->ecalEnergy() ;
#else
                        double ecalEnergy = elec->correctedEcalEnergy() ;
#endif
                        double err2 = 0.0;
                        if (elec->isEB()) {
                                err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                                err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                                err2 += 1.00e-02*1.00e-02;
                        } else if (elec->isEE()) {
                                err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                                err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                                err2 += 1.94e-03*1.94e-03;
                        }
                        perr = ecalEnergy * sqrt(err2);
         }

         double scaleFactor_el = 1.0;

         int xbin = x_elpTaxis->FindBin(elec->pt()); int ybin = y_eletaaxis->FindBin(fabs(elec->eta()));
         if(elec->pt()>minPt && elec->pt()<maxPt ){  scaleFactor_el = el_corr->GetBinContent(xbin,ybin);  }

         double pterr = scaleFactor_el*(perr*elec->pt()/elec->p());

         return pterr;

}
*/

double HelperFunction::pterr(TLorentzVector ph){

         if(debug_) cout<<"perr for pf photon"<<endl;

         double perr = PFEnergyResolution().getEnergyResolutionEm(ph.E(), ph.Eta());

         double pterr = perr*ph.Pt()/ph.P();

         return pterr;
}

//
// const member functions
//

//
// static member functions
//
#endif
