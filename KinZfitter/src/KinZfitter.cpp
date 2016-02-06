/*************************************************************************
 *  Authors:   Tongguang Cheng
 *************************************************************************/
#ifndef KinZfitter_cpp
#define KinZfitter_cpp

/// KinFitter header
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
#include "KinZfitter/HelperFunction/interface/HelperFunction.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

///----------------------------------------------------------------------------------------------
/// KinZfitter::KinZfitter - constructor/
///----------------------------------------------------------------------------------------------

KinZfitter::KinZfitter(bool isData)
{    

     PDFName_ = "GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8";

     debug_ = false;

     if(debug_) std::cout << "KinZfitter. The debug flag is ON with "<<PDFName_<< std::endl;
	
     /// Initialise HelperFunction
     helperFunc_ = new HelperFunction();
     isCorrPTerr_ = true; 
     isData_ = isData; 

}


void KinZfitter::Setup(std::vector< reco::Candidate* > selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons){

     // reset everything for each event
     idsZ1_.clear(); idsZ2_.clear();      
     idsFsrZ1_.clear(); idsFsrZ2_.clear();

     p4sZ1_.clear(); p4sZ2_.clear(); p4sZ1ph_.clear(); p4sZ2ph_.clear();
     p4sZ1REFIT_.clear(); p4sZ2REFIT_.clear(); p4sZ1phREFIT_.clear(); p4sZ2phREFIT_.clear();
     
     pTerrsZ1_.clear(); pTerrsZ2_.clear(); pTerrsZ1ph_.clear(); pTerrsZ2ph_.clear();
     pTerrsZ1REFIT_.clear(); pTerrsZ2REFIT_.clear(); pTerrsZ1phREFIT_.clear(); pTerrsZ2phREFIT_.clear();

     initZs(selectedLeptons, selectedFsrPhotons);

     if(debug_) cout<<"list ids"<<endl;
     if(debug_) cout<<"IDs[0] "<<idsZ1_[0]<<" IDs[1] "<<idsZ1_[1]<<" IDs[2] "<<idsZ2_[0]<<" IDs[3] "<<idsZ2_[1]<<endl;

     fs_=""; 
     if(abs(idsZ1_[0])==11 && abs(idsZ2_[0])==11) fs_="4e";
     if(abs(idsZ1_[0])==13 && abs(idsZ2_[0])==13) fs_="4mu";   
     if(abs(idsZ1_[0])==11 && abs(idsZ2_[0])==13) fs_="2e2mu";
     if(abs(idsZ1_[0])==13 && abs(idsZ2_[0])==11) fs_="2mu2e";

     if(debug_) cout<<"fs is "<<fs_<<endl;

     /////////////
     edm::FileInPath pdfFileWithFullPath("KinZfitter/KinZfitter/ParamZ1/dummy.txt");     
     string paramZ1_dummy = pdfFileWithFullPath.fullPath(); 
    
     TString paramZ1 = TString( paramZ1_dummy.substr(0,paramZ1_dummy.length() - 9));

     paramZ1+=PDFName_;
     paramZ1+="_";
     paramZ1+=+fs_;
     paramZ1+=".txt";

     if(debug_) cout<<"paramZ1 in "<<paramZ1<<endl;

     std::ifstream input(paramZ1);
     std::string line;
     while (!input.eof() && std::getline(input,line))
      {
         std::istringstream iss(line);
         string p; double val;
         if(iss >> p >> val) {
          if(p=="sg")  { sgVal_ = val; }//cout<<"sg is "<<sgVal_<<endl;}
          if(p=="a" )  { aVal_ = val;  }//cout<<"a is  "<<aVal_<<endl;}
          if(p=="n" )  { nVal_ = val;  }//cout<<"n is  "<<nVal_<<endl;}
          if(p=="f")   { fVal_ = val;  }//cout<<"f is  "<<fVal_<<endl;}

          if(p=="mean" )  { meanVal_ = val;  }//cout<<"mean is  "<<meanVal_<<endl;}
          if(p=="sigma" )  { sigmaVal_ = val;  }//cout<<"sigma is  "<<sigmaVal_<<endl;}
          if(p=="f1")   { f1Val_ = val;  }//cout<<"f1 is  "<<f1Val_<<endl;}
         }
      }

}



///----------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------

void KinZfitter::initZs(std::vector< reco::Candidate* > selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons){

        if(debug_) cout<<"init leptons"<<endl;

        for(unsigned int il = 0; il<selectedLeptons.size(); il++)
         {
            double pTerr = 0; TLorentzVector p4;

            reco::Candidate * c = selectedLeptons[il];              
            pTerr = helperFunc_->pterr(c ,  isData_);
            p4.SetPxPyPzE(c->px(),c->py(),c->pz(),c->energy());  
            int pdgId = c->pdgId();

            if(debug_) cout<<"pdg id "<<pdgId<<endl;

            if(il<2){
              idsZ1_.push_back(pdgId);
              pTerrsZ1_.push_back(pTerr);
              p4sZ1_.push_back(p4);            

            }
            else{

              idsZ2_.push_back(pdgId);
              pTerrsZ2_.push_back(pTerr);
              p4sZ2_.push_back(p4);
            }

         }

        if(debug_) cout<<"init fsr photons"<<endl;

        for(unsigned int ifsr = 0; ifsr<4; ifsr++)
         {

            TLorentzVector p4 = selectedFsrPhotons[ifsr];
            if(selectedFsrPhotons[ifsr].Pt()==0) continue;

            if(debug_) cout<<"ifsr "<<ifsr<<endl;

            double pTerr = 0;

            pTerr = helperFunc_->pterr(p4); //,isData_);

            if(debug_) cout<<" pt err is "<<pTerr<<endl;

            if(ifsr<2){

                if(debug_) cout<<"for fsr Z1 photon"<<endl;

                pTerrsZ1ph_.push_back(pTerr);
                p4sZ1ph_.push_back(p4);
                idsFsrZ1_.push_back(idsZ1_[ifsr]);

              }
            else{

                if(debug_) cout<<"for fsr Z2 photon"<<endl;

                pTerrsZ2ph_.push_back(pTerr);
                p4sZ2ph_.push_back(p4);
                idsFsrZ2_.push_back(idsZ2_[ifsr-2]);

            }

         }

         if(debug_) cout<<"p4sZ1ph_ "<<p4sZ1ph_.size()<<" p4sZ2ph_ "<<p4sZ2ph_.size()<<endl;
  
}

void KinZfitter::SetZ1Result(double l1, double l2, double lph1, double lph2)
{

  if(debug_) cout<<"start set Z1 result"<<endl;

  // pT scale after refitting w.r.t. reco pT
  lZ1_l1_ = l1; lZ1_l2_ = l2;
  if(debug_) cout<<"l1 "<<l1<<" l2 "<<l2<<endl;
  lZ1_ph1_ = lph1; lZ1_ph2_ = lph2;

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];
  TLorentzVector Z2_1 = p4sZ2_[0]; TLorentzVector Z2_2 = p4sZ2_[1];

  TLorentzVector Z1_1_True(0,0,0,0);
  Z1_1_True.SetPtEtaPhiM(lZ1_l1_*Z1_1.Pt(),Z1_1.Eta(),Z1_1.Phi(),Z1_1.M());
  TLorentzVector Z1_2_True(0,0,0,0);
  Z1_2_True.SetPtEtaPhiM(lZ1_l2_*Z1_2.Pt(),Z1_2.Eta(),Z1_2.Phi(),Z1_2.M());

  p4sZ1REFIT_.push_back(Z1_1_True); p4sZ1REFIT_.push_back(Z1_2_True);

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){

      TLorentzVector Z1ph = p4sZ1ph_[ifsr1];
      TLorentzVector Z1phTrue(0,0,0,0);
  
      double l = 1.0;
      if(ifsr1==0) l = lZ1_ph1_; if(ifsr1==1) l = lZ1_ph2_;
  
      Z1phTrue.SetPtEtaPhiM(l*Z1ph.Pt(),Z1ph.Eta(),Z1ph.Phi(),Z1ph.M());
  
      p4sZ1phREFIT_.push_back(Z1phTrue);
  
  }

  // since it is Z1 refit result, Z2 kinematics keep at what it is 
  p4sZ2REFIT_.push_back(p4sZ2_[0]); p4sZ2REFIT_.push_back(p4sZ2_[1]);
  pTerrsZ2REFIT_.push_back(pTerrsZ2_[0]); pTerrsZ2REFIT_.push_back(pTerrsZ2_[1]);

  for(unsigned int iz2 = 0; iz2 < p4sZ2ph_.size(); iz2++){

     p4sZ2phREFIT_.push_back(p4sZ2ph_[iz2]);
     pTerrsZ2phREFIT_.push_back(pTerrsZ2ph_[iz2]);

  }

  if(debug_) cout<<"end set Z1 result"<<endl;

}

double KinZfitter::GetRefitM4l()
{

  vector<TLorentzVector> p4s = GetRefitP4s();

  TLorentzVector pH(0,0,0,0); 
  for(unsigned int i = 0; i< p4s.size(); i++){
     pH = pH + p4s[i];
  }

  return pH.M();

}

double KinZfitter::GetRefitMZ1()
{

  vector<TLorentzVector> p4s = GetRefitP4s();

  TLorentzVector pZ1(0,0,0,0);

  pZ1 = p4s[0] + p4s[1];

  return pZ1.M();

}

double KinZfitter::GetRefitM4lErr()
{

  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1REFIT_[0]);p4s.push_back(p4sZ1REFIT_[1]);
  p4s.push_back(p4sZ2REFIT_[0]);p4s.push_back(p4sZ2REFIT_[1]);

  // patch when MINIUT FAILS
  if(pTerrsZ1REFIT_[0]==0||pTerrsZ1REFIT_[1]==0)
   return GetM4lErr();

  pTErrs.push_back(pTerrsZ1REFIT_[0]); pTErrs.push_back(pTerrsZ1REFIT_[1]);
  pTErrs.push_back(pTerrsZ2REFIT_[0]); pTErrs.push_back(pTerrsZ2REFIT_[1]);

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1phREFIT_.size(); ifsr1++){

      p4s.push_back(p4sZ1phREFIT_[ifsr1]);
      pTErrs.push_back(pTerrsZ1phREFIT_[ifsr1]);

  }

  for(unsigned int ifsr2 = 0; ifsr2<p4sZ2phREFIT_.size(); ifsr2++){

      p4s.push_back(p4sZ2phREFIT_[ifsr2]);
      pTErrs.push_back(pTerrsZ2phREFIT_[ifsr2]);

  }

  return helperFunc_->masserror(p4s,pTErrs);

}

double KinZfitter::GetRefitM4lErrFullCov()
{


  vector<TLorentzVector> Lp4s = GetRefitP4s();
  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1REFIT_[0]);p4s.push_back(p4sZ1REFIT_[1]);
  pTErrs.push_back(pTerrsZ1REFIT_[0]); pTErrs.push_back(pTerrsZ1REFIT_[1]);

  // patch when MINUIT FAILS
  if(pTerrsZ1REFIT_[0]==0||pTerrsZ1REFIT_[1]==0)

   return GetM4lErr();

  if(p4sZ1phREFIT_.size()>=1){
   p4s.push_back(p4sZ1phREFIT_[0]); pTErrs.push_back(pTerrsZ1phREFIT_[0]);
  } 
  if(p4sZ1phREFIT_.size()==2){
   p4s.push_back(p4sZ1phREFIT_[1]); pTErrs.push_back(pTerrsZ1phREFIT_[1]);
  }

  p4s.push_back(p4sZ2REFIT_[0]);p4s.push_back(p4sZ2REFIT_[1]);
  pTErrs.push_back(pTerrsZ2REFIT_[0]); pTErrs.push_back(pTerrsZ2REFIT_[1]);

  if(p4sZ2phREFIT_.size()>=1){
   p4s.push_back(p4sZ2phREFIT_[0]); pTErrs.push_back(pTerrsZ2phREFIT_[0]);
  }
  if(p4sZ2phREFIT_.size()==2){
   p4s.push_back(p4sZ2phREFIT_[1]); pTErrs.push_back(pTerrsZ2phREFIT_[1]);
  }

  double errorUncorr = helperFunc_->masserror(p4s,pTErrs);

  vector<double> pTErrs1; vector<double> pTErrs2;
  for(unsigned int i = 0; i<pTErrs.size(); i++){
     if(i==0) pTErrs1.push_back(pTErrs[i]);
     else pTErrs1.push_back(0.0);
  }
  for(unsigned int i = 0; i<pTErrs.size(); i++){
     if(i==1) pTErrs2.push_back(pTErrs[i]);
     else pTErrs2.push_back(0.0);
  }

  double error1 = helperFunc_->masserror(p4s,pTErrs1);
  double error2 = helperFunc_->masserror(p4s,pTErrs2);

  double errorph1 = 0.0; double errorph2 = 0.0; 
  if(p4sZ2phREFIT_.size()>=1){ 
 
   vector<double> pTErrsph1;
   for(unsigned int i = 0; i<pTErrs.size(); i++){
     if(i==2) pTErrsph1.push_back(pTErrs[i]);
     else pTErrsph1.push_back(0.0);
   }

   errorph1 = helperFunc_->masserror(p4s,pTErrsph1);

  }

  if(p4sZ2phREFIT_.size()>=2){

   vector<double> pTErrsph2;
   for(unsigned int i = 0; i<pTErrs.size(); i++){
     if(i==3) pTErrsph2.push_back(pTErrs[i]);
     else pTErrsph2.push_back(0.0);
   }
  
   errorph2 = helperFunc_->masserror(p4s,pTErrsph2);  

  }

  if(debug_) cout<<"error1 "<<error1<<" error2 "<<error2<<endl;

  ////
  // covariance matrix

  double delta12 = error1*error2*covMatrixZ1_(0,1)/sqrt(covMatrixZ1_(0,0)*covMatrixZ1_(1,1));
  double delta1ph1 = 0.0; double delta1ph2 = 0.0;
  double delta2ph1 = 0.0; double delta2ph2 = 0.0;
  double deltaph1ph2 = 0.0;

  if(p4sZ1phREFIT_.size()>=1){
     delta1ph1 = error1*errorph1*covMatrixZ1_(0,2)/sqrt(covMatrixZ1_(0,0)*covMatrixZ1_(2,2));
     delta2ph1 = error2*errorph1*covMatrixZ1_(1,2)/sqrt(covMatrixZ1_(1,1)*covMatrixZ1_(2,2));
  }

  if(p4sZ1phREFIT_.size()>=2){
     delta1ph2 = error1*errorph2*covMatrixZ1_(0,3)/sqrt(covMatrixZ1_(0,0)*covMatrixZ1_(3,3));
     delta2ph2 = error2*errorph2*covMatrixZ1_(1,3)/sqrt(covMatrixZ1_(1,1)*covMatrixZ1_(3,3));
     delta1ph2 = errorph1*errorph2*covMatrixZ1_(2,3)/sqrt(covMatrixZ1_(2,2)*covMatrixZ1_(3,3));
  }


  double correlation = delta12+delta1ph1+delta1ph2+delta2ph1+delta2ph2+deltaph1ph2;

  double err = sqrt(errorUncorr*errorUncorr+correlation);

  return err;
  
}

double KinZfitter::GetM4lErr()
{
  
  vector<TLorentzVector> p4s;
  vector<double> pTErrs;
  
  p4s.push_back(p4sZ1_[0]);p4s.push_back(p4sZ1_[1]);
  p4s.push_back(p4sZ2_[0]);p4s.push_back(p4sZ2_[1]);
  
  pTErrs.push_back(pTerrsZ1_[0]); pTErrs.push_back(pTerrsZ1_[1]);
  pTErrs.push_back(pTerrsZ2_[0]); pTErrs.push_back(pTerrsZ2_[1]);
  
  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){
      
      p4s.push_back(p4sZ1ph_[ifsr1]);
      pTErrs.push_back(pTerrsZ1ph_[ifsr1]);
  
  }
  
  for(unsigned int ifsr2 = 0; ifsr2<p4sZ2ph_.size(); ifsr2++){
      
      p4s.push_back(p4sZ2ph_[ifsr2]);
      pTErrs.push_back(pTerrsZ2ph_[ifsr2]);
  
  }
  
  return helperFunc_->masserror(p4s,pTErrs);

}

double KinZfitter::GetMZ1Err()
{

  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1_[0]);p4s.push_back(p4sZ1_[1]);
  pTErrs.push_back(pTerrsZ1_[0]); pTErrs.push_back(pTerrsZ1_[1]);

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){

      p4s.push_back(p4sZ1ph_[ifsr1]);
      pTErrs.push_back(pTerrsZ1ph_[ifsr1]);

  }

  return helperFunc_->masserror(p4s,pTErrs);

}


vector<TLorentzVector> KinZfitter::GetRefitP4s()
{

  TLorentzVector Z1_1 = p4sZ1REFIT_[0]; TLorentzVector Z1_2 = p4sZ1REFIT_[1];
  TLorentzVector Z2_1 = p4sZ2REFIT_[0]; TLorentzVector Z2_2 = p4sZ2REFIT_[1];

  /// fsr photons

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1phREFIT_.size(); ifsr1++){

      int id_fsr1 = idsFsrZ1_[ifsr1];
      TLorentzVector Z1ph = p4sZ1phREFIT_[ifsr1];

      if(id_fsr1==idsZ1_[0]) Z1_1 = Z1_1 + Z1ph;
      if(id_fsr1==idsZ1_[1]) Z1_2 = Z1_2 + Z1ph;

  }

  for(unsigned int ifsr2 = 0; ifsr2<p4sZ2phREFIT_.size(); ifsr2++){

      int id_fsr2 = idsFsrZ2_[ifsr2];
      TLorentzVector Z2ph = p4sZ2phREFIT_[ifsr2];

      if(id_fsr2==idsZ2_[0]) Z2_1 = Z2_1 + Z2ph;
      if(id_fsr2==idsZ2_[1]) Z2_2 = Z2_2 + Z2ph;

  }

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); p4s.push_back(Z1_2);
  p4s.push_back(Z2_1); p4s.push_back(Z2_2);

  return p4s;

}

vector<TLorentzVector> KinZfitter::GetP4s()
{

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];
  TLorentzVector Z2_1 = p4sZ2_[0]; TLorentzVector Z2_2 = p4sZ2_[1];

  /// fsr photons

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){

      int id_fsr1 = idsFsrZ1_[ifsr1];
      TLorentzVector Z1ph = p4sZ1ph_[ifsr1];

      if(id_fsr1==idsZ1_[0]) Z1_1 = Z1_1 + Z1ph;
      if(id_fsr1==idsZ1_[1]) Z1_2 = Z1_2 + Z1ph;

  }

  for(unsigned int ifsr2 = 0; ifsr2<p4sZ2ph_.size(); ifsr2++){

      int id_fsr2 = idsFsrZ2_[ifsr2];
      TLorentzVector Z2ph = p4sZ2ph_[ifsr2];

      if(id_fsr2==idsZ2_[0]) Z2_1 = Z2_1 + Z2ph;
      if(id_fsr2==idsZ2_[1]) Z2_2 = Z2_2 + Z2ph;

  }

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); 
  p4s.push_back(Z1_2);
  p4s.push_back(Z2_1); 
  p4s.push_back(Z2_2);

  return p4s;

}

void KinZfitter::KinRefitZ1()
{
  double l1,l2,lph1,lph2;
  l1 = 1.0; l2 = 1.0; lph1 = 1.0; lph2 = 1.0; 

  PerZ1Likelihood(l1,l2,lph1,lph2);
  if(debug_) cout<<"l1 "<<l1<<"; l2 "<<l2<<" lph1 "<<lph1<<" lph2 "<<lph2<<endl;
  SetZ1Result(l1,l2,lph1,lph2);
  if(debug_) cout<<"Z1 refit done"<<endl;
}


/*
void KinZfitter::KinRefitZZ(std::vector< pat::Muon > selectedMuons, std::vector< pat::Electron > selectedElectrons, std::vector< pat::PFParticle > selectedZ1FsrPhotons, std::vector< pat::PFParticle > selectedZ2FsrPhotons)
}
*/

int KinZfitter::PerZ1Likelihood(double & l1, double & l2, double & lph1, double & lph2)
{

    l1= 1.0; l2 = 1.0;
    lph1 = 1.0; lph2 = 1.0;

    if(debug_) cout<<"start Z1 refit"<<endl;

    TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];

    double RECOpT1 = Z1_1.Pt(); double RECOpT2 = Z1_2.Pt();
    double pTerrZ1_1 = pTerrsZ1_[0]; double pTerrZ1_2 = pTerrsZ1_[1];

    if(debug_)cout<<"pT1 "<<RECOpT1<<" pTerrZ1_1 "<<pTerrZ1_1<<endl;
    if(debug_)cout<<"pT2 "<<RECOpT2<<" pTerrZ1_2 "<<pTerrZ1_2<<endl;

    //////////////

    TLorentzVector Z1_ph1, Z1_ph2;
    double pTerrZ1_ph1, pTerrZ1_ph2;
    double RECOpTph1, RECOpTph2;

    TLorentzVector nullFourVector(0, 0, 0, 0);
    Z1_ph1=nullFourVector; Z1_ph2=nullFourVector;
    RECOpTph1 = 0; RECOpTph2 = 0;
    pTerrZ1_ph1 = 0; pTerrZ1_ph2 = 0;

    if(p4sZ1ph_.size()>=1){

      Z1_ph1 = p4sZ1ph_[0]; pTerrZ1_ph1 = pTerrsZ1ph_[0];
      RECOpTph1 = Z1_ph1.Pt();
      if(debug_) cout<<"put in Z1 fsr photon 1 pT "<<RECOpTph1<<" pT err "<<pTerrZ1_ph1<<endl; 
    }
    if(p4sZ1ph_.size()==2){
      //if(debug_) cout<<"put in Z1 fsr photon 2"<<endl;
      Z1_ph2 = p4sZ1ph_[1]; pTerrZ1_ph2 = pTerrsZ1ph_[1];
      RECOpTph2 = Z1_ph2.Pt();     
    }

    RooRealVar* pT1RECO = new RooRealVar("pT1RECO","pT1RECO", RECOpT1, 5, 500);
    RooRealVar* pT2RECO = new RooRealVar("pT2RECO","pT2RECO", RECOpT2, 5, 500);
   
    double RECOpT1min = max(5.0, RECOpT1-2*pTerrZ1_1);
    double RECOpT2min = max(5.0, RECOpT2-2*pTerrZ1_2);

    RooRealVar* pTph1RECO = new RooRealVar("pTph1RECO","pTph1RECO", RECOpTph1, 5, 500);
    RooRealVar* pTph2RECO = new RooRealVar("pTph2RECO","pTph2RECO", RECOpTph2, 5, 500);

    double RECOpTph1min = max(0.5, RECOpTph1-2*pTerrZ1_ph1);
    double RECOpTph2min = max(0.5, RECOpTph2-2*pTerrZ1_ph2);

    // observables pT1,2,ph1,ph2
    RooRealVar* pT1 = new RooRealVar("pT1", "pT1FIT", RECOpT1, RECOpT1min, RECOpT1+2*pTerrZ1_1 );
    RooRealVar* pT2 = new RooRealVar("pT2", "pT2FIT", RECOpT2, RECOpT2min, RECOpT2+2*pTerrZ1_2 );

    RooRealVar* m1 = new RooRealVar("m1","m1", Z1_1.M());
    RooRealVar* m2 = new RooRealVar("m2","m2", Z1_2.M());

    if(debug_) cout<<"m1 "<<m1->getVal()<<" m2 "<<m2->getVal()<<endl;

    double Vtheta1, Vphi1, Vtheta2, Vphi2;
    Vtheta1 = (Z1_1).Theta(); Vtheta2 = (Z1_2).Theta();
    Vphi1 = (Z1_1).Phi(); Vphi2 = (Z1_2).Phi();

    RooRealVar* theta1 = new RooRealVar("theta1","theta1",Vtheta1);
    RooRealVar* phi1   = new RooRealVar("phi1","phi1",Vphi1);
    RooRealVar* theta2 = new RooRealVar("theta2","theta2",Vtheta2);
    RooRealVar* phi2   = new RooRealVar("phi2","phi2",Vphi2);

    // dot product to calculate (p1+p2+ph1+ph2).M()
    RooFormulaVar E1("E1","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)",
                          RooArgList(*pT1,*theta1,*m1));
    RooFormulaVar E2("E2","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)",
                          RooArgList(*pT2,*theta2,*m2));
    if(debug_) cout<<"E1 "<<E1.getVal()<<"; E2 "<<E2.getVal()<<endl;

    /////

    RooRealVar* pTph1 = new RooRealVar("pTph1", "pTph1FIT", RECOpTph1, RECOpTph1min, RECOpTph1+2*pTerrZ1_ph1 );
    RooRealVar* pTph2 = new RooRealVar("pTph2", "pTph2FIT", RECOpTph2, RECOpTph2min, RECOpTph2+2*pTerrZ1_ph2 );

    double Vthetaph1, Vphiph1, Vthetaph2, Vphiph2;
    Vthetaph1 = (Z1_ph1).Theta(); Vthetaph2 = (Z1_ph2).Theta();
    Vphiph1 = (Z1_ph1).Phi(); Vphiph2 = (Z1_ph2).Phi();

    RooRealVar* thetaph1 = new RooRealVar("thetaph1","thetaph1",Vthetaph1);
    RooRealVar* phiph1   = new RooRealVar("phiph1","phiph1",Vphiph1);
    RooRealVar* thetaph2 = new RooRealVar("thetaph2","thetaph2",Vthetaph2);
    RooRealVar* phiph2   = new RooRealVar("phiph2","phi2",Vphiph2);

    RooFormulaVar Eph1("Eph1","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))", 
                              RooArgList(*pTph1,*thetaph1));
    RooFormulaVar Eph2("Eph2","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))", 
                              RooArgList(*pTph2,*thetaph2));

    //// dot products of 4-vectors

    // 3-vector DOT
    RooFormulaVar* p1v3D2 = new RooFormulaVar("p1v3D2",
         "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3)))+(TMath::Cos(@4-@5)))",
         RooArgList(*pT1,*pT2,*theta1,*theta2,*phi1,*phi2));    
    if(debug_) cout<<"p1 DOT p2 is "<<p1v3D2->getVal()<<endl;
    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar p1D2("p1D2","@0*@1-@2",RooArgList(E1,E2,*p1v3D2));

    //lep DOT fsrPhoton1

    // 3-vector DOT
    RooFormulaVar* p1v3Dph1 = new RooFormulaVar("p1v3Dph1",
         "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
         RooArgList(*pT1,*pTph1,*theta1,*thetaph1,*phi1,*phiph1));

    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar p1Dph1("p1Dph1","@0*@1-@2",RooArgList(E1,Eph1,*p1v3Dph1));

    // 3-vector DOT
    RooFormulaVar* p2v3Dph1 = new RooFormulaVar("p2v3Dph1",
         "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
         RooArgList(*pT2,*pTph1,*theta2,*thetaph1,*phi2,*phiph1));
    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar p2Dph1("p2Dph1","@0*@1-@2",RooArgList(E2,Eph1,*p2v3Dph1));

    // lep DOT fsrPhoton2 

    // 3-vector DOT
    RooFormulaVar* p1v3Dph2 = new RooFormulaVar("p1v3Dph2",
         "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
         RooArgList(*pT1,*pTph2,*theta1,*thetaph2,*phi1,*phiph2));

    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar p1Dph2("p1Dph2","@0*@1-@2",RooArgList(E1,Eph2,*p1v3Dph2));

    // 3-vector DOT
    RooFormulaVar* p2v3Dph2 = new RooFormulaVar("p2v3Dph2",
         "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
         RooArgList(*pT2,*pTph2,*theta2,*thetaph2,*phi2,*phiph2));
    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar p2Dph2("p2Dph2","@0*@1-@2",RooArgList(E2,Eph2,*p2v3Dph2));

    // fsrPhoton1 DOT fsrPhoton2

    // 3-vector DOT
    RooFormulaVar* ph1v3Dph2 = new RooFormulaVar("ph1v3Dph2",
         "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
         RooArgList(*pTph1,*pTph2,*thetaph1,*thetaph2,*phiph1,*phiph2));    
    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar ph1Dph2("ph1Dph2","@0*@1-@2",RooArgList(Eph1,Eph2,*ph1v3Dph2));

    // mZ1

    RooFormulaVar* mZ1;
    mZ1 = new RooFormulaVar("mZ1","TMath::Sqrt(2*@0+@1*@1+@2*@2)",RooArgList(p1D2,*m1,*m2));
    if(p4sZ1ph_.size()==1)
      mZ1 = new RooFormulaVar("mZ1","TMath::Sqrt(2*@0+2*@1+2*@2+@3*@3+@4*@4)",
                                    RooArgList(p1D2, p1Dph1, p2Dph1, *m1,*m2));
    if(p4sZ1ph_.size()==2)
      mZ1 = new RooFormulaVar("mZ1","TMath::Sqrt(2*@0+2*@1+2*@2+2*@3+2*@4+2*@5+@6*@6+@7*@7)",
                              RooArgList(p1D2,p1Dph1,p2Dph1,p1Dph2,p2Dph2,ph1Dph2, *m1,*m2));

    if(debug_) cout<<"mZ1 is "<<mZ1->getVal()<<endl;

    // pTerrs, 1,2,ph1,ph2
    RooRealVar sigmaZ1_1("sigmaZ1_1", "sigmaZ1_1", pTerrZ1_1);
    RooRealVar sigmaZ1_2("sigmaZ1_2", "sigmaZ1_2", pTerrZ1_2);

    RooRealVar sigmaZ1_ph1("sigmaZ1_ph1", "sigmaZ1_ph1", pTerrZ1_ph1);
    RooRealVar sigmaZ1_ph2("sigmaZ1_ph2", "sigmaZ1_ph2", pTerrZ1_ph2);

    // resolution for decay products
    RooGaussian gauss1("gauss1","gaussian PDF", *pT1RECO, *pT1, sigmaZ1_1);
    RooGaussian gauss2("gauss2","gaussian PDF", *pT2RECO, *pT2, sigmaZ1_2);

    RooGaussian gaussph1("gaussph1","gaussian PDF", *pTph1RECO, *pTph1, sigmaZ1_ph1);
    RooGaussian gaussph2("gaussph2","gaussian PDF", *pTph2RECO, *pTph2, sigmaZ1_ph2);

    RooRealVar bwMean("bwMean", "m_{Z^{0}}", 91.187);
    RooRealVar bwGamma("bwGamma", "#Gamma", 2.5);

    RooRealVar sg("sg", "sg", sgVal_);
    RooRealVar a("a", "a", aVal_);
    RooRealVar n("n", "n", nVal_);

    RooCBShape CB("CB","CB",*mZ1,bwMean,sg,a,n);
    RooRealVar f("f","f", fVal_);

    RooRealVar mean("mean","mean",meanVal_);
    RooRealVar sigma("sigma","sigma",sigmaVal_);
    RooRealVar f1("f1","f1",f1Val_);

    RooGenericPdf RelBW("RelBW","1/( pow(mZ1*mZ1-bwMean*bwMean,2)+pow(mZ1,4)*pow(bwGamma/bwMean,2) )", RooArgSet(*mZ1,bwMean,bwGamma) );

    RooAddPdf RelBWxCB("RelBWxCB","RelBWxCB", RelBW, CB, f);
    RooGaussian gauss("gauss","gauss",*mZ1,mean,sigma);
    RooAddPdf RelBWxCBxgauss("RelBWxCBxgauss","RelBWxCBxgauss", RelBWxCB, gauss, f1);

    RooProdPdf *PDFRelBWxCBxgauss;
    PDFRelBWxCBxgauss = new RooProdPdf("PDFRelBWxCBxgauss","PDFRelBWxCBxgauss", 
                                     RooArgList(gauss1, gauss2, RelBWxCBxgauss) );
    if(p4sZ1ph_.size()==1)    
      PDFRelBWxCBxgauss = new RooProdPdf("PDFRelBWxCBxgauss","PDFRelBWxCBxgauss", 
                                     RooArgList(gauss1, gauss2, gaussph1, RelBWxCBxgauss) );
    if(p4sZ1ph_.size()==2)
      PDFRelBWxCBxgauss = new RooProdPdf("PDFRelBWxCBxgauss","PDFRelBWxCBxgauss", 
                                     RooArgList(gauss1, gauss2, gaussph1, gaussph2, RelBWxCBxgauss) );

    // observable set
    RooArgSet *rastmp;
      rastmp = new RooArgSet(*pT1RECO,*pT2RECO);
    if(p4sZ1ph_.size()==1)
      rastmp = new RooArgSet(*pT1RECO,*pT2RECO,*pTph1RECO);
    if(p4sZ1ph_.size()>=2)
      rastmp = new RooArgSet(*pT1RECO,*pT2RECO,*pTph1RECO,*pTph2RECO);

    RooDataSet* pTs = new RooDataSet("pTs","pTs", *rastmp);
    pTs->add(*rastmp); 

    //RooAbsReal* nll;
    //nll = PDFRelBWxCBxgauss->createNLL(*pTs);
    //RooMinuit(*nll).migrad();

    RooFitResult* r = PDFRelBWxCBxgauss->fitTo(*pTs,RooFit::Save(),RooFit::PrintLevel(-1));
    const TMatrixDSym& covMatrix = r->covarianceMatrix();
   
    const RooArgList& finalPars = r->floatParsFinal();
    for (int i=0 ; i<finalPars.getSize(); i++){
        TString name = TString(((RooRealVar*)finalPars.at(i))->GetName());

        if(debug_) cout<<"name list of RooRealVar for covariance matrix "<<name<<endl;

    }

    int size = covMatrix.GetNcols();
    //TMatrixDSym covMatrixTest_(size);
    covMatrixZ1_.ResizeTo(size,size);
    covMatrixZ1_ = covMatrix;   

    if(debug_) cout<<"save the covariance matrix"<<endl;
    
    l1 = pT1->getVal()/RECOpT1; l2 = pT2->getVal()/RECOpT2;
    double pTerrZ1REFIT1 = pT1->getError(); double pTerrZ1REFIT2 = pT2->getError();

    pTerrsZ1REFIT_.push_back(pTerrZ1REFIT1);
    pTerrsZ1REFIT_.push_back(pTerrZ1REFIT2);

    if(p4sZ1ph_.size()>=1){

      if(debug_) cout<<"set refit result for Z1 fsr photon 1"<<endl;

      lph1 = pTph1->getVal()/RECOpTph1;
      double pTerrZ1phREFIT1 = pTph1->getError();
      if(debug_) cout<<"scale "<<lph1<<" pterr "<<pTerrZ1phREFIT1<<endl;  
   
      pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT1);

    } 
    if(p4sZ1ph_.size()==2){

      lph2 = pTph2->getVal()/RECOpTph2;
      double pTerrZ1phREFIT2 = pTph2->getError();
      pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT2);

    }

    //delete nll;
    delete r;
    delete mZ1;
    delete pT1; delete pT2; delete pTph1; delete pTph2;
    delete pT1RECO; delete pT2RECO; delete pTph1RECO; delete pTph2RECO;
    delete ph1v3Dph2; delete p1v3Dph1; delete p2v3Dph1; delete p1v3Dph2; delete p2v3Dph2;
    delete PDFRelBWxCBxgauss;
    delete pTs;
    delete rastmp;

    if(debug_) cout<<"end Z1 refit"<<endl;

    return 0;

}

/*
TMatrixDSym GetRefitZ1BigCov(){

  


}
*/
     
/*
RooFormulaVar KinZfitter::p1DOTp2(RooRealVar pT1, RooRealVar theta1, RooRealVar phi1, RooRealVar m1, TString index1, RooRealVar pT2, RooRealVar theta2, RooRealVar phi2, RooRealVar m2, TString index2)
{

    RooFormulaVar cosDeltaPhi("cosDeltaPhi"+index1+index2,"TMath::Cos(@0-@1)",RooArgList(phi1,phi2));

    RooFormulaVar cosTheta1("cosTheta1"+index1,"TMath::Cos(@0)",RooArgList(theta1));
    RooFormulaVar sinTheta1("sinTheta1"+index1,"TMath::Sin(@0)",RooArgList(theta1));
    RooFormulaVar cosTheta2("cosTheta2"+index2,"TMath::Cos(@0)",RooArgList(theta2));
    RooFormulaVar sinTheta2("sinTheta2"+index2,"TMath::Sin(@0)",RooArgList(theta2));

    RooFormulaVar E1("E"+index1,"TMath::Sqrt((@0*@0)/(@1*@1)+@2*@2)", RooArgList(pT1,sinTheta1,m1));
    RooFormulaVar E2("E"+index2,"TMath::Sqrt((@0*@0)/(@1*@1)+@2*@2)", RooArgList(pT2,sinTheta2,m2));
    // 3-vector DOT
    RooFormulaVar p1Dp2("p"+index1+"Dp"+index2,"@0*@1*(@2+(@3*@4)/(@5*@6))",RooArgList(pT1,pT2,cosDeltaPhi,cosTheta1,cosTheta2,sinTheta1,sinTheta2));
    // 4-vector DOT metric 1 -1 -1 -1
    RooFormulaVar p1Dv4p2("p"+index1+"DOTp"+index2,"@0*@1-@2",RooArgList(E1,E2,p1Dp2));


    return p1Dv4p2;
}
*/






#endif
