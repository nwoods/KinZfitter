# KinZfitter
Using Z mass kinematic constraint(s) to refit lepton momenta in H-to-ZZ*

- Code structure

KinZfitter : the class read the inputs from leptons and fsr photons that form the Higgs Candidate,, do the refitting, and get the refitted results
HelperFunction : the class that read lepton/photon pT errors by accessing the pat:Electron, pat:Muon and pat:PFPartticle and also provide the function to help calculate mass4l error (including the covariance matrix)

To include the refit in your analyzer:

0.Check out package

  cd $CMSSW_BASE/src
  
  git clone https://github.com/tocheng/KinZfitter.git
  
  scram b -j 8 

In your main analyzer:


0. Add the package into your BuildFile.xml
   <use   name="KinZfitter/KinZfitter"/>

1.include the head file

  #include "KinZfitter/KinZfitter/interface/KinZfitter.h"

2.Declare and then initialize the KinZfitter class when initializing your analyzer

    KinZfitter *kinZfitter;
    kinZfitter = new KinZfitter(isData);
    //(In data, (isData=true). In mc (isData=false))

3.Prepare inputs after Higgs candidate is formed:

  leptons: 
   Suppose Lep_Z1_1,Lep_Z1_2, Lep_Z2_1,Lep_Z2_2 are pat:Electron or pat:Muon from Z1 and Z2 decay.
   
   In your analyzer, do:

     vector<reco::Candidate *> selectedLeptons;
     reco::Candidate *cZ1_1 = dynamic_cast<reco::Candidate* >(&Lep_Z1_1);
     reco::Candidate *cZ1_2 = dynamic_cast<reco::Candidate* >(&Lep_Z1_2);
     reco::Candidate *cZ2_1 = dynamic_cast<reco::Candidate* >(&Lep_Z2_1);
     reco::Candidate *cZ2_2 = dynamic_cast<reco::Candidate* >(&Lep_Z2_2);
     selectedLeptons.push_back(cZ1_1);
     selectedLeptons.push_back(cZ1_2);
     selectedLeptons.push_back(cZ2_1);
     selectedLeptons.push_back(cZ2_2);

  fsrPhotons :
    Supporse find pat::PFParticle fsrPhoton which is the fsr photon.
    To fill the array if the photon is associated to a certain lepton from Z1 or Z2 decay, 
    do something like:

     pat::PFParticle[4] selectedFsrPhotonsArray;
     if(associateLeptonZ1_1) selectedFsrPhotonsArray[0] = fsrPhoton;
     if(associateLeptonZ1_2) selectedFsrPhotonsArray[1] = fsrPhoton;
     if(associateLeptonZ2_1) selectedFsrPhotonsArray[2] = fsrPhoton;
     if(associateLeptonZ2_2) selectedFsrPhotonsArray[3] = fsrPhoton;
 
    Keep the order the same as vector<reco::Candidate *> selectedLeptons

4.Setup, refit and get the refitted results:

   In your analyzer, do

      kinZfitter->Setup(selectedLeptons, selectedFsrPhotonsArray);
      kinZfitter->KinRefitZ1();
      
      // refit mass4l
      double mass4lREFIT = kinZfitter->GetRefitM4l();
      // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
      vector < TLorentzVector > p4 = kinZfitter->GetRefitP4s(); 

  There is a function called GetRefitM4lErr() which calculates mass4l error after refitting 
  assuming that all lepton momenta are UNcorrelated 
  which is good approximation for reco lepton momenta BUT UNTRUE for refitted lepton momenta.
  This function will be extended to calculated the mass4l using full covariance matrix soon.





