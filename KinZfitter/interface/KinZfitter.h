/*************************************************************************
*  Authors:   Tongguang CHeng(IHEP, Beijing) Hualin Mei(UF)
*************************************************************************/
#ifndef KinZfitter_h
#define KinZfitter_h

// C++ includes
#include <iostream>
#include <complex>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
// ROOT includes
#include "TString.h"
#include "TLorentzVector.h"

// CMSSW related pT error calculator
#include "KinZfitter/HelperFunction/interface/HelperFunction.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// ROOFIT

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooCBShape.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"

#include "RooFitResult.h"

// fit result covariance matrix
#include <TMatrixDSym.h>

#include <iostream>
#include <map>

// helper function calculate lepton/pfphoton (un)corrected pT error
class HelperFunction;

using namespace std;


class KinZfitter {
public:
	
        KinZfitter(bool isData);

	/// Kinematic fit of lepton momenta
        /// HelperFunction class to calcluate per lepton(+photon) pT error
        void Setup(std::vector< reco::Candidate* > selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons);

        ///
        void KinRefitZ1();

        int  PerZ1Likelihood(double & l1, double & l2, double & lph1, double & lph2);
        void SetZ1Result(double l1, double l2, double lph1, double lph2);

        // result wrappers
        double GetRefitM4l();
        double GetM4l();
        double GetRefitMZ1();

        double GetMZ1Err();
        double GetRefitM4lErr();
        double GetM4lErr();
        double GetRefitM4lErrFullCov();

        // cov matrix change for spherical coordinate to Cartisean coordinates

        void SetZ1BigCov();
        void SetZ2BigCov();

/*
        TMatrixDSym GetRefitZ1BigCov();        

        // cov matrix when refitting Z2
        //void SetBigCovZZ(); cov matrix when 4e/4mu final state and both Z needs to be refitted

        TMatrixDSym GetRefitZZBigCov();
        //  
        TMatrixDSym GetRefitZZSFBigCov();
        TMatrixDSym GetRefitZZOFBigCov();
*/

        std::vector<TLorentzVector> GetRefitP4s();
        std::vector<TLorentzVector> GetP4s();

        ////////////////

        //Need to check in deep whether this function is doable 
        //RooFormulaVar p1DOTp2(RooRealVar pT1, RooRealVar theta1, RooRealVar phi1, RooRealVar m1, TString index1, RooRealVar pT2, RooRealVar theta2, RooRealVar phi2, RooRealVar m2, TString index2);

private:

        /// True mZ/mZ1 shape, final states
        TString PDFName_, fs_;      
	
	/// debug flag
	bool debug_;
       
        /// whether use correction for pT error
        bool isCorrPTerr_; 	
        /// whether use data or mc correction
        bool isData_;

        /// HelperFunction class to calcluate per lepton(+photon) pT error
        HelperFunction * helperFunc_;

        void initZs(std::vector< reco::Candidate* > selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhoton);

        /// lepton ids for Z1 Z2
        std::vector<int> idsZ1_, idsZ2_;
        /// lepton ids that fsr photon associated to
        std::vector<int> idsFsrZ1_, idsFsrZ2_;
        /// (Four) TLorentzVectors that form the Higgs Candidate 
        std::vector<TLorentzVector> p4sZ1_, p4sZ2_, p4sZ1ph_, p4sZ2ph_;
        std::vector<TLorentzVector> p4sZ1REFIT_, p4sZ2REFIT_, p4sZ1phREFIT_, p4sZ2phREFIT_;

        /// pTerr vector
        std::vector<double> pTerrsZ1_, pTerrsZ2_, pTerrsZ1ph_, pTerrsZ2ph_;
        std::vector<double> pTerrsZ1REFIT_, pTerrsZ2REFIT_, pTerrsZ1phREFIT_, pTerrsZ2phREFIT_;

        // covariance matrix 
        // what directly coming from Refit
        TMatrixDSym covMatrixZ1_, covMatrixZ2_, covMatrixZZ_;
        // covariance matrix in the Cartesian coordinates
        TMatrixDSym bigCovMatrix_;

        // refit energy scale with respect to reco pT
        double lZ1_l1_, lZ1_l2_, lZ2_l1_, lZ2_l2_;
        double lZ1_ph1_, lZ1_ph2_, lZ2_ph1_, lZ2_ph2_;
        // True mZ1 shape parameters
        double sgVal_, aVal_, nVal_, fVal_, meanVal_, sigmaVal_, f1Val_;

};

#endif
