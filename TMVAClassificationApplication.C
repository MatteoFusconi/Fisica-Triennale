/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <THnSparse.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void TMVAClassificationApplication(Float_t ptmin, Float_t ptmax, TString myMethodList = "") 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader

   Float_t massK0S, tImpParBach, tImpParV0, bachelorPt, CtK0S, cosPAK0S, CosThetaStar, signd0, bachelorP, nSigmaTOFpr, nSigmaTPCpr, nSigmaTPCpi, nSigmaTPCka, bachTPCmom, nSigmaTOFpi, nSigmaTOFka, nSigmapr, dcaV0;
   reader->AddVariable("massK0S", &massK0S);
   reader->AddVariable("tImpParBach", &tImpParBach);
   reader->AddVariable("tImpParV0", &tImpParV0);
   reader->AddVariable("CtK0S := DecayLengthK0S*0.497/v0P", &CtK0S);
   reader->AddVariable("cosPAK0S", &cosPAK0S);
   //reader->AddVariable("CosThetaStar", &CosThetaStar);
   //reader->AddVariable( "nSigmapr := nSigmaTOFpr > -900 ? sqrt(nSigmaTOFpr*nSigmaTOFpr + nSigmaTPCpr*nSigmaTPCpr) : nSigmaTPCpr", &nSigmapr);
   //reader->AddVariable("signd0", &signd0);
   //reader->AddVariable("dcaV0", &dcaV0);
   reader->AddVariable("nSigmaTOFpr", &nSigmaTOFpr);
   reader->AddVariable("nSigmaTOFpi", &nSigmaTOFpi);
   reader->AddVariable("nSigmaTOFka", &nSigmaTOFka);
   reader->AddVariable("nSigmaTPCpr", &nSigmaTPCpr);
   reader->AddVariable("nSigmaTPCpi", &nSigmaTPCpi);
   reader->AddVariable("nSigmaTPCka", &nSigmaTPCka);

   // Spectator variables declared in the training have to be added to the reader, too

   Float_t massLc2K0Sp, LcPt, massLc2Lambdapi, massLambda, massLambdaBar, V0positivePt, V0negativePt, dcaV0pos, dcaV0neg, v0Pt, V0positiveEta, bachelorEta, centrality;
   reader->AddSpectator("massLc2K0Sp", &massLc2K0Sp);   
   reader->AddSpectator("LcPt", &LcPt);
   //reader->AddSpectator("massLambda", &massLambda);
   //reader->AddSpectator("massLambdaBar", &massLambdaBar);
   //reader->AddSpectator("cosPAK0S", &cosPAK0S);
   reader->AddSpectator("V0positivePt", &V0positivePt);
   reader->AddSpectator("V0negativePt", &V0negativePt);
   //reader->AddSpectator("dcaV0pos", &dcaV0pos);
   //reader->AddSpectator("dcaV0neg", &dcaV0neg);
   reader->AddSpectator("v0Pt", &v0Pt);   
   reader->AddSpectator("dcaV0", &dcaV0);
   //reader->AddSpectator("V0positiveEta", &V0positiveEta);
   reader->AddSpectator("bachelorEta", &bachelorEta);
   reader->AddSpectator("centrality", &centrality);

   // --- Book the MVA methods

   TString dir    = "dataset/weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
	 TString weightfile;
	 weightfile = dir + prefix + TString("_") + TString(it->first) + Form("_20220504_%0.0f_%0.0f_11", ptmin, ptmax) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   UInt_t nbin = 100;
   TH2F   *histBDTVsInvMass(0);
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt_prompt(0), *histBdt_bfd(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           {
     histBDTVsInvMass     = new TH2F( "MVA_BDT_vs_InvMass", "MVA_BDT_vs_InvMass; BDT; m_{inv}(pK^{0}_{S})[GeV/#it{c}^{2}]", 10000, -1, 1, 1000, 2.05, 2.55);
     histBdt_prompt       = new TH1F( "MVA_BDT_prompt", "MVA_BDT", 1000, -1.0, 1.0 );
     histBdt_bfd          = new TH1F( "MVA_BDT_bfd", "MVA_BDT", 1000, -1.0, 1.0 );
   }
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   std::cout << "--- Select signal sample" << std::endl;
   TChain* theTree = new TChain("treeList_0_24_0_24_Sgn");
   theTree->AddFile("AnalysisResults_EventMixing_Template.root");
   //theTree->AddFile("../treeData/3520_LHC2016_deghjop/AnalysisResults.root");
   //theTree->AddFile("../treeData/3522_LHC2016_kl/AnalysisResults.root");
   //theTree->AddFile("../treeData/3521_LHC2017_cefhijklmor/AnalysisResults.root");
   //theTree->AddFile("../treeData/3523_LHC2018_bdefghijklmnop/AnalysisResults.root");
   //theTree->AddFile("../treeMC/3002_LHC20I3_P82016/AnalysisResults.root");
   //theTree->AddFile("../treeMC/3003_LHC20I3_P82017/AnalysisResults.root");
   //theTree->AddFile("../treeMC/3004_LHC20I3_P82018/AnalysisResults.root");
   
   //theTree->AddFile("4068_LHC2016_deghjop/AnalysisResults.root");
   //theTree->AddFile("4069_LHC2016_kl/AnalysisResults.root");
   //theTree->AddFile("4070_LHC2017_cefhijklmor/AnalysisResults.root");
   //theTree->AddFile("4071_LHC2018_bdefghijklmnop/AnalysisResults.root");
   
   Float_t massGamma, combinedProtonProb, LcEta, V0negativeEta, TPCProtonProb, TOFProtonProb, LcP, v0P, V0positiveP, V0negativeP, v0Eta, DecayLengthLc, DecayLengthK0S, bachCode, k0SCode, alphaArm, ptArm, weightPtFlat, weightFONLL5overLHC13d3, weightFONLL5overLHC13d3Lc, weightNch, NtrkRaw, NtrkCorr, NtrkAll, origin, SigmacPt, CosThetaStarSoftPi, deltaM;
   
   theTree->SetBranchAddress( "massLc2K0Sp", &massLc2K0Sp);
   //theTree->SetBranchAddress( "alphaArm", &alphaArm);
   theTree->SetBranchAddress( "massK0S", &massK0S);
   //theTree->SetBranchAddress( "massLambda", &massLambda);
   //theTree->SetBranchAddress( "massLambdaBar", &massLambdaBar);
   theTree->SetBranchAddress( "cosPAK0S", &cosPAK0S);
   theTree->SetBranchAddress( "dcaV0", &dcaV0);
   theTree->SetBranchAddress( "tImpParBach", &tImpParBach);
   theTree->SetBranchAddress( "tImpParV0", &tImpParV0);
   theTree->SetBranchAddress( "nSigmaTPCpr", &nSigmaTPCpr);
   theTree->SetBranchAddress( "nSigmaTPCpi", &nSigmaTPCpi);
   theTree->SetBranchAddress( "nSigmaTPCka", &nSigmaTPCka);
   theTree->SetBranchAddress( "nSigmaTOFpr", &nSigmaTOFpr);
   theTree->SetBranchAddress( "nSigmaTOFpi", &nSigmaTOFpi);
   theTree->SetBranchAddress( "nSigmaTOFka", &nSigmaTOFka);
   theTree->SetBranchAddress( "bachelorPt", &bachelorPt);
   theTree->SetBranchAddress( "V0positivePt", &V0positivePt);
   theTree->SetBranchAddress( "V0negativePt", &V0negativePt);
   //theTree->SetBranchAddress( "dcaV0pos", &dcaV0pos);
   //theTree->SetBranchAddress( "dcaV0neg", &dcaV0neg);
   theTree->SetBranchAddress( "v0Pt", &v0Pt);
   //theTree->SetBranchAddress( "SigmacPt", &SigmacPt);
   theTree->SetBranchAddress( "LcPt", &LcPt);
   //theTree->SetBranchAddress( "combinedProtonProb", &combinedProtonProb);
   //theTree->SetBranchAddress( "V0positiveEta", &V0positiveEta);
   theTree->SetBranchAddress( "bachelorEta", &bachelorEta);
   theTree->SetBranchAddress( "v0P", &v0P);
   theTree->SetBranchAddress( "DecayLengthK0S", &DecayLengthK0S);
   //theTree->SetBranchAddress( "NtrkRaw", &NtrkRaw);
   //theTree->SetBranchAddress( "NtrkCorr", &NtrkCorr);
   //theTree->SetBranchAddress( "NtrkAll", &NtrkAll);
   //theTree->SetBranchAddress( "ptArm", &ptArm);
   theTree->SetBranchAddress( "CosThetaStar", &CosThetaStar);
   //theTree->SetBranchAddress( "signd0", &signd0);
   theTree->SetBranchAddress( "centrality", &centrality);
   //theTree->SetBranchAddress( "origin", &origin);
   //theTree->SetBranchAddress( "deltaM", &deltaM);					 
   //theTree->SetBranchAddress( "CosThetaStarSoftPi", &CosThetaStarSoftPi);
   theTree->SetBranchAddress( "bachelorP", &bachelorP);
  
   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   Float_t massLc2K0Sp_old = 0;
   Bool_t breplica = 0;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   Double_t mPr_PDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
   Double_t mK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
   Double_t mLcPDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
   
   Float_t Ebachelor, EbachelorT, Ev0T, Ev0, E2, Et2, P2, Pt2, massBkg, cosThetaInvMass, massT2, massBkgPt, theta1, theta2, y1, y2, deltaY;
   TRandom *rnd = new TRandom(65539);

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
   //for (Long64_t ievt=0; ievt<5000000;ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      CtK0S = DecayLengthK0S*0.497/v0P;
      nSigmapr = nSigmaTOFpr > -900 ? sqrt(nSigmaTOFpr*nSigmaTOFpr + nSigmaTPCpr*nSigmaTPCpr) : nSigmaTPCpr;
      
      if(!( (LcPt < ptmax) && (LcPt > ptmin) )) continue;
      
      // --- Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["BDT"          ])   {
	if (origin == 4) histBdt_prompt->Fill( reader->EvaluateMVA( "BDT method"           ) );
   	else if (origin == 5) histBdt_bfd->Fill( reader->EvaluateMVA( "BDT method"           ) );
	histBDTVsInvMass->Fill( reader->EvaluateMVA("BDT method"), massLc2K0Sp );
      }
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // --- Write histograms
   TFile *target  = new TFile( Form("TMVAApp_BDT_SigmacPt_20220329_%0.0f_%0.0f_RotationalBackground.root", ptmin, ptmax),"RECREATE" );
   //TFile *target  = new TFile( Form("TMVAApp_BDT_SigmacPt_20220504_%0.0f_%0.0f.root", ptmin, ptmax),"RECREATE" );

   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["BDT"          ])   {
     histBdt_prompt->Write();
     histBdt_bfd->Write();
     histBDTVsInvMass->Write();
   }
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
   target->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
  
   delete reader;
    
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
