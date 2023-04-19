/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers.
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables.
/// The methods to be used can be switched on and off by means of booleans, or
/// via the prompt command, for example:
///
///     root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)
///
/// (note that the backslashes are mandatory)
/// If no method given, a default set of classifiers is used.
/// The output file "TMVA.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// Launch the GUI via the command:
///
///     root -l ./TMVAGui.C
///
/// You can also compile and run the example with the following commands
///
///     make
///     ./TMVAClassification <Methods>
///
/// where: `<Methods> = "method1 method2"` are the TMVA classifier names
/// example:
///
///     ./TMVAClassification Fisher LikelihoodPCA BDT
///
/// If no method given, a default set is of classifiers is used
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAClassification
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
//#include "TMVA/AliRDHFCutsLctoV0.h"

int TMVAClassification(Float_t ptmin = 0, Float_t ptmax = 1, TString suffix = "", TString myMethodList = "")
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

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
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   } 

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TDatime date;
   Int_t year = date.GetYear();
   Int_t month = date.GetMonth();
   Int_t day = date.GetDay();
   
   TString outfileName( Form("TMVA_Lc_%s_%d%02d%02d_ptBin_%.0f_%.0f_11.root", suffix.Data(), year, month, day, ptmin, ptmax ));
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
					       "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );
   
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

   //test di correlazione
   dataloader->AddVariable( "massK0S", 'F' );
   dataloader->AddVariable( "tImpParBach", 'F' );
   dataloader->AddVariable( "tImpParV0", 'F' );
   dataloader->AddVariable( "CtK0S := DecayLengthK0S*0.497/v0P", 'F');
   dataloader->AddVariable( "cosPAK0S", 'F');
   //dataloader->AddVariable( "CosThetaStar", 'F');
   //dataloader->AddVariable( "nSigmapr := nSigmaTOFpr > -900 ? sqrt(nSigmaTOFpr*nSigmaTOFpr + nSigmaTPCpr*nSigmaTPCpr) : nSigmaTPCpr", 'F' );
   //dataloader->AddVariable( "signd0", 'F');
   dataloader->AddVariable( "nSigmaTOFpr", 'F' );
   dataloader->AddVariable( "nSigmaTOFpi", 'F' );
   dataloader->AddVariable( "nSigmaTOFka", 'F' );
   dataloader->AddVariable( "nSigmaTPCpr", 'F' );
   dataloader->AddVariable( "nSigmaTPCpi", 'F' );
   dataloader->AddVariable( "nSigmaTPCka", 'F' );  
   //dataloader->AddVariable( "massLc2K0Sp", 'F' ); 
   //dataloader->AddVariable( "LcPt", 'F' ); 
   //dataloader->AddVariable( "dcaV0", 'F' ); 
   //dataloader->AddVariable( "v0Pt", 'F' );
   //dataloader->AddVariable( "combinedProtonProb", 'F' ); 
   //dataloader->AddVariable( "SigmacPt", 'F');
   //dataloader->AddVariable( "CosThetaStarSoftPi", 'F');
   //dataloader->AddVariable( "deltaM", 'F');

   
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   dataloader->AddSpectator( "massLc2K0Sp", "mass Lc-->K0S+p", "units", 'F' );
   dataloader->AddSpectator( "LcPt",  "Lc Pt", "units", 'F' );
   //dataloader->AddSpectator( "massLc2Lambdapi",  "mass Lc -->L(1520)+pi", "units", 'F' );
   //dataloader->AddSpectator( "massLambda", "mass V0 = Lambda", "units", 'F' );
   //dataloader->AddSpectator( "massLambdaBar", "mass V0 = LambdaBar", "units", 'F' );
   //dataloader->AddSpectator( "cosPAK0S", "cosPointingAngle K0S", "units", 'F' );
   dataloader->AddSpectator( "V0positivePt", "V0 positive Pt", "units", 'F' );
   dataloader->AddSpectator( "V0negativePt", "V0 negative Pt", "units", 'F' );
   //dataloader->AddSpectator( "dcaV0pos", "dca V0 positive", "units", 'F' );
   //dataloader->AddSpectator( "dcaV0neg", "dca V0 negative", "units", 'F' );
   dataloader->AddSpectator( "v0Pt", "K0S Pt", "units", 'F');	
   dataloader->AddSpectator( "dcaV0", "dca V0", "units", 'F');	
   //dataloader->AddSpectator( "V0positiveEta", "V0pos Eta", "units", 'F');	
   dataloader->AddSpectator( "bachelorEta", "eta bachelor", "units", 'F');	
   dataloader->AddSpectator( "centrality", "centrality", "units", 'F');
    
   // Read training and test data
   TString fnameSgn1, fnameSgn2, fnameSgn3;
   //fnameSgn1 = "../treeMC/2687_LHC20I3_2016/AnalysisResults.root";
   //fnameSgn2 = "../treeMC/2688_LHC20I3_2017/AnalysisResults.root";
   //fnameSgn3 = "../treeMC/2689_LHC20I3_2018/AnalysisResults.root";
   fnameSgn1 = "3002_LHC20I3_P82016/AnalysisResults.root";
   fnameSgn2 = "3003_LHC20I3_P82017/AnalysisResults.root";
   fnameSgn3 = "3004_LHC20I3_P82018/AnalysisResults.root";
   
   if (gSystem->AccessPathName( fnameSgn1 ) || gSystem->AccessPathName( fnameSgn2 ) || gSystem->AccessPathName( fnameSgn3 )){  // file does not exist in local directory 
     Printf("Signal File for training does not exist, check please, and retry. Now I'll return..."); return 0; 
   }
   
   TString fnameBkg1, fnameBkg2, fnameBkg3, fnameBkg4;
   //fnameBkg1 = "../treeData/3520_LHC2016_deghjop/AnalysisResults.root";
   //fnameBkg2 = "../treeData/3522_LHC2016_kl/AnalysisResults.root";
   //fnameBkg3 = "../treeData/3521_LHC2017_cefhijklmor/AnalysisResults.root";
   //fnameBkg4 = "../treeData/3523_LHC2018_bdefghijklmnop/AnalysisResults.root";
   fnameBkg1 = "4068_LHC2016_deghjop/AnalysisResults.root";
   fnameBkg2 = "4069_LHC2016_kl/AnalysisResults.root";
   fnameBkg3 = "4070_LHC2017_cefhijklmor/AnalysisResults.root";
   fnameBkg4 = "4071_LHC2018_bdefghijklmnop/AnalysisResults.root";
   
   if (gSystem->AccessPathName( fnameBkg1 ) || gSystem->AccessPathName( fnameBkg2 ) || gSystem->AccessPathName( fnameBkg3 ) || gSystem->AccessPathName( fnameBkg4 )){  // file does not exist in local directory
     Printf("Signal File for training does not exist, check please, and retry. Now I'll return..."); return 0;
   }

   TFile *inputSgn1 = TFile::Open( fnameSgn1 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputSgn1->GetName() << std::endl;
   TFile *inputSgn2 = TFile::Open( fnameSgn2 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputSgn2->GetName() << std::endl;
   TFile *inputSgn3 = TFile::Open( fnameSgn3 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputSgn3->GetName() << std::endl;
   TFile *inputBkg1 = TFile::Open( fnameBkg1 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputBkg1->GetName() << std::endl;
   TFile *inputBkg2 = TFile::Open( fnameBkg2 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputBkg2->GetName() << std::endl;
   TFile *inputBkg3 = TFile::Open( fnameBkg3 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputBkg3->GetName() << std::endl;
   TFile *inputBkg4 = TFile::Open( fnameBkg4 );
   std::cout << "--- TMVAClassification       : Using input file: " << inputBkg4->GetName() << std::endl;   
   // --- Register the training and test trees
   TTree *signal1      = (TTree*)inputSgn1->Get("treeList_0_24_0_24_Sgn");
   TTree *signal2      = (TTree*)inputSgn2->Get("treeList_0_24_0_24_Sgn");
   TTree *signal3      = (TTree*)inputSgn3->Get("treeList_0_24_0_24_Sgn");
   TTree *background1  = (TTree*)inputBkg1->Get("treeList_0_24_0_24_Sgn");
   TTree *background2  = (TTree*)inputBkg2->Get("treeList_0_24_0_24_Sgn");
   TTree *background3  = (TTree*)inputBkg3->Get("treeList_0_24_0_24_Sgn");
   TTree *background4  = (TTree*)inputBkg4->Get("treeList_0_24_0_24_Sgn");
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signal1,      signalWeight     );
   dataloader->AddSignalTree    ( signal2,      signalWeight     );
   dataloader->AddSignalTree    ( signal3,      signalWeight     );
   dataloader->AddBackgroundTree( background1,  backgroundWeight );
   dataloader->AddBackgroundTree( background2,  backgroundWeight );
   dataloader->AddBackgroundTree( background3,  backgroundWeight );
   dataloader->AddBackgroundTree( background4,  backgroundWeight );

   Int_t nSigmaSideBands = 3;
   Float_t sideBands = 0;
   if(ptmin==0 && ptmax==1) {
     sideBands = nSigmaSideBands * 0.0076;
   }
   else if(ptmin==1 && ptmax==2) {
     sideBands = nSigmaSideBands * 0.0076;
   }
   else if(ptmin==2 && ptmax==3) {
     sideBands = nSigmaSideBands * 0.0076;
   }
   else if(ptmin==2 && ptmax==4) {
     sideBands = nSigmaSideBands * 0.0077;
   }
   else if(ptmin==3 && ptmax==4) {
     sideBands = nSigmaSideBands * 0.0079;
   }
   else if(ptmin==4 && ptmax==5) {
     sideBands = nSigmaSideBands * 0.0084;
   }
   else if(ptmin==4 && ptmax==6) {
     sideBands = nSigmaSideBands * 0.0078;
   }
   else if(ptmin==5 && ptmax==6) {
     sideBands = nSigmaSideBands * 0.0090;
   }
   else if(ptmin==6 && ptmax==7) {
     sideBands = nSigmaSideBands * 0.0096;
   }
   else if(ptmin==6 && ptmax==8) {
     sideBands = nSigmaSideBands * 0.0098;
   }
   else if(ptmin==7 && ptmax==8) {
     sideBands = nSigmaSideBands * 0.0101;
   }
   else if(ptmin==8 && ptmax==10) {
     sideBands = nSigmaSideBands * 0.0108;
   }
   else if(ptmin==8 && ptmax==12) {
     sideBands = nSigmaSideBands * 0.0111;
   }   
   else if(ptmin==10 && ptmax==12) {
     sideBands = nSigmaSideBands * 0.0118;
   }
   else if(ptmin==12 && ptmax==24) {
     sideBands = nSigmaSideBands * 0.0136;
   }

   // Apply additional cuts on the signal and background samples (can be different)
   
   TCut mycuts, mycutb;
   Int_t nTrainingEventsSgn, nTrainingEventsBkg, nTestingEventsBkg;
   Float_t isFromC = 4;
   
   mycuts = Form("LcPt < %f && LcPt > %f && origin == %f", ptmax, ptmin, isFromC); //only prompt
   mycutb = Form("LcPt < %f && LcPt > %f && abs(massLc2K0Sp - 2.286) > %f", ptmax, ptmin, sideBands);
   nTrainingEventsSgn = TMath::Min((signal1->GetEntries(mycuts)+signal2->GetEntries(mycuts)+signal3->GetEntries(mycuts)) * 0.5, 500000.);
   nTrainingEventsBkg = TMath::Min((background1->GetEntries(mycutb)+background2->GetEntries(mycutb)+background3->GetEntries(mycutb)+background4->GetEntries(mycutb)) * 0.6, 500000.);
   nTestingEventsBkg = TMath::Min((background1->GetEntries(mycutb)+background2->GetEntries(mycutb)+background3->GetEntries(mycutb)+background4->GetEntries(mycutb)) * 0.4, 500000.);
   
   
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
   					   Form("nTrain_Signal=%0d:nTest_Signal=%0d:nTrain_Background=%0d:nTest_Background=%0d:SplitMode=Random:NormMode=NumEvents:!V", nTrainingEventsSgn, nTrainingEventsSgn, nTrainingEventsBkg, nTestingEventsBkg) );

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher",
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators


   // Multi-architecture DNN implementation.
   if (Use["DNN_CPU"] or Use["DNN_GPU"]) {
      // General layout.
      TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

      // Training strategies.
      TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
      TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training0 + "|" + training1 + "|" + training2;

      // General Options.
      TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                          "WeightInitialization=XAVIERUNIFORM");
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

      // Cuda implementation.
      if (Use["DNN_GPU"]) {
         TString gpuOptions = dnnOptions + ":Architecture=GPU";
         factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_GPU", gpuOptions);
      }
      // Multi-core CPU implementation.
      if (Use["DNN_CPU"]) {
         TString cpuOptions = dnnOptions + ":Architecture=CPU";
         factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", cpuOptions);
      }
   }

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
 
   //Default Boosted Decision Trees
   TString methodTitle;
   if (Use["BDT"]){  // Adaptive Boost
     methodTitle = Form("BDT_Default_%.0f_%.0f", ptmin, ptmax);  
     Printf("Default ON: methodTitle will be %s", methodTitle.Data());
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
			  "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   }
   
   if (Use["BDTB"]) // Bagging
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory
   //
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAt001","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------


   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   
   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   //if (!gROOT->IsBatch()) TMVAGui( outfileName );
   
   return 0;

}
   
int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   Float_t ptmin = 0, ptmax = 1;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return TMVAClassification(ptmin, ptmax, methodList,"");
}
