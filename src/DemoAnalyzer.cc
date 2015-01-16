// -*- C++ -*-
//
// Package:    DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/src/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Pesaresi
//         Created:  Tue Oct 15 18:32:56 BST 2013
// $Id$
//
// Modified by: Alexander Morton
//


// system include files
#include <memory>
#include <cstdlib>
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/Common/interface/DetSet.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"

#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2D.h"

#include <stdlib.h>
#include <time.h>

const bool PIXEL_DIGI_LOOP ( true );
// const bool PIXEL_DIGI_LOOP ( false );

const bool TRACKING_PARTICLE_LOOP ( true );
// const bool TRACKING_PARTICLE_LOOP ( false );

const bool STUB_LOOPS ( true );
// const bool STUB_LOOPS ( false );

// const bool EXCLUDE_INTERACTING_PRTS ( true );
const bool EXCLUDE_INTERACTING_PRTS ( false );

// const bool WRITE_TEXTFILE ( true );
const bool WRITE_TEXTFILE ( false );

const bool WRITE_TEXTFILE_EXC_2GEV ( true );
// const bool WRITE_TEXTFILE_EXC_2GEV ( false );

// const bool TEXTFILE_OFFSET( true );
const bool TEXTFILE_OFFSET ( false );

const int PARTICLE_ID (221);           // Type of particles to be written to text file, #13 for muons, #11 for electrons, #221 postive pions
const int NUM_STUBS_OFFSET ( 5 );     // Max number of stubs to offset during offsetting
const int MAX_OFFSET_AMOUNT ( 200 );  // In microns

const std::string lFileName ("/home/hep/adm10/BDT4PT/source_data/Outfile_Pions_2-250GeV.txt");

class DemoAnalyzer : public edm::EDAnalyzer {

   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


//       virtual void beginRun(edm::Run const&, edm::EventSetup const&);
//       virtual void endRun(edm::Run const&, edm::EventSetup const&);
//       virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
//       virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      // ----------member data ---------------------------

      edm::ESHandle<TrackerGeometry> GeometryHandle;
      edm::ESHandle< StackedTrackerGeometry >  stackedGeometryHandle;

      const StackedTrackerGeometry*            theStackedGeometry;

      std::ofstream myFile; // Output Particle Data to textfile for track fitting

//    VIVA PLOTS

  // Efficiency
      TH1F * h_PrimaryParticlePt;
      TH1F * h_FivePlusStubsPrts;     // Frequency of particles with 5 or more stubs
      TH1F * h_LessThanFiveStubsPrts;     // Frequency of particles with less than stubs
      TH1F * h_FrequencyOfNumStubs;      // Frequency of number of accepted stubs produced by Tracking Particles
      TH1F * h_FivePlusStubsFraction; // Proportion of partcles (pT) that create 5 or more stubs
      TH1F * h_numberStubs;
      TH1F * h_numberStubsFrac;

      TH1F * h_TrkPrtPt_Frac;

      TH1F * h_NumTrkPrtsFivePlusAccStubs;
      TH1F * h_NumTrkPrtsFivePlusRejStubs;

//    Not implemented yet
      TH1F * h_AccStubFrequency;      // Frequency of Accepted Stubs(pT)
      TH1F * h_AccGenStubFrequency;      // Frequency of Genuine Accepted Stubs(pT)
      TH1F * h_AllStubsFrequency;     // Frequency of Stubs (pT)
      TH1F * h_AccStubsOutOfAllStubs; // Proportion of Accepted Stubs compared to All Stubs (pT)

      TH1F * h_AccStubEta;
      TH1F * h_AccStubPhi;
      TH2F * h_AccStubEtaPhi;
      TH2F * h_AccStubEtaPt;
      TH2F * h_AccStubPhiPt;

      TH1F * h_RejStubEta;
      TH1F * h_RejStubPhi;
      TH2F * h_RejStubEtaPhi;
      TH2F * h_RejStubEtaPt;
      TH2F * h_RejStubPhiPt;

  // Likelihood

  // Higgs

      // Higgs->4thatcreate4acceptedstubs
      // proportion that are gen?

  // Track Fitting - See BD4PT PDF folder

// Old Plots

      TH2D * h_etaphi;  // Pixel Digis eta/phi
      TH1F * h_pTdist;  // pT dist of TrkPrts that create stubs
      TH1F * h_pTdistAccGenStub; // pT dist of TrkPrts that create Acc Gen Stubs
      TH1F * h_pTdistRejGenStub; // pT dist of TrkPrts that create Rej Gen Stubs
      TH1F * h_TrkPrtPt;  // pT dist of all TrkPrts

      TH1F * h_TrkPrtAccAllStubsPt; // pT dist of Tracking Particles that produced  Accepted Stubs (for each stub)
      TH1F * h_TrkPrtAccGenStubsPt; // pT dist of Tracking Particles that produced Accepted Genuine Stubs (for each stub)
      TH1F * h_TrkPrtAccComStubsPt; // pT dist of Tracking Particles that produced Accepted Combinatorial Stubs (for each stub)
      TH1F * h_TrkPrtAccUnkStubsPt; // pT dist of Tracking Particles that produced Accepted Unknown Stubs (for each stub)
      TH1F * h_TrkPrtRejAllStubsPt; // pT dist of Tracking Particles that produced Rejected Stubs (for each stub)
      TH1F * h_TrkPrtRejGenStubsPt; // pT dist of Tracking Particles that produced Rejected Genuine Stubs (for each stub)
      TH1F * h_TrkPrtRejComStubsPt; // pT dist of Tracking Particles that produced Rejected Combinatorial Stubs (for each stub)
      TH1F * h_TrkPrtRejUnkStubsPt; // pT dist of Tracking Particles that produced Rejected Unknown Stubs (for each stub)

      TH1F * h_TrkPrtAccAllPt; // pT dist of Tracking Particles that produced Accepted Stubs (for each TrkPrt)
      TH1F * h_TrkPrtAccGenPt; // pT dist of Tracking Particles that produced Accepted Genuine Stubs (for each TrkPrt)
      TH1F * h_TrkPrtAccComPt; // pT dist of Tracking Particles that produced Accepted Combinatorial Stubs (for each TrkPrt)
      TH1F * h_TrkPrtAccUnkPt; // pT dist of Tracking Particles that produced Accepted Unknown Stubs (for each TrkPrt)
      TH1F * h_TrkPrtRejAllPt; // pT dist of Tracking Particles that produced Rejected Stubs (for each TrkPrt)
      TH1F * h_TrkPrtRejGenPt; // pT dist of Tracking Particles that produced Rejected Genuine Stubs (for each TrkPrt)
      TH1F * h_TrkPrtRejComPt; // pT dist of Tracking Particles that produced Rejected Combinatorial Stubs (for each TrkPrt)
      TH1F * h_TrkPrtRejUnkPt; // pT dist of Tracking Particles that produced Rejected Unknown Stubs (for each TrkPrt)


      TH1F * h_pTdistAccGenStub_layer[6]; // pT dist of TrkPrts that create Acc Gen Stubs for each of the Barrel Layers
      TH1F * h_pTdistRejGenStub_layer[6]; // pT dist of TrkPrts that create Rej Gen Stubs for each of the Barrel Layers

      TH1F * h_pTdistAccGenStub_disk[5]; // pT dist of TrkPrts that create Rej Gen Stubs for each of the Disk Layers
      TH1F * h_pTdistRejGenStub_disk[5]; // pT dist of TrkPrts that create Rej Gen Stubs for each of the Disk Layers
//
      TH2D * h_stub_phi; // Acc Gen Stub pT as a function of phi
      TH2D * h_stub_eta; // Acc Gen Stub pT as a function of eta
      TH2D * h_stub_etaphi; // Acc Gen Stub Eta/Phi Distribution
      TH2D * h_stub_xy; // Acc Gen Stub pT xy dist
      TH2D * h_stub_zr; // Acc Gen Stub zr dist

      TH2D * h_rejstub_phi; // Acc Rej Stub pT as a function of phi
      TH2D * h_rejstub_eta; // Acc Gen Stub pT as a function of eta
      TH2D * h_rejstub_etaphi; // Acc Gen Stub Eta/Phi Distribution
      TH2D * h_rejstub_xy; // Acc Gen Stub pT xy dist
      TH2D * h_rejstub_zr; // Acc Gen Stub zr dist
//
      TH1F * h_stub_eta_distb; // Number of Acc Gen Stubs over eta
      TH1F * h_stub_phi_distb; // Number of Acc Gen Stubs over phi
      TH1F * h_stub_perp_barrel_distb; // Number of Acc Gen Stubs over perp (barrel)
      TH1F * h_stub_perp_endcap_distb; // Number of Acc Gen Stubs over perp (disk)
      TH1F * h_stub_pT_distb; // Number of Acc Gen Stubs over pT
      TH1F * h_stub_pT_layer_distb[6]; // Number of Acc Gen Stubs over pT for each layer
      TH1F * h_stub_pT_disk_distb[5]; // Number of Acc Gen Stubs over pT for each disk

      TH1F * h_rejstub_eta_distb; // Number of Rej Gen Stubs over eta
      TH1F * h_rejstub_phi_distb; // Number of Rej Gen Stubs over phi
      TH1F * h_rejstub_perp_barrel_distb; // Number of Rej Gen Stubs over perp (barrel)
      TH1F * h_rejstub_perp_endcap_distb; // Number of Rej Gen Stubs over perp (disk)
      TH1F * h_rejstub_pT_distb; // Number of Rej Gen Stubs over pT
      TH1F * h_rejstub_pT_layer_distb[6]; // Number of Rej Gen Stubs over pT for each layer
      TH1F * h_rejstub_pT_disk_distb[5]; // Number of Rej Gen Stubs over pT for each disk
//
      TH1F * h_stubdiff_pT_distb;       // Fraction of Acc/Rej Gen Stubs pT distribution
      TH1F * h_accTrkPrt_fraction;      // Fraction of Acc Gen Trk Prts/Total Genuine Stub Creating Tracking Particles, pT dist

      TH1F * h_stubdiff_pT_layer_distb[6]; // Fraction of Acc/Rej Gen Stubs pT distribution for each layer
      TH1F * h_accTrkPrt_layer_fraction[6];  // Fraction of Acc Gen Trk Prts/Total Genuine Stub Creating Tracking Particles, pT dist, for each layer

      TH1F * h_stubdiff_pT_disk_distb[5]; // Fraction of Acc/Rej Gen Stubs pT distribution for each disk
      TH1F * h_accTrkPrt_disk_fraction[5];  // Fraction of Acc Gen Trk Prts/Total Genuine Stub Creating Tracking Particles, pT dist, for each disk

      TH1F * h_AccStubPrts_TotalTrkPrts_fraction; // Proportion of Acc Stubs created from all TrkPrts at certain pTs
      TH1F * h_RejStubPrts_TotalTrkPrts_fraction; // Proportion of Rej Stubs created from all TrkPrts at certain pTs
      TH1F * h_GenAccStubPrts_TotalTrkPrts_fraction; // Proportion of Gen Acc Stubs created from all TrkPrts at certain pTs
      TH1F * h_GenAccStubPrts_AccStubPrts_fraction; // Proportion of Acc Stubs created from all Acc Stubs at certain pTs

      TH1F * h_pTMuonsFromHiggs;                // pT Dist of the (4) muons produced from Higgs decay
      TH1F * h_pTElectronsFromHiggs;            // pT Dist of the (4) electrons produced from Higgs decay
      TH1F * h_pTElectronsMuonsFromHiggs;       // pT Dist of the (2) muons and (2) electrons produced from Higgs decay

      TH1F * h_pTMuonsFromHiggsGenStubs;          // pT Dist of the (4) muons produced from Higgs decay that produce Gen Acc Stubs
      TH1F * h_pTElectronsFromHiggsGenStubs;      // pT Dist of the (4) electrons produced from Higgs decay that produce Gen Acc Stubs
      TH1F * h_pTElectronsMuonsFromHiggsGenStubs; // pT Dist of the (2) muons and (2) electrons produced from Higgs decay that produce Gen Acc Stubs

      int MuonPrtCounter;         // Counter for number of Muons in file
      int ElectronPrtCounter;     // Counter for number of Electrons in file
      int TauPrtCounter;          // Counter for number of Taus in file
      int ZPrtCounter;            // Counter for number of Z's in file
      int WPrtCounter;            // Counter for number of W's in file
      int ChargedPionPrtCounter;  // Counter for number of Charged Pions in file
      int NeutralPionPrtCounter;  // Counter for number of Pions in file
      int PhotonPrtCounter;       // Counter for number of Photons in file
      int HiggsPrtCounter;        // Counter for number of Higgs in file
      int OtherPrtCounter;        // Counter for number of other particles in file
      int lPrtCounter;

      int GenAccStubCreatingPrtCounter; // Number of TrkPrts that create GenAccStubs
      int StubCreatingPrtCounter;       // Number of TrkPrts that create GenStubs
      int l4MuCounter;                  // Number of clean H->4mu stubs created
      int l4eCounter;                   // Number of clean H->4e stubs created
      int l2e2MuCounter;                // Number of clean H->2e2mu stubs created

      int lCounterGenStub;
      int lCounterCombStub;
      int lCounterUnkStub;

// Mariza's Bit?

///////////////////////////////////////////////////////////////////////
//stubs vs real particles per energy bin
      TH1F * h_stubTrkPrtdiff_distb;    // Mariza's Plot?

//counters for genuine/non genuine


///////////////////////////////////////////////////////////////////////


};

//
// constants, enums   
// and typedefsprocess.mix.input.nbPileupEvents.averageNumber = cms.double(140)
//
//
// static data member definitions
//

//
// constructors and destructorintel atom
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<edm::RandomNumberGenerator> rng;
  if(!rng.isAvailable()) {
    throw cms::Exception("Configuration")
      << "DemoAnalyzer requires the RandomNumberGeneratorService,\n"
         "which is not present in the configuration file. You must add\n"
         "the service in the configuration file or remove the modules that\n"
         "require it.\n";
  }
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

//   double lFlagFirstStubOnlyAccStubPt (0.0);
//   double lFlagFirstStubOnlyRejStubPt (0.0);

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
// Handles and Labels
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------

  // Geometry setup
  iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
  theStackedGeometry = stackedGeometryHandle.product();
  iSetup.get<TrackerDigiGeometryRecord>().get(GeometryHandle);

// Pixel Digi Handles
   edm::Handle<DetSetVector<PixelDigi> > theDigis;
   iEvent.getByLabel("simSiPixelDigis", theDigis);

// Tracking Particle Handles

  edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
  iEvent.getByLabel("mix", "MergedTrackTruth", TruthTrackContainer);

//  Stub Handles

  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubRejectedHandle;
  iEvent.getByLabel( "TTStubsFromPixelDigis", "StubAccepted", TTStubHandle );
  iEvent.getByLabel( "TTStubsFromPixelDigis", "StubRejected", TTStubRejectedHandle );

// MC Truth for Stubs Handles and iEvent.getByLabels

  edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > >                MCTruthTTStubHandle;   // TTStubAssociation Map Handle for Accepted Stubs
  iEvent.getByLabel( "TTStubAssociatorFromPixelDigis", "StubAccepted", MCTruthTTStubHandle );

  edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > >                MCTruthTTStubRejectedHandle; // TTStubAssociation Map Handle for Rejected Stubs
  iEvent.getByLabel( "TTStubAssociatorFromPixelDigis", "StubRejected", MCTruthTTStubRejectedHandle );

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
// Write Stub info to textfile for trackfitting
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
  int lFlagFirstStubTextFileLoop (0);
  int lNumerOffsetStubs (0);            // Counts how many stubs have been offset

  int lTempPrtCountTxtLoop(0);  // Counter to be incremented after each particle in the event has looped over, to correct the next one.

  // WRITE STUB TO TEXT FILE LOOP

  if ( WRITE_TEXTFILE )
  {
    for ( std::vector<TrackingParticle>::const_iterator lTrkPrtIt = TruthTrackContainer->begin(); lTrkPrtIt != TruthTrackContainer->end(); lTrkPrtIt++)
    {
      edm::Ptr<TrackingParticle> tempTPPtr(TruthTrackContainer, lTempPrtCountTxtLoop++);  // Create temporary TrkPrt Pointer, first arguement states we are using the TrkPrt Truth Container, second arguement states which particle in the container in the current event we are looking at.

      const double lPrtPt(tempTPPtr->pt()), lPrtEta(tempTPPtr->eta()), lPrtPhi(tempTPPtr->phi());  // Grab TrkPrt Info
      const int lPrtId(tempTPPtr->pdgId());

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > > theStubs = MCTruthTTStubHandle->findTTStubRefs(tempTPPtr);

//         std::cout << "Number of stubs: " << theStubs.size() << std::endl;
// --------------------------------------------------------------------------------------------------------------------------------------------------
    // Accepted stubs associated with Tracking Particle
// --------------------------------------------------------------------------------------------------------------------------------------------------
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > >::const_iterator lItTrkAccStubs = theStubs.begin();
      for ( ; lItTrkAccStubs != theStubs.end(); ++lItTrkAccStubs )
      {

        if ( MCTruthTTStubHandle->isGenuine(*lItTrkAccStubs) )
        {

//           std::cout << "Genuine Stub (Textfile Output Loop)" << std::endl;

          std::stringstream lStr;
          // Grab associated TrkPt Parameters

          if ( lPrtPt < 2.0 && WRITE_TEXTFILE_EXC_2GEV ) { break; }

          GlobalPoint stubPosition = theStackedGeometry->findGlobalPosition(&(**lItTrkAccStubs)); // Grab the position of the accepted stub

          double lStubRho(stubPosition.mag()), lStubEta(stubPosition.eta()), lStubPhi (stubPosition.phi()), lStubZ(stubPosition.z()); // Stub position info


          DetId lStubDetId ((*lItTrkAccStubs)->getDetId());
          uint32_t lStubId (lStubDetId.rawId());

          if (TEXTFILE_OFFSET && lPrtId == abs(PARTICLE_ID))                                        // If recorded stubs are to be offset, they are offset within this loop
          {

            edm::Service<edm::RandomNumberGenerator> rng;
            CLHEP::HepRandomEngine& engine = rng->getEngine();

            double randomNumber1 = (int((CLHEP::RandFlat::shoot(&engine, 0.0, 100.0)))%2);  // Generate random number to randomly decide which stub to offset

//             std::cout << "randomNumber1: " << randomNumber1 << std::endl;

            if ( lNumerOffsetStubs < NUM_STUBS_OFFSET && randomNumber1 == 1 )              // If haven't offset max number of stubs yet and is stub to be offset
            {
              double randomNumber2 = CLHEP::RandFlat::shoot(&engine, -100.0, 100.0);  // Generate random number to determine the severity of the offset
              double randomNumber3 = CLHEP::RandFlat::shoot(&engine, -100.0, 100.0);  // Generate random number to determine the severity of the offset
              double randomNumber4 = CLHEP::RandFlat::shoot(&engine, -100.0, 100.0);  // Generate random number to determine the severity of the offset
//             std::cout << "randomNumber2: " << randomNumber2 << std::endl;
//             std::cout << "randomNumber2: " << randomNumber3 << std::endl;
//             std::cout << "randomNumber2: " << randomNumber4 << std::endl;

              ++lNumerOffsetStubs;                                    // Increment counter for registering how many stubs have been offset so far
//               std::cout << "Stub Offset" << std::endl;
//               std::cout << lRandomBinaryInt << ": random binary discriminant used for random offseting." << std::endl;
//               std::cout << "Number of OffSettedStubs/NUM_STUBS_OFFSET: " << lNumerOffsetStubs << "/" << NUM_STUBS_OFFSET << std::endl;

              if ( lStubZ < 1300.0 )                                  // Offset in Barrel
              {
                lStubEta += cos(randomNumber2*MAX_OFFSET_AMOUNT*0.01) + sin(randomNumber3*MAX_OFFSET_AMOUNT*0.01);
              }
              else if ( lStubZ >= 1300.0 )                                  // Offset in Endcaps
              {
                lStubRho += MAX_OFFSET_AMOUNT*0.0001*randomNumber4;
              }
            }
//             else { std::cout << "Maximum number of stubs to be offset reached" << std::endl;}
          }
//           std::cout << "Particle ID: " <<  lPrtId << std::endl;
//           std::cout << "lFlagFirstStubTextFileLoop: " << lFlagFirstStubTextFileLoop << std::endl;

          if ( lPrtId == abs(PARTICLE_ID) )
          {
//           std::cout << "abs(PARTICLE_ID): " << abs(PARTICLE_ID) << std::endl;
            if (lFlagFirstStubTextFileLoop == 0)
            {
             std::cout << "Particle ID (Write to file 1st time): " <<  lPrtId << std::endl;
              lStr << std::hex << lStubId; // Insert Stub Id as hex to lStr
              myFile << " T " << lPrtPt << " " << lPrtEta << " " << lPrtPhi << " A " << lStr.str() << " " << lStubRho << " " << lStubEta << " " << lStubPhi;  // Write info of TrkPrt that created Acc Gen Stub, and first stub pos to file
              ++lFlagFirstStubTextFileLoop;
            }
            else
            {
              std::cout << "Particle ID (Write to file subsequent hits): " <<  lPrtId << std::endl;
              lStr << std::hex << lStubId; // Insert Stub Id as hex to lStr
              myFile << " A " << lStr.str() << " " << lStubRho << " " << lStubEta << " " << lStubPhi;
            }
          }
        }
      }
    }
    if ( lFlagFirstStubTextFileLoop != 0 )
    {
      myFile << "\n"; // After all stubs have been looped over for this particle, new line
    }
  }

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
// Pixel Digis loops - prints the eta/phi dist of all the pixel digi hits
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------

  if (PIXEL_DIGI_LOOP)
  {
    edm::DetSetVector<PixelDigi>::const_iterator module=theDigis->begin();
    for ( ; module!=theDigis->end(); module++ ) 
    {
      DetId tkId(module->id);
      if(tkId.subdetId()==1)
      {
        PXBDetId pxbId(tkId);
        if(pxbId.layer()==10)
        {
         edm::DetSet<PixelDigi>::const_iterator digi=module->data.begin();
         for ( ; digi!=module->data.end(); digi++ )
         {
           const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*>(GeometryHandle->idToDet(module->id));
           const PixelTopology* pixelTopol = dynamic_cast<const PixelTopology*>( &(pixelDet->specificTopology()));

           MeasurementPoint mp0( digi->row()+0.5, digi->column()+0.5 );
           GlobalPoint pos =  (GeometryHandle->idToDet(module->id))->surface().toGlobal(pixelTopol->localPosition(mp0));

           h_etaphi->Fill(pos.eta(),pos.phi());

           h_etaphi->GetXaxis()->SetTitle("#eta");
           h_etaphi->GetYaxis()->SetTitle("#phi");
           h_etaphi->GetXaxis()->CenterTitle();
           h_etaphi->GetYaxis()->CenterTitle();
         }
        }
      }
    }
  }

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
// Tracking Particles Loop
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------



//   double pTHiggsMuons[4] = {0.0};     // For Higgs Analysis, pT array for the daughter muons
//   double pTHiggsElectrons[4] = {0.0}; // For Higgs Analysis, pT array for the daughter electrons

//   double pTHiggsMuonsAccGenStubs[4] = {0.0};     // For Higgs Analysis, pT array for the daughter muons that create a genuine accepted stub
//   double pTHiggsElectronsAccGenStubs[4] = {0.0}; // For Higgs Analysis, pT array for the daughter electrons that create a genuine accepted stub

  int lFlagFirstAccStubOnly(0);     // Counter 
  int lFlagFirstAccGenStubOnly(0);  
  int lFlagFirstAccComStubOnly(0);  
  int lFlagFirstAccUnkStubOnly(0);  
  int lFlagFirstRejStubOnly(0);     
  int lFlagFirstStubRejGenOnly(0);  
  int lFlagFirstStubRejComOnly(0);  
  int lFlagFirstStubRejUnkOnly(0);  

  int lMuonCounterEvent(0);         // Counter for how many muons are present in the event
  int lElectronCounterEvent(0);     // Counter for how many electrons are present in the event
  int lAntiMuonCounterEvent(0);     // Counter for how many anti-muons are present in the event
  int lAntiElectronCounterEvent(0); // Counter for how many positrons are present in the event

  int lTempPrtCount(0);  // Counter to be incremented after each particle in the event has looped over, to correct the next one.

// --------------------------------------------------------------------------------------------------------------------------------------------------
  // Loop over all Tracking Particles in file
// --------------------------------------------------------------------------------------------------------------------------------------------------

  if ( TRACKING_PARTICLE_LOOP )
  {
    for ( std::vector<TrackingParticle>::const_iterator lTrkPrtIt = TruthTrackContainer->begin(); lTrkPrtIt != TruthTrackContainer->end(); lTrkPrtIt++)
    {
      const double PrtPt (lTrkPrtIt->pt()); // Grab TrkPrt pT
      const int PrtId (lTrkPrtIt->pdgId()); // Grab TrkPt Id

      h_TrkPrtPt->Fill(PrtPt);
      h_TrkPrtPt->GetXaxis()->SetTitle("Tracking Particle p_{T}");
      h_TrkPrtPt->GetXaxis()->CenterTitle();

      if (lTempPrtCount == 0)
      {
        h_PrimaryParticlePt->Fill(PrtPt);
      }

  //     std::cout << "PiD: " << lTrkPrtIt->pdgId() << " TrackParticlePt: " << PrtPt << "\t Charge: " << (lTrkPrtIt->threeCharge())/3 << "\t Eta: " << lTrkPrtIt->eta() << "\t Phi: " << lTrkPrtIt->phi() << std::endl;
      if ( PrtId == 13 )                          // If a muon
      {
        ++MuonPrtCounter;                         // Add one to the Muon Counter
        ++lMuonCounterEvent;                      // Add one to the Muon Counter (used for this event only)
//        std::cout << "Muon Pt: " << PrtPt << std::endl;
      }
      else if ( PrtId == -13 )                    // If an anti-muon
      {
        ++lAntiMuonCounterEvent;
//        std::cout << "Anti-muon Pt: " << PrtPt << std::endl;
      }
      else if ( PrtId == 11 )                     // If an electron
      {
        ++ElectronPrtCounter;
        ++lElectronCounterEvent;
//        std::cout << "Electron Pt: " << PrtPt << std::endl;
      }
      else if ( PrtId == -11 )                     // If a positron
      {
        ++lAntiElectronCounterEvent;
//         std::cout << "Electron Pt: " << PrtPt << std::endl;
      }
      else if ( PrtId == 15 || PrtId == -15 )                     // If a tauon or anti-tauon
      {
        ++TauPrtCounter;
      }
      else if ( PrtId == 23 )                                     // If a Z
      {
        ++ZPrtCounter;
      }
      else if ( PrtId == 221 || PrtId == -221 )                   // If Pi(+/-)
      {
        ++ChargedPionPrtCounter;
//         std::cout << "Charged Pion Pt: " << PrtPt << std::endl;
      }
      else if ( PrtId == 111  )                                   // If Pi(0)
      {
        ++NeutralPionPrtCounter;
      }
      else if ( PrtId == 22  )                                    // If a Photon
      {
        ++PhotonPrtCounter;
      }
      else if ( PrtId == 25  )                                    // If a Higgs
      {
        ++HiggsPrtCounter;
      }
      else                                                        // If a different particle
      {
        ++OtherPrtCounter;
      }

      ++lPrtCounter;

      edm::Ptr<TrackingParticle> tempTPPtr(TruthTrackContainer, lTempPrtCount++);  // Create temporary TrkPrt Pointer, first arguement states we are using the TrkPrt Truth Container, second arguement states which particle in the container in the current event we are looking at.
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > > theStubs = MCTruthTTStubHandle->findTTStubRefs(tempTPPtr);
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > > theRejStubs = MCTruthTTStubRejectedHandle->findTTStubRefs(tempTPPtr);

//         std::cout << "Number of stubs: " << theStubs.size() << std::endl;

        if ( theStubs.size() >= 5 ) 
        {
          h_FivePlusStubsPrts->Fill(PrtPt);
        }
        else if ( theStubs.size() < 5 )
        {
          h_LessThanFiveStubsPrts->Fill(PrtPt);
        }

        if ( theStubs.size() >= 5 && lTempPrtCount == 1)
        {
          h_NumTrkPrtsFivePlusAccStubs->Fill(PrtPt);
        }

        else if ( theRejStubs.size() >= 5 && lTempPrtCount == 1)
        {
          h_NumTrkPrtsFivePlusRejStubs->Fill(PrtPt);
        }

        h_FrequencyOfNumStubs->Fill(theStubs.size());

// --------------------------------------------------------------------------------------------------------------------------------------------------
    // Accepted stubs associated with Tracking Particle
// --------------------------------------------------------------------------------------------------------------------------------------------------
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > >::const_iterator lItTrkAccStubs = theStubs.begin();
      for ( ; lItTrkAccStubs != theStubs.end(); ++lItTrkAccStubs )
      {

        if ( MCTruthTTStubHandle->isGenuine(*lItTrkAccStubs) )
        {
//             if (lTempPrtCount == 0)    // If The first tracking particle in the event (ie. the original particle) - non-pileup
//             {
              if (lFlagFirstAccStubOnly == 0)  // If the first accepted stub for a Trk Prt, fill up the relevant histograms
              {
                h_TrkPrtAccAllPt->Fill(PrtPt);
                h_TrkPrtAccAllPt->GetXaxis()->SetTitle("Accepted stub creating Tracking Particle p_{T}");
                h_TrkPrtAccAllPt->GetXaxis()->CenterTitle(); ++lFlagFirstAccStubOnly;
              }
              h_TrkPrtAccGenPt->Fill(PrtPt);
              h_TrkPrtAccGenPt->GetXaxis()->SetTitle("Genuine Accepted stub creating Tracking Initial Particle p_{T}");
              h_TrkPrtAccGenPt->GetXaxis()->CenterTitle();

              if ( !EXCLUDE_INTERACTING_PRTS )
              {
                ++lFlagFirstAccGenStubOnly; // If Excluding Interacting Particles in non-pileup scenerio, set flag so that only the first tracking particle that creates a stub is read.
              }
  
            h_TrkPrtAccAllStubsPt->Fill(PrtPt);
            h_TrkPrtAccGenStubsPt->Fill(PrtPt);
  
//             }
  
//           std::cout << "Accepted Stub from Tracking Particle is Genuine" << std::endl;
        }
        if ( MCTruthTTStubHandle->isCombinatoric(*lItTrkAccStubs) )
        {
/*            if (lTempPrtCount == 0)
            {*/
              if (lFlagFirstAccStubOnly ==0)
              {
                h_TrkPrtAccAllPt->Fill(PrtPt);
                h_TrkPrtAccAllPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
                h_TrkPrtAccAllPt->GetXaxis()->CenterTitle(); ++lFlagFirstAccStubOnly;
              }
              h_TrkPrtAccComPt->Fill(PrtPt);
              h_TrkPrtAccComPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
              h_TrkPrtAccComPt->GetXaxis()->CenterTitle();
              if ( !EXCLUDE_INTERACTING_PRTS )
              {
                ++lFlagFirstAccComStubOnly; // If Excluding Interacting Particles in non-pileup scenerio, set flag so that only the first tracking particle that creates a stub is read.
              }
//             }
            h_TrkPrtAccAllStubsPt->Fill(PrtPt);
            h_TrkPrtAccComStubsPt->Fill(PrtPt);
  //         std::cout << "Accepted Stub from Tracking Particle is Combinatoric" << std::endl;
        }
        if ( MCTruthTTStubHandle->isUnknown(*lItTrkAccStubs) )
        {
  //         if ( lTrkPrtIt->pdgId() == 13 && MUONS_ONLY )
  //         {
/*            if (lTempPrtCount == 0)
            {*/
              if (lFlagFirstAccStubOnly ==0)
              {
                h_TrkPrtAccAllPt->Fill(PrtPt);
                h_TrkPrtAccAllPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
                h_TrkPrtAccAllPt->GetXaxis()->CenterTitle(); ++lFlagFirstAccStubOnly;
              }
              h_TrkPrtAccUnkPt->Fill(PrtPt);
              h_TrkPrtAccUnkPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
              h_TrkPrtAccUnkPt->GetXaxis()->CenterTitle();
              if ( !EXCLUDE_INTERACTING_PRTS )
              {
                ++lFlagFirstAccUnkStubOnly; // If Excluding Interacting Particles in non-pileup scenerio, set flag so that only the first tracking particle that creates a stub is read.
              }
//             }
            h_TrkPrtAccAllStubsPt->Fill(PrtPt);
            h_TrkPrtAccUnkStubsPt->Fill(PrtPt);

        }
  //         std::cout << "Accepted Stub from Tracking Particle is Unknown" << std::endl;
      }
      // Rejected stubs associated with Tracking Particle

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > >::const_iterator lItTrkRejStubs = theRejStubs.begin();
      for ( ; lItTrkRejStubs != theRejStubs.end(); ++lItTrkRejStubs )
      {
        if ( MCTruthTTStubRejectedHandle->isGenuine(*lItTrkRejStubs) )
        {
/*            if (lTempPrtCount == 0)
            {*/
              if (lFlagFirstRejStubOnly ==0)
              {
                h_TrkPrtRejAllPt->Fill(PrtPt);
                h_TrkPrtRejAllPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
                h_TrkPrtRejAllPt->GetXaxis()->CenterTitle(); ++lFlagFirstRejStubOnly;
              }
              h_TrkPrtRejGenPt->Fill(PrtPt);
              h_TrkPrtRejGenPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
              h_TrkPrtRejGenPt->GetXaxis()->CenterTitle();
              if ( !EXCLUDE_INTERACTING_PRTS )
              {
                ++lFlagFirstStubRejGenOnly; // If Excluding Interacting Particles in non-pileup scenerio, set flag so that only the first tracking particle that creates a stub is read.
              }
//             }
            h_TrkPrtRejAllStubsPt->Fill(PrtPt);
            h_TrkPrtRejGenStubsPt->Fill(PrtPt);
  //         std::cout << "Rejected Stub from Tracking Particle is Genuine" << std::endl;
        }
        if ( MCTruthTTStubRejectedHandle->isCombinatoric(*lItTrkRejStubs) )
        {
/*            if (lTempPrtCount == 0)
            {*/
              if (lFlagFirstRejStubOnly ==0)
              {
                h_TrkPrtRejAllPt->Fill(PrtPt);
                h_TrkPrtRejAllPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
                h_TrkPrtRejAllPt->GetXaxis()->CenterTitle(); ++lFlagFirstRejStubOnly;
              }
              h_TrkPrtRejComPt->Fill(PrtPt);
              h_TrkPrtRejComPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
              h_TrkPrtRejComPt->GetXaxis()->CenterTitle();
              if ( !EXCLUDE_INTERACTING_PRTS )
              {
                ++lFlagFirstStubRejComOnly; // If Excluding Interacting Particles in non-pileup scenerio, set flag so that only the first tracking particle that creates a stub is read.
              }
//             }
            h_TrkPrtRejAllStubsPt->Fill(PrtPt);
            h_TrkPrtRejComStubsPt->Fill(PrtPt);
  //         std::cout << "Rejected Stub from Tracking Particle is Combinatoric" << std::endl;
        }
        if ( MCTruthTTStubRejectedHandle->isUnknown(*lItTrkRejStubs) )
        {
/*            if (lTempPrtCount == 0)
            {*/
              if (lFlagFirstRejStubOnly ==0)
              {
                h_TrkPrtRejAllPt->Fill(PrtPt);
                h_TrkPrtRejAllPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
                h_TrkPrtRejAllPt->GetXaxis()->CenterTitle(); ++lFlagFirstRejStubOnly;
              }
              h_TrkPrtRejUnkPt->Fill(PrtPt);
              h_TrkPrtRejUnkPt->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
              h_TrkPrtRejUnkPt->GetXaxis()->CenterTitle();
              if ( !EXCLUDE_INTERACTING_PRTS )
              {
                ++lFlagFirstStubRejUnkOnly; // If Excluding Interacting Particles in non-pileup scenerio, set flag so that only the first tracking particle that creates a stub is read.
              }
//             }
            h_TrkPrtRejAllStubsPt->Fill(PrtPt);
            h_TrkPrtRejUnkStubsPt->Fill(PrtPt);
          }

  //         std::cout << "Rejected Stub from Tracking Particle is Unknown" << std::endl;
      }
    }

  }
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
// Stub Loops
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
  if ( STUB_LOOPS )
  {
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
  // Accepted Stubs
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------

    int lFlagFirstStubOnly(0);
  
    for ( edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator lIterDSV=TTStubHandle->begin(); lIterDSV!=TTStubHandle->end(); lIterDSV++ )
    {
      DetId thisStackedDetId = lIterDSV->id();
      edmNew::DetSet< TTStub< Ref_PixelDigi_ > > theStubs = (*TTStubHandle)[ thisStackedDetId ];
  
      for ( edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator lIterTTStub = theStubs.begin(); lIterTTStub != theStubs.end(); lIterTTStub++)
      {
        StackedTrackerDetId stubDetId = lIterTTStub->getDetId();
  // 
        uint iLayer = stubDetId.iLayer();   // Returns which layer in the Barrel the stub was formed in - returns 999999 if endcap
  //       uint iRing = stubDetId.iRing();     // Returns which ring in the Endcap the stub was formed in  - returns 999999 if barrel
        uint iDisk = stubDetId.iDisk();     // Returns which disk in the Endcap the stub was formed in  - returns 999999 if barrel
  //       uint iPhi = stubDetId.iPhi();
  //       uint iZ = stubDetId.iZ();
  
  //       std::cout << "iLayer = " << iLayer << std::endl;
  
  
        GlobalPoint stubPosition = theStackedGeometry->findGlobalPosition(&(*lIterTTStub)); // Grab the position of the accepted stub
  
  //       std::cout << "StubPos (r,eta,phi): " << stubPosition.perp() << "," << stubPosition.eta() << "," << stubPosition.phi() << std::endl;
  
        edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > theRef = edmNew::makeRefTo(TTStubHandle, lIterTTStub );
  
  //       std::cout<<"STUB :: "<<lIterDSV->id()<<" : " <<MCTruthTTStubHandle->isGenuine(theRef)<<" : " <<MCTruthTTStubHandle->isCombinatoric(theRef)<<" : "
  //                <<MCTruthTTStubHandle->isUnknown(theRef) <<std::endl;
  
        edm::Ptr< TrackingParticle > lParticlePtr (MCTruthTTStubHandle->findTrackingParticlePtr(theRef)); // Get tracking particle associated with stub (only for acc gen stubs)
  
        if (MCTruthTTStubHandle->isGenuine(theRef) == 1)  // If a genuine stub
        {
          const double lTrueAccPt = lParticlePtr->pt();   // Grab associated TrkPt pT
  
          if (lFlagFirstStubOnly == 0 ) // If first stub the event has created ...
          {
            ++StubCreatingPrtCounter;
            ++GenAccStubCreatingPrtCounter;
            h_pTdist->Fill(lTrueAccPt);                                                                   // Fill histogram of Genuine Stub pT
            h_pTdist->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
            h_pTdist->GetXaxis()->CenterTitle();
            h_pTdistAccGenStub->Fill(lTrueAccPt);                                                         // Fill histogram of Acc Gen Stub pT
            h_pTdistAccGenStub->GetXaxis()->SetTitle("Genuine Accepted Tracking Initial Particle p_{T}");
            h_pTdistAccGenStub->GetXaxis()->CenterTitle();
            ++lFlagFirstStubOnly; // Set flag to 1, to prevent refilling histograms with same Trk Prt and to output data to file in correct format
          }
  
          if ( iLayer != 999999 ) // If in the Barrel
          {
            h_pTdistAccGenStub_layer[iLayer-1]->Fill(lTrueAccPt);
            h_pTdistAccGenStub_layer[iLayer-1]->GetXaxis()->SetTitle("Accepted Tracking Initial Particle p_{T}");
            h_pTdistAccGenStub_layer[iLayer-1]->GetXaxis()->CenterTitle();
          }
  
          if ( iDisk != 999999 )  // If in one of the Endcaps
          {
            h_pTdistAccGenStub_disk[iDisk-1]->Fill(lTrueAccPt);
            h_pTdistAccGenStub_disk[iDisk-1]->GetXaxis()->SetTitle("Accepted Tracking Initial Particle p_{T}");
            h_pTdistAccGenStub_disk[iDisk-1]->GetXaxis()->CenterTitle();
          }
  
          h_stub_phi->Fill(stubPosition.phi(), lTrueAccPt);                 // Fill histogram with phi of Gen Acc stub and pT of TrkPrt that created the stub
          h_stub_phi->GetXaxis()->SetTitle("#phi");
          h_stub_phi->GetYaxis()->SetTitle("stub p_{T}");
          h_stub_phi->GetXaxis()->CenterTitle();
          h_stub_phi->GetYaxis()->CenterTitle();
  
          h_stub_eta->Fill(stubPosition.eta(), lTrueAccPt);
          h_stub_eta->GetXaxis()->SetTitle("#eta");
          h_stub_eta->GetYaxis()->SetTitle("stub p_{T}");
          h_stub_eta->GetXaxis()->CenterTitle();
          h_stub_eta->GetYaxis()->CenterTitle();
  
        h_stub_etaphi->Fill(stubPosition.eta(),stubPosition.phi());         // Fill histogram with eta of Gen Acc stub and pT of TrkPrt that created the stub
          h_stub_etaphi->GetXaxis()->SetTitle("#eta");
          h_stub_etaphi->GetYaxis()->SetTitle("#phi");
          h_stub_etaphi->GetXaxis()->CenterTitle();
          h_stub_etaphi->GetYaxis()->CenterTitle();
  
          h_stub_xy->Fill(stubPosition.x()*10.0, stubPosition.y()*10.0);    // Fill histogram with x and y pos of Gen Acc stub
          h_stub_xy->GetXaxis()->SetTitle("x [mm]");
          h_stub_xy->GetYaxis()->SetTitle("y [mm]");
          h_stub_xy->GetXaxis()->CenterTitle();
          h_stub_xy->GetYaxis()->CenterTitle();
  
          h_stub_zr->Fill(stubPosition.z()*10.0, stubPosition.perp()*10.0); // Fill histogram with z and perp pos of Gen Acc stub
          h_stub_zr->GetXaxis()->SetTitle("z [mm]");
          h_stub_zr->GetYaxis()->SetTitle("r [mm]");
          h_stub_zr->GetXaxis()->CenterTitle();
          h_stub_zr->GetYaxis()->CenterTitle();
  
          // TH1F's - distributions
  
          h_stub_eta_distb->Fill(stubPosition.eta());                       // Fill histo with Eta distribution of Gen Acc Stubs
          h_stub_eta_distb->GetXaxis()->SetTitle("#eta");
          h_stub_eta_distb->GetXaxis()->CenterTitle();
  
          h_stub_phi_distb->Fill(stubPosition.phi());                       // Fill histo with Phi distribution of Gen Acc Stubs
          h_stub_phi_distb->GetXaxis()->SetTitle("phi");
          h_stub_phi_distb->GetXaxis()->CenterTitle();
  
          if ( (stubPosition.z()) <= double(130.0) )                        // If in the barrel region ...
          {
            h_stub_perp_barrel_distb->Fill(stubPosition.perp()*10.0);       // Fill histo with perp distribution of Gen Acc Stubs in the barrel
            h_stub_perp_barrel_distb->GetXaxis()->SetTitle("r [mm]");
            h_stub_perp_barrel_distb->GetXaxis()->CenterTitle();
          }
          else if ( (stubPosition.z()) > double(130.0) )                    // If in the endcap region ...
          {
            h_stub_perp_endcap_distb->Fill(stubPosition.perp()*10.0);       // Fill histo with perp distribution of Gen Acc Stubs in the endcaps
            h_stub_perp_endcap_distb->GetXaxis()->SetTitle("r [mm]");
            h_stub_perp_endcap_distb->GetXaxis()->CenterTitle();
          }
  
          h_stub_pT_distb->Fill(lTrueAccPt);                                // Fill histo with pT distribution of Acc Gen Stubs
          h_stub_pT_distb->GetXaxis()->SetTitle("stub p_{T}");
          h_stub_pT_distb->GetXaxis()->CenterTitle();
  
          if ( iLayer != 999999 )                                                   // If in the barrel ...
          {
            h_stub_pT_layer_distb[iLayer-1]->Fill(lTrueAccPt);                      // Fill histo with pT Gen Acc Stubs in said layer
            h_stub_pT_layer_distb[iLayer-1]->GetXaxis()->SetTitle("Stub p_{T}");
            h_stub_pT_layer_distb[iLayer-1]->GetXaxis()->CenterTitle();
          }
  
          if ( iDisk != 999999 )                                                    // If in the endcaps ...
          {
            h_stub_pT_disk_distb[iDisk-1]->Fill(lTrueAccPt);
            h_stub_pT_disk_distb[iDisk-1]->GetXaxis()->SetTitle("Stub p_{T}");      // Fill histo with pT Gen Acc Stubs in said layer
            h_stub_pT_disk_distb[iDisk-1]->GetXaxis()->CenterTitle();
          }
  
  //       std::cout << "Accepted Stub is Genuine" << std::endl;
          lCounterGenStub++;
        }
  
  
        else if (MCTruthTTStubHandle->isCombinatoric(theRef) == 1 )
        {
  //         std::cout << "Accepted Stub is Combinatoric" << std::endl;
  	lCounterCombStub++;
        }
  
        else if (MCTruthTTStubHandle->isUnknown(theRef) == 1 )
        {
  //         std::cout << "Accepted Stub is Unknown" << std::endl;
  	lCounterUnkStub++;
        }
  
      }
    }

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
// Rejected Stubs loops
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------

    uint lRejFlagFirstStubOnly (0);
  
    for ( edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator lIterDSV=TTStubRejectedHandle->begin(); lIterDSV!=TTStubRejectedHandle->end(); lIterDSV++ )
    {
      DetId thisStackedDetId = lIterDSV->id();
      edmNew::DetSet< TTStub< Ref_PixelDigi_ > > theRejStubs = (*TTStubRejectedHandle)[ thisStackedDetId ];
  
      for ( edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator lIterRejTTStub = theRejStubs.begin(); lIterRejTTStub != theRejStubs.end(); lIterRejTTStub++)
      {
  // 
        StackedTrackerDetId RejStubDetId = lIterRejTTStub->getDetId();
  // 
        uint iLayer = RejStubDetId.iLayer();   // Barrel - returns 999999 if endcap
        uint iDisk = RejStubDetId.iDisk();     // Endcap - returns 999999 if barrel
  
/*        double RejRoughStubPt = theStackedGeometry->findRoughPt(3.80,&(*lIterRejTTStub));    // Can use aprox of B field = 3.8T due to homogenity of field in Barrel
        if (RejRoughStubPt>10000.0) {RejRoughStubPt=9999.99;}*/
  //       std::cout << "RejectedRoughStubPt: " << RejRoughStubPt << std::endl;
        GlobalPoint RejStubPosition = theStackedGeometry->findGlobalPosition(&(*lIterRejTTStub));
  
        edm::Ref< edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > theRejRef = edmNew::makeRefTo(TTStubRejectedHandle, lIterRejTTStub );
  
        edm::Ptr< TrackingParticle > lParticlePtr (MCTruthTTStubRejectedHandle->findTrackingParticlePtr(theRejRef));
  
        if (MCTruthTTStubRejectedHandle->isGenuine(theRejRef) == 1)
        {
          double lTrueRejPt = lParticlePtr->pt();
  
          if (lRejFlagFirstStubOnly == 0 )
          {
            h_pTdist->Fill(lTrueRejPt);
            h_pTdist->GetXaxis()->SetTitle("Tracking Initial Particle p_{T}");
            h_pTdist->GetXaxis()->CenterTitle();
            h_pTdistRejGenStub->Fill(lTrueRejPt);
            h_pTdistRejGenStub->GetXaxis()->SetTitle("Rejected Tracking Initial Particle p_{T}");
            h_pTdistRejGenStub->GetXaxis()->CenterTitle();
  
          if (lFlagFirstStubOnly == 0)
            {
              ++StubCreatingPrtCounter;
            }
          lRejFlagFirstStubOnly++;
          }
  
  
          if ( iLayer != 999999 )
          {
            h_pTdistRejGenStub_layer[iLayer-1]->Fill(lTrueRejPt);
            h_pTdistRejGenStub_layer[iLayer-1]->GetXaxis()->SetTitle("Rejected Tracking Initial Particle p_{T}");
            h_pTdistRejGenStub_layer[iLayer-1]->GetXaxis()->CenterTitle();
          }
  
          if ( iDisk != 999999 )
          {
            h_pTdistRejGenStub_disk[iDisk-1]->Fill(lTrueRejPt);
            h_pTdistRejGenStub_disk[iDisk-1]->GetXaxis()->SetTitle("Rejected Tracking Initial Particle p_{T}");
            h_pTdistRejGenStub_disk[iDisk-1]->GetXaxis()->CenterTitle();
          }
  
          h_rejstub_phi->Fill(RejStubPosition.phi(), lTrueRejPt);
          h_rejstub_phi->GetXaxis()->SetTitle("#phi");
          h_rejstub_phi->GetYaxis()->SetTitle("stub p_{T}");
          h_rejstub_phi->GetXaxis()->CenterTitle();
          h_rejstub_phi->GetYaxis()->CenterTitle();
  
          h_rejstub_eta->Fill(RejStubPosition.eta(), lTrueRejPt);
          h_rejstub_eta->GetXaxis()->SetTitle("#eta");
          h_rejstub_eta->GetYaxis()->SetTitle("stub p_{T}");
          h_rejstub_eta->GetXaxis()->CenterTitle();
          h_rejstub_eta->GetYaxis()->CenterTitle();
  
          h_rejstub_etaphi->Fill(RejStubPosition.eta(),RejStubPosition.phi());
          h_rejstub_etaphi->GetXaxis()->SetTitle("#eta");
          h_rejstub_etaphi->GetYaxis()->SetTitle("#phi");
          h_rejstub_etaphi->GetXaxis()->CenterTitle();
          h_rejstub_etaphi->GetYaxis()->CenterTitle();
  
          h_rejstub_xy->Fill(RejStubPosition.x()*10.0, RejStubPosition.y()*10.0);
          h_rejstub_xy->GetXaxis()->SetTitle("x [mm]");
          h_rejstub_xy->GetYaxis()->SetTitle("y [mm]");
          h_rejstub_xy->GetXaxis()->CenterTitle();
          h_rejstub_xy->GetYaxis()->CenterTitle();
  
          h_rejstub_zr->Fill(RejStubPosition.z()*10.0, RejStubPosition.perp()*10.0);
          h_rejstub_zr->GetXaxis()->SetTitle("z [mm]");
          h_rejstub_zr->GetYaxis()->SetTitle("r [mm]");
          h_rejstub_zr->GetXaxis()->CenterTitle();
          h_rejstub_zr->GetYaxis()->CenterTitle();
  
          h_rejstub_eta_distb->Fill(RejStubPosition.eta());
          h_rejstub_eta_distb->GetXaxis()->SetTitle("#eta");
          h_rejstub_eta_distb->GetXaxis()->CenterTitle();
  
          h_rejstub_phi_distb->Fill(RejStubPosition.phi());
          h_rejstub_phi_distb->GetXaxis()->SetTitle("phi");
          h_rejstub_phi_distb->GetXaxis()->CenterTitle();
  
          if ( (RejStubPosition.z()) <= double(125.0) )
          {
            h_rejstub_perp_barrel_distb->Fill(RejStubPosition.perp()*10.0);
            h_rejstub_perp_barrel_distb->GetXaxis()->SetTitle("r [mm]");
            h_rejstub_perp_barrel_distb->GetXaxis()->CenterTitle();
          }
          else if ( (RejStubPosition.z()) > double(125.0) )
  
          {
            h_rejstub_perp_endcap_distb->Fill(RejStubPosition.perp()*10.0);
            h_rejstub_perp_endcap_distb->GetXaxis()->SetTitle("r [mm]");
            h_rejstub_perp_endcap_distb->GetXaxis()->CenterTitle();
          }
  
          h_rejstub_pT_distb->Fill(lTrueRejPt);
          h_rejstub_pT_distb->GetXaxis()->SetTitle("stub p_{T}");
          h_rejstub_pT_distb->GetXaxis()->CenterTitle();
  
          if ( iLayer != 999999 )
          {
            h_rejstub_pT_layer_distb[iLayer-1]->Fill(lTrueRejPt);
            h_rejstub_pT_layer_distb[iLayer-1]->GetXaxis()->SetTitle("Stub p_{T}");
            h_rejstub_pT_layer_distb[iLayer-1]->GetXaxis()->CenterTitle();
          }
  
          if ( iDisk != 999999 )
          {
            h_rejstub_pT_disk_distb[iDisk-1]->Fill(lTrueRejPt);
            h_rejstub_pT_disk_distb[iDisk-1]->GetXaxis()->SetTitle("Stub p_{T}");
            h_rejstub_pT_disk_distb[iDisk-1]->GetXaxis()->CenterTitle();
          }
  
  //         std::cout << "Rejected Stub is Genuine" << std::endl;
  //         CounterRejGenPrt++;//Mariza's Gen Rej Stub Counter
        }
  
        else if (MCTruthTTStubRejectedHandle->isCombinatoric(theRejRef) == 1 )
        {
  //         std::cout << "Rejected Stub is Combinatoric" << std::endl;
  // 	CounterRejNonGenPrt++;// Mariza's Non Gen Rej Stub Counter
        }
  
        else if (MCTruthTTStubRejectedHandle->isUnknown(theRejRef) == 1 )
        {
  //         std::cout << "Rejected Stub is Unknown" << std::endl;
  // 	CounterRejNonGenPrt++;// Mariza's Non Gen Rej Stub Counter
        }
       }
      }
    }

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
}
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------



// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
  if ( WRITE_TEXTFILE )
  {
    myFile.open(lFileName.c_str());
    std::cout << "myFile.open ..." << std::endl;
  }
  edm::Service<TFileService> fs;

  h_PrimaryParticlePt = fs->make<TH1F>("h_PrimaryParticlePt", "pT Dist of Primary Tracking Particles in event", 600, 0.0, 600.0);
  h_FivePlusStubsPrts = fs->make<TH1F>("h_FivePlusStubsPrts", "pT Distribution of Tracking Particles that create 5 accepted stubs or more", 600, 0.0, 600.0);
  h_LessThanFiveStubsPrts = fs->make<TH1F>("h_LessThanFiveStubsPrts", "pT Distribution of Tracking Particles that create less than 5 accepted stubs", 600, 0.0, 600.0);
  h_FivePlusStubsFraction = fs->make<TH1F>("h_FivePlusStubsFraction", "Proportion of Tracking Particles that create five of more stubs as a function of pT", 600, 0.0, 600.0);
  h_FrequencyOfNumStubs = fs->make<TH1F>("h_FrequencyOfNumStubs", "Frequency of number of accepted stubs", 20, 0.0, 20.0);
  h_numberStubs = fs->make<TH1F>("h_numberStubs", "Frequency of the types of accepted stubs produced", 3, 0.5, 3.5);
  h_numberStubsFrac = fs->make<TH1F>("h_numberStubsFrac", "Fraction of the types of accepted stubs produced", 3, 0.5, 3.5);
  h_TrkPrtPt_Frac = fs->make<TH1F>("h_TrkPrtPt_Frac", "Fraction of pT Distribution of all the Tracking Particles", 600, 0.0, 600.0);

  h_TrkPrtPt_Frac->GetXaxis()->SetTitle("Tracking Particle p_{T}");


  h_numberStubs->GetXaxis()->SetBinLabel( 1 , "Genuine Stubs" );
  h_numberStubs->GetXaxis()->SetBinLabel( 2 , "Combinatorial Stubs" );
  h_numberStubs->GetXaxis()->SetBinLabel( 3 , "Unknown Stubs" );

  h_numberStubsFrac->GetXaxis()->SetBinLabel( 1 , "Genuine Stubs" );
  h_numberStubsFrac->GetXaxis()->SetBinLabel( 2 , "Combinatorial Stubs" );
  h_numberStubsFrac->GetXaxis()->SetBinLabel( 3 , "Unknown Stubs" );

  h_NumTrkPrtsFivePlusAccStubs = fs->make<TH1F>("h_NumTrkPrtsFivePlusAccStubs", "pT Dist of TrkPrts that create at least five acc stubs", 600, 0.0, 600.0);
  h_NumTrkPrtsFivePlusRejStubs = fs->make<TH1F>("h_NumTrkPrtsFivePlusRejStubs", "pT Dist of TrkPrts that create at least five rej stubs", 600, 0.0, 600.0);

  h_etaphi = fs->make<TH2D>("h_etaphi","Pixel Digis #eta vs. #phi",50,-2.5,2.5,50,-3.14,3.14);
  h_pTdist = fs->make<TH1F>("h_pTdist", "pT Distribution of Tracking Particles that create stubs", 600, 0.0, 600.0);
  h_pTdistAccGenStub = fs->make<TH1F>("h_pTdistAccGenStub", "pT Distribution of Accepted Genuine Tracking Particles", 600, 0.0, 600.0);
  h_pTdistRejGenStub = fs->make<TH1F>("h_pTdistRejGenStub", "pT Distribution of Rejected Tracking Particles", 600, 0.0, 600.0);

  h_TrkPrtPt = fs->make<TH1F>("h_TrkPrtPt", "pT Distribution of all the Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtAccAllStubsPt = fs->make<TH1F>("h_TrkPrtAccAllStubsPt", "pT Distribution of Accepted Stubs", 600, 0.0, 600.0);
  h_TrkPrtAccGenStubsPt = fs->make<TH1F>("h_TrkPrtAccGenStubsPt", "pT Distribution of Genuine Accepted Stubs", 600, 0.0, 600.0);
  h_TrkPrtAccComStubsPt = fs->make<TH1F>("h_TrkPrtAccComStubsPt", "pT Distribution of Combinatorial Accepted Stubs", 600, 0.0, 600.0);
  h_TrkPrtAccUnkStubsPt = fs->make<TH1F>("h_TrkPrtAccUnkStubsPt", "pT Distribution of Unknown Accepted Stubs", 600, 0.0, 600.0);
  h_TrkPrtRejAllStubsPt = fs->make<TH1F>("h_TrkPrtRejAllStubsPt", "pT Distribution of Rejected Stubs", 600, 0.0, 600.0);
  h_TrkPrtRejGenStubsPt = fs->make<TH1F>("h_TrkPrtRejGenStubsPt", "pT Distribution of Genuine Rejected Stubs", 600, 0.0, 600.0); 
  h_TrkPrtRejComStubsPt = fs->make<TH1F>("h_TrkPrtRejComStubsPt", "pT Distribution of Combinatorial Rejected Stubs", 600, 0.0, 600.0);
  h_TrkPrtRejUnkStubsPt = fs->make<TH1F>("h_TrkPrtRejUnkStubsPt", "pT Distribution of Unknown Rejected Stubs", 600, 0.0, 600.0);

  h_TrkPrtAccAllPt = fs->make<TH1F>("h_TrkPrtAccAllPt", "pT Distribution of Accepted Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtAccGenPt = fs->make<TH1F>("h_TrkPrtAccGenPt", "pT Distribution of Genuine Accepted Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtAccComPt = fs->make<TH1F>("h_TrkPrtAccComPt", "pT Distribution of Combinatorial Accepted Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtAccUnkPt = fs->make<TH1F>("h_TrkPrtAccUnkPt", "pT Distribution of Unknown Accepted Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtRejAllPt = fs->make<TH1F>("h_TrkPrtRejAllPt", "pT Distribution of Rejected Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtRejGenPt = fs->make<TH1F>("h_TrkPrtRejGenPt", "pT Distribution of Genuine Rejected Tracking Particles", 600, 0.0, 600.0); 
  h_TrkPrtRejComPt = fs->make<TH1F>("h_TrkPrtRejComPt", "pT Distribution of Combinatorial Rejected Tracking Particles", 600, 0.0, 600.0);
  h_TrkPrtRejUnkPt = fs->make<TH1F>("h_TrkPrtRejUnkPt", "pT Distribution of Unknown Rejected Tracking Particles", 600, 0.0, 600.0);

  h_AccStubPrts_TotalTrkPrts_fraction = fs->make<TH1F>("h_AccStubPrts_TotalTrkPrts_fraction", "Fraction of Accepted Tracking Particles/Total Tracking Particles", 600, 0.0, 600.0);
  h_RejStubPrts_TotalTrkPrts_fraction = fs->make<TH1F>("h_RejStubPrts_TotalTrkPrts_fraction", "Fraction of Rejected Tracking Particles/Total Tracking Particles", 600, 0.0, 600.0);
  h_GenAccStubPrts_TotalTrkPrts_fraction = fs->make<TH1F>("h_GenAccStubPrts_TotalTrkPrts_fraction", "Fraction of Genuine Accepted Tracking Particles/Total Tracking Particles", 600, 0.0, 600.0);
  h_GenAccStubPrts_AccStubPrts_fraction = fs->make<TH1F>("h_GenAccStubPrts_AccStubPrts_fraction", "Fraction of Genuine Accepted Tracking Particles/Accepted Tracking Particles", 600, 0.0, 600.0);

  h_AccStubPrts_TotalTrkPrts_fraction->GetXaxis()->SetTitle("Tracking Particle p_{T}");
  h_RejStubPrts_TotalTrkPrts_fraction->GetXaxis()->SetTitle("Tracking Particle p_{T}");
  h_GenAccStubPrts_TotalTrkPrts_fraction->GetXaxis()->SetTitle("Tracking Particle p_{T}");
  h_GenAccStubPrts_AccStubPrts_fraction->GetXaxis()->SetTitle("Tracking Particle p_{T}");


  h_stub_phi = fs->make<TH2D>("h_stub_phi", "Accepted Stub p_{T} vs. #phi", 100, -3.14, 3.14, 1200, 0.0, 600.0);
  h_stub_eta = fs->make<TH2D>("h_stub_eta", "Accepted Stub p_{T} vs. #eta", 100, -2.5, 2.5, 1200, 0.0, 600.0);
  h_stub_etaphi = fs->make<TH2D>("h_stub_etaphi", "Accepted Stub #eta vs. #phi", 50, -2.5, 2.5, 100, -3.14, 3.14);
  h_stub_xy = fs->make<TH2D>("h_stub_xy", "Accepted Stub x vs. y", 2400, -1200.0, 1200.0, 2400, -1200.0, 1200.0);
  h_stub_zr = fs->make<TH2D>("h_stub_zr", "Accepted Stub z vs. r", 5000, 0, 2500.0, 1200, 0, 1200.0);

  h_rejstub_phi = fs->make<TH2D>("h_rejstub_phi", "Rejected Stub p_{T} vs. #phi", 100, -3.14, 3.14, 1200, 0.0, 600.0);
  h_rejstub_eta = fs->make<TH2D>("h_rejstub_eta", "Rejected Stub p_{T} vs. #eta", 100, -2.5, 2.5, 1200, 0.0, 600.0);
  h_rejstub_etaphi = fs->make<TH2D>("h_rejstub_etaphi", "Rejected Stub #eta vs. #phi", 50, -2.5, 2.5, 100, -3.14, 3.14);
  h_rejstub_xy = fs->make<TH2D>("h_rejstub_xy", "Rejected Stub x vs. y", 2400, -1200.0, 1200.0, 2400, -1200.0, 1200.0);
  h_rejstub_zr = fs->make<TH2D>("h_rejstub_zr", "Rejected Stub z vs. r", 5000, 0, 2500.0, 1200, 0.0, 1200.0);

  h_stub_eta_distb = fs->make<TH1F>("h_stub_eta_distb", "Accepted Stub #eta Distribution", 100, -2.50, 2.50);
  h_stub_phi_distb = fs->make<TH1F>("h_stub_phi_distb", "Accepted Stub #phi Distribution", 100, -3.14, 3.14);
  h_stub_perp_barrel_distb = fs->make<TH1F>("h_stub_perp_barrel_distb", "Accepted Stub Radial (perp) Distribution (Barrel)", 1200, 0.0, 1200.0);
  h_stub_perp_endcap_distb = fs->make<TH1F>("h_stub_perp_endcap_distb", "Accepted Stub Radial (perp) Distribution (Endcap)", 1200, 0.0, 1200.0);
  h_stub_pT_distb = fs->make<TH1F>("h_stub_pT_distb", "Accepted Stub p_{T} Distribution", 600, 0.0, 600.0);

  h_rejstub_eta_distb  = fs->make<TH1F>("h_rejstub_eta_distb", "Rejected Stub #eta Distribution", 100, -2.50, 2.50);
  h_rejstub_phi_distb = fs->make<TH1F>("h_rejstub_phi_distb", "Rejected Stub #phi Distribution", 100, -3.14, 3.14);
  h_rejstub_perp_barrel_distb  = fs->make<TH1F>("h_rejstub_perp_barrel_distb", "Rejected Stub Radial (perp) Distribution (Barrel)", 1200, 0.0, 1200.0);
  h_rejstub_perp_endcap_distb  = fs->make<TH1F>("h_rejstub_perp_endcap_distb", "Rejected Stub Radial (perp) Distribution (Endcap)", 1200, 0.0, 1200.0);
  h_rejstub_pT_distb = fs->make<TH1F>("h_rejstub_pT_distb", "Rejected Stub p_{T} Distribution", 600, 0.0, 600.0);

//   h_stubdiff_eta_distb  = fs->make<TH1F>("h_stubdiff_eta_distb", "Difference in Acc'd/Rej'd Stub #eta Distribution", 100, -2.50, 2.50);
//   h_stubdiff_phi_distb = fs->make<TH1F>("h_stubdiff_phi_distb", "Difference in Acc'd/Rej'd Stub #phi Distribution", 100, -3.14, 3.14);
//   h_stubdiff_perp_distb  = fs->make<TH1F>("h_stubdiff_perp_distb", "Difference in Acc'd/Rej'd Stub Radial (perp) Distribution", 1200, 0.0, 1200.0);

  h_stubdiff_pT_distb = fs->make<TH1F>("h_stubdiff_pT_distb", "Fraction of Acc'd/Rej Stubs p_{T} Distribution", 600, 0.0, 600.0);
  h_accTrkPrt_fraction = fs->make<TH1F>("h_accTrkPrt_fraction", "Fraction of Acc'd TrkParticles p_{T} Distribution", 600, 0.0, 600.0);

  for ( uint lIt=0; lIt !=6; ++lIt )
  {
    std::stringstream lSStr;
    lSStr << (lIt+1);
    h_pTdistAccGenStub_layer[lIt] = fs->make<TH1F>(("h_pTdistAccGenStub_layer"+lSStr.str()).c_str(), ("pT Distribution of Accepted Tracking Particles for Layer"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_pTdistRejGenStub_layer[lIt] = fs->make<TH1F>(("h_pTdistRejGenStub_layer"+lSStr.str()).c_str(), ("pT Distribution of Rejected Tracking Particles for Layer"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_stub_pT_layer_distb[lIt] = fs->make<TH1F>(("h_stub_pT_layer"+lSStr.str()+"_distb").c_str(), ("Accepted Stub p_{T} Distribution for Layer"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_rejstub_pT_layer_distb[lIt] = fs->make<TH1F>(("h_rejstub_pT_layer"+lSStr.str()+"_distb").c_str(), ("Accepted Stub p_{T} Distribution for Layer"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_stubdiff_pT_layer_distb[lIt] = fs->make<TH1F>(("h_stubdiff_pT_layer"+lSStr.str()+"_distb").c_str(), ("Fraction of Acc'd/Rej Stubs p_{T} Distribution for Layer"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_accTrkPrt_layer_fraction[lIt] = fs->make<TH1F>(("h_accTrkPrt_layer"+lSStr.str()+"_fraction").c_str(), ("Fraction of Acc'd TrkParticles p_{T} Distribution for Layer"+lSStr.str()).c_str(), 600, 0.0, 600.0);
  }

  for ( uint lIt=0; lIt !=5; ++lIt )
  {
    std::stringstream lSStr;
    lSStr << (lIt+1);
    h_pTdistAccGenStub_disk[lIt] = fs->make<TH1F>(("h_pTdistAccGenStub_disk"+lSStr.str()).c_str(), ("pT Distribution of Accepted Tracking Particles for Disk"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_pTdistRejGenStub_disk[lIt] = fs->make<TH1F>(("h_pTdistRejGenStub_disk"+lSStr.str()).c_str(), ("pT Distribution of Rejected Tracking Particles for Disk"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_stub_pT_disk_distb[lIt] = fs->make<TH1F>(("h_stub_pT_disk"+lSStr.str()+"_distb").c_str(), ("Accepted Stub p_{T} Distribution for Disk"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_rejstub_pT_disk_distb[lIt] = fs->make<TH1F>(("h_rejstub_pT_disk"+lSStr.str()+"_distb").c_str(), ("Accepted Stub p_{T} Distribution for Disk"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_stubdiff_pT_disk_distb[lIt] = fs->make<TH1F>(("h_stubdiff_pT_disk"+lSStr.str()+"_distb").c_str(), ("Fraction of Acc'd/Rej Stubs p_{T} Distribution for Disk"+lSStr.str()).c_str(), 600, 0.0, 600.0);
    h_accTrkPrt_disk_fraction[lIt] = fs->make<TH1F>(("h_accTrkPrt_disk"+lSStr.str()+"_fraction").c_str(), ("Fraction of Acc'd TrkParticles p_{T} Distribution for Disk"+lSStr.str()).c_str(), 600, 0.0, 600.0);
  }

//////////////////////////////////////////////////////////////////////
  h_stubTrkPrtdiff_distb = fs->make<TH1F>("h_stubTrkPrtdiff_distb", "Number of stubs normalised against tracking particles p_{T} distribution", 600, 0.0, 600.0);
//////////////////////////////////////////////////////////////////////

  MuonPrtCounter = 0;
  ElectronPrtCounter = 0;
  TauPrtCounter = 0;
  ZPrtCounter = 0;
  WPrtCounter = 0;
  ChargedPionPrtCounter = 0;
  NeutralPionPrtCounter = 0;
  PhotonPrtCounter = 0;
  OtherPrtCounter = 0;

  lPrtCounter = 0;

  l4MuCounter = 0;
  l4eCounter = 0;
  l2e2MuCounter = 0;
  StubCreatingPrtCounter = 0;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
  if ( WRITE_TEXTFILE )
  {
    myFile.close();
    std::cout << "myFile.close()" << std::endl;
  }

  const int lTotalStubs (lCounterGenStub+lCounterCombStub+lCounterUnkStub);


  h_numberStubs->SetBinContent(1, lCounterGenStub);
  h_numberStubs->SetBinContent(2, lCounterCombStub);
  h_numberStubs->SetBinContent(3, lCounterUnkStub);

  h_numberStubsFrac->SetBinContent(1, double(double(lCounterGenStub)/double(lTotalStubs)));
  h_numberStubsFrac->SetBinContent(2, double(double(lCounterCombStub)/double(lTotalStubs)));
  h_numberStubsFrac->SetBinContent(3, double(double(lCounterUnkStub)/double(lTotalStubs)));

  const uint lCounterMax (600);
//   int lTotalAccPrts (0);

//   for ( uint lTotAccCounter = 1; lTotAccCounter != (lCounterMax+1); ++lTotAccCounter )
//   {
//     lTotalAccPrts += lAccTrkPrtClone->GetBinContent(lTotAccCounter);
// //     std::cout << "lTotalAccStubs: " << lTotalAccStubs << std::endl;
//   }

  // Get bin contents
  for ( uint lCounter = 1; lCounter != (lCounterMax+1); ++lCounter )
  {
//     std::cout << "lCounter: " << lCounter << std::endl;

    int lFivePlusStubs = h_FivePlusStubsPrts->GetBinContent(lCounter);
    int lLessFiveStubs = h_LessThanFiveStubsPrts->GetBinContent(lCounter);
    double lFiveStubsFraction(double(lFivePlusStubs)/double(lLessFiveStubs+lFivePlusStubs));
    if (lLessFiveStubs == 0) { lFiveStubsFraction = 0.0; }
    h_FivePlusStubsFraction->Fill(lCounter, lFiveStubsFraction);
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
//  Associated Stubs from TrkPrt loop data
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------

//     int lTrkPrtAllAccStubsPt = h_TrkPrtAccAllStubsPt->GetBinContent(lCounter); //  pT of each accepted stub
//     int lTrkPrtAllRejStubsPt = h_TrkPrtRejAllStubsPt->GetBinContent(lCounter); //  pT of each rejected st

    int lTrkPrtAllAccStubs = h_TrkPrtAccAllPt->GetBinContent(lCounter); //  pT of each accepted stub
    int lTrkPrtAllRejStubs = h_TrkPrtRejAllPt->GetBinContent(lCounter); //  pT of each rejected stub

    int lTrkPrtGenAccStubs = h_TrkPrtAccGenPt->GetBinContent(lCounter); //  pT of each Genuine stub
    int lTrkPrtGenComStubs = h_TrkPrtAccComPt->GetBinContent(lCounter); //  pT of each Comb stub
    int lTrkPrtRejAccStubs = h_TrkPrtRejGenPt->GetBinContent(lCounter); //  pT of each Rej Gen stub
    int lTrkPrtRejComStubs = h_TrkPrtRejComPt->GetBinContent(lCounter); //  pT of each Rej Comb stub
//

//     int lPrimaryTrkPrtPt = h_PrimaryParticlePt->GetBinContent(lCounter);

    int lTrkPrtFivePlusAccStubs = h_NumTrkPrtsFivePlusAccStubs->GetBinContent(lCounter);
//     int lTrkPrtFivePlusRejStubs = h_NumTrkPrtsFivePlusRejStubs->GetBinContent(lCounter);

//
    int lTotalTrkPrts = h_TrkPrtPt->GetBinContent(lCounter);            //  pT of particle
//     int lTotalTrkPrtsStubs = lTrkPrtAllAccStubsPt+lTrkPrtAllRejStubsPt;

    int lTotalStubs = lTrkPrtGenAccStubs+lTrkPrtGenComStubs+lTrkPrtRejAccStubs+lTrkPrtRejComStubs;
    int lTotalAccStubs = lTrkPrtGenAccStubs+lTrkPrtGenComStubs;

//     std::cout << "lTrkPrtAllAccStubs: " << lTrkPrtAllAccStubs << std::endl;
//     std::cout << "lTrkPrtGenAccStubs: " << lTrkPrtGenAccStubs << std::endl;

//     std::cout << "lTrkPrtAllAccStubs: " << lTrkPrtAllAccStubs << std::endl;
//     std::cout << "lTrkPrtAllRejStubs: " << lTrkPrtAllRejStubs << std::endl;
//     std::cout << "lTrkPrtGenAccStubs: " << lTrkPrtGenAccStubs << std::endl;
//     std::cout << "lTotalTrkPrts: " << lTotalTrkPrts << std::endl;

    double l_pT_frac (double(lTotalTrkPrts)/double(lPrtCounter));

    double h_acc_frac (double(lTrkPrtFivePlusAccStubs)/double(lTotalTrkPrts));
    double h_rej_frac (double(lTrkPrtAllRejStubs)/double(lTotalTrkPrts));
    double h_gen_frac (double(lTrkPrtGenAccStubs)/double(lTotalStubs));
    double h_gen_acc_prop (double(lTrkPrtGenAccStubs)/double(lTotalAccStubs));

//     if ( lPrimaryTrkPrtPt == 0 )
//     {
//       h_acc_frac = 0.0;
//       h_rej_frac = 0.0;
//     }

//     std::cout << "l_pT_frac: " << l_pT_frac << std::endl;

    if ( lTotalTrkPrts == 0 || lTrkPrtAllAccStubs ==0 )
    {
      h_acc_frac = 0.0;
      h_rej_frac = 0.0;
      h_gen_frac = 0.0;
    }

    if ( lTrkPrtAllAccStubs == 0)
    {
      h_gen_acc_prop = 0.0;
    }



    h_TrkPrtPt_Frac->Fill(lCounter, l_pT_frac);
    h_AccStubPrts_TotalTrkPrts_fraction->Fill(lCounter,h_acc_frac);
    h_RejStubPrts_TotalTrkPrts_fraction->Fill(lCounter,h_rej_frac);
    h_GenAccStubPrts_TotalTrkPrts_fraction->Fill(lCounter,h_gen_frac);
    h_GenAccStubPrts_AccStubPrts_fraction->Fill(lCounter,h_gen_acc_prop);

// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
//  Associated TrkPrts from stub loop data
// --------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------
    int lAccTrkPrts = h_pTdistAccGenStub->GetBinContent(lCounter);
    int lAccStubs = h_stub_pT_distb->GetBinContent(lCounter);
    int lRejTrkPrts = h_pTdistRejGenStub->GetBinContent(lCounter);
    int lRejStubs = h_rejstub_pT_distb->GetBinContent(lCounter);

//     int lTrkPrts

    int lAccTrkPrtsLayer[6];
    int lAccStubsLayer[6];
    int lRejTrkPrtsLayer[6];
    int lRejStubsLayer[6];

    int lAccTrkPrtsDisk[5];
    int lAccStubsDisk[5];
    int lRejTrkPrtsDisk[5];
    int lRejStubsDisk[5];

//     std::cout << __FILE__ << ":" << __LINE__ << std::endl;

    for ( int lIt = 0; lIt != 6; ++lIt )
    {
       lAccTrkPrtsLayer[lIt] = h_pTdistAccGenStub_layer[lIt]->GetBinContent(lCounter);
       lAccStubsLayer[lIt] = h_stub_pT_layer_distb[lIt]->GetBinContent(lCounter);
       lRejTrkPrtsLayer[lIt] = h_pTdistRejGenStub_layer[lIt]->GetBinContent(lCounter);
       lRejStubsLayer[lIt] = h_rejstub_pT_layer_distb[lIt]->GetBinContent(lCounter);
    }
    for ( int lIt = 0; lIt != 5; ++lIt )
    {
       lAccTrkPrtsDisk[lIt] = h_pTdistAccGenStub_disk[lIt]->GetBinContent(lCounter);
       lAccStubsDisk[lIt] = h_stub_pT_disk_distb[lIt]->GetBinContent(lCounter);
       lRejTrkPrtsDisk[lIt] = h_pTdistRejGenStub_disk[lIt]->GetBinContent(lCounter);
       lRejStubsDisk[lIt] = h_rejstub_pT_disk_distb[lIt]->GetBinContent(lCounter);
    }

//     std::cout << "lAccStubs/lRejStubs: " << lAccStubs << "/" << lRejStubs << std::endl;
   // Perform Calc
    int lTotalGenStubs = lAccStubs+lRejStubs;  // Total Genuine Stubs - surely should be all stubs? Combinatoric and Unknown also?
    int lTotalPrts = lAccTrkPrts+lRejTrkPrts;
    int lTotalGenStubsLayer[6];
    int lTotalPrtsLayer[6];
    int lTotalGenStubsDisk[5];
    int lTotalPrtsDisk[5];


    for ( int lIt = 0; lIt != 6; ++lIt )
    {
      lTotalGenStubsLayer[lIt] = lAccStubsLayer[lIt]+lRejStubsLayer[lIt];
      lTotalPrtsLayer[lIt] = lAccTrkPrtsLayer[lIt]+lRejTrkPrtsLayer[lIt];
    }
    for ( int lIt = 0; lIt != 5; ++lIt )
    {
      lTotalGenStubsDisk[lIt] = lAccStubsDisk[lIt]+lRejStubsDisk[lIt];
      lTotalPrtsDisk[lIt] = lAccTrkPrtsDisk[lIt]+lRejTrkPrtsDisk[lIt];
    }
//     std::cout << "lTotalGenStubs: " << lTotalGenStubs << std::endl;
    double lGenEffSq = (double(lAccStubs))/(double(lTotalGenStubs));
    double lWeightSquared = ((double(lAccTrkPrts))/(double(lTotalPrts)));
//////////////////////////////////////////////////////////////////////
    double lStubsVsPrts=0.0;
    if(lTotalPrts!=0){lStubsVsPrts = (double(lTotalGenStubs))/(double(lTotalPrts));}
//////////////////////////////////////////////////////////////////////
    double lGenEffSqLayer[6] = {0.0};
    double lWeightSquaredLayer[6] = {0.0};

    double lGenEffSqDisk[5] = {0.0};
    double lWeightSquaredDisk[5] = {0.0};

    for ( int lIt = 0; lIt != 6; ++lIt )
    {
      lGenEffSqLayer[lIt] = (double(lAccStubsLayer[lIt]))/(double(lTotalGenStubsLayer[lIt]));
      lWeightSquaredLayer[lIt] = ((double(lAccTrkPrtsLayer[lIt]))/(double(lTotalPrtsLayer[lIt])));
    }

    for ( int lIt = 0; lIt != 5; ++lIt )
    {
      lGenEffSqDisk[lIt] = (double(lAccStubsDisk[lIt]))/(double(lTotalGenStubsDisk[lIt]));
      lWeightSquaredDisk[lIt] = ((double(lAccTrkPrtsDisk[lIt]))/(double(lTotalPrtsDisk[lIt])));
    }

//     std::cout << "lWeightSquared: " << lWeightSquared << std::endl;
//     std::cout << "lWeight: " << lWeight << std::endl;
//     std::cout << "lEff: " << lEff << std::endl;
    // Print contents to new graph

    if ( lAccStubs == 0 )
    {
      lGenEffSq = 0;
    }

    if ( lAccTrkPrts == 0 )
    {
      lWeightSquared = 0;
    }

    for ( int lIt = 0; lIt != 6; ++lIt )
    {
      if ( lAccStubsLayer[lIt] == 0 )
      {
        lGenEffSqLayer[lIt] = 0;
      }

      if ( lAccTrkPrtsLayer[lIt] == 0 )
      {
        lWeightSquaredLayer[lIt] = 0;
      }
    }

    for ( int lIt = 0; lIt != 5; ++lIt )
    {
      if ( lAccStubsDisk[lIt] == 0 )
      {
        lGenEffSqDisk[lIt] = 0;
      }

      if ( lAccTrkPrtsDisk[lIt] == 0 )
      {
        lWeightSquaredDisk[lIt] = 0;
      }
    }

    h_stubdiff_pT_distb->Fill( lCounter, lGenEffSq );
    h_accTrkPrt_fraction->Fill( lCounter, lWeightSquared);

//////////////////////////////////////////////////////////////////////
    h_stubTrkPrtdiff_distb->Fill( lCounter, lStubsVsPrts);
///////////////////////////////////////////////////////////////////

    for ( int lIt = 0; lIt != 6; ++lIt )
    {
      h_stubdiff_pT_layer_distb[lIt]->Fill( lCounter, lGenEffSqLayer[lIt]);
      h_accTrkPrt_layer_fraction[lIt]->Fill( lCounter, lWeightSquaredLayer[lIt]);
    }

    for ( int lIt = 0; lIt != 5; ++lIt )
    {
      h_stubdiff_pT_disk_distb[lIt]->Fill( lCounter, lGenEffSqDisk[lIt]);
      h_accTrkPrt_disk_fraction[lIt]->Fill( lCounter, lWeightSquaredDisk[lIt]);
    }
  }
  h_stubdiff_pT_distb->GetXaxis()->SetTitle("stub p_{T}");
  h_stubdiff_pT_distb->GetXaxis()->CenterTitle();
  h_stubdiff_pT_distb->GetYaxis()->SetTitle("Acc'd Stubs/Total Stubs");
  h_stubdiff_pT_distb->GetYaxis()->CenterTitle();
//   h_stubdiff_pT_distb->Rebin(2);

  h_accTrkPrt_fraction->GetXaxis()->SetTitle("Tracking Particle p_{T}");
  h_accTrkPrt_fraction->GetXaxis()->CenterTitle();
  h_accTrkPrt_fraction->GetYaxis()->SetTitle("Acc'd Tracking Particle/Total Tracking Particles");
  h_accTrkPrt_fraction->GetYaxis()->CenterTitle();
//   h_accTrkPrt_fraction->Rebin(2);

///////////////////////////////////////////////////////////////////
  h_stubTrkPrtdiff_distb->GetXaxis()->SetTitle("p_{T}");
  h_stubTrkPrtdiff_distb->GetXaxis()->CenterTitle();
  h_stubTrkPrtdiff_distb->GetYaxis()->SetTitle("Number of stubs/Tracking Particles");
  h_stubTrkPrtdiff_distb->GetYaxis()->CenterTitle();
///////////////////////////////////////////////////////////////////

  std::cout << MuonPrtCounter << " muons present." << std::endl;
  std::cout << ElectronPrtCounter << " electrons present." << std::endl;
  std::cout << TauPrtCounter << " tauons present." << std::endl;
  std::cout << ZPrtCounter << " Z's present." << std::endl;
  std::cout << WPrtCounter << " W's present." << std::endl;
  std::cout << ChargedPionPrtCounter << " pi{+}'s present." << std::endl;
  std::cout << NeutralPionPrtCounter << " pi{0}'s present." << std::endl;
  std::cout << PhotonPrtCounter << " photons present." << std::endl;
  std::cout << HiggsPrtCounter << " SM/MSUSY Higgs present" << std::endl;
  std::cout << OtherPrtCounter << " Other particles present." << std::endl;

//   if (HIGGS_ANALYSIS)
//   {
//     std::cout << "Found " << l4MuCounter << " cases where only four muons have interacted with the tracker." << std::endl;
//     std::cout << "Found " << l4eCounter << " cases where only four electrons have interacted with the tracker." << std::endl;
//     std::cout << "Found " << l2e2MuCounter << " cases where only four electrons have interacted with the tracker." << std::endl;
// 
//   }

  std::cout << GenAccStubCreatingPrtCounter << " TrkPrts create genuine accepted stubs" << std::endl;
  std::cout << StubCreatingPrtCounter << " TrkPrts create genuine (accepted and rejected) stubs." << std::endl;
//   myFile.close();
  //std::cout << "myFile.close()" << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
