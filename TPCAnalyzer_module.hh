////////////////////////////////////////////////////////////////////////
// Class:       TPCAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        TPCAnalyzer_module.cc
//
// Generated at Mon Mar 15 04:43:39 2021 by Marina Bravo using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"


#include "TTree.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"

#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>

#define fXFidCut 200
#define fYFidCut 200
#define fZFidCut1 5
#define fZFidCut2 495

#define fDefaulNeutrinoID 99999

namespace test {
  class TPCAnalyzer
;
}


class test::TPCAnalyzer : public art::EDAnalyzer {
public:
  explicit TPCAnalyzer
(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCAnalyzer(TPCAnalyzer const&) = delete;
  TPCAnalyzer(TPCAnalyzer&&) = delete;
  TPCAnalyzer & operator=(TPCAnalyzer const&) = delete;
  TPCAnalyzer & operator=(TPCAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  void resetVars();
  int VertexToDriftTick(double vt, double vx);

  std::string fMCTruthLabel;
  std::string fSimEnergyDepositLabel;
  std::string fSimEnergyDepositInstanceLabel;
  std::string fSimEnergyDepositLabelOut;
  std::string fSimEnergyDepositInstanceLabelOut;
  std::string fRawDigitLabel;
  std::string fRecobWireLabel;
  std::string fHitLabel;
  std::string fSpacePointLabel;
  std::string fVertexLabel;
  bool fSaveTruth;
  bool fSaveSimED;
  bool fSaveSimEDOut;
  bool fSaveWaveforms;
  bool fSaveWires;
  bool fSaveHits;
  bool fSaveSpacePoints;
  bool fSaveVertex;
  bool fCreateTPCMap;
  bool fApplyFiducialCut;
  bool fApplyVertexSCE;

  TTree* fTree;
  int fEventID, fRunID, fSubRunID;

  //True variables
  std::vector<int> fTruePrimariesPDG;
  std::vector<double> fTruePrimariesE;
  double fTrueVx;
  double fTrueVy;
  double fTrueVz;
  double fTrueVt;
  int fTrueVU;
  int fTrueVV;
  int fTrueVC;
  int fTrueVTimeTick;
  double fTrueVEnergy;

  //True SimEnergyDeposits
  std::vector<double> fEnDepE;
  std::vector<double> fEnDepX;
  std::vector<double> fEnDepY;
  std::vector<double> fEnDepZ;
  std::vector<double> fEnDepT;

  //True SimEnergyDeposits Out
  std::vector<double> fEnDepEOut;
  std::vector<double> fEnDepXOut;
  std::vector<double> fEnDepYOut;
  std::vector<double> fEnDepZOut;
  std::vector<double> fEnDepTOut;

  //Hit variables
  std::vector<double> fHitsPeakTime;
  std::vector<double> fHitsIntegral;
  std::vector<double> fHitsChannel;
  std::vector<double> fHitsChi2;
  std::vector<double> fHitsNDF;

  //Space Point Variables
  std::vector<double> fSpacePointX;
  std::vector<double> fSpacePointY;
  std::vector<double> fSpacePointZ;
  std::vector<double> fSpacePointIntegral;

  //Waveforms
  std::vector<std::vector<double>> fRawChannelADC;
  std::vector<int> fRawChannelID;
  std::vector<double> fRawChannelPedestal;

  //Recob Wires
  unsigned int fNROIs;
  std::vector<std::vector<float>> fWireADC;
  std::vector<unsigned int> fWireID;
  std::vector<int> fWireStampTime;

  //Reconstructed vertex
  double fRecoVx;
  double fRecoVy;
  double fRecoVz;
  int fRecoVU;
  int fRecoVV;
  int fRecoVC;
  int fRecoVTimeTick;

  int fNAnalyzedEvents;

  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  unsigned int fNChannels;
  unsigned int fReadoutWindow;
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;
};





void test::TPCAnalyzer::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fTree=tfs->make<TTree>("AnaTPCTree", "Analysis Output Tree");
  fTree->Branch("RunID", &fRunID, "RunID/I");
  fTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");

  if(fSaveTruth){
    fTree->Branch("TruePrimariesPDG", &fTruePrimariesPDG);
    fTree->Branch("TruePrimariesE", &fTruePrimariesE);
    fTree->Branch("TrueVx", &fTrueVx, "TrueVx/D");
    fTree->Branch("TrueVy", &fTrueVy, "TrueVy/D");
    fTree->Branch("TrueVz", &fTrueVz, "TrueVz/D");
    fTree->Branch("TrueVt", &fTrueVt, "TrueVt/D");
    fTree->Branch("TrueVU", &fTrueVU, "TrueVU/I");
    fTree->Branch("TrueVV", &fTrueVV, "TrueVV/I");
    fTree->Branch("TrueVC", &fTrueVC, "TrueVC/I");
    fTree->Branch("TrueVTimeTick", &fTrueVTimeTick, "TrueVC/I");
    fTree->Branch("TrueVEnergy", &fTrueVEnergy, "TrueVEnergy/D");
  }

  if(fSaveSimED){
    fTree->Branch("EnDepE", &fEnDepE);
    fTree->Branch("EnDepX", &fEnDepX);
    fTree->Branch("EnDepY", &fEnDepY);
    fTree->Branch("EnDepZ", &fEnDepZ);
    fTree->Branch("EnDepT", &fEnDepT);
  }

  if(fSaveSimEDOut){
    fTree->Branch("EnDepEOut", &fEnDepEOut);
    fTree->Branch("EnDepXOut", &fEnDepXOut);
    fTree->Branch("EnDepYOut", &fEnDepYOut);
    fTree->Branch("EnDepZOut", &fEnDepZOut);
    fTree->Branch("EnDepTOut", &fEnDepTOut);
  }

  if(fSaveWaveforms){
    fTree->Branch("RawChannelADC", &fRawChannelADC);
    //fTree->Branch("RawChannelID", &fRawChannelID);
    //fTree->Branch("RawChannelPedestal", &fRawChannelPedestal);
  }

  if(fSaveWires){
    fTree->Branch("NROIs", &fNROIs);
    fTree->Branch("WireID", &fWireID);
    fTree->Branch("WireStampTime", &fWireStampTime);
    fTree->Branch("WireADC", &fWireADC);
    //fTree->Branch("RawChannelID", &fRawChannelID);
    //fTree->Branch("RawChannelPedestal", &fRawChannelPedestal);
  }

  if(fSaveWires){
    fTree->Branch("RawChannelADC", &fRawChannelADC);
    //fTree->Branch("RawChannelID", &fRawChannelID);
    //fTree->Branch("RawChannelPedestal", &fRawChannelPedestal);
  }

  if(fSaveHits){
    fTree->Branch("HitsIntegral", &fHitsIntegral);
    fTree->Branch("HitsPeakTime", &fHitsPeakTime);
    fTree->Branch("HitsChannel", &fHitsChannel);
    fTree->Branch("HitsChi2", &fHitsChi2);
    fTree->Branch("HitsNDF", &fHitsNDF);
  }

  if(fSaveSpacePoints){
    fTree->Branch("SpacePointX", &fSpacePointX);
    fTree->Branch("SpacePointY", &fSpacePointY);
    fTree->Branch("SpacePointZ", &fSpacePointZ);
    fTree->Branch("SpacePointIntegral", &fSpacePointIntegral);
  }

  if(fSaveVertex){
    fTree->Branch("RecoVx", &fRecoVx, "RecoVx/D");
    fTree->Branch("RecoVy", &fRecoVy, "RecoVy/D");
    fTree->Branch("RecoVz", &fRecoVz, "RecoVz/D");
    fTree->Branch("RecoVU", &fRecoVU, "RecoVU/I");
    fTree->Branch("RecoVV", &fRecoVV, "RecoVV/I");
    fTree->Branch("RecoVC", &fRecoVC, "RecoVC/I");
    fTree->Branch("RecoVTimeTick", &fRecoVTimeTick, "RecoVC/I");
  }

  fNAnalyzedEvents=0;
}

void test::TPCAnalyzer::endJob(){

  if(fCreateTPCMap){
    std::ofstream fileout("TPCMapping.txt");
    std::ofstream fileoutXYZ("TPCMappingXYZ.txt");
    if(fileout.is_open()){
        double xyz_start[3], xyz_end[3];;
        for(unsigned int ch=0; ch<fGeom->Nchannels(); ch++){
          std::vector<geo::WireID> wireV = fGeom->ChannelToWire(ch);
          for(size_t w=0; w<wireV.size(); w++){
            fGeom->WireEndPoints(wireV[w], xyz_start, xyz_end);
            fileout<<ch<<" "<<wireV[w].Plane<<" "<<wireV[w].TPC<<std::endl;
            fileoutXYZ<<ch<<" "<<wireV[w].Plane<<" "<<wireV[w].TPC<<" ";
            fileoutXYZ<<xyz_start[0]<<" "<<xyz_start[1]<<" "<<xyz_start[2]<<" ";
            fileoutXYZ<<xyz_end[0]<<" "<<xyz_end[1]<<" "<<xyz_end[2]<<std::endl;
            std::cout<<ch<<":"<<w<<" ID="<<wireV[w].Wire<<" Plane="<<wireV[w].Plane<<" TPC="<<wireV[w].TPC<<std::endl;
          }
        }
    }
    fileout.close();
    fileoutXYZ.close();
  }

}

DEFINE_ART_MODULE(test::TPCAnalyzer)
