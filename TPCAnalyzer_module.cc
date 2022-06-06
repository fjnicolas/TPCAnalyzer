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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"


#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"

#include "lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h"
#include "lardata/RecoBaseProxy/ProxyBase/withAssociated.h"
#include "lardata/RecoBaseProxy/ProxyBase/withParallelData.h"
#include "lardata/RecoBaseProxy/ProxyBase/withZeroOrOne.h"
#include "lardata/RecoBaseProxy/ProxyBase/getCollection.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/SpacePointSolver/TripletFinder.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/ParticleInventoryService.h"


#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "larevt/SpaceCharge/SpaceCharge.h"

#include "sbndcode/RecoUtils/RecoUtils.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"


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
#define fZFidCut1 0
#define fZFidCut2 500

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
  std::string fRawDigitLabel;
  std::string fRecobWireLabel;
  std::string fHitLabel;
  std::string fSpacePointLabel;
  bool fSaveTruth;
  bool fSaveSimED;
  bool fSaveWaveforms;
  bool fSaveWires;
  bool fSaveHits;
  bool fSaveSpacePoints;
  bool fCreateTPCMap;
  bool fApplyFiducialCut;

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

  //Hit variables
  std::vector<double> fHitsPeakTime;
  std::vector<double> fHitsIntegral;
  std::vector<double> fHitsChannel;

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

  int fNAnalyzedEvents;

  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  unsigned int fNChannels;
  unsigned int fReadoutWindow;
  double fTriggerOffsetTPC;
  double fTickPeriodTPC;
  double fDriftVelocity;
  double fWirePlanePosition;
};


test::TPCAnalyzer::TPCAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCTruthLabel( p.get<std::string>("MCTruthLabel", "generator") ),
  fSimEnergyDepositLabel( p.get<std::string>("SimEnergyDepositLabel", "ionandscint") ),
  fSimEnergyDepositInstanceLabel( p.get<std::string>("SimEnergyDepositInstanceLabel", "priorSCE") ),
  fRawDigitLabel( p.get<std::string>("RawDigitLabel", "daq") ),
  fRecobWireLabel( p.get<std::string>("RecobWireLabel", "caldata") ),
  fHitLabel( p.get<std::string>("HitLabel", "gaushit") ),
  fSpacePointLabel( p.get<std::string>("SpacePointLabel", "pandora") ),
  fSaveTruth( p.get<bool>("SaveTruth", "true") ),
  fSaveSimED( p.get<bool>("SaveSimED", "true") ),
  fSaveWaveforms( p.get<bool>("SaveWaveforms", "false") ),
  fSaveWires( p.get<bool>("SaveWires", "false") ),
  fSaveHits( p.get<bool>("SaveHits", "true") ),
  fSaveSpacePoints( p.get<bool>("SaveSpacePoints", "false") ),
  fCreateTPCMap( p.get<bool>("CreateTPCMap", "false") ),
  fApplyFiducialCut( p.get<bool>("ApplyFiducialCut", "true") ),
  fNChannels(fGeom->Nchannels())
  // More initializers here.
{

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

  fTriggerOffsetTPC = clockData.TriggerOffsetTPC(); //in us
  fTickPeriodTPC = clockData.TPCClock().TickPeriod(); //in us
  fReadoutWindow = detProp.ReadOutWindowSize();
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  fWirePlanePosition = std::abs( fGeom->Plane(1).GetCenter()[0] );


  std::cout<<"  - Read TPC clocks...  ReadOutWindowSize: "<<fReadoutWindow<<"  TriggerOffsetTPC: "<<fTriggerOffsetTPC;
  std::cout<<"  TickPeriodTPC: "<<fTickPeriodTPC<<std::endl;
  std::cout<<"  - Drift Velocity: "<<fDriftVelocity<<"  WirePlanePosition: "<<fWirePlanePosition<<std::endl;

  fRawChannelADC.resize(fNChannels, std::vector<double>(0));
  for(size_t k=0; k<fNChannels; k++) fRawChannelADC[k].reserve(fReadoutWindow);
}


void test::TPCAnalyzer::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  std::cout<<"Running TPCAnalyzer---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<"\n";
  //auto const fClockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  //auto const fDetProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, fClockData);

  //............................Event General Info
  fNAnalyzedEvents++;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fEventID = e.id().event();
  //Reset tree variables
  resetVars();


  //............................Read Truth Objects
  if(fSaveTruth){
    art::Handle<std::vector<simb::MCTruth>> mctruths;
    e.getByLabel(fMCTruthLabel, mctruths);
    std::cout<<" --- Saving MCTruth\n";

    for (auto const& truth : *mctruths) {

      for (int p = 0; p < truth.NParticles(); p++){
        simb::MCParticle const& TruePart = truth.GetParticle(p);

        if( TruePart.StatusCode()==1 ){
          fTruePrimariesPDG.push_back( TruePart.PdgCode() );
          fTruePrimariesE.push_back( TruePart.E() );
        }

        if( TruePart.Mother()==-1 && ( abs(TruePart.PdgCode())==12 || abs(TruePart.PdgCode())==14 ) ){
          fTrueVx=TruePart.EndX();
          fTrueVy=TruePart.EndY();
          fTrueVz=TruePart.EndZ();
          fTrueVt=TruePart.T();
          fTrueVEnergy=TruePart.E();

          if(fApplyFiducialCut && std::abs(fTrueVx)<fXFidCut && std::abs(fTrueVy)<fYFidCut && fTrueVz>fZFidCut1 && fTrueVz<fZFidCut2){
            const double p[3]={fTrueVx, fTrueVy, fTrueVz};
            std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(p))<<" "<<fGeom->FindTPCAtPosition(p).TPC<<std::endl;


            if( fGeom->HasTPC(fGeom->FindTPCAtPosition(p)) ){
              unsigned int tpcID=fGeom->FindTPCAtPosition(p).TPC;
              fTrueVU=fGeom->NearestChannel(p, 0, tpcID, 0);
              fTrueVV=fGeom->NearestChannel(p, 1, tpcID, 0);
              fTrueVC=fGeom->NearestChannel(p, 2, tpcID, 0);
              fTrueVTimeTick=VertexToDriftTick(fTrueVt, fTrueVx);
            }
          }
          else{
            fTrueVU=-1;
            fTrueVV=-1;
            fTrueVC=-1;
            fTrueVTimeTick=-1;
          }

          std::cout<<"  -- Vertex: "<<fTrueVx<<" "<<fTrueVy<<" "<<fTrueVz<<std::endl;
          std::cout<<"   - VertexWire: "<<fTrueVU<<" "<<fTrueVV<<" "<<fTrueVC<<std::endl;
        }
      }
    }
  }


  //............................Read SimEnergyDeposits
  if(fSaveSimED){
    art::Handle<std::vector<sim::SimEnergyDeposit> > SimEDHandle;
    e.getByLabel(fSimEnergyDepositLabel, fSimEnergyDepositInstanceLabel, SimEDHandle);
    std::cout<<"  ---- Reading SimEnergyDeposition from handle: "<<SimEDHandle.provenance()->moduleLabel();
    std::cout<<":"<<SimEDHandle.provenance()->productInstanceName()<<" ----\n";


    for (auto const& SimED : *SimEDHandle){
      fEnDepE.push_back(SimED.Energy());
      fEnDepX.push_back(SimED.MidPointX());
      fEnDepY.push_back(SimED.MidPointY());
      fEnDepZ.push_back(SimED.MidPointZ());
      fEnDepT.push_back( (SimED.StartT()+SimED.EndT())/2. );
    }
  }


  //............................Read Waveforms
  if(fSaveWaveforms){
    art::Handle<std::vector<raw::RawDigit>> eventRawDigits;
    std::vector<art::Ptr<raw::RawDigit>> eventRawDigitsVect;
    std::cout<<" --- Saving RawDigits\n";

    e.getByLabel(fRawDigitLabel, eventRawDigits);
    art::fill_ptr_vector(eventRawDigitsVect, eventRawDigits);

    //size_t NWaveforms=eventRawDigits->size();

    for (const art::Ptr<raw::RawDigit> &RD: eventRawDigitsVect){
      //std::cout<<i<<std::endl;
      //Store channel ID and wvf pedestal
      //fRawChannelID[i] = RD->Channel();
      //fRawChannelPedestal[i] = RD->GetPedestal();
      //Save ADC values

      for (long unsigned int j = 0; j<RD->ADCs().size(); j++){
        fRawChannelADC[RD->Channel()].push_back(RD->ADCs().at(j));
      }
    }
  }


  //............................Read Wires
  if(fSaveWires){
    art::Handle<std::vector<recob::Wire>> eventWires;
    std::vector<art::Ptr<recob::Wire>> eventWiresVect;
    std::cout<<" --- Saving Wires\n";

    e.getByLabel(fRecobWireLabel, eventWires);
    art::fill_ptr_vector(eventWiresVect, eventWires);

    for (const art::Ptr<recob::Wire> &W: eventWiresVect){

      const recob::Wire::RegionsOfInterest_t&       WireROIVect = W->SignalROI();
      unsigned int fWireChannel = W->Channel();

      for(const auto& ROI : WireROIVect.get_ranges()) {
        const std::vector<float>& WirePulse = ROI.data();
        //std::cout<<<<":"<<signal.size()<<":"<<range.begin_index()<<" ";
        fNROIs++;
        fWireID.push_back(fWireChannel);
        fWireStampTime.push_back(ROI.begin_index());
        fWireADC.push_back(WirePulse);
      }

    }
  }


  //............................Read Hits
  if(fSaveHits || fSaveSpacePoints){
    art::Handle<std::vector<recob::Hit>> hitsHandle;
    std::vector<art::Ptr<recob::Hit>> hitsVect;
    std::cout<<" --- Saving recob::Hit\n";
    e.getByLabel(fHitLabel, hitsHandle);
    art::fill_ptr_vector(hitsVect, hitsHandle);

    if(fSaveHits){
      for (const art::Ptr<recob::Hit> &hit: hitsVect){
        fHitsPeakTime.push_back(hit->PeakTime());
        fHitsIntegral.push_back(hit->Integral());
        fHitsChannel.push_back(hit->Channel());
      }
    }

    if(fSaveSpacePoints){
      art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
      std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
      std::cout<<" --- Saving recob::SpacePoints\n";

      e.getByLabel(fSpacePointLabel, eventSpacePoints);
      art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);

      art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);

      for (const art::Ptr<recob::SpacePoint> &SP: eventSpacePointsVect){

        std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());

        if (SPHit.at(0)->WireID().Plane==2){
          fSpacePointX.push_back(SP->position().X());
          fSpacePointY.push_back(SP->position().Y());
          fSpacePointZ.push_back(SP->position().Z());
          fSpacePointIntegral.push_back(SPHit.at(0)->Integral());
        }

      }
    }
  }


  fTree->Fill();
}

int test::TPCAnalyzer::VertexToDriftTick(double vt, double vx){
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}


void test::TPCAnalyzer::resetVars()
{
  if(fSaveTruth){
    fTruePrimariesPDG.clear();
    fTruePrimariesE.clear();
    fTrueVx=-1e3;
    fTrueVy=-1e3;
    fTrueVz=-1e3;
    fTrueVt=-1e3;
    fTrueVU=-1;
    fTrueVV=-1;
    fTrueVC=-1;
    fTrueVTimeTick=-1;
    fTrueVEnergy=-1e3;
  }

  if(fSaveSimED){
    fEnDepE.clear();
    fEnDepX.clear();
    fEnDepY.clear();
    fEnDepZ.clear();
    fEnDepT.clear();
  }

  if(fSaveWaveforms){
    fRawChannelID.clear();
    fRawChannelID.resize(fNChannels, -1);
    fRawChannelADC.clear();
    fRawChannelADC.resize(fNChannels, std::vector<double>(0));
    //for(size_t k=0; k<fNChannels; k++){fRawChannelADC[k].reserve(fReadoutWindow);}
    fRawChannelPedestal.resize(fNChannels, -1e3);
  }

  if(fSaveWires){
    fNROIs=0;
    fWireID.clear();
    fWireStampTime.clear();
    fWireADC.clear();
    //fWireADC.reserve(fNChannels);
    //fWireID.reserve(fNChannels);
  }

  if(fSaveHits){
    fHitsIntegral.clear();
    fHitsPeakTime.clear();
    fHitsChannel.clear();
  }

  if(fSaveSpacePoints){
    fSpacePointX.clear();
    fSpacePointY.clear();
    fSpacePointZ.clear();
    fSpacePointIntegral.clear();
  }

}



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
  }

  if(fSaveSpacePoints){
    fTree->Branch("SpacePointX", &fSpacePointX);
    fTree->Branch("SpacePointY", &fSpacePointY);
    fTree->Branch("SpacePointZ", &fSpacePointZ);
    fTree->Branch("SpacePointIntegral", &fSpacePointIntegral);
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
