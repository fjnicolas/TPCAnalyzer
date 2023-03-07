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
#include "lardataobj/RecoBase/Slice.h"
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
  std::string fRawDigitLabel;
  std::string fRecobWireLabel;
  std::string fHitLabel;
  std::string fSliceLabel;
  std::string fTrackLabel;
  std::string fSpacePointLabel;
  std::string fVertexLabel;
  bool fSaveTruth;
  bool fSaveSimED;
  bool fSaveWaveforms;
  bool fSaveWires;
  bool fSaveHits;
  bool fSaveSpacePoints;
  bool fSaveVertex;
  bool fSaveTrack;
  bool fCreateTPCMap;
  bool fApplyFiducialCut;
  int fMinHitsSlice;

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

  //Hit variables from track
  std::vector<double> fTrkHitsPeakTime;
  std::vector<double> fTrkHitsIntegral;
  std::vector<double> fTrkHitsChannel;

  int fnSlices;

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


test::TPCAnalyzer::TPCAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCTruthLabel( p.get<std::string>("MCTruthLabel", "generator") ),
  fSimEnergyDepositLabel( p.get<std::string>("SimEnergyDepositLabel", "ionandscint") ),
  fSimEnergyDepositInstanceLabel( p.get<std::string>("SimEnergyDepositInstanceLabel", "priorSCE") ),
  fRawDigitLabel( p.get<std::string>("RawDigitLabel", "daq") ),
  fRecobWireLabel( p.get<std::string>("RecobWireLabel", "caldata") ),
  fHitLabel( p.get<std::string>("HitLabel", "gaushit") ),
  fSliceLabel( p.get<std::string>("SliceLabel", "pandora") ),
  fTrackLabel( p.get<std::string>("TrackLabel", "pandoraTrack") ),
  fSpacePointLabel( p.get<std::string>("SpacePointLabel", "pandora") ),
  fVertexLabel( p.get<std::string>("VertexLabel", "pandora") ),
  fSaveTruth( p.get<bool>("SaveTruth", "true") ),
  fSaveSimED( p.get<bool>("SaveSimED", "true") ),
  fSaveWaveforms( p.get<bool>("SaveWaveforms", "false") ),
  fSaveWires( p.get<bool>("SaveWires", "false") ),
  fSaveHits( p.get<bool>("SaveHits", "true") ),
  fSaveSpacePoints( p.get<bool>("SaveSpacePoints", "false") ),
  fSaveVertex( p.get<bool>("SaveVertex", "true") ),
  fSaveTrack( p.get<bool>("SaveTrack", "true") ),
  fCreateTPCMap( p.get<bool>("CreateTPCMap", "false") ),
  fApplyFiducialCut( p.get<bool>("ApplyFiducialCut", "true") ),
  fMinHitsSlice( p.get<int>("MinHitsSlice", 25) ),
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

  ::art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel("gaushit", hitHandle);
  // Loop over slices
  std::vector< art::Ptr<recob::Hit> > hhVect;
  art::fill_ptr_vector(hhVect, hitHandle);
  std::cout<<" Number of total hits: "<<hhVect.size()<<std::endl;

  //............................Read Recob Slice
  ::art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceLabel, sliceHandle);
  //............................Read PFPs
  ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fSliceLabel, pfpHandle);
  //............................Read Recob Tracks
  ::art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  //Vector for recob hits
  std::vector<art::Ptr<recob::Hit>> hitVect;
  //Vector for recob PFParticles
  std::vector<art::Ptr<recob::PFParticle>> pfpVect;
  //Slice association for vertex
  art::FindManyP<recob::Hit> slice_hit_assns (sliceHandle, e, fSliceLabel);
  //Slice association for PFParticles
  art::FindManyP<recob::PFParticle> slice_pfp_assns (sliceHandle, e, fSliceLabel);
  // Loop over slices
  std::vector< art::Ptr<recob::Slice> > sliceVect;
  art::fill_ptr_vector(sliceVect, sliceHandle);
  std::cout<<" Number of slices: "<<sliceVect.size()<<std::endl;
  fnSlices = sliceVect.size();
  for(auto & slice:sliceVect){
    //----- get all hits associated to the slice
    hitVect = slice_hit_assns.at(slice.key());

    int nHits = hitVect.size();
    if(nHits<fMinHitsSlice){
      std::cout<<"   Number of hits in the slice "<<slice.key()<<": "<<hitVect.size()<<" ...skipping"<<std::endl;
      continue;
    }
    else{
      std::cout<<"   Number of hits in the slice "<<slice.key()<<": "<<hitVect.size()<<" ...saving"<<std::endl;
    }

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


    //............................Read Recob Hits
    if(fSaveHits){
      std::cout<<"    --- Saving recob::Hit\n";
      for (const art::Ptr<recob::Hit> &hit: hitVect){
        fHitsPeakTime.push_back(hit->PeakTime());
        fHitsIntegral.push_back(hit->Integral());
        fHitsChannel.push_back(hit->Channel());
      }
    }


    //----- get the space points
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

    //Get Vertex Association
    art::FindManyP<recob::Vertex> vertexAssoc (pfpHandle, e, fVertexLabel);

    //----- get the pfparticles for vertex
    if(fSaveVertex){
      std::cout<<" --- Saving Reconstructed Vertex\n";
      pfpVect = slice_pfp_assns.at(slice.key());
      size_t neutrinoID = fDefaulNeutrinoID;
      std::cout<<"   *** PFParticle size:"<<pfpVect.size()<<std::endl;
      for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
        std::cout<<"PFParticlePDG:"<<pfp->PdgCode()<<" Primary="<<pfp->IsPrimary()<<std::endl;
        if( !( pfp->IsPrimary() && ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ) ) continue;
        neutrinoID = pfp->Self();
      }
      std::cout<<"    ** NeutrinoID:"<<neutrinoID<<"\n\n";

      if(neutrinoID != fDefaulNeutrinoID){


        //PFParticle loop
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){

          //Get PFParticle Vertex
          std::vector< art::Ptr<recob::Vertex> > vertexVec = vertexAssoc.at(pfp.key());

          //PFParticle's vertex loop
          std::cout<<"     PFParticle: "<<pfp->Self()<<std::endl;
          for(const art::Ptr<recob::Vertex> &ver : vertexVec){
            double xyz_vertex[3];
            ver->XYZ(xyz_vertex);
            std::vector<double> xyz_vec(std::begin(xyz_vertex), std::end(xyz_vertex));
            //fpfpVertexPosition.push_back( xyz_vec );
            double chi2=ver->chi2(), chi2ndof=ver->chi2PerNdof();
            std::cout<<"  --VERTEX  ID="<<ver->ID()<<"  x,y,z="<<xyz_vertex[0]<<","<<xyz_vertex[1]<<","<<xyz_vertex[2];
            std::cout<<" Chi2="<<chi2<<" Chi2/DoF="<<chi2ndof<<" Status:"<<ver->status()<<",\n";
          }

          //Fill neutrino vertex
          if(pfp->Self()==neutrinoID){
            double xyz_vertex[3];
            vertexVec[0]->XYZ(xyz_vertex);
            std::cout<<"  Filling neutrino vertex...\n";

            fRecoVx= xyz_vertex[0];
            fRecoVy= xyz_vertex[1];
            fRecoVz= xyz_vertex[2];

            if(fApplyFiducialCut && std::abs(fRecoVx)<fXFidCut && std::abs(fRecoVy)<fYFidCut && fRecoVz>fZFidCut1 && fRecoVz<fZFidCut2){
              //const double p[3]={fVx, fTrueVy, fTrueVz};
              std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_vertex))<<" "
                <<fGeom->FindTPCAtPosition(xyz_vertex).TPC<<std::endl;

              if( fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_vertex)) ){
                unsigned int tpcID=fGeom->FindTPCAtPosition(xyz_vertex).TPC;
                fRecoVU=fGeom->NearestChannel(xyz_vertex, 0, tpcID, 0);
                fRecoVV=fGeom->NearestChannel(xyz_vertex, 1, tpcID, 0);
                fRecoVC=fGeom->NearestChannel(xyz_vertex, 2, tpcID, 0);
                fRecoVTimeTick=VertexToDriftTick(fTrueVt, fRecoVx);
              }
            }
            else{
              fRecoVU=-1;
              fRecoVV=-1;
              fRecoVC=-1;
              fRecoVTimeTick=-1;
            }
          }
        }
      }

    }


    art::FindManyP<recob::Track> pfp_track_assoc (pfpHandle, e, fTrackLabel);
    std::vector< art::Ptr<recob::Track> > trackVect;
    art::FindManyP<recob::Hit> track_hit_assoc (trackHandle, e, fTrackLabel);
    std::vector< art::Ptr<recob::Hit> > hitTrackVec;

    //----- get the pfparticles for vertex
    if(fSaveTrack){
      std::cout<<" --- Saving Longest Reco Track\n";
      pfpVect = slice_pfp_assns.at(slice.key());


      std::cout<<"   *** PFParticle size:"<<pfpVect.size()<<std::endl;

      //PFParticle loop
      double longestTrackLength=-1;
      art::Ptr<recob::Track> longest_trk;
      for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
        trackVect = pfp_track_assoc.at(pfp.key());

        std::cout<<"     PFParticle: "<<pfp->Self()<<" "<<pfp->PdgCode()<<" "<<trackVect.size();
        if(trackVect.size()!=0)std::cout<<" "<<trackVect[0]->Length()<<std::endl;
        else std::cout<<std::endl;
        if(trackVect.size()!=0){
          for(const art::Ptr<recob::Track> &trk : trackVect){
            std::cout<<" In track: "<<trk.key()<<std::endl;
            if(trk->Length()>longestTrackLength){
              longestTrackLength = trk->Length();
              longest_trk = trk;
            }
            hitTrackVec = track_hit_assoc.at(trk.key());
            std::cout<<" NHits: "<<hitTrackVec.size()<<std::endl;
          }
        }
      }


      if(longest_trk){
        std::cout<<" THE LONGEST LENGTH IS: "<<longest_trk->Length()<<std::endl;
        hitTrackVec = track_hit_assoc.at(longest_trk.key());
        std::cout<<"    --- Saving associated recob::Hit N="<<hitTrackVec.size()<<"\n";
        for (const art::Ptr<recob::Hit> &hit: hitTrackVec){
          fTrkHitsPeakTime.push_back(hit->PeakTime());
          fTrkHitsIntegral.push_back(hit->Integral());
          fTrkHitsChannel.push_back(hit->Channel());
        }
      }
    }

    fTree->Fill();
  }


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
    fnSlices = 0;
  }

  if(fSaveTrack){
    fTrkHitsIntegral.clear();
    fTrkHitsPeakTime.clear();
    fTrkHitsChannel.clear();
  }

  if(fSaveSpacePoints){
    fSpacePointX.clear();
    fSpacePointY.clear();
    fSpacePointZ.clear();
    fSpacePointIntegral.clear();
  }

  if(fSaveVertex){
    fRecoVx=-1e3;
    fRecoVy=-1e3;
    fRecoVz=-1e3;
    fRecoVU=-1;
    fRecoVV=-1;
    fRecoVC=-1;
    fRecoVTimeTick=-1;
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
    fTree->Branch("nSlices", &fnSlices, "nSlices/I");
  }

  if(fSaveTrack){
    fTree->Branch("TrkHitsIntegral", &fTrkHitsIntegral);
    fTree->Branch("TrkHitsPeakTime", &fTrkHitsPeakTime);
    fTree->Branch("TrkHitsChannel", &fTrkHitsChannel);
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

}

DEFINE_ART_MODULE(test::TPCAnalyzer)
