////////////////////////////////////////////////////////////////////////
// Class:       TPCAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        TPCAnalyzer_module.cc
//
// Generated at Mon Mar 15 04:43:39 2021 by Marina Bravo using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "sbndcode/TPCAnalyzer/TPCAnalyzer_module.hh"


// Constructor
test::TPCAnalyzer::TPCAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCTruthLabel( p.get<std::string>("MCTruthLabel", "generator") ),
  fSimEnergyDepositLabel( p.get<std::string>("SimEnergyDepositLabel", "ionandscint") ),
  fSimEnergyDepositInstanceLabel( p.get<std::string>("SimEnergyDepositInstanceLabel", "priorSCE") ),
  fRawDigitLabel( p.get<std::string>("RawDigitLabel", "daq") ),
  fRecobWireLabel( p.get<std::string>("RecobWireLabel", "caldata") ),
  fHitLabel( p.get<std::string>("HitLabel", "gaushit") ),
  fSpacePointLabel( p.get<std::string>("SpacePointLabel", "pandora") ),
  fVertexLabel( p.get<std::string>("VertexLabel", "pandora") ),
  fSaveTruth( p.get<bool>("SaveTruth", "true") ),
  fSaveSimED( p.get<bool>("SaveSimED", "true") ),
  fSaveWaveforms( p.get<bool>("SaveWaveforms", "false") ),
  fSaveWires( p.get<bool>("SaveWires", "false") ),
  fSaveHits( p.get<bool>("SaveHits", "true") ),
  fSaveSpacePoints( p.get<bool>("SaveSpacePoints", "false") ),
  fSaveVertex( p.get<bool>("SaveVertex", "false") ),
  fCreateTPCMap( p.get<bool>("CreateTPCMap", "false") ),
  fApplyFiducialCut( p.get<bool>("ApplyFiducialCut", "true") ),
  fApplyVertexSCE( p.get<bool>("ApplyVertexSCE", "true") ),
  fNChannels(fGeom->Nchannels())
  // More initializers here.
{

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

  fTriggerOffsetTPC = clockData.TriggerOffsetTPC(); //in us
  fTickPeriodTPC = clockData.TPCClock().TickPeriod(); //in us
  fReadoutWindow = detProp.ReadOutWindowSize();
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  constexpr geo::TPCID tpcid{0, 0};
  fWirePlanePosition = std::abs( fGeom->Plane(geo::PlaneID{tpcid, 1}).GetCenter().X() );


  std::cout<<"  - Read TPC clocks...  ReadOutWindowSize: "<<fReadoutWindow<<"  TriggerOffsetTPC: "<<fTriggerOffsetTPC;
  std::cout<<"  TickPeriodTPC: "<<fTickPeriodTPC<<std::endl;
  std::cout<<"  - Drift Velocity: "<<fDriftVelocity<<"  WirePlanePosition: "<<fWirePlanePosition<<std::endl;

  fRawChannelADC.resize(fNChannels, std::vector<double>(0));
  for(size_t k=0; k<fNChannels; k++) fRawChannelADC[k].reserve(fReadoutWindow);
}


// Main function
void test::TPCAnalyzer::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  auto fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  std::cout<<"Running TPCAnalyzer---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<"\n";

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

          if(fApplyVertexSCE){

            auto vOffset = fSCE->GetPosOffsets({fTrueVx, fTrueVy, fTrueVz});
            fTrueVx -= vOffset.X();
            fTrueVy += vOffset.Y();
            fTrueVz += vOffset.Z();

            std::cout<<"Applying SCE: "<<vOffset.X()<<" "<<vOffset.Y()<<" "<<vOffset.Z()<<std::endl;
            std::cout<<fSCE->EnableSimSpatialSCE()<<" "<<fSCE->EnableSimEfieldSCE()<<" "<<fSCE->EnableCorrSCE()<<" "<<fSCE->EnableCalSpatialSCE()<<std::endl;
          }
          
          if(fApplyFiducialCut && std::abs(fTrueVx)<fXFidCut && std::abs(fTrueVy)<fYFidCut && fTrueVz>fZFidCut1 && fTrueVz<fZFidCut2){
            geo::Point_t po ={fTrueVx, fTrueVy, fTrueVz};

            std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(po))<<" "<<fGeom->FindTPCAtPosition(po).TPC<<std::endl;


            if( fGeom->HasTPC(fGeom->FindTPCAtPosition(po)) ){
              
              unsigned int tpcID=fGeom->FindTPCAtPosition(po).TPC;
              geo::PlaneID plane(0, tpcID, 0);
              //fTrueVU=fGeom->NearestChannel(po, 0, tpcID, 0);
              fTrueVU=fGeom->NearestChannel(po, geo::PlaneID(0, tpcID, 0));
              fTrueVV=fGeom->NearestChannel(po, geo::PlaneID(0, tpcID, 1));
              fTrueVC=fGeom->NearestChannel(po, geo::PlaneID(0, tpcID, 2));
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
        fHitsChi2.push_back(hit->GoodnessOfFit());
        fHitsNDF.push_back(hit->DegreesOfFreedom());
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

  //............................Read Reco Vertex
  if(fSaveVertex){
    art::Handle< std::vector<recob::PFParticle> > eventPFParticle;
    std::vector< art::Ptr<recob::PFParticle> > pfparticleVect;
    if(e.getByLabel(fVertexLabel, eventPFParticle))
      art::fill_ptr_vector(pfparticleVect, eventPFParticle);
    std::cout<<" --- Saving Reconstructed Vertex\n";

    size_t neutrinoID = fDefaulNeutrinoID;
    std::cout<<"   *** PFParticle size:"<<pfparticleVect.size()<<std::endl;
    for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
      std::cout<<"PFParticlePDG:"<<pfp->PdgCode()<<" Primary="<<pfp->IsPrimary()<<std::endl;
      if( !( pfp->IsPrimary() && ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ) ) continue;
      neutrinoID = pfp->Self();
    }
    std::cout<<"    ** NeutrinoID:"<<neutrinoID<<"\n\n";

    if(neutrinoID != fDefaulNeutrinoID){
      //Get Vertex Association
      art::FindManyP<recob::Vertex> vertexAssoc (pfparticleVect, e, fVertexLabel);

      //PFParticle loop
      for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){

        //Get PFParticle Vertex
        std::vector< art::Ptr<recob::Vertex> > vertexVec = vertexAssoc.at(pfp.key());

        /*//PFParticle's track loop
        if(pfpTracks.empty()) {
          continue;
        }*/

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
            geo::Point_t xyz_vertexP{fRecoVx, fRecoVy, fRecoVz};
            std::cout<<"HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_vertexP))<<" "
              <<fGeom->FindTPCAtPosition(xyz_vertexP).TPC<<std::endl;

            if( fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_vertexP)) ){
              unsigned int tpcID=fGeom->FindTPCAtPosition(xyz_vertexP).TPC;
              
              fRecoVU=fGeom->NearestChannel(xyz_vertexP, geo::PlaneID(0, tpcID, 0) );
              fRecoVV=fGeom->NearestChannel(xyz_vertexP, geo::PlaneID(0, tpcID, 1) );
              fRecoVC=fGeom->NearestChannel(xyz_vertexP, geo::PlaneID(0, tpcID, 2) );
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
    fHitsChi2.clear();
    fHitsNDF.clear();
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