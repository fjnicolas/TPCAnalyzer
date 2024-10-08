////////////////////////////////////////////////////////////////////////
// Class:       TPCAnalyzer
// Plugin Type: analyzer (art v3_05_01)
// File:        TPCAnalyzer_module.cc
////////////////////////////////////////////////////////////////////////

#include "sbndcode/TPCAnalyzer/TPCAnalyzer_module.hh"

// Constructor
test::TPCAnalyzer::TPCAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCTruthLabel( p.get<std::string>("MCTruthLabel", "generator") ),
  fMCLabel( p.get<std::string>("MCLabel", "largeant") ),
  fSimEnergyDepositLabel( p.get<std::string>("SimEnergyDepositLabel", "ionandscint") ),
  fSimEnergyDepositInstanceLabel( p.get<std::string>("SimEnergyDepositInstanceLabel", "priorSCE") ),
  fSimEnergyDepositLabelOut( p.get<std::string>("SimEnergyDepositLabelOut", "ionandscintout") ),
  fSimEnergyDepositInstanceLabelOut( p.get<std::string>("SimEnergyDepositInstanceLabelOut", "") ),
  fSimChannelLabel( p.get<std::string>("SimChannelLabel") ),
  fRawDigitLabel( p.get<std::string>("RawDigitLabel", "daq") ),
  fRecobWireLabel( p.get<std::string>("RecobWireLabel", "caldata") ),
  fHitLabel( p.get<std::string>("HitLabel", "gaushit") ),
  fReco2Label( p.get<std::string>("Reco2Label", "pandora") ),
  fTrackLabel( p.get<std::string>("TrackLabel", "pandoraTrack") ),
  fClusterLabel( p.get<std::string>("ClusterLabel", "pandora") ),
  fSpacePointLabel( p.get<std::string>("SpacePointLabel", "pandora") ),
  fVertexLabel( p.get<std::string>("VertexLabel", "pandora") ),
  fCalorimetryLabel( p.get<std::string>("CalorimetryLabel", "pandoraCalo") ),
  fParticleIDLabel( p.get<std::string>("ParticleIDLabel", "pandoraPid") ),
  fSaveReco2( p.get<bool>("SaveReco2", "false") ),
  fSaveTruth( p.get<bool>("SaveTruth", "true") ),
  fSaveSimED( p.get<bool>("SaveSimED", "true") ),
  fSaveSimEDOut( p.get<bool>("SaveSimEDOut", "false") ),
  fSaveWaveforms( p.get<bool>("SaveWaveforms", "false") ),
  fSaveWires( p.get<bool>("SaveWires", "false") ),
  fSaveHits( p.get<bool>("SaveHits", "true") ),
  fSaveSpacePoints( p.get<bool>("SaveSpacePoints", "false") ),
  fSaveVertex( p.get<bool>("SaveVertex", "true") ),
  fCreateTPCMap( p.get<bool>("CreateTPCMap", "false") ),
  fApplyFiducialCut( p.get<bool>("ApplyFiducialCut", "true") ),
  fApplyVertexSCE( p.get<bool>("ApplyVertexSCE", "true") ),
  fUseSlices( p.get<bool>("UseSlices", "true") ),
  fUseSimChannels( p.get<bool>("UseSimChannels", "false") ),
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

  fCos60 = std::cos(  60 * (M_PI / 180.0) );
  fSin60 = std::sin(  60 * (M_PI / 180.0) );
  fWirePitch = 0.3;

  std::cout<<"Cos and sin 60 "<<fCos60<<" "<<fSin60<<std::endl;

  std::cout<<"  - Read TPC clocks...  ReadOutWindowSize: "<<fReadoutWindow<<"  TriggerOffsetTPC: "<<fTriggerOffsetTPC;
  std::cout<<"  TickPeriodTPC: "<<fTickPeriodTPC<<std::endl;
  std::cout<<"  - Drift Velocity: "<<fDriftVelocity<<"  WirePlanePosition: "<<fWirePlanePosition<<std::endl;

  fRawChannelADC.resize(fNChannels, std::vector<double>(0));
  for(size_t k=0; k<fNChannels; k++) fRawChannelADC[k].reserve(fReadoutWindow);
}

// Fill Hits function
void test::TPCAnalyzer::FillHits(int clusterId, std::vector<art::Ptr<recob::Hit>> hitVect, std::map<int, art::Ptr<recob::SpacePoint>> hitToSpacePointMap){
  for (const art::Ptr<recob::Hit> &hit: hitVect){
    fHitsView.push_back(hit->View());
    fHitsPeakTime.push_back(hit->PeakTime());
    fHitsIntegral.push_back(hit->Integral());
    fHitsSummedADC.push_back(hit->SummedADC());
    fHitsChannel.push_back(hit->Channel());
    fHitsAmplitude.push_back(hit->PeakAmplitude());
    fHitsRMS.push_back(hit->RMS());
    fHitsStartT.push_back(hit->StartTick());
    fHitsEndT.push_back(hit->EndTick());
    fHitsWidth.push_back( std::abs(hit->StartTick()-hit->EndTick()) );
    fHitsChi2.push_back(hit->GoodnessOfFit());
    fHitsNDF.push_back(hit->DegreesOfFreedom());
    fHitsClusterID.push_back(clusterId);
    if(fSaveSpacePoints){
      if(hitToSpacePointMap.find(hit.key())!=hitToSpacePointMap.end()){
        fHitsX.push_back(hitToSpacePointMap[hit.key()]->XYZ()[0]);
        fHitsY.push_back(hitToSpacePointMap[hit.key()]->XYZ()[1]);
        fHitsZ.push_back(hitToSpacePointMap[hit.key()]->XYZ()[2]);
      }
      else{
        fHitsX.push_back(-999);
        fHitsY.push_back(-999);
        fHitsZ.push_back(-999);
      }
    }
    else{
      fHitsX.push_back(-999);
      fHitsY.push_back(-999);
      fHitsZ.push_back(-999);
    }
  }
}

// Fill reco2 function
void test::TPCAnalyzer::FillReco2(art::Event const& e, std::vector<art::Ptr<recob::PFParticle>> pfpVect, std::map<int, art::Ptr<recob::SpacePoint>> hitToSpacePointMap){

    resetRecoVars();
    
    //Read PFPs
    ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    e.getByLabel(fReco2Label, pfpHandle);
    //Read Recob Tracks
    ::art::Handle<std::vector<recob::Track>> trackHandle;
    e.getByLabel(fTrackLabel, trackHandle);
    //Read Recob Cluster
    ::art::Handle<std::vector<recob::Cluster>> clusterHandle;
    e.getByLabel(fClusterLabel, clusterHandle);

    //Vector for recob hits
    std::vector<art::Ptr<recob::Hit>> hitVect;
    //PF to cluster
    art::FindManyP<recob::Cluster> pfp_cluster_assns (pfpHandle, e, fReco2Label);
    //PF to track
    art::FindManyP<recob::Track> pfp_track_assns (pfpHandle, e, fTrackLabel);
    //Cluster to hit
    art::FindManyP<recob::Hit> cluster_hit_assns (clusterHandle, e, fReco2Label);
    //Track to hit
    art::FindManyP<recob::Hit> track_hit_assns (trackHandle, e, fTrackLabel);
    //PF to vertex
    art::FindManyP<recob::Vertex> pfp_vertex_assns(pfpHandle, e, fReco2Label);
    // Track to calorimetry
    art::FindManyP<anab::Calorimetry> track_to_calo_assns(trackHandle, e, fCalorimetryLabel);
    // Track to PID
    art::FindManyP<anab::ParticleID> track_to_pid_assns(trackHandle, e, fParticleIDLabel);

    //PFParticle loop -- GetPrimary
    bool isNeutrino = false;
    for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){

      std::cout<<"   ** PFParticle: "<<pfp->Self()<<"      PDG:"<<pfp->PdgCode()<<"  Primary="<<pfp->IsPrimary()<<" Mother="<<pfp->Parent()<<std::endl;

      // Save reconstructed neutrino vertex
      if(  pfp->IsPrimary() && ( !fUseSlices || ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 ) ) ){
        std::cout<<"    This is a reconstructed netrino!\n";
        isNeutrino=true;
        //Get PFParticle Vertex
        std::vector< art::Ptr<recob::Vertex> > vertexVec = pfp_vertex_assns.at(pfp.key());
        for(const art::Ptr<recob::Vertex> &ver : vertexVec){
          geo::Point_t xyz_vertex = ver->position();
          double chi2=ver->chi2(), chi2ndof=ver->chi2PerNdof();
          std::cout<<"    --VERTEX  ID="<<ver->ID()<<"  x,y,z="<<xyz_vertex.X()<<","<<xyz_vertex.Y()<<","<<xyz_vertex.Z();
          std::cout<<" Chi2="<<chi2<<" Chi2/DoF="<<chi2ndof<<" Status:"<<ver->status()<<",\n";

          fRecoVx= xyz_vertex.X();
          fRecoVy= xyz_vertex.Y();
          fRecoVz= xyz_vertex.Z();

          std::cout<<"  Reco vertex in FV "<<PointInFV(fRecoVx, fRecoVy, fRecoVz)<<std::endl;

          if(fApplyFiducialCut && PointInFV(fRecoVx, fRecoVy, fRecoVz)){
            //const double p[3]={fVx, fTrueVy, fTrueVz};
            std::cout<<"     HasTPC: "<<fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_vertex))<<" "
              <<fGeom->FindTPCAtPosition(xyz_vertex).TPC<<std::endl;

            if( fGeom->HasTPC(fGeom->FindTPCAtPosition(xyz_vertex)) ){
              unsigned int tpcID=fGeom->FindTPCAtPosition(xyz_vertex).TPC;

              fRecoVU=fGeom->NearestChannel(xyz_vertex, geo::PlaneID(0, tpcID, 0));
              fRecoVV=fGeom->NearestChannel(xyz_vertex, geo::PlaneID(0, tpcID, 1));
              fRecoVC=fGeom->NearestChannel(xyz_vertex, geo::PlaneID(0, tpcID, 2));
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

      //Read cluster and store hits
      if(fUseSlices){
        std::vector<art::Ptr<recob::Cluster>> cluster_v = pfp_cluster_assns.at(pfp.key());
        for(size_t i=0; i<cluster_v.size(); i++){
          std::vector<art::Ptr<recob::Hit>> hitVect = cluster_hit_assns.at(cluster_v[i].key());
          std::cout<<"  ClusterID="<<cluster_v[i]->ID()<<" Hits: "<<hitVect.size()<<std::endl;
          FillHits(pfp->Self(), hitVect, hitToSpacePointMap);
        }
      }
     


      //Read the tracks and store the PFParticle start/end points
      std::vector<art::Ptr<recob::Track>> track_v = pfp_track_assns.at(pfp.key());
      for(size_t i=0; i<track_v.size(); i++){
        std::cout<<"     * Track number "<<i<<std::endl;
        std::cout<<"   "<<track_v[i]->Vertex()<<std::endl;
        std::cout<<"   "<<track_v[i]->End()<<std::endl;
        
        std::vector<double> start {track_v[i]->Vertex().X(), track_v[i]->Vertex().Y(), track_v[i]->Vertex().Z()};
        std::vector<double> end {track_v[i]->End().X(), track_v[i]->End().Y(), track_v[i]->End().Z()};
        fPFTrackStart.push_back(start);
        fPFTrackEnd.push_back(end);
        fPFPDGCode.push_back(pfp->PdgCode());

        if(!fUseSlices){
          std::vector<art::Ptr<recob::Hit>> hitVect = track_hit_assns.at(track_v[i].key());
          std::cout<<"  Hits: "<<hitVect.size()<<std::endl;
          FillHits(pfp->Self(), hitVect, hitToSpacePointMap);
        }
      
        // --- Get the associated calorimetry and PID objects
        std::vector<art::Ptr<anab::Calorimetry>> caloV = track_to_calo_assns.at(track_v[i].key());
        std::vector<art::Ptr<anab::ParticleID>> pidV = track_to_pid_assns.at(track_v[i].key());

        for (size_t j = 0; j < caloV.size(); ++j) {
          anab::Calorimetry calo = *caloV[j];

          // Collection plane
          std::cout<<"  Calorimetry object in plane "<<calo.PlaneID().Plane<<std::endl;
          std::cout<<"  Kinetic Energy: "<<calo.KineticEnergy()<<std::endl;
          
        } // end of calorimetry loop



      }
    }

    // if save reco1, save one slice per entry 
    // one entry per slice
    if(isNeutrino || !fUseSlices){
      fTree->Fill();
    }
    

}

// Main function
void test::TPCAnalyzer::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  auto fSCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  std::cout<<"\n -------------------------------------- \n";
  std::cout<<"Running TPCAnalyzer---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<"\n";

  //............................Event General Info
  fNAnalyzedEvents++;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fEventID = e.id().event();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

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
          if( TruePart.PdgCode()==2212 ) fIntNProtons++;
          if( TruePart.PdgCode()==2112 ) fIntNNeutrons++;
          if( TruePart.PdgCode()==111 ) fIntNPi0++;
          if( TruePart.PdgCode()==211 ) fIntNPip++;
          if( TruePart.PdgCode()==-211 ) fIntNPim++;
          if( TruePart.PdgCode()==13 ) fIntNMuonP++;
          if( TruePart.PdgCode()==-13 ) fIntNMuonM++;
          if( TruePart.PdgCode()==11 ) fIntNElectronP++;
          if( TruePart.PdgCode()==-11 ) fIntNElectronM++;
          if( TruePart.PdgCode()==3122 ) fIntNLambda++;
        }

        if( TruePart.Mother()==-1 && ( abs(TruePart.PdgCode())==12 || abs(TruePart.PdgCode())==14 ) ){
          fTrueVx=TruePart.EndX();
          fTrueVy=TruePart.EndY();
          fTrueVz=TruePart.EndZ();
          fTrueVt=TruePart.T();
          fTrueVEnergy=TruePart.E();

          if(fApplyVertexSCE){

            auto vOffset = fSCE->GetPosOffsets({fTrueVx, fTrueVy, fTrueVz});
            if(fTrueVx<0)
              fTrueVx = fTrueVx - vOffset.X();
            else
              fTrueVx = fTrueVx + vOffset.X();
            fTrueVy = fTrueVy + vOffset.Y();
            fTrueVz = fTrueVz + vOffset.Z();

            std::cout<<"Applying the SCE correction: "<<vOffset.X()<<" "<<vOffset.Y()<<" "<<vOffset.Z()<<std::endl;
            std::cout<<fSCE->EnableSimSpatialSCE()<<" "<<fSCE->EnableSimEfieldSCE()<<" "<<fSCE->EnableCorrSCE()<<" "<<fSCE->EnableCalSpatialSCE()<<std::endl;
          }
          
          fIntMode = truth.GetNeutrino().Mode();
          fIntCCNC = truth.GetNeutrino().CCNC();
          fIntInFV = PointInFV(fTrueVx, fTrueVy, fTrueVz);

          if(fApplyFiducialCut && PointInFV(fTrueVx, fTrueVy, fTrueVz)){
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


    // --- Fill the true information
    // --- Fill the true lambda variables
    // Get the handles
    art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthVect;
    e.getByLabel(fMCTruthLabel, mctruthHandle);
    art::fill_ptr_vector(mctruthVect, mctruthHandle);
      //............................Read MCParticles
    art::Handle< std::vector<simb::MCParticle> > mcparticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcpVect;
    e.getByLabel(fMCLabel, mcparticleHandle);
    art::fill_ptr_vector(mcpVect, mcparticleHandle);

    bool fFillLambdaTrue = true;
    if(fFillLambdaTrue){
      LambdaTruthManager lambdaMgr(mctruthVect, mcpVect);
      if(lambdaMgr.HasLambdaVDecayed()){
        fLambdaPionPDir.push_back(lambdaMgr.PionMomentumDirection().X());
        fLambdaPionPDir.push_back(lambdaMgr.PionMomentumDirection().Y());
        fLambdaPionPDir.push_back(lambdaMgr.PionMomentumDirection().Z());
        fLambdaProtonPDir.push_back(lambdaMgr.ProtonMomentumDirection().X());
        fLambdaProtonPDir.push_back(lambdaMgr.ProtonMomentumDirection().Y());
        fLambdaProtonPDir.push_back(lambdaMgr.ProtonMomentumDirection().Z());
      }
    }

  }

 
  //............................Read SimEnergyDeposits
  if(fSaveSimED){
    
    std::cout<<"  ---- Reading SimEnergyDeposition ----\n";

    if(!fUseSimChannels){

      art::Handle<std::vector<sim::SimEnergyDeposit> > SimEDHandle;
      e.getByLabel(fSimEnergyDepositLabel, fSimEnergyDepositInstanceLabel, SimEDHandle);
      std::cout<<"  ---- Reading SimEnergyDeposition from handle: "<<SimEDHandle.provenance()->moduleLabel();
      std::cout<<":"<<SimEDHandle.provenance()->productInstanceName()<<" ----\n";

      if(SimEDHandle.isValid()){
        std::cout<<"  ---- SimEnergyDeposition handle is valid ----\n";

        for (auto const& SimED : *SimEDHandle){
          float y = SimED.MidPointY();
          float z = SimED.MidPointZ();
          fEnDepE.push_back(SimED.Energy());
          fEnDepX.push_back(SimED.MidPointX());
          fEnDepY.push_back(y);
          fEnDepZ.push_back(z);
          fEnDepU.push_back( (z*fCos60-y*fSin60)/fWirePitch );
          fEnDepV.push_back( (z*fCos60+y*fSin60)/fWirePitch );
          fEnDepC.push_back( z/fWirePitch );
          fEnDepT.push_back( (SimED.StartT()+SimED.EndT())/2. );
          fEnDepPDG.push_back( SimED.PdgCode() );
        }

      }
      else{
        std::cout<<"  ---- SimEnergyDeposition handle is NOT valid ----\n";
      }

    }
    else{

      // Use SimChannels
      std::cout<<"  ---- Using SimChannels ---- " << fSimChannelLabel << "\n";
      std::unordered_map<int, std::vector<sim::IDE>> SimIDEMap;
      art::Handle< std::vector<sim::SimChannel> > SimChannelHandle;
      e.getByLabel(fSimChannelLabel, SimChannelHandle);
      std::cout << "  ---- Reading SimChannels from handle: " << SimChannelHandle.provenance()->moduleLabel();


      if(SimChannelHandle.isValid()){
        
      
        for(auto const &sc: *SimChannelHandle){ 
          
          // TDC-IDE map: for a given SimChannel, get a map (time tick in TDC, vector of ID)
          const auto & tdcidemap = sc.TDCIDEMap();
          
          // Loop over all of the tdc IDE map objects
          for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
            
            // Second is a vector of IDEs
            const std::vector<sim::IDE> idevec = (*mapitr).second;
            
            // Loop over all of the IDEs in a given simchannel
            for(size_t iv = 0; iv < idevec.size(); ++iv){
              
              if(SimIDEMap.find(idevec[iv].trackID)==SimIDEMap.end()){
                SimIDEMap[idevec[iv].trackID]=std::vector<sim::IDE>();
                SimIDEMap[idevec[iv].trackID].push_back(idevec[iv]);
              }
              else
                SimIDEMap[idevec[iv].trackID].push_back(idevec[iv]);
            }
          }
        }

        for(auto &simide: SimIDEMap){
          for(size_t k=0; k<simide.second.size(); k++){
            float y = simide.second[k].y;
            float z = simide.second[k].z;
            int pdg = pi_serv->TrackIdToParticle_P( std::abs( simide.second[k].trackID ) )->PdgCode();
            int mcpTime = pi_serv->TrackIdToParticle_P( std::abs( simide.second[k].trackID ) )->T();
            fEnDepE.push_back(simide.second[k].energy/3);
            fEnDepX.push_back(simide.second[k].x);
            fEnDepY.push_back(y);
            fEnDepZ.push_back(z);
            fEnDepU.push_back( (z*fCos60-y*fSin60)/fWirePitch );
            fEnDepV.push_back( (z*fCos60+y*fSin60)/fWirePitch );
            fEnDepC.push_back( z/fWirePitch );
            fEnDepT.push_back( mcpTime );
            fEnDepPDG.push_back( pdg );
          }
        }
      }
      else{
        std::cout<<"  ---- SimChannel handle is NOT valid ----\n";
      }

    }
  }


  //............................Read SimEnergyDeposits out volume
  if(fSaveSimEDOut){
    art::Handle<std::vector<sim::SimEnergyDeposit> > SimEDHandle;
    e.getByLabel(fSimEnergyDepositLabelOut, fSimEnergyDepositInstanceLabelOut, SimEDHandle);
    std::cout<<"  ---- Reading SimEnergyDeposition from handle: "<<SimEDHandle.provenance()->moduleLabel();
    std::cout<<":"<<SimEDHandle.provenance()->productInstanceName()<<" ----\n";

    if(SimEDHandle.isValid()){
      std::cout<<"  ---- SimEnergyDeposition handle is valid ----\n";
    

      for (auto const& SimED : *SimEDHandle){
        fEnDepEOut.push_back(SimED.Energy());
        fEnDepXOut.push_back(SimED.MidPointX());
        fEnDepYOut.push_back(SimED.MidPointY());
        fEnDepZOut.push_back(SimED.MidPointZ());
        fEnDepTOut.push_back( (SimED.StartT()+SimED.EndT())/2. );
      }
    
    }
    else{
      std::cout<<"  ---- SimEnergyDeposition handle is NOT valid ----\n";
    }

  }


  //............................Read Waveforms
  if(fSaveWaveforms){
    art::Handle<std::vector<raw::RawDigit>> eventRawDigits;
    std::vector<art::Ptr<raw::RawDigit>> eventRawDigitsVect;
    std::cout<<" --- Saving RawDigits\n";

    e.getByLabel(fRawDigitLabel, eventRawDigits);
    art::fill_ptr_vector(eventRawDigitsVect, eventRawDigits);

    for (const art::Ptr<raw::RawDigit> &RD: eventRawDigitsVect){

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
    resetWireVars();
    art::Handle<std::vector<recob::Wire>> eventWires;
    std::vector<art::Ptr<recob::Wire>> eventWiresVect;
    std::cout<<" --- Saving Wires\n";

    e.getByLabel(fRecobWireLabel, eventWires);
    art::fill_ptr_vector(eventWiresVect, eventWires);

    std::cout<<"  Number of wires: "<<eventWiresVect.size()<<std::endl;
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


  //............................Read Hits (Reco1 version)
  if(fSaveReco2==false){

    if(fSaveHits){
      art::Handle<std::vector<recob::Hit>> hitsHandle;
      std::vector<art::Ptr<recob::Hit>> hitsVect;
      std::cout<<" --- Saving recob::Hit\n";
      e.getByLabel(fHitLabel, hitsHandle);
      art::fill_ptr_vector(hitsVect, hitsHandle);

      for (const art::Ptr<recob::Hit> &hit: hitsVect){
        fHitsPeakTime.push_back(hit->PeakTime());
        fHitsIntegral.push_back(hit->Integral());
        fHitsSummedADC.push_back(hit->SummedADC());
        fHitsChannel.push_back(hit->Channel());
        fHitsAmplitude.push_back(hit->PeakAmplitude());
        fHitsRMS.push_back(hit->RMS());
        fHitsStartT.push_back(hit->StartTick());
        fHitsEndT.push_back(hit->EndTick());
        fHitsWidth.push_back( std::abs(hit->StartTick()-hit->EndTick()) );
        fHitsChi2.push_back(hit->GoodnessOfFit());
        fHitsNDF.push_back(hit->DegreesOfFreedom());
        fHitsClusterID.push_back(0);
      }
    }
    

    if(fSaveSpacePoints){
      art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
      std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
      std::cout<<"  --- Saving recob::SpacePoints\n";

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

  //............................Read Hits (Reco2 version)
  if( fSaveReco2 ){

    std::cout<<" --- Saving reco2\n";
    
    // map to store the space points associated to hits
    std::map<int, art::Ptr<recob::SpacePoint>> hitToSpacePointMap;

    if(fSaveSpacePoints){
      art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
      std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
      std::cout<<"      ** Creating hit-space point map\n";

      e.getByLabel(fSpacePointLabel, eventSpacePoints);
      art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);

      art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);

      for (const art::Ptr<recob::SpacePoint> &SP: eventSpacePointsVect){

        std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());

        //std::cout<<"  SpacePoint: "<<SP->ID()<<" "<<SPHit.at(0).key()<<std::endl;
        hitToSpacePointMap[SPHit.at(0).key()] = SP;
  
      }
    }

    // Necessary handles
    //Read Recob Slice
    ::art::Handle<std::vector<recob::Slice>> sliceHandle;
    e.getByLabel(fReco2Label, sliceHandle);
    
    //Vector for recob PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfpVect;


    if(fUseSlices){
      //Vector for recob Slices
      std::vector<art::Ptr<recob::Slice>> sliceVect;
      //Slice to PFParticles association
      art::FindManyP<recob::PFParticle> slice_pfp_assns (sliceHandle, e, fReco2Label);
    
      // Loop over slices
      art::fill_ptr_vector(sliceVect, sliceHandle);
      fNSlices = sliceVect.size();
      std::cout<<" Number of slices to analyze: "<<fNSlices<<std::endl;
    
      for(auto & slice:sliceVect){
          // Get the slices PFPs
          pfpVect = slice_pfp_assns.at(slice.key());
          size_t slice_ix = std::distance(sliceVect.begin(), std::find(sliceVect.begin(), sliceVect.end(), slice));
          std::cout<<std::endl<<slice_ix<<"  -- Slice ID="<<slice->ID()<<" NPFPs="<<pfpVect.size()<<std::endl;
          
          // Fill reco2
          FillReco2(e, pfpVect, hitToSpacePointMap);
      }//end slice loop
    }
    else{
      //Read Recob PFParticles
      ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
      e.getByLabel(fReco2Label, pfpHandle);
      art::fill_ptr_vector(pfpVect, pfpHandle);
      std::cout<<" Number of PFParticles to analyze: "<<pfpVect.size()<<std::endl;
      FillReco2(e, pfpVect, hitToSpacePointMap);
    }

  } //end SaveReco2 block

  // if save reco1, save one event per entry
  if(fSaveReco2==false){
    fTree->Fill();
  }
  
}


int test::TPCAnalyzer::VertexToDriftTick(double vt, double vx){
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}


bool test::TPCAnalyzer::PointInFV(double x, double y, double z){
  return ( std::abs(x)>fXFidCut1 && std::abs(x)<fXFidCut2 && std::abs(y)<fYFidCut && z>fZFidCut1 && z<fZFidCut2 );
}

void test::TPCAnalyzer::resetTrueVars(){
  
  if(fSaveTruth){
    fTruePrimariesPDG.clear();
    fTruePrimariesE.clear();
    fTruePrimariesStartP.clear();
    fTrueVx=-1e3;
    fTrueVy=-1e3;
    fTrueVz=-1e3;
    fTrueVt=-1e3;
    fTrueVU=-1;
    fTrueVV=-1;
    fTrueVC=-1;
    fTrueVTimeTick=-1;
    fTrueVEnergy=-1e3;

    fIntMode = -1;
    fIntCCNC = -1;
    fIntNProtons = 0;
    fIntNNeutrons = 0;
    fIntNPi0 = 0;
    fIntNPip = 0;
    fIntNPim = 0;
    fIntNMuonP = 0;
    fIntNMuonM = 0;
    fIntNElectronP = 0;
    fIntNElectronM = 0;
    fIntNLambda = 0;

    fLambdaProtonPDir.clear();
    fLambdaPionPDir.clear();
  }

  if(fSaveSimED){
    fEnDepE.clear();
    fEnDepX.clear();
    fEnDepY.clear();
    fEnDepZ.clear();
    fEnDepU.clear();
    fEnDepV.clear();
    fEnDepC.clear();
    fEnDepT.clear();
    fEnDepPDG.clear();
  }

  if(fSaveSimEDOut){
    fEnDepEOut.clear();
    fEnDepXOut.clear();
    fEnDepYOut.clear();
    fEnDepZOut.clear();
    fEnDepTOut.clear();
  }
}

void test::TPCAnalyzer::resetSimVars(){
  if(fSaveWaveforms){
    fRawChannelID.clear();
    fRawChannelID.resize(fNChannels, -1);
    fRawChannelADC.clear();
    fRawChannelADC.resize(fNChannels, std::vector<double>(0));
    //for(size_t k=0; k<fNChannels; k++){fRawChannelADC[k].reserve(fReadoutWindow);}
    fRawChannelPedestal.resize(fNChannels, -1e3);
  }
}

void test::TPCAnalyzer::resetWireVars()
{
  fNROIs=0;
  fWireID.clear();
  fWireStampTime.clear();
  fWireADC.clear();
}

void test::TPCAnalyzer::resetRecoVars()
{

  if(fSaveHits){
    fHitsView.clear();
    fHitsIntegral.clear();
    fHitsSummedADC.clear();
    fHitsPeakTime.clear();
    fHitsChannel.clear();
    fHitsAmplitude.clear();
    fHitsRMS.clear();
    fHitsStartT.clear();
    fHitsEndT.clear();
    fHitsWidth.clear();
    fHitsChi2.clear();
    fHitsNDF.clear();
    fHitsClusterID.clear();
    fHitsX.clear();
    fHitsY.clear();
    fHitsZ.clear();
  }

  fNSlices = 0;

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

  if(fSaveReco2){
    fPFTrackStart.clear();
    fPFTrackEnd.clear();
    fPFPDGCode.clear();
  }
}

void test::TPCAnalyzer::resetVars()
{
  resetTrueVars();
  resetSimVars();
  resetRecoVars();
}
