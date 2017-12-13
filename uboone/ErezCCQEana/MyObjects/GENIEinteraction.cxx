#ifndef GENIEinteraction_CXX
#define GENIEinteraction_CXX

#include "GENIEinteraction.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GENIEinteraction::GENIEinteraction (Int_t frun, Int_t fsubrun, Int_t fevent, Int_t fmcevent_id):
run(frun),
subrun(fsubrun),
event(fevent),
mcevent_id(fmcevent_id)
{
    fastest_proton_momentum = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetKinematics (Float_t fQ2, Float_t fXb, Float_t fW, Float_t fy, Int_t fccnc, Int_t fmode){
    Q2 = fQ2;
    Xb = fXb;
    W = fW;
    y = fy;
    ccnc = fccnc;
    mode = fmode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetMomentumTransfer (){
    // July-2017
    // should only be called after GENIEinteraction::SetNuMomentum() and GENIEinteraction::SetLeptonMomentum()
    q = nu - muon;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::AddPrimary ( Int_t fpdg                  // pdg code
                                   ,TLorentzVector fmomentum    // 4-momentum
                                   ,Int_t fstatus_code          // status code
                                   ,Int_t fmother               // mother
                                   ,std::string fprocess        // process
){
    
    momentum = fmomentum;
     // interaction neutrino
    switch (fpdg) {
            
        case 12: // ν(e)
            Nnu ++;
            Nnu_e ++;
            break;
            
        case 14: // ν(µ)
            Nnu ++;
            Nnu_mu ++;
            break;
            
        default:
            break;
    }
    
    if (fstatus_code==1) { // status code 0 particles are unstable or do not exit the nucleus and are thus irrelevant
        
        pdg.push_back(fpdg);
        mother.push_back(fmother);
        process.push_back(fprocess);
        status_code.push_back(fstatus_code);
        
        Nprimaries++;
        
        switch (pdg.back()) {
                
            case 13: // µ-
                Nmu++;
                Nmu_minus++;
                break;
                
            case -13: // µ+
                Nmu++;
                Nmu_plus++;
                break;
                
            case 2212: // p
                p3vect.push_back( momentum.Vect() ) ;
                Np++;
                break;
                
                
            case 2112: // n
                n3vect.push_back( momentum.Vect() ) ;
                Nn++;
                break;
                
            case 211: // π+
                Npi++;
                Npi_plus++;
                break;
                
            case -211: // π-
                Npi++;
                Npi_minus++;
                break;
                
            case 111: // π0
                Npi++;
                Npi_0++;
                break;
                
            case 11: // e-
                Ne_minus++;
                Nel++;
                break;
                
            case -11: // e+
                Ne_plus++;
                Nel++;
                break;
                
            case 22: // photon
                Ngamma++;
                break;
                
            default:
                break;
        }
    }
    return true;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::AddTrack(PandoraNuTrack ftrack){
    
    tracks.push_back(ftrack);
    switch (ftrack.GetMCpdgCode()) {
            // muon
        case 13:
            Is_mu_TrackReconstructed = true;
            muonTrack = ftrack;
            if (muonTrack.IsTrackInFV()) Is_mu_TrackInFV=true;
            break;
            
            // proton
        case 2212:
            Is_p_TrackReconstructed = true;
            // the "proton track" is the reconstructed track of the fastest proton
            // which reduces to the only reconstructed proton track,
            // in the case of a single reconstructed proton track
            // (1mu-1p and reconstructed CC1p0pi events)
            if ( ftrack.GetTruthMomentum().P() > fastest_proton_momentum
                &&
                ftrack.IsTrackInFV() ){
                
                fastest_proton_momentum = ftrack.GetTruthMomentum().P();
                protonTrack = ftrack;
                Is_p_TrackInFV=true;
            }
            break;
            
        default:
            break;
    }
    if (Is_mu_TrackReconstructed && Is_p_TrackReconstructed){
        IsVertexReconstructed = true;
        
        if ( Is_mu_TrackInFV && Is_p_TrackInFV ){
            IsVertexInFV = true;
        }
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::SortNucleons(){
    
    
    for (auto i: sort_by_momentum_magnitude( p3vect )){
        protons.push_back( TLorentzVector( p3vect.at(i) , sqrt( p3vect.at(i).Mag2() + 0.938*0.938 ) ) );
    }
    
    for (auto i: sort_by_momentum_magnitude( n3vect )){
        neutrons.push_back( TLorentzVector( n3vect.at(i) , sqrt( n3vect.at(i).Mag2() + 0.939*0.939 ) ) );
    }
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OåOooo........oooOO0OOooo......
vector<size_t> GENIEinteraction::sort_by_momentum_magnitude(const vector<TVector3> &v) {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1].Mag() > v[i2].Mag();});
    return idx;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::ComputePmissPrec(){
    
    if(Np>0){
        
        Plead = protons.at(0);
        
        // The leading proton was a neutron before the reaction!
        n_miss = Plead - q;
        
        for (auto proton: protons) {
            Pcm += proton;
        }
        Pcm  -= q;
        Prec  = Pcm - n_miss;
        
        
        // for SRC
        theta_pq    = (Plead.P()>0 && q.P()>0) ? TMath::RadToDeg() * Plead.Vect().Angle(q.Vect()) : -9999;
        p_over_q    = Plead.P()/q.P();
        
        // A(ν,𝝁p) missing mass M²(miss) = (q + (Mp+Mn) - Plead)² , all 4-vectors
        Mmiss       = (q + TLorentzVector( TVector3() , 0.938+0.939 ) - Plead).Mag();
        
    }
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetVertexPosition (TVector3 fpos){
    vertex_position = fpos;
    CheckContainement ();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetReco_mu_p_distance(){
    // the reconstructed distance between the muon and the proton tracks
    // only relevant for genie interactions in which both a muon and at least one proton were reconsturcted
    // for multiple reconstructed proton tracks, take the leading (fastest) proton
    // which is protons[0], after the protons were sorted
    // by GENIEinteraction::SortNucleons()
    if ( IsVertexReconstructed == true ){
        reco_mu_p_distance = muonTrack.ClosestDistanceToOtherTrack( protonTrack ) ;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::IncludesTrack(Int_t ftrackID) const{
    // ask if the track (ftrackID) is in the list of tracks belonging to this genie interaction
    for (auto t:tracks){
        if (t.GetTrackID() == ftrackID) return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetTruthTopology(){
    
    // determine the value of truth-definition flags
    
    Int_t Np_200MeVc = 0 , Nn_200MeVc = 0;
    for (auto Pp : p3vect){
        if (Pp.Mag()>=0.2) Np_200MeVc++;
    }
    for (auto Pn : n3vect){
        if (Pn.Mag()>=0.2) Nn_200MeVc++;
    }

    
    // IsCC_Np_200MeVc
    if ( ccnc==0 && Nmu>=1 && Np_200MeVc>=1 ){
        IsCC_Np_200MeVc = true;
        
        
        // IsCC_1p_200MeVc
        if ( Np_200MeVc==1 ){
            IsCC_1p_200MeVc = true;
            
            // IsCC_1p_200MeVc_0pi
            // which means that the final state includes
            // 1 protons with momentum >= 200 MeV/c
            // 0 neutrons with momentum >= 200 MeV/c
            // 0 pions
            // 0 electrons
            // 0 photons
            if ( Npi==0 && Nn_200MeVc==0 && Nel==0 && Ngamma==0 ){
                IsCC_1p_200MeVc_0pi = true;
                
                // if this genie is a CC 1p 0pi, all the tracks in this interactions should be flagged accordingly
                for (auto & t : tracks){
                    t.SetCC_1p_200MeVc_0pi();
                }
            }
        }
    }
    
    
    // IsCCQE : is QE or not - genie's "mode" flag: QE=0
    if ( mode==0 && ccnc==0 ){
        SetIsCCQE( true );
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetReconstructedTopology(){
    
    // determine the value of reconstructed-definition flags
    
    // tag vertices
    // ------------
    // contained        : vertex contained in active volume of the detector (256 x 233 x 1037 cm^3)
    // µ-reconstructed  : the muon track was reconstructed
    // p-reconstructed  : at least one proton was reconstructed
    // µp-reconstructed : both the muon and at least one proton were reconstructed
    // µp               : a muon and only one p were reconstructed, and nothing else was reconstructed
    
    SetIsInActiveVolume( CheckContainement() );
    
    // IsVertexReconstructed is determined in GENIEinteraction::AddTrack()
    // together with Is_mu_TrackReconstructed and Is_p_TrackReconstructed
    if (IsVertexReconstructed){
        SetReco_mu_p_distance();
        
        // 1mu-1p is a topology-based definitino:
        // a vertex in which only two tracks were reconstructed in the f.s.,
        // one is a muon and the other a proton
        if (tracks.size()==2){
            SetIs1mu1p( true );
        }
    }
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::Print(bool DoPrintTracks) const{

    cout << "\033[31m" << "GENIE interaction " << mcevent_id << "\n~~~~~~~~~~~~~~~~~~~~~ "<< "\033[0m" << endl;
    SHOW(mcevent_id);
    SHOWTVector3(vertex_position);
    if (!IsVertexContained) Printf("vertex not contained!");
    
    // std::vector-s of all primary paritcles
    SHOWstdVector( pdg );
    SHOWstdVector( mother );
    SHOWstdVector( process );
    
    SHOWTLorentzVector( nu );
    SHOWTLorentzVector( muon );
    SHOWTLorentzVector( q );
    if(!protons.empty()){
        cout << "\033[35m" << protons.size() << " protons:" << endl;
        for (auto proton:protons) {
            SHOWTLorentzVector( proton );
        }
        cout << "\033[31m";
    }
    if(!neutrons.empty()){
        cout << "\033[35m" << neutrons.size() << " neutrons:" << endl;
        for (auto neutron:neutrons) {
            SHOWTLorentzVector( neutron );
        }
        cout << "\033[31m";
    }
    SHOW2( Xb , Q2 );
    SHOW3( theta_pq , p_over_q , Mmiss );
    SHOW( Nprimaries );
    SHOW3( Np , Nn , Npi );
    SHOW3( Nmu , Nel , Ngamma );
    SHOW2( ccnc , mode );
    SHOW( IsVertexContained );
    SHOW( IsCCQE );
    SHOW3( Is_mu_TrackReconstructed , Is_p_TrackReconstructed , IsVertexReconstructed );
    SHOW3( Is_mu_TrackInFV , Is_p_TrackInFV , IsVertexInFV );
    SHOW2( Is1mu1p, IsCC_1p_200MeVc_0pi );
    
    if ( IsCC_1p_200MeVc_0pi ) {
        SHOW2( muonTrack.IsTrackContainedSoft() , protonTrack.IsTrackContainedSoft() );
    }

    if(DoPrintTracks && !tracks.empty()){
        cout << "\033[33m" << "--------------\n"
        << tracks.size() << " pandoraNu tracks in GENIE interaction "
        << mcevent_id << "\n--------------" <<  "\033[37m" << endl;
        SHOWstdVector( trackID );
        for (auto t: tracks) {
            t.Print(true);
        }
    }
    
    if (Is_mu_TrackReconstructed){
        cout  << "\033[33m" << "muon track:"  << "\033[30m" << endl;
        muonTrack.Print();
    } else {
        cout  << "\033[33m" << "no muon track reconstructed"  << "\033[30m" << endl;
    }
    if (Is_p_TrackReconstructed){
        cout  << "\033[33m" << "proton track:"  << "\033[30m" << endl;
        protonTrack.Print();
    } else {
        cout  << "\033[33m" << "no proton track reconstructed"  << "\033[30m" << endl;
    }
    if (IsVertexReconstructed) SHOW(reco_mu_p_distance);

}


#endif
