#ifndef TRUNCMEAN_CXX
#define TRUNCMEAN_CXX

#include "TruncMean.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TruncMean::TruncMean(int fdebug, float frad):
debug(fdebug),
rad(frad)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TruncMean::CalcTruncMeanProfile(const std::vector< float > & rr_v,
                                     const std::vector< float > & dq_v,
                                     std::vector< float > & dq_trunc_v,
                                     const float & 	nsigma
                                     ) {
    
    // how many points to sample
    int Nneighbor = (int)(rad * 3 * 2);
    
    dq_trunc_v.clear();
    dq_trunc_v.reserve( rr_v.size() );
    
    int Nmax = dq_v.size()-1;
    
    for (size_t n=0; n < dq_v.size(); n++) {
        
        // current residual range
        float rr = rr_v.at(n);
        
        int nmin = n - Nneighbor;
        int nmax = n + Nneighbor;
        
        if (nmin < 0) nmin = 0;
        if (nmax > Nmax) nmax = Nmax;
        
        // vector for local dq values
        std::vector<float> dq_local_v;
        
        for (int i=nmin; i < nmax; i++) {
            
            float dr = rr - rr_v[i];
            if (dr < 0) dr *= -1;
            
            if (dr > rad) continue;
            
            dq_local_v.push_back( dq_v[i] );
            
        }// for all ticks we want to scan
        
        if (dq_local_v.size() == 0) {
            dq_trunc_v.push_back( dq_v.at(n) );
            continue;
        }
        
        // calculate median and rms
        float median = Median(dq_local_v);
        float rms    = RMS(dq_local_v);
        
        float truncated_dq = 0.;
        int npts = 0;
        for (auto const& dq : dq_local_v) {
            if ( ( dq < (median+rms * nsigma) ) && ( dq > (median-rms * nsigma) ) ){
                truncated_dq += dq;
                npts += 1;
            }
        }
        
        dq_trunc_v.push_back( truncated_dq / npts );
    }// for all values
    
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TruncMean::CalcTruncMean(const std::vector<float>& resrange,
                              const std::vector<float>& dQ,
                              std::vector<float>& dQtruncated) {
    //int debug_level = 0;
    //Debug( debug_level , "TruncMean::CalcTruncMean()");
    
    // how many neighboring points to sample for the truncation
    int Nneighbor = (int)(rad * 3 * 2);
    //Debug( debug_level , "Nneighbor: %",Nneighbor);

    dQtruncated.clear();
    dQtruncated.reserve( resrange.size() );
    int Nmax = dQ.size()-1;
    //Debug(debug_level,"Nmax:%",Nmax);

    for (int n=0; n < int(dQ.size()); n++) {
        
        // current residual range
        float c_resrange = resrange.at(n);
        
        int nmin = std::max(n - Nneighbor , 0);
        int nmax = std::min(n + Nneighbor , Nmax);
        //Debug(debug_level,"n:%, sampling between % < n < % for the truncation",n, nmin,nmax);
        
        // vector for local dQ values
        std::vector<float> dQlocal;
        for (int i=nmin; i < nmax; i++) {
            float dr = fabs(c_resrange - resrange[i]);
            //Debug(debug_level,"resrange distance dr: % (radius=%)",dr,rad);
            if (dr < rad) {
                dQlocal.push_back( dQ[i] );
                //Debug(debug_level,"added % into dQlocal",dQ[i]);
            }
            else{
                //Debug(debug_level,"did not add % into dQlocal",dQ[i]);
            }
        }// for all ticks we want to scan
        
        //Debug(debug_level,"dQlocal.size():%",dQlocal.size());
        //        if (dQlocal.size() == 0) {
        if (dQlocal.size() <= 1) { // this case dQlocal.size()==1 makes the rms diverge, and we want to avoid it
            dQtruncated.push_back( dQ.at(n) );
            continue;
        }
        
        // calculate median and rms
        float median = Median(dQlocal);
        float rms    = RMS(dQlocal);
        //Debug( debug_level , "median=% rms=%",median,rms);
        
        float c_dQlocal = 0.;
        float c_dQtruncated = 0.;
        int npts = 0;
        for (auto const& dq : dQlocal) {
            c_dQlocal += dq;
            if (      ( dq < (median+rms) )
                &&    ( dq > (median-rms) ) ){
                c_dQtruncated += dq;
                npts += 1;
            }
        }
        dQtruncated.push_back( c_dQtruncated / npts );
        //Debug( debug_level , "pushed % into dQtruncated",dQtruncated.back());

    }// for all values
    
    return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float TruncMean::Median(const std::vector<float>& v){
    
    if (v.size() == 1) return v[0];
    
    std::vector<float> vcpy = v;
    
    std::sort(vcpy.begin(), vcpy.end());
    
    float median = vcpy[ vcpy.size() / 2 ];
    
    return median;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float TruncMean::RMS(const std::vector<float>& v){
    
    float avg = 0.;
    for (auto const& val : v) avg += val;
    avg /= v.size();
    float rms = 0.;
    for (auto const& val : v) rms += (val-avg)*(val-avg);
    rms = sqrt( rms / ( v.size() -  1 ) );
    
    return rms;
}

#endif
