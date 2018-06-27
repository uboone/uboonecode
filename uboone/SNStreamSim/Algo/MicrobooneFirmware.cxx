#ifndef MICROBOONEFIRMWARE_CXX
#define MICROBOONEFIRMWARE_CXX

#include "MicrobooneFirmware.h"
#include <limits>
#include <cstddef>

namespace compress {
  

  MicrobooneFirmware::MicrobooneFirmware(fhicl::ParameterSet const& pset)
    : CompressionAlgoBase()
  {

    _thresh = std::vector<int>(3,0.);
    _pol    = std::vector<int>(3,0);
    std::vector<std::vector<int> > tmp(3,std::vector<int>(2,0));
    _buffer = tmp;
    
    std::vector<int> fCompressThresholds   = pset.get< std::vector<int> > ("CompressThresholds");
    std::vector<int> fPolarity             = pset.get< std::vector<int> > ("Polarity");
    _maxADC                                = pset.get< int >              ("MaxADC");
    std::vector<int> fPlaneBuffers         = pset.get< std::vector<int> > ("PlaneBuffers");
    _block                                 = pset.get< int >              ("BlockSize");
    _deltaB                                = pset.get< int >              ("BaselineThreshold");
    _deltaV                                = pset.get< int >              ("VarianceThreshold");
    bool fDebug                            = pset.get< bool >             ("Debug");

    // the input for _deltaB is in ADC but the value used in the algorithm is in ADC^2:
    _deltaB *= _deltaB;

    if (fPlaneBuffers.size() != 6) {
      throw std::runtime_error("ERROR in MicrobooneFirmware: incorrect size of input for plane-buffer values.");
    }
    
    SetUVYplaneBuffer(fPlaneBuffers[0],fPlaneBuffers[1],fPlaneBuffers[2],
		      fPlaneBuffers[3],fPlaneBuffers[4],fPlaneBuffers[5]);

    if (fPolarity.size() != 3) {
      throw std::runtime_error("ERROR in MicrobooneFirmware: incorrect size of input for Polarity values.");
    }

    SetPolarity(fPolarity[0], fPolarity[1], fPolarity[2]);

    if (fCompressThresholds.size() != 3) {
      throw std::runtime_error("ERROR in MicrobooneFirmware: incorrect size of input for CompressThresholds values.");
    }

    SetCompressThresh(fCompressThresholds[0], fCompressThresholds[1], fCompressThresholds[2]);

    SetDebug(fDebug);

  }

  void MicrobooneFirmware::SetUVYplaneBuffer(int upre, int upost, int vpre, int vpost, int ypre, int ypost){
    
    _buffer[0][0]  = upre;
    _buffer[0][1] = upost;
    _buffer[1][0]  = vpre;
    _buffer[1][1] = vpost;
    _buffer[2][0]  = ypre;
    _buffer[2][1] = ypost;
    return;
  }
  
  
  void MicrobooneFirmware::ApplyCompression(const std::vector<short> &waveform, const int mode, const unsigned int ch){
    
    // iterator to the beginning and end of the waveform
    _begin = waveform.begin();
    _end   = waveform.end();
    
    std::vector<short> outputwf;
    _baselines.clear();
    _variances.clear();
    
    unsigned int _baseline[3];
    unsigned int _variance[3];
    _baseline[0] = std::numeric_limits<unsigned int>::max();
    _baseline[1] = std::numeric_limits<unsigned int>::max();
    _baseline[2] = std::numeric_limits<unsigned int>::max();
    _variance[0] = std::numeric_limits<unsigned int>::max();
    _variance[1] = std::numeric_limits<unsigned int>::max();
    _variance[2] = std::numeric_limits<unsigned int>::max();
    
    // Tick @ which waveform starts
    int start = 0;

    // start & end tick for each Region Of Interest (ROI) saved
    tick s;
    tick e;
    std::pair<tick,tick> thisRange; 
    
    _pl = mode;

    if (_debug) { std::cout << "\t algo operating on wf of size " << waveform.size() << std::endl; }

    for (size_t n = 0; n < waveform.size(); n++) {

      _thisTick = _begin + n;

      // have we reached the end of a segment?
      if (n % (_block - 1) == 0){

	if (_debug) { std::cout << "\t\t reached end of segment @ tick " << n << std::endl; }

	// is this the 1st block? if so calculate mean and variance
	if ( (n + 1 - _block) == 0){
	  int baseline = 0;
	  int var  = 0;
	  int diff = 0;
	  tick t = _thisTick - _block + 1;
	  int mm = 0;
	  for (; t < _thisTick + 1; t++){
	    baseline += *t;
	    mm += 1;
	  }
	  baseline /= _block;
	  t = _thisTick - _block + 1;
	  for (; t < _thisTick + 1; t++){
	    diff = ( (*t) - baseline ) * ( (*t) - baseline );
	    if (diff < _maxADC) var += diff;
	    else var += _maxADC;
	  }
	  var = var >> 6;
	  if (_debug){
	    std::cout << "\t\t B for block 1 : " << baseline << std::endl;
	    std::cout << "\t\t V for block 1 : " << var      << std::endl;
	  }
	  _baseline[1] = baseline;
	  _variance[1] = var;

	  baseline = 0;
	  var      = 0;
	  diff     = 0;
	  t = _thisTick + 1;
	  for (; t < _thisTick + _block + 1; t++)
	    baseline += *t;
	  baseline = baseline >> 6;
	  t = _thisTick + 1;
	  int nn = 0;
	  for (; t < _thisTick + _block + 1; t++){
	    diff = ( (*t) - baseline ) * ( (*t) - baseline );
	    if (diff < _maxADC) var += diff;
	    else var += _maxADC;
	    nn += 1;
	  }
	  var = var >> 6;
	  if (_debug){
	    std::cout << "\t\t B for block 2 : " << baseline << std::endl;
	    std::cout << "\t\t V for block 2 : " << var      << std::endl;
	  }
	  _baseline[2] = baseline;
	  _variance[2] = var;
	}// 1st block updating

	// always compute baseline and variance for next block, if it exists
	if ( (n + _block) < waveform.size() ){
	  int baseline = 0;
	  int var      = 0;
	  int diff     = 0;
	  tick t = _thisTick + 1;
	  for (; t < _thisTick + _block + 1; t++)
	    baseline += *t;
	  baseline = baseline >> 6;
	  t = _thisTick + 1;
	  for (; t < _thisTick + _block + 1; t++){
	    diff = ( (*t) - baseline ) * ( (*t) - baseline );
	    if (diff < _maxADC) var += diff;
	    else var += _maxADC;
	  }
	  var = var >> 6;

	  if (_debug){
	    std::cout << "\t\t B for block 3 : " << baseline << std::endl;
	    std::cout << "\t\t V for block 3 : " << var      << std::endl;
	  }

	  // update baselines and variances

	  // shift baselines and variances by 1 to accomodate
	  // for the new value from the last block
	  _baseline[0] = _baseline[1];
	  _variance[0] = _variance[1];
	  _baseline[1] = _baseline[2];
	  _variance[1] = _variance[2];
	  // add the newly calculated value
	  // to the last element
	  _baseline[2] = baseline;
	  _variance[2] = var;

	}// if we are updating the NEXT block

	if (_debug){
	  std::cout << "Baseline. Block 1: " << _baseline[0] << "\tBlock 2: " << _baseline[1] << "\tBlock 3: " << _baseline[2] << std::endl;
	  std::cout << "Variance. Block 1: " << _variance[0] << "\tBlock 2: " << _variance[1] << "\tBlock 3: " << _variance[2] << std::endl;
	}


	// Determine if the 3 blocks are quiet enough to update the baseline
	if ( ( (_baseline[2] - _baseline[1]) * (_baseline[2] - _baseline[1]) < _deltaB ) && 
	     ( (_baseline[2] - _baseline[0]) * (_baseline[2] - _baseline[0]) < _deltaB ) && 
	     ( (_baseline[1] - _baseline[0]) * (_baseline[1] - _baseline[0]) < _deltaB ) &&
	     ( (_variance[2] - _variance[1]) * (_variance[2] - _variance[1]) < _deltaV ) &&
	     ( (_variance[2] - _variance[0]) * (_variance[2] - _variance[0]) < _deltaV ) &&
	     ( (_variance[1] - _variance[0]) * (_variance[1] - _variance[0]) < _deltaV ) ){
	  _baselineMap[ch] = _baseline[1];
	  if (_debug) std::cout << "Baseline updated to value " << _baselineMap[ch] << std::endl;
	}
      }// if we hit the end of a new block

      _v1 = _variance[0];
      _v2 = _variance[1];
      _v3 = _variance[2];
      _b1 = _baseline[0];
      _b2 = _baseline[1];
      _b3 = _baseline[2];

      
      if ( (_baselineMap.find(ch) != _baselineMap.end()) ){
	
	double base = _baselineMap[ch];
	
	// reset maxima
	_max = 0;

	// Then go through the 3 blocks again trying to find a waveform to save
	
	// save will keep track of tick at which waveform goes above threshold
	// == 0 if not -> use as logic method to decide if to push back or not
	_save = 0;
	// also keep track of each new sub-waveform to push back to output waveform
	outputwf.clear();
	
	// loop over central block to search for above-threshold regions
	

	//for (t = _thisTick + _block; t < _thisTick + 2 * _block; t++){
	
	tick t = _thisTick;
	
	double thisADC = *t;
	if (thisADC-base > _max) { _max = thisADC-base; }
	
	if ( PassThreshold(thisADC, base) ){
	  if (_verbose) { std::cout << "+ "; }
	  _save = 1;
	  // yay -> active
	  // if start == 0 it means it's a new pulse! (previous tick was quiet)
	  // keep track of maxima
	  if ( start == 0 ){
	    start = int(t-_begin);
	    // also, since we just started...add "backwards ticks" to account for padding
	    if ( (t-_buffer[mode][0]) > _begin) { s = t-_buffer[mode][0]; }
	    else { s = _begin; }
	    if (_verbose) { std::cout << "found start-tick " << s-_begin << std::endl; }
	  }
	  // add bin content to temporary output waveform
	  outputwf.push_back(*t);
	}
	
	else{
	  // we are in a sub-threshold region.
	  // 2 possibilities:
	  // 1) we were in a sub-threshold region at the previous tick -> then just carry on
	  // 2) we were in an active region in the previous tick -> we just "finished" this mini-waveform.
	  //    then Complete padding and save to output
	  if ( start > 0 ){
	    // finish padding
	    if ( (t+_buffer[mode][1]) < _end) { e = t+_buffer[mode][1]; }
	    else { e = _end; }
	    // push back waveform and time-tick
	    if (_verbose) {
	      std::cout << std::endl;
	      std::cout << "saving [" << s-_begin << ", " << e-_begin << "]" << std::endl;
	    }
	    // if the beginning is before the end of the previous section -> just expand
	    if (_timeRange.size() > 0){
	      if (s < _timeRange.back().second){
		// this new range starts before the last one ends -> edit the last one
		_timeRange.back().second = e;
	      }
	      else{
		thisRange = std::make_pair(s,e);
		_timeRange.push_back(thisRange);
	      }
	    }
	    else{
	      thisRange = std::make_pair(s,e);
	      _timeRange.push_back(thisRange);
	    }
	    _save = 0;
	    start = 0;
	  }
	}
      }

    }// for all ticks

    return;
  }


  bool MicrobooneFirmware::PassThreshold(double thisADC, double base){

    if (_pol[_pl] == 0){ //unipolar setting set at command line

        //if positive threshold
	  if (_thresh[_pl] >= 0){
          if (thisADC > base + _thresh[_pl])
    	return true;
       }

	// if negative threshold
        else{
          if (thisADC < base + _thresh[_pl])
    	return true;
        }

	  }

    else { //bipolar setting set at command line
	  if  (thisADC >= base + _thresh[_pl]) {
	    return true;
	  }

	  if (thisADC <= base - _thresh[_pl]) {
	      return true;
	    }
	  }

    return false;
  }
  
}
#endif
