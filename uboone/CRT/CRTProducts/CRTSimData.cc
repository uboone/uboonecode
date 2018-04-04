#include "uboone/CRT/CRTProducts/CRTSimData.hh"

namespace crt{

  CRTSimData::CRTSimData(): fChannel(0), fT0(0), fT1(0){
  }
  CRTSimData::CRTSimData(uint32_t channel, uint32_t t0, 
    uint32_t t1, uint32_t adc):
    fChannel(channel),
    fT0(t0),
    fT1(t1),
    fADC(adc) {
    }
  CRTSimData::~CRTSimData(){
  }
  uint32_t CRTSimData::Channel(){
    return this->fChannel;
  }
  uint32_t CRTSimData::T0(){
    return this->fT0;
  }
  uint32_t CRTSimData::T1(){
    return this->fT1;
  }
  uint32_t CRTSimData::ADC(){
    return this->fADC;
  }
}
