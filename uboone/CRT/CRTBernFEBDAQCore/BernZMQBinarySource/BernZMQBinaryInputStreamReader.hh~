#include "canvas/Persistency/Provenance/Timestamp.h"

#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "bernfebdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"

#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <memory>

#include "fhiclcpp/ParameterSet.h"

namespace bernfebdaq{
  typedef std::map< uint64_t,
		    std::unique_ptr<std::vector<artdaq::Fragment> > > FragMap_t;
  class BernZMQBinaryInputStreamReader;
}

class bernfebdaq::BernZMQBinaryInputStreamReader {

public:
  BernZMQBinaryInputStreamReader( fhicl::ParameterSet const& );
  BernZMQBinaryInputStreamReader( fhicl::ParameterSet const& p, std::ifstream& is)
    : BernZMQBinaryInputStreamReader(p)
  { SetInputStream(is); }

  void SetInputStream (std::ifstream& is ) { fInputStream = &is; }

  art::Timestamp ReadUntilSpecialEvent(FragMap_t &);
  
private:

  std::ifstream*      fInputStream;

  typedef std::vector< std::vector<BernZMQEvent> > FEBBuffer_t;
  typedef std::unordered_map< uint32_t,uint32_t >  FEBStuckEventMap_t;

  //configs
  size_t       fNFEBBuffers;
  int          fVerbosity;            
  unsigned int fMaxTSDiff;
  unsigned int fMaxStuckEvents;
  uint32_t     fPullTimeLatency; //ms
  uint32_t     fFragsPerSecond;

  uint32_t fNanosecondsPerFragment;
  
  std::unordered_map< uint64_t, FEBBuffer_t  >         fFEBBuffers;
  std::unordered_map< uint64_t, uint32_t  >            fPullLastEvent;
  std::unordered_map< uint64_t, uint32_t  >            fLastTS1Ref;
  std::unordered_map< uint64_t, size_t  >              fFEBCurrentBuffer;
  std::unordered_map< uint64_t, FEBStuckEventMap_t  >  fFEBStuckEventMaps;

  //from the slf machines
  struct mytimeb_t
  {
    long int time;                /* Seconds since epoch, as from `time'.  */
    unsigned short int millitm; /* Additional milliseconds.  */
    short int timezone;         /* Minutes west of GMT.  */
    short int dstflag;          /* Nonzero if Daylight Savings Time used.  */
  };

  art::Timestamp GetNextPullTime();
  art::Timestamp ConvertMyTimebType( struct mytimeb_t const& );
  
  art::Timestamp fLastPullTime;
  
  art::Timestamp SpecialEventTimestamp( BernZMQEvent const& );

  size_t fEventsPerPullCounter;

  void InitializeBuffer(uint16_t mac);
  bool KeepCurrentEvent(std::vector<BernZMQEvent>&,BernZMQEvent const&,bool&);

  void ClearNextBuffer(BernZMQEvent const& zmq_ev);
  void ClearNextBufferAndInsert(BernZMQEvent const& zmq_ev);

  void FillFragments(std::vector<BernZMQEvent> const&,FragMap_t &);

};
