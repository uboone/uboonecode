#ifndef CRTFileManager_h
#define CRTFileManager_h
////////////////////////////////////////////////////////////////////////
// Class:       CRTFileManager
// Plugin Type: service (art v2_05_01)
// File:        CRTFileManager.h
//
// Generated at Thu Jan  3 14:09:49 2019 by Herbert Greenlee using cetskelgen
// from cetlib version v1_21_00.
//
// Purpose:  This service manages multiple instances of open gallery files
//           for CRT swizzled data.  It can perform the following tasks.
//
//           1.  Locating swizzled CRT files that match an art event time stamp.
//
//           2.  Opening swizzled CRT files (using xrootd), and positioning them
//               at the correct position for matching.
//
//           Swizzled CRT files are maintained as open files accross multiple
//           art modules and art events.
//
// FCL paramters:
//
// Debug      - Debug flag.
// MaxFiles   - Maximum number of open gallery files to cache.
// ubversion_CRTHits - Swizzled CRT version.
// ubversion_CRTHits_top - Swizzled CRT version for top panel CRT hits (stream 1).
//                         Optional, default same as ubversion_CRTHits.
// CRTMetadataDir - Directory containing extracted CRT metadata.  If defined,
//                  CRT metadata will be looked up in this directory instead
//                  of being queried from SAM.
// 
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <list>
#include <map>
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "gallery/Event.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "uboone/CRT/CRTProducts/CRTHit.hh"

namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}

namespace crt {
  class CRTFileManager;
}


class crt::CRTFileManager {
public:

  // Nested struct with information about CRT binary files.

  struct CRTFileInfo
  {
    std::string fFileName;                 // Name of CRT binary file.
    boost::posix_time::ptime fStartTime;   // Start time of CRT binary file.
    boost::posix_time::ptime fEndTime;     // End time of CRT binary file.
    std::vector<std::string> fSwizzled;    // Names of CRT swizzled files.
  };

  explicit CRTFileManager(fhicl::ParameterSet const & p, art::ActivityRegistry & areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Prefetch extracted CRT metadata into metadata cache.

  void prefetch(const std::vector<boost::posix_time::ptime>& posix_times) const;

  // Find matching swizzled CRT files based on art time stamp.
  // Use GPS time from DAQHeaderTimeUBooNE data product.

  std::vector<std::string> findMatchingCRTFiles(art::Timestamp event_time) const;

  // Open a gallery file.
  // Return value is a reference to a gallery event that is owned by this service.
  // Use function isValid to test whether returned file was opened successfully.
  // In case of a successful open, the newly opened file is added to the list of 
  // cached open gallery files.
  // The oldest open gallery file may be closed if the number of open gallery
  // files would exceed the maximum cache size.

  gallery::Event& openFile(std::string file_name);

  // Rewind gallery file.

  bool rewind(gallery::Event&) const;

  // Get next nonempty filtered and sorted CRT hit collection starting from 
  // current gallery event.
  // Return empty collection if end of file is reached.

  void get_crt_hits(gallery::Event&, std::vector<crt::CRTHit>& output_hits) const;

private:

  // FCL parameters.

  bool fDebug;               // Debug mode.
  unsigned int fMaxFiles;    // Maximum number of open files to cache.
  std::string fCRTHitLabel;  // CRT hit module label.
  std::string fCRTVersion;   // Swizzled CRT version (ub_project.version).
  std::string fCRTVersionTop;// Swizzled CRT version for top CRT panels (ub_project.version).
  std::string fCRTMetadataDir; // Directoroy containing extracted CRT metadata.

  // This data structure contains information that is known from sam metadata about
  // CRT binary and swizzled files.
  // The map key is the name of a CRT swizzled file.
  // The map value contains information about the corresponding CRT binary file.

  mutable std::map<std::string, CRTFileInfo> fCRTFiles;

  // Open gallery file cache.
  // The list is sorted in order of decreasing age, most recently opened file at end.

  std::list<std::pair<std::string, std::unique_ptr<gallery::Event> > > fCRTEvents;
};

DECLARE_ART_SERVICE(crt::CRTFileManager, LEGACY)
#endif /* CRTFileManager_h */
