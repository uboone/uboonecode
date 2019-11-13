/**
 * \file TruncMean.h
 *
 * \ingroup 3DMichel
 *
 * \brief Class def header for a class TruncMean
 *
 * @author david caratelli
 * @author erez cohen
 */

/** \addtogroup 3DMichel
 
 @{*/
#ifndef TRUNCMEAN_H
#define TRUNCMEAN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Rtypes.h"
using namespace std;


/**
 \class TruncMean
 User defined class TruncMean ... these comments are used to generate
 doxygen documentation!
 */

class TruncMean{
    
public:
    
    /// constructors
    TruncMean () = default;
    TruncMean (int fdebug, float frad=10.0);
    
    
    /**
     @brief Given residual range and dq vectors return truncated local dq.
     Input vectors are assumed to be match pair-wise
     (nth entry in rr_v corresponds to nth entry in dq_v vector).
     Input rr_v values are also assumed to be ordered: monotonically increasing
     or decreasing.
     For every dq value a truncated linear dq value is calculated as follows:
     0) all dq values within a rr range set by the class variable _rad are selected.
     1) the median and rms of these values is calculated.
     2) the subset of local dq values within the range [median-rms, median+rms] is selected.
     3) the resulting local truncated dq is the average of this truncated subset.
     */
    void        CalcTruncMean (const std::vector<float>& rr_v,
                               const std::vector<float>& dq_v,
                               std::vector<float>& dq_trunc_v);
    
    void CalcTruncMeanProfile (const std::vector< float > & rr_v,
                               const std::vector< float > & dq_v,
                               std::vector< float > & dq_trunc_v,
                               const float & nsigma = 1
                               );
    
    
    
    void              SetRad (float frad) {rad = frad;};
    void            SetDebug (int fdebug) {debug = fdebug;};
    
    float             Median (const std::vector<float>& v);
    float                RMS (const std::vector<float>& v);
    
    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    void Debug(Int_t verobosity_level, const char* format) // base function
    {
        if ( debug < verobosity_level ) return;
        std::cout << format << std::endl;
    }
    template<typename T, typename... Targs>
    void Debug(Int_t verobosity_level, const char* format, T value, Targs... Fargs) // recursive variadic function
    {
        if ( debug < verobosity_level ) return;
        for ( ; *format != '\0'; format++ ) {
            if ( *format == '%' ) {
                std::cout << value << " ";
                Debug(verobosity_level, format+1, Fargs...); // recursive call
                return;
            }
            std::cout << *format;
        }
        std::cout << endl;
    }
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -

    
    
    int     debug=0;
    double  rad=10.;
    
};

#endif
/** @} */ // end of doxygen group

