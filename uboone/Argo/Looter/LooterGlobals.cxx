#include <string>
#include "LooterGlobals.h"

static std::string gLooterOutput = "";
std::string& looterOutput()
{
  return gLooterOutput;
}


