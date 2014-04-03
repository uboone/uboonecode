#include <stdlib.h>
#include <stdio.h>

#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <iostream>
int Error; 

extern void InitGui(); // Initializer for GUI
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
#include "art/Framework/Core/RootDictionaryManager.h"
#include "art/Framework/Art/InitRootHandlers.h"

// Initialize the ROOT system
TROOT root("Rint", "The ROOT Interactive Interface", initfuncs);

int main(int argc, char **argv)
{
  std::cout << "Loading libraries." << std::endl;
  art::RootDictionaryManager rdm;
  art::completeRootHandlers();
  std::cout << "...done" << std::endl;
  
  TRint *theApp = new TRint("ROOT example", &argc, argv, NULL, 0, 0);

//  gStyle->SetPalette(1,0);
  
  // Run interactive interface
  theApp->Run();
  
 
  return(0);
}

