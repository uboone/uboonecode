#include "uboone/CRT/CRTAuxFunctions.hh"
#include "cetlib/exception.h"

int crt::auxfunctions::getFEBN(uint64_t febid){
  
  int FEBN = (febid  - 413240800282 + 26)%100; //FEB ID                                                                                                          
  return FEBN;
}


void crt::auxfunctions::FillPos(std::string filePos,   std::map <int, std::vector<double> >& sensor_pos){ //key = FEB*100+ch                
  
  std::cout<<"Reading SiPMs position information"<<std::endl;

  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(filePos,fname);
  
  std::ifstream in;
  in.open(fname.c_str());
  std::cout <<"File "<< fname << "open? "<<in.is_open()<<std::endl;



  if(in.is_open()){
    std::cout<<"File open: "<<fname.c_str()<<std::endl;
  }

  while (!in.eof()){

    double x,y,z;
    int id, pl, subpl, ori;
    std::vector<double> pos;
    
    in>>id>>x>>y>>z>>pl>>subpl>>ori;
    
    pos.push_back(x);
    pos.push_back(y);
    pos.push_back(z);
    pos.push_back(pl);
    pos.push_back(subpl);
    pos.push_back(ori);
    
    sensor_pos[id] = pos;
    
  }
}

std::vector<double> crt::auxfunctions::getPos(int ID,  std::map <int, std::vector<double> >& sensor_pos){
  
  auto it2 = sensor_pos.find(ID);
  std::vector<double> pos = (*it2).second;

  return pos;
}

void crt::auxfunctions::FillFEBDel(std::string fileFEBDel,   std::map <int,double> & FEBDel){ //key = FEB                                                
  
  std::cout<<"Reading FEB delay information"<<std::endl;
 
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fileFEBDel,fname);

  std::ifstream in;
  in.open(fname.c_str());
  std::cout <<"File open? "<<in.is_open()<<std::endl;

  if(in.is_open()){
    std::cout<<"File open: "<<fname.c_str()<<std::endl;
  }

  while (!in.eof()){

    double del;
    int id;

    in>>id>>del;
    FEBDel[id] = del;
    
  }
}

double crt::auxfunctions::getFEBDel(int ID,  std::map <int,double >& FEBDel){

  auto it2 = FEBDel.find(ID);
  double del = (*it2).second;

  return del;
}


void crt::auxfunctions::FillGain(std::string fileGain,   std::map <int, std::pair<double,double> >& sensor_gain){ //key = FEB*100+ch            
  
  std::cout<<"Reading SiPMs gain information"<<std::endl;

  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fileGain,fname);
  
  std::ifstream in;
  in.open(fname.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;
  if(in.is_open()){
    std::cout<<"File open: "<<fname.c_str()<<std::endl;
  }
  
  while (!in.eof()){
    
    double g, g_err;
    int id;
    std::pair<double,double> gain;
    
    in>>id>>g>>g_err;
    gain = std::make_pair(g,g_err);
    sensor_gain[id] = gain;
      
  }
}


std::pair<double,double> crt::auxfunctions::getGain(int ID,  std::map <int, std::pair<double,double> >& sensor_gain){

  auto it2 = sensor_gain.find(ID);
  std::pair<double,double> gain = (*it2).second;
  return gain;
}


std::vector<double> crt::auxfunctions::inter_X(double S1, std::vector<double>& pS1, double S2, std::vector<double>& pS2, double L){
  
  std::vector<double> interpos;

  if(S1<0) S1 = 0.000000001;
  if(S2<0) S2 = 0.000000001;
  
  
  //ret=STRIPW/2.*atan(log(1.*sR/sL)); 
  //double ret = (L/2)*atan(log(1.*S2/S1)); //empiric. CHECK with MC
  
  double x =  ((S1*pS1[0] + S2*pS2[0])/(S1+S2)); //barycenter
  double y =  ((S1*pS1[1] + S2*pS2[1])/(S1+S2)); //barycenter
  double z =  ((S1*pS1[2] + S2*pS2[2])/(S1+S2)); //barycenter
  
  double pl = (pS1[3] + pS2[3])/2;
  double ly = (pS1[4] + pS2[4])/2;
  double ori = (pS1[5] + pS2[5])/2;
  double inter = -1;

  
  if(pS1[0] != pS2[0]){ inter=1;}
  // x = pS1[0] + ret;}  //interpol x
  if(pS1[1] != pS2[1]){ inter=2;}
  // y = pS1[1] +ret;}  //interpol y
  if(pS1[2] != pS2[2]){ inter=3;}
  // z = pS1[2] + ret;}  //interpol z
  
  interpos.push_back(x);
  interpos.push_back(y);
  interpos.push_back(z);
  interpos.push_back(pl);
  interpos.push_back(ly);
  interpos.push_back(ori);
  interpos.push_back(inter);

  return interpos;
}

double crt::auxfunctions::inter_X_error(double S1, double S2, double L){

  if(S1<0) S1 = 0.000000001;
  if(S2<0) S2 = 0.000000001;
  
  double dy = 0.344677; //pol1
  double d0y = -1.92045;

  double x =  ((S1*0 + S2*L)/(S1+S2)); //barycenter
  x =  x + dy * x + d0y; //barycenter Norm pol1

  double x_error =  1.92380e+00+1.47186e-02*x-5.29446e-03*x*x;
  
  if(x_error<1.45) x_error=1.45;
  if(x_error>1.99) x_error=1.99;

  return x_error;
}


double crt::auxfunctions::getTcorr(std::vector<double>& inpos1, std::vector<double>& inpos2 ,double T){                                                              
  double L = 0;                                                                                                                                              
  if(inpos1[3]==0 || inpos1[3]==3){//bottom&top (XZ)                                                                                                        
    if(inpos1[6]==1){//correction in z                                                                                                                 
      
      L = abs(inpos1[2]-inpos2[2]);                                                                                                                          
    }                                                                                                                                                            
    if(inpos1[6]==3){//correction in x                                                                                                                       
      L = abs(inpos1[0]-inpos2[0]);                                                                                                                          
    }                                                                                                                                                        
  }                                                                                                                                                            
  if(inpos1[3]==1 || inpos1[3]==2){//FT&pipe (YZ)                                                                                                               
    if(inpos1[6]==2){//correction in z                                                                                                                       
      L = abs(inpos1[2]-inpos2[2]);                                                                                                                          
    }                                                                                                                                                            
    if(inpos1[6]==3){//correction in y                                                                                                                       
      L = abs(inpos1[1]-inpos2[1]);                                                                                                                          
    }                                                                                                                                                            
  }                                                                                                                                                         
  double Tcorr = T - (L*6.2/100);                                                                                                                            
  return Tcorr; 
}

void crt::auxfunctions::FillPartTop(std::string fileTop, int mac_buffer[3][100]){
  
  std::cout<<"Reading Part Mac of Top information"<<std::endl;

  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fileTop,fname);
  
  std::ifstream in;
  in.open(fname.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;
  if(in.is_open()){
    std::cout<<"File open: "<<fname.c_str()<<std::endl;
  }
  int counter=0;
  int febnr1;
  int febnr2;
  int febnr3;
  while (!in.eof()){
    
    in>>febnr1>>febnr2>>febnr3;
    mac_buffer[0][counter]=febnr1;
    mac_buffer[1][counter]=febnr2;
    mac_buffer[2][counter]=febnr3;
    counter++;
      
  }
  mac_buffer[0][99]=counter;
  mac_buffer[1][99]=counter;
  mac_buffer[2][99]=counter;
}
