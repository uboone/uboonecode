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

void crt::auxfunctions::Init_TS0_corr(std::string filename,crt::TS0_CORRECTION correctionpoints[50] , uint32_t start_s,uint32_t end_s){
  for(int i=0;i<50;i++){
    correctionpoints[i].sec=0;
    correctionpoints[i].offset=0;
    correctionpoints[i].scale=0;
  }
  
  uint32_t sec;
  double offset;
  double scale;
  int counter=0;
  
  std::cout<<"Reading GPS correstion for ts0 information"<<std::endl;

  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(filename,fname);
  
  std::ifstream in;
  in.open(fname.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;
  if(in.is_open()){
    std::cout<<"File open: "<<fname.c_str()<<std::endl;
  }
    while (!in.eof()) {
      in>>sec>>offset>>scale;
      if(sec>=(start_s-10*3600) && sec<(end_s+10*3600)){
        correctionpoints[counter].sec=sec;
        correctionpoints[counter].offset=offset;
        correctionpoints[counter].scale=scale;
        counter++;
      }
    }      
    in.close();
  std::cout << "Found: " << counter << " points for the GPS correction (matched trigger)." << std::endl;
  for(int i=0; i<counter;i++){
    if(correctionpoints[i].sec>=(start_s) && correctionpoints[i].sec<(end_s)) printf("sec: %d, offset: %lf, scale: %lf\n",correctionpoints[i].sec, correctionpoints[i].offset, correctionpoints[i].scale);
  }
}

void crt::auxfunctions::Init_mspoll_delay(std::string file_FEB_MS_delay, double Ms_delay[200]){
   std::cout<<"Reading MS delay for polls information"<<std::endl;

  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(file_FEB_MS_delay,fname);
  
  std::ifstream in;
  in.open(fname.c_str());

  std::cout <<"File open? "<<in.is_open()<<std::endl;
  if(in.is_open()){
    std::cout<<"File open: "<<fname.c_str()<<std::endl;
  }
  else
   {
    std::cout << "Error opening poll ms delay file";
   }
    int febnr;
    double factor;
    while (!in.eof()) {
      in>>febnr>>factor;
      Ms_delay[febnr]=factor;
      //std::cout << "FEB: " << febnr << " with factor: " << factor << std::endl;
    }      
    in.close();
   
}

double crt::auxfunctions::CRT_Only_Offset(uint32_t sec){
   double offset=0;
  offset=(3.43953e+09+634.757*(sec-1496e6));                                          //f1
  if(offset>=0 && offset<1e9) return offset;
  offset=(2.60674e+09+717.633*(sec-1496e6)+1.02278e-05*(sec-1496e6)*(sec-1496e6));    //f2
  if(offset>=0 && offset<1e9) return offset;
  offset=(1.63414e+09+741.45*(sec-1496e6)+1.57055e-05*(sec-1496e6)*(sec-1496e6));     //f3
  if(offset>=0 && offset<1e9) return offset;
  offset=(1.52973e+09+773.749*(sec-1496e6)+2.14305e-05*(sec-1496e6)*(sec-1496e6));    //f4
  if(offset>=0 && offset<1e9) return offset;
  offset=(5.14068e+08+746.236*(sec-1496e6)+1.18375e-05*(sec-1496e6)*(sec-1496e6));    //f5
  if(offset>=0 && offset<1e9) return offset;
  offset=(-4.85346e+08+743.584*(sec-1496e6)+1.57855e-05*(sec-1496e6)*(sec-1496e6));   //f6
  if(offset>=0 && offset<1e9) return offset;
  offset=(-1.54728e+09+803.891*(sec-1496e6)+7.86822e-07*(sec-1496e6)*(sec-1496e6));   //f7
  if(offset>=0 && offset<1e9) return offset;
  offset=(-2.49471e+09+755.792*(sec-1496e6)+1.09478e-05*(sec-1496e6)*(sec-1496e6));   //f8
  if(offset>=0 && offset<1e9) return offset;
  offset=(-3.34555e+09+692.4*(sec-1496e6)+1.80768e-05*(sec-1496e6)*(sec-1496e6));     //f9
  if(offset>=0 && offset<1e9) return offset;
  offset=(-1.27155e+09+833.209*(sec-1500e6)+1.43991e-05*(sec-1500e6)*(sec-1500e6));   //f10
  if(offset>=0 && offset<1e9) return offset;
  offset=(-2.32269e+09+871.886*(sec-1500e6)+7.55479e-06*(sec-1500e6)*(sec-1500e6));   //f11
  if(offset>=0 && offset<1e9) return offset;
  offset=(-2.83569e+09+617.373*(sec-1500e6)+4.10254e-05*(sec-1500e6)*(sec-1500e6));   //f12
  if(offset>=0 && offset<1e9) return offset;
  offset=(-3.19383e+09+945.459*(sec-1500e6)+0*(sec-1500e6)*(sec-1500e6));             //f13
  if(offset>=0 && offset<1e9) return offset;
  offset=(-4.14647e+09+910.346*(sec-1500e6)+5.5141e-06*(sec-1500e6)*(sec-1500e6));    //f14
  if(offset>=0 && offset<1e9) return offset;
  offset=(-7.83914e+09+1020.95*(sec-1500e6)+-1.00822e-06*(sec-1500e6)*(sec-1500e6));  //f15
  if(offset>=0 && offset<1e9) return offset;
  offset=(-9.35939e+09+1108.93*(sec-1500e6)-4.25095e-06*(sec-1500e6)*(sec-1500e6));   //f16
  if(offset>=0 && offset<1e9) return offset;
  offset=(-8.50541e+09+728.814*(sec-1500e6)+1.52889e-05*(sec-1500e6)*(sec-1500e6));   //f17
  if(offset>=0 && offset<1e9) return offset;
  offset=(-1.09097e+10+982.475*(sec-1500e6)+3.88011e-06*(sec-1500e6)*(sec-1500e6));   //f18
  if(offset>=0 && offset<1e9) return offset;
  offset=(-4.78612e+09-228.385*(sec-1500e6)+5.53339e-05*(sec-1500e6)*(sec-1500e6));   //f19
  if(offset>=0 && offset<1e9) return offset;
  
  return -1;
}











