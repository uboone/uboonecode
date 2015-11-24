////////////////////////////////////////////////////////////////////////
//
//  UBDatabaseUtil_plugin
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <iostream>
#include <fstream>

// LArSoft includes
#include "UBDatabaseUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

//-----------------------------------------------
util::UBDatabaseUtil::UBDatabaseUtil(fhicl::ParameterSet const& pset, art::ActivityRegistry & /* reg */)
{
  conn = NULL;
  this->reconfigure(pset);
  fChannelMap.clear();
  fChannelReverseMap.clear();
}

//------------------------------------------------
util::UBDatabaseUtil::~UBDatabaseUtil() 
{

}

//----------------------------------------------
int util::UBDatabaseUtil::Connect(int conn_wait)
{
  if(!fShouldConnect)
    return -1;
  
  if(conn_wait)
    sleep(conn_wait);
    
  conn = PQconnectdb(connection_str);
  if (PQstatus(conn) == CONNECTION_BAD) {
    mf::LogWarning("UBDatabaseUtil") << "Connection to database failed, "<<PQerrorMessage(conn)<<"\n";
    if( ( strstr(PQerrorMessage(conn),"remaining connection slots are reserved")!=NULL || 
	  strstr(PQerrorMessage(conn),"sorry, too many clients already")!=NULL )
	&& conn_wait<20 ) {
      conn_wait+=2;
      mf::LogWarning("UBDatabaseUtil") << "retrying connection after " << conn_wait << " seconds \n";
      return this->Connect(conn_wait);
    }
    if(fToughErrorTreatment)
      throw cet::exception("DataBaseUtil") << " DB connection failed\n";
   
  } else {
    LOG_DEBUG("UBDatabaseUtil")<<"Connected OK\n";
    return 1;
  }
  return -1;
}


int util::UBDatabaseUtil::DisConnect()
{
  if(!fShouldConnect)
    return -1;
  //close connection
  LOG_DEBUG("UBDatabaseUtil")<<"Closing Connection \n";
  PQfinish(conn);
  return 1;
}



//------------------------------------------------
void util::UBDatabaseUtil::reconfigure(fhicl::ParameterSet const& pset)
{
  fDBHostName            = pset.get< std::string >("DBHostName"          );
  fDBName      		 = pset.get< std::string >("DBName"     );
  fDBUser 		 = pset.get< std::string >("DBUser");
  fTableName 		 = pset.get< std::string >("TableName");
  fPort      		 = pset.get< int >("Port"     );
  fPassword 		 = " ";
  fToughErrorTreatment   = pset.get< bool >("ToughErrorTreatment");
  fShouldConnect   	 = pset.get< bool >("ShouldConnect");
  
  // constructor decides if initialized value is a path or an environment variable
  std::string passfname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(pset.get< std::string >("PassFileName"), passfname);
    
  std::ifstream in(passfname.c_str());
  if(in.is_open())   {
    std::getline(in, fPassword);
    in.close();
  }
   
  sprintf(connection_str,"host=%s dbname=%s user=%s port=%d password=%s",fDBHostName.c_str(),fDBName.c_str(),fDBUser.c_str(),fPort,fPassword.c_str());
   
  return;
}




int util::UBDatabaseUtil::SelectSingleFieldByQuery(std::vector<std::string> &value,const char * query)
{
  PGresult *result;
  char * string_val;
  
  if(this->Connect()==-1)  {
    if(fShouldConnect)
      mf::LogWarning("UBDatabaseUtil")<< "DB Connection error \n";
    else
      mf::LogInfo("UBDatabaseUtil")<< "Not connecting to DB by choice. \n";
    return -1;
  }
    
  result = PQexec(conn, query);

  if (!result) {
    mf::LogInfo("UBDatabaseUtil")<< "PQexec command failed, no error code\n";
    return -1;
  } 
  else if(PQresultStatus(result)!=PGRES_TUPLES_OK) {
    if(PQresultStatus(result)==PGRES_COMMAND_OK) 
      LOG_DEBUG("UBDatabaseUtil")<<"Command executed OK, "<< PQcmdTuples(result) <<" rows affected\n";
    else
      mf::LogWarning("UBDatabaseUtil")<<"Command failed with code "
				    <<PQresStatus(PQresultStatus(result)) <<", error message "
				    <<PQresultErrorMessage(result)<<"\n";
		 
    PQclear(result);	
    this->DisConnect();
    return -1;
  }
  else {
    //  mf::LogInfo("UBDatabaseUtil")<<"Query may have returned data\n";
    //  mf::LogInfo("UBDatabaseUtil")<<"Number of rows returned: "<<PQntuples(result)
    //   <<", fields: "<<PQnfields(result)<<" \n";
    
    if(PQntuples(result)>=1){
      for(int i=0;i<PQntuples(result);i++)
	{
	  string_val=PQgetvalue(result,i,0);
	  value.push_back(string_val);
	  LOG_DEBUG("UBDatabaseUtil")<<" extracted value: "<<value[i] << "\n";
	}
      PQclear(result);
      this->DisConnect();
      return 0;
    }
    else {
      mf::LogWarning("UBDatabaseUtil")<<"wrong number of rows returned:"<<PQntuples(result)<<"\n";
      PQclear(result);
      this->DisConnect();
      return -1;
    }	
  }
  
  
}



int util::UBDatabaseUtil::GetTemperatureFromDB(int run,double &temp_real)
{
  std::vector<std::string> retvalue;
  char cond[30];
  sprintf(cond,"run = %d",run);
  int err=SelectFieldByName(retvalue,"temp",cond,fTableName.c_str());
  
  if(err!=-1 && retvalue.size()==1){
    char * endstr;
    temp_real=std::strtod(retvalue[0].c_str(),&endstr); 
    return 0; 
  }

  return -1;
 
  
}




int util::UBDatabaseUtil::GetEfieldValuesFromDB(int run,std::vector<double> &efield)
{
  
  std::vector<std::string> retvalue;
    
  char query[200];
  sprintf(query,"SELECT EFbet FROM EField,%s WHERE Efield.FID = %s.FID AND run = %d ORDER BY planegap",fTableName.c_str(),fTableName.c_str(),run);
  int err=SelectSingleFieldByQuery(retvalue,query);
  
  if(err!=-1 && retvalue.size()>=1){
    efield.clear();    //clear value before setting new values
    for(unsigned int i=0;i<retvalue.size();i++) {
      char * endstr;
      efield.push_back(std::strtod(retvalue[i].c_str(),&endstr)); 
    }
    return 0; 
  }

  return -1;
    
}
      


int util::UBDatabaseUtil::SelectFieldByName(std::vector<std::string> &value,
					  const char * field,
					  const char * condition,
					  const char * table)	{
  
  char query[100];
  sprintf(query,"SELECT %s FROM %s WHERE %s",field, table, condition);
  
  return SelectSingleFieldByQuery(value,query);
    
}
      







int util::UBDatabaseUtil::GetLifetimeFromDB(int run,double &lftime_real) {

  //  char query[100];
  //  sprintf(query,"SELECT tau FROM argoneut_test WHERE run = %d",run);

  std::vector<std::string> retvalue;
  char cond[30];
  sprintf(cond,"run = %d",run);
  int err=SelectFieldByName(retvalue,"tau",cond,fTableName.c_str());
  
  if(err!=-1 && retvalue.size()==1){
    char * endstr;
    lftime_real=std::strtod(retvalue[0].c_str(),&endstr); 
    return 0; 
  }


  return -1;

}

int util::UBDatabaseUtil::GetTriggerOffsetFromDB(int run,double &T0_real) {

  //  char query[100];
  //  sprintf(query,"SELECT tau FROM argoneut_test WHERE run = %d",run);

  std::vector<std::string> retvalue;
  char cond[30];
  sprintf(cond,"run = %d",run);
  int err=SelectFieldByName(retvalue,"T0",cond,fTableName.c_str());
  
  if(err!=-1 && retvalue.size()==1){
    char * endstr;
    T0_real=std::strtod(retvalue[0].c_str(),&endstr); 
    return 0; 
  }


  return -1;

}


int util::UBDatabaseUtil::GetPOTFromDB(int run,long double &POT) {

  //  char query[100];
  //  sprintf(query,"SELECT tau FROM argoneut_test WHERE run = %d",run);

  std::vector<std::string> retvalue;
  char cond[30];
  sprintf(cond,"run = %d",run);
  int err=SelectFieldByName(retvalue,"pot",cond,fTableName.c_str());
  
  if(err!=-1 && retvalue.size()==1){
    char * endstr;
    POT=std::strtold(retvalue[0].c_str(),&endstr); 
    return 0; 
  }


  return -1;

}

namespace util {

  void UBDatabaseUtil::LoadUBChannelMap( int data_taking_timestamp, int  swizzling_timestamp) {

    if ( fChannelMap.size()>0 ) {
      // Use prevously grabbed data to avoid repeated call to database.
      // Also this avoids inglorious segfault. 
      return;
    }

    //Place the query http://dbdata0vm.fnal.gov:8188/QE/uboone/query?F=get_map_double_sec&a=-1
    const char *url = "http://dbdata0vm.fnal.gov:8188/QE/uboone/query?F=get_map_double_sec&a=-1";
    const char *uagent = "";
    int *error = 0;
    Dataset mapping = getData(url, uagent, error);
    Tuple row = getFirstTuple(mapping);
    int eightthousandtwohundredfiftysix = getNfields(row);
    mf::LogError("WHAT THE HECK")<< eightthousandtwohundredfiftysix;


    if ( conn==NULL )
      Connect( 0 );

    if(PQstatus(conn)!=CONNECTION_OK) {
      mf::LogError("") << __PRETTY_FUNCTION__ << ": Couldn't open connection to postgresql interface"  << PQdb(conn) <<":"<<PQhost(conn);
      PQfinish(conn);
      throw art::Exception( art::errors::FileReadError )
        << "Failed to get channel map from DB."<< std::endl;
    }
    
    fChannelMap.clear();
    fChannelReverseMap.clear();

    PGresult *res  = PQexec(conn, "BEGIN");    
    if (PQresultStatus(res) != PGRES_COMMAND_OK) { 
      mf::LogError("")<< "postgresql BEGIN failed";
      PQclear(res);
      PQfinish(conn);
      throw art::Exception( art::errors::FileReadError )
        << "postgresql BEGIN failed." << std::endl;
    }
    
    // Jason St. John's updated call to versioned database.
    // get_map_double_sec (data_taking_timestamp int DEFAULT now() , 
    //                     swizzling_timestamp   int DEFAULT now() )
    // Returns rows of: crate, slot, fem_channel, larsoft_channel 
    // Both arguments are optional, or can be passed their default of now(), or can be passed an explicit timestamp:
    // Example: "SELECT get_map_double_sec(1438430400);"
    PQclear(res);

    char dbquery[1] = {' '};  //I hate C++ strong typing and string handling so very, very much. 
    sprintf(dbquery, "SELECT get_map_double_sec(%i,%i);", data_taking_timestamp, swizzling_timestamp);
    res = PQexec(conn, dbquery); 

    if ((!res) || (PQresultStatus(res) != PGRES_TUPLES_OK) || (PQntuples(res) < 1))
      {
	mf::LogError("")<< "SELECT command did not return tuples properly. \n" << PQresultErrorMessage(res) << "Number rows: "<< PQntuples(res);
        PQclear(res);
        PQfinish(conn);
        throw art::Exception( art::errors::FileReadError )
          << "postgresql SELECT failed." << std::endl;
      } 

    int num_records=PQntuples(res);            //One record per channel, ideally.

    for (int i=0;i<num_records;i++) {
      std::string tup = PQgetvalue(res, i, 0); // (crate,slot,FEMch,larsoft_chan) format
      tup = tup.substr(1,tup.length()-2);      // Strip initial & final parentheses.
      std::vector<std::string> fields;
      split(tup, ',', fields);                 // Explode substrings into vector with comma delimiters. 
      
      int crate_id     = atoi( fields[0].c_str() );
      int slot         = atoi( fields[1].c_str() );
      int boardChan    = atoi( fields[2].c_str() );
      int larsoft_chan = atoi( fields[3].c_str() );
      
      UBDaqID daq_id(crate_id,slot,boardChan);
      std::pair<UBDaqID, UBLArSoftCh_t> p(daq_id,larsoft_chan);

      if ( fChannelMap.find(daq_id) != fChannelMap.end() ){
	std::cout << __PRETTY_FUNCTION__ << ": ";
        std::cout << "Multiple entries!" << std::endl;
        mf::LogWarning("")<< "Multiple DB entries for same (crate,card,channel). "<<std::endl
			  << "Redefining (crate,card,channel)=>id link ("
			  << daq_id.crate<<", "<< daq_id.card<<", "<< daq_id.channel<<")=>"
			  << fChannelMap.find(daq_id)->second;
      }
             
      fChannelMap.insert( p );
      fChannelReverseMap.insert( std::pair< UBLArSoftCh_t, UBDaqID >( larsoft_chan, daq_id ) );
    }
  }// end of LoadUBChannelMap

  UBChannelMap_t UBDatabaseUtil::GetUBChannelMap( int data_taking_timestamp, int swizzling_timestamp ) {
    LoadUBChannelMap( data_taking_timestamp, swizzling_timestamp );
    return fChannelMap;
  }

  UBChannelReverseMap_t UBDatabaseUtil::GetUBChannelReverseMap( int data_taking_timestamp, int swizzling_timestamp ) {
    LoadUBChannelMap( data_taking_timestamp, swizzling_timestamp );
    return fChannelReverseMap;
  }

  // Handy, typical string-splitting-to-vector function. 
  // I hate C++ strong typing and string handling so very, very much. 
  std::vector<std::string> & UBDatabaseUtil::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
  }
 
 
}



namespace util{
 
  DEFINE_ART_SERVICE(UBDatabaseUtil)

} // namespace util
 
  
