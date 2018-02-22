#include "RegionInputer.h"

RegionInputer::RegionInputer(){
	
};

RegionInputer::~RegionInputer(){

}

void RegionInputer::getCO2(RegionData *rd){
	//netcdf error
	NcError err(NcError::silent_nonfatal);

	string filename = md->reginputdir +"co2.nc";
	NcFile co2File(filename.c_str(), NcFile::ReadOnly);
 	if(!co2File.is_valid()){
 		string msg = filename+" is not valid";
 		char* msgc = const_cast< char* > ( msg.c_str());
 		throw Exception(msgc, I_NCFILE_NOT_EXIST);
 	}
 	
 	NcDim* yrD = co2File.get_dim("YEAR");
 	if(!yrD->is_valid()){
  		throw Exception("YEAR Dimension is not Valid in RegionInputer::getCO2", I_NCDIM_NOT_EXIST);
 	}
 	 	
 	NcVar* co2yrV = co2File.get_var("YEAR");
 	NcVar* co2V = co2File.get_var("CO2");
 	if(co2yrV==NULL || co2V==NULL){
 	   throw Exception("Cannot get CO2 in RegionInputer::getCO2 ", I_NCVAR_NOT_EXIST);
 	}

	NcBool nb1 = co2yrV->get(&rd->co2year[0], MAX_CO2_DRV_YR);
	NcBool nb2 = co2V->get(&rd->co2[0], MAX_CO2_DRV_YR);
	if(!nb1 || !nb2){
	 throw Exception("problem in reading CO2 in RegionInputer::getCO2 ", I_NCVAR_GET_ERROR);
	}

	co2File.close();
}

void RegionInputer::setModelData(ModelData* mdp){
   	md = mdp;
};

