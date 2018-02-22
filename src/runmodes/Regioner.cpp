#include "Regioner.h"

Regioner::Regioner(){
	
};

Regioner::~Regioner(){
	
};


void Regioner::init(string& controlfile){

	try{

		configin.controlfile = controlfile;
		//Input and output initialization		
 		configin.ctrl4regnrun(&md);

 		md.myid = 0;
 		md.numprocs = 1;

 		md.checking4run();
 
 		fd.useseverity = md.useseverity;
 		
 		md.consoledebug = false;

 		//create a list of cohort id, each process should run through all cohorts in the list
 		createCohorList4Run();

 		//region-level input
 		rin.setModelData(&md);
 		rin.getCO2(&rd);

 		//grid-level input
 		gin.setModelData(&md);
 		gin.init();

 		//cohort-level input
 		cin.setModelData(&md);
		cin.init();	
		
 		//inputers
 		if(md.initmode==2){
	 		sitein.initSiteinFile(md.initialfile);
	 		runcht.setSiteinInputer(&sitein);
 		} else if(md.initmode==3){
 		 	if(md.runeq){
 		 		cout <<"cannot set initmode as restart for equlibrium run  \n";
 		 		cout <<"reset to 'lookup'\n";
 		 		md.initmode=1;
 		 	} else {
 		 		resin.init(md.initialfile);
 		 		runcht.setRestartInputer(&resin);
 		 	}
 		} else {
 			md.initmode = 1;
 		}

 		runcht.setGridInputer(&gin);
 		runcht.setCohortInputer(&cin);
 		 
 		//output variables, outputers and initialization
 		   //1) multiple-cohort/yearly output
 		string stage("?");
 		MAX_OREGN_YR = 0;
		if(md.runeq){
			stage="-eq";
 		} else {
 			if(md.runsp){
 				stage="-sp";
 				MAX_OREGN_YR = MAX_SP_YR;
 			}

 			if(md.runtr){
 				stage="-tr";
 				MAX_OREGN_YR = MAX_TR_YR;
 			}

			runcht.setRegionalOutputer(&rout);
			rout.setOutData(&regnod);
			rout.init(md.outputdir, md.myid, stage, MAX_OREGN_YR);
 			runcht.cht.setRegnOutData(&regnod);

 		}
 		resout.setRestartOutData(&resod);
		resout.init(md.outputdir, stage, md.numprocs, md.myid);
 		runcht.cht.setRestartOutData(&resod);
 		runcht.setRestartOutputer(&resout);

 		//error output
		errout.init(md.outputdir, md.myid, stage);

 		//set up data (inputs and processes) connecntion (initialization)
		bd.setEnvData(&ed);
		rgrid.setEnvData(&ed);
 		rgrid.setRegionData(&rd);
 		rgrid.setGridData(&gd);

 		runcht.cht.setTime(&timer);
 		runcht.cht.setProcessData(&ed, &bd, &fd);
 		runcht.cht.setModelData(&md);
 		runcht.cht.setInputData(&rd, &gd, &cd);
 		runcht.cht.setAtmData(&rgrid);

 		//ONE cohort initialization
  		runcht.cht.init(); //after set everying
 
	}catch (Exception &exception){
  		cout <<"problem in initialize in Regioner::init\n";
  		exception.mesg();
  		exit(-1);
	}
    
};


void Regioner::run(){
	
	int prevcruid =-1;
	
	//error initialization
	int errcount = 0;
	errout.errorid = 0;

	list<int>::iterator jj ; 
	for ( jj=runchtlist.begin() ; jj!=runchtlist.end(); jj++){
		int chtid = *jj;
		
		errout.chtid = chtid;
		
		//get the eqchtid, spchtid/trchtid, restart-id, and cruid
		int eqcid = 0;  //the record order in the input files, NOT the cohort ID (chtid)
		int cid = 0;    //the record order in the input files, NOT the cohort ID (chtid)
		int rescid = 0; //the record order in the input files, NOT the cohort ID (chtid)

		try {

			//for regional run, only one of the following can be true;
			if(md.runeq){
				cd.eqchtid = chtid;
				eqcid=cin.getEqRecID(cd.eqchtid);  //needed for cruid searching
				cid=eqcid;                         //
			}

			if(md.runsp){
				cd.spchtid = chtid;
				cid=cin.getSpRecID(cd.spchtid);

				cin.getEqchtid5SpFile(cd.eqchtid, cid);   
				eqcid=cin.getEqRecID(cd.eqchtid);
			
				cd.reschtid = cd.eqchtid;
				if (md.initmode==3) rescid = resin.getRecordId(cd.reschtid);
			}
		
			if(md.runtr){
				cd.trchtid = chtid;			
				cid=cin.getTrRecID(cd.trchtid);
				cin.getSpchtid5TrFile(cd.spchtid, cid); 

				int spcid=cin.getSpRecID(cd.spchtid);
				cin.getEqchtid5SpFile(cd.eqchtid, spcid);
				eqcid=cin.getEqRecID(cd.eqchtid);
			
				cd.reschtid = cd.spchtid;
				if (md.initmode==3) rescid = resin.getRecordId(cd.reschtid);
			}
			
			cin.getCRUID(cd.cruid, eqcid); //(eq/sp/tr)cid: starting from ZERO

		} catch (Exception &exception){
			errout.errorid = -1;
			errout.outputVariables(errcount);
			errcount+=1;

			if(md.consoledebug){
				cout <<"problem in setting IDs in Regioner::run\n";
				exception.mesg();
			}

		}

		if(cd.cruid>=0 && cid>=0 && eqcid>=0 && rescid>=0){
   			int error = 0;

   			//grid-level data for a cohort
			try {
				int gid = gin.getGridRecID(cd.cruid);
				gd.gid=gid;
    			if(cd.cruid!=prevcruid){
    				gin.getGridData(&gd, gid);   

    				error = rgrid.reinit(gid); //reinit for a new grid
    				if (error!=0) {

    					if(md.consoledebug){
    						cout <<"problem in grid data in Regioner::run\n";
    					}

    		    		errout.errorid = -3;
    					errout.outputVariables(errcount);
    					errcount+=1;

    					continue;     //jump over to next cohort, due to grid-data error
    				}

					prevcruid =gid;
    			}
	    	} catch (Exception &exception){
				exception.mesg();
				if(md.consoledebug){
					cout <<"problem in reinitializing grid in Regioner::run\n";
				}

	    		errout.errorid = -3;
				errout.outputVariables(errcount);
				errcount+=1;

				continue;     //jump over to next cohort, due to grid-data error

 	    	}
 
	    	//cohort-level data for a cohort
			runcht.jcalifilein=true;  // for reading Jcalinput.txt, the default is true (must be done before re-initiation)
			runcht.ccdriverout=false;  //don't change to true for regioner
				
			error = runcht.reinit(cid, eqcid, rescid); //reinit for a new cohort

			//run a cohort
   			try {
  				if (error!=0) {
   					cout<<"Error for reinitializing cohort: "<<chtid<<" - SKIPPED! \n";

   					errout.errorid = -5;
    				errout.outputVariables(errcount);
    				errcount+=1;

					continue;     //jump over to next cohort, due to cohort reinit error

  				} else {
    				if (md.consoledebug)
    					cout<<"cohort: "<<chtid<<" @ "<<md.runstages
							<<" - running! \n";
    				runcht.run();
    			}

    		} catch (Exception &exception){
				exception.mesg();
    			if(md.consoledebug){
    					cout <<"problem in running cohort in Regioner::run\n";
    			}

    			rout.missingValues(MAX_OREGN_YR, runcht.cht.cohortcount);

    			errout.errorid = -4;
    			errout.outputVariables(errcount);
    			errcount+=1;

				continue;     //jump over to next cohort, due to run cohort error

    		}

		} else { // end of cruid >=0 && other IDs>=0
			cout<<"No grid exists for cohort: "<<chtid<<" - SKIPPED! \n";

			rout.missingValues(MAX_OREGN_YR, runcht.cht.cohortcount);

			errout.errorid = -2;
			errout.outputVariables(errcount);
			errcount+=1;
		} // end of cruid >=0 && errout.errorid ==0

		runcht.cht.cohortcount++;
		 
	}// end of cohort loop
	
};

void Regioner::createCohorList4Run(){
	// read in a list of cohorts to run

	//netcdf error
	NcError err(NcError::silent_nonfatal);

	//open file and check if valid
	string filename = md.runchtfile;
	NcFile runFile(filename.c_str(), NcFile::ReadOnly);
 	if(!runFile.is_valid()){
 		string msg = filename+" is not valid";
 		char* msgc = const_cast< char* > ( msg.c_str());
 		throw Exception(msgc, I_NCFILE_NOT_EXIST);
 	}
 	
 	NcDim* chtD = runFile.get_dim("CHTID");
 	if(!chtD->is_valid()){
 		throw Exception("CHT Dimension is no Valid in createCohortList4Run", I_NCDIM_NOT_EXIST);
 	}
 	
 	NcVar* chtV = runFile.get_var("CHTID");
 	if(chtV==NULL){
 	   throw Exception("Cannot get CHTID in createCohortList4Run ", I_NCVAR_NOT_EXIST);
 	}

 	int numcht = chtD->size();
 	
	int chtid  = -1;
	int chtid0 = -1;
	int chtidx = -1;
	for (int i=0; i<numcht; i++){
		chtV->set_cur(i);
   		chtV->get(&chtid, 1);
   		runchtlist.push_back(chtid);
	   	
	   	if (i==0) chtid0=chtid;
	   	if (i==numcht-1) chtidx=chtid;
   	}

	cout <<md.casename << ": " <<numcht <<"  cohorts to be run @" <<md.runstages<< "\n";
	cout <<"   from:  " <<chtid0<<"  to:  " <<chtidx <<"\n";
   
};



