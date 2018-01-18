//this makes a XB cluster-level rootfile

#include "xb_root_writer.h"
unsigned int default_beam_out = 81; 

//this one just writes an XBc rootfile
void XB::rwrite(char* f_root_out, std::vector<XB::clusterZ> &event_klZ ){
	xb_ball the_cb;
	//event
	unsigned int Xbevnt=0;
	unsigned int Xbcmult=0;
	//clusters in event
	unsigned int Xbcn[162],Xbci[162];
	float Xbcth[162],Xbcsume[162];
	//crystals in cluster
	//unsigned int Xbcii[162][162];

	unsigned int default_beam_out = 81; //the default beam in
	//sort the xb cluster data by event id before writing
	std::sort( event_klZ.begin(), event_klZ.end(), evnt_id_comparison );

	TFile *fout = TFile::Open(f_root_out, "recreate"); //output file
	TTree *hxbc = new TTree("hxbc", "An XB_cluster TTree"); //the XB cluster tree
	//branches
	hxbc->Branch("Xbevnt",&Xbevnt); //event number
	hxbc->Branch("Xbcmult",&Xbcmult); //number of clusters in an event
	hxbc->Branch("Xbcn",&Xbcn,"Xbcn[Xbcmult]/i"); //number of crystals in a cluster
	hxbc->Branch("Xbci",&Xbci,"Xbci[Xbcmult]/i"); //index of the centroid
	hxbc->Branch("Xbcth",&Xbcth,"Xbcth[Xbcmult]/F"); //centroid theta ( inclination towards beam line )
	hxbc->Branch("Xbcsume",&Xbcsume,"Xbcsume[Xbcmult]/F"); //gamma sum energy of cluster 
	//hxbc->Branch("Xbcii",&Xbcii,"Xbcii[Xbcmult][162]/i");//indices of paricipating crystals, organised clusterwise

	//start eventloop
	for( unsigned int i=0; i < event_klZ.size(); ++i ){
		Xbevnt=event_klZ[i].evnt;	
		Xbcmult=event_klZ[i].n;
		if (Xbcmult>0){ //are there any clusters?
			for ( unsigned int k=0; k < event_klZ[i].n; ++k ){//loop over the clusters
				Xbcn[k]=event_klZ[i].clusters[k].n;
				Xbci[k]=event_klZ[i].clusters[k].centroid_id;
				Xbcth[k]=angular_distance( the_cb.at( default_beam_out ).altitude, the_cb.at( default_beam_out ).azimuth, 
						event_klZ[i].clusters[k].c_altitude, event_klZ[i].clusters[k].c_azimuth );
				Xbcsume[k]=event_klZ[i].clusters[k].sum_e;
/*				for ( unsigned int cr=0; cr < event_klZ[i].clusters[k].n; cr++ ){//loop over crystals in cluster
					Xbcii[k][cr]=event_klZ[i].clusters[k].crys.at(cr);
				} 
*/			}
		}
	hxbc->Fill();
	}//end eventloop
	hxbc->Write();
	fout->Close();	
}

// this one will stitch
void XB::rwrite(char* f_root_out, char* stch_r_file, std::vector<XB::clusterZ> &event_klZ, bool v_flag){
	xb_ball the_cb;
        unsigned int Xbcmult=0;
        unsigned int Xbcn[162],Xbci[162];
        float Xbcth[162],Xbcsume[162];
	//crystals in cluster
	//unsigned int Xbcii[162][162]; //not very useful atm, so I leave it out

	unsigned int default_beam_out = 81; //the default beam in

	TFile stitch_file( stch_r_file );
	if( stitch_file.IsZombie() ){
                throw XB::error( "File error!", "XB::rwrite" );
        }

	//relating to the root file we stitch the XB data to
	//get the name of the TTree
	TIter nextkey( stitch_file.GetListOfKeys() );
        TKey *key;
        TObject *obj;
        TTree *data_tree;
        while( (key = (TKey*)nextkey()) ){
                obj = key->ReadObj();
                if( obj->IsA()->InheritsFrom( TTree::Class() ) ){
                        data_tree = (TTree*)obj;
                        break;
                }
        }

	//we need the event id from the tree
	int Evnt;
	TBranch *b_Evnt;

	data_tree->SetBranchAddress("Evnt",&Evnt, &b_Evnt);	

	//clone the tree into a new root file	
	TFile* fout = new TFile( f_root_out, "recreate" );
	TTree *newtr = data_tree->CloneTree(0);
	//add new branches
        TBranch *bxbcmult = newtr->Branch("Xbcmult",&Xbcmult); //number of clusters in an event
        TBranch *bxbcn = newtr->Branch("Xbcn",&Xbcn,"Xbcn[Xbcmult]/i"); //number of crystals in a cluster
        TBranch *bxbci = newtr->Branch("Xbci",&Xbci,"Xbci[Xbcmult]/i"); //index of the centroid
        TBranch *bxbth = newtr->Branch("Xbcth",&Xbcth,"Xbcth[Xbcmult]/F"); //centroid altitude
        TBranch *bxbcsume = newtr->Branch("Xbcsume",&Xbcsume,"Xbcsume[Xbcmult]/F"); //gamma sum energy of cluster
	//TBranch *bxbcii = newtr->Branch("Xbcii",&Xbcii,"Xbcii[Xbcmult][162]/i");//indices of paricipating crystals, organised clusterwise

	//number of events
	int nevents=data_tree->GetEntries();

	//guess the ordering
	int ord=0;
	srand(nevents);
	for ( int chk=0;chk<10;chk++){  
		int ri=rand()%(nevents-1);	
		data_tree->GetEvent( ri );unsigned int ev1=Evnt;
		data_tree->GetEvent( ri+1 );unsigned int ev2=Evnt;	
		if ( ev1 > ev2 ) ord--;
		else ord++;
	}

	if ( ord==-10 && v_flag ) puts( "Rootfile has reverse ordering by event id " );
	else if ( ord==10 && v_flag ) puts( "Rootfile is sorted by event id " );
	else if ( v_flag ) puts ( "Rootfile is not ordered, this will be slow "); 

	//sort the xb cluster data by event id (it's probably sorted already) 
	std::sort( event_klZ.begin(), event_klZ.end(), evnt_id_comparison );
	//reverse if root file is sorted oppositely (shouldn't be the case, but anyway...)
	if ( ord==-10 ) std::reverse( event_klZ.begin(), event_klZ.end() );

	//iterator
	unsigned int cl_i = 0; 
	unsigned int v_cl_i=cl_i; 

	//loop over events 
	for (int jentry=0;jentry<nevents;jentry++){
		data_tree->GetEvent(jentry);	
		if ( v_flag && (float)(jentry/25000.)==(int)(jentry/25000.) ) printf( "%i %% done \n", (int)(100.*jentry/nevents) );  
		bool match = false;
		//find the match in the cluster data
		while ( !match && cl_i >= 0 && cl_i < event_klZ.size() ) {
			if ( event_klZ[cl_i].evnt ==  Evnt  ) match = true;
			else cl_i++;
		}//end of search over cluster events
		//reset counter to last entry found
		if ( match&&abs(ord)==10 ) v_cl_i=cl_i;
		else if ( abs(ord)==10 ) cl_i=v_cl_i;
		else cl_i=0; 
                Xbcmult=event_klZ[cl_i].n;
                if (Xbcmult>0){ //are there any clusters?
                        for( unsigned int k=0; k < event_klZ[cl_i].n; ++k ){//loop over the clusters
                                Xbcn[k]=event_klZ[cl_i].clusters[k].n;
                                Xbci[k]=event_klZ[cl_i].clusters[k].centroid_id;
                                Xbcth[k]=angular_distance( the_cb.at( default_beam_out ).altitude, the_cb.at( default_beam_out ).azimuth,
                                                event_klZ[cl_i].clusters[k].c_altitude, event_klZ[cl_i].clusters[k].c_azimuth );
                                Xbcsume[k]=event_klZ[cl_i].clusters[k].sum_e;
/*				for ( unsigned int cr=0; cr < event_klZ[cl_i].clusters[k].n; cr++ ){//loop over all crystals in cluster
					Xbcii[k][cr]=event_klZ[cl_i].clusters[k].crys.at(cr);
				} 
*/			}
                }
		//Fill the new branches
		newtr->Fill();
	}//end of eventloop (stitch file)
	newtr->Write();
	delete newtr;
	fout->Close();
}

//the comparison utility
//for the event holder structure
bool evnt_id_comparison( const XB::event_holder &one, const XB::event_holder &two ){
        return one.evnt < two.evnt;
}

