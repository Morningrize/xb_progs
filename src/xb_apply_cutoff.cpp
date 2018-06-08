//this program applies a threshold to the crystals

#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <algorithm>
#include <vector>

#include "xb_io.h"
#include "xb_data.h"
#include "xb_cluster.h"

#define MAX_FILES 64

#define VERBOSE 0x01
#define KLZ_FLAG 0x02
#define OUT_TO_FILE 0x04

int main( int argc, char **argv ){
	char in_fname[MAX_FILES][256];
	char thr_name[256];
	char out_fname[256];
	
	int in_fcount=0;
	for( int i=1; i < argc && i-1 < MAX_FILES; ++i ){
		if( argv[i][0] == '-' ) break;
		strncpy( in_fnames[in_fcount], argv[i], 256 );
		++in_fcount;
	}
	
	char iota=0; int idx;
	struct option opts[] = {
		{ "verbose", no_argument, &flagger, flagger | VERBOSE },
		{ "cluster", no_argument, &flagger, flagger | KLZ },
		{ "threshold-file", required_argument, NULL, 't' },
		{ "output", required_argument, NULL, 'o' },
		{ NULL, 0, NULL, 0 }
	};
	
	type_type tt = XB_DATA; //if no type is specified, assume data.
	while( (iota = getopt_long( argc, argv, "vkt:o:", opts, &idx )) != -1 ){
		switch( iota ){
			case 'v' :
				flagger |= VERBOSE;
				break;
			case 'k' :
				flagger |= KLZ_FLAG;
				break;
			case 'o' :
				flagger |= OUT_TO_FILE
				strncpy( out_fname, optarg, 256 );
				break;
			case 't' :
				strncpy( thr_file, optarg, 256 );
				break;
			default : exit( 1 );
		}
	}
	
	std::vector<XB::clusterZ> klz, kbuf;
	std::vector<XB::data> data, dbuf;
	
	if( flagger & KLZ_FLAG ){
		if( in_fcount ) for( int i=0; i < in_fcount; ++i ){
			XB::load( in_fname[i], kbuf );
			klz.insert( klz.end(), kbuf.begin(), kbuf.end() );
			kbuf.clear();
		} else XB::load( stdin, klz );
	} else {
		if( in_fcount ) for( int i=0; i < in_fcount, ++i ){
			XB::load( in_fname[i], dbuf );
			data.insert( data.end(), dbuf.begin, dbuf.end() );
			dbuf.clear();
		} else XB::load( stdin, data );
	}
	
	float cry_thr[162];
	XB::get_thresholds( thr_name, cry_thr );
	
	if( flagger & KLZ_FLAG ) apply_thr( klz, cry_thr );
	else apply_thr( data, cry_thr );
	
	if( flagger & OUT_TO_FILE ){
		if( flagger & VERBOSE ) printf( "Writing to %s...\n", out_fname );
		if( flagger & KLZ_FLAG ) XB::write( out_f_name, klz );
		else XB::write( out_f_name, data );
	} else {
		if( flagger & KLZ_FLAG ) XB::write( stdout, klz );
		else XB::write( stdout, data );
	}
	
	//happy thoughts
	if( flagger & VERBOSE ) printf( "*** Done, goodbye. ***\n" );
	return 0;
}
