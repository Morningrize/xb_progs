//this program applies a threshold to the crystals

#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <limits.h>
#include <algorithm>
#include <vector>

#include "xb_io.h"
#include "xb_data.h"
#include "xb_cluster.h"

#define MAX_FILES 64

#define VERBOSE 0x01
#define KLZ_FLAG 0x02
#define OUT_TO_FILE 0x04

int get_thresholds( char *thr_name, float *cry_thr );
void apply_thr( std::vector<XB::data> &data, float *cry_thr );
void apply_thr( std::vector<XB::clusterZ> &klz, float *cry_thr );

int main( int argc, char **argv ){
	char in_fname[MAX_FILES][256];
	char thr_name[256] = "/dev/null";
	char out_fname[256];
	int flagger=0;
	
	int in_fcount=0;
	for( int i=1; i < argc && i-1 < MAX_FILES; ++i ){
		if( argv[i][0] == '-' ) break;
		strncpy( in_fname[in_fcount], argv[i], 256 );
		++in_fcount;
	}
	
	char iota=0; int idx;
	struct option opts[] = {
		{ "verbose", no_argument, &flagger, flagger | VERBOSE },
		{ "cluster", no_argument, &flagger, flagger | KLZ_FLAG },
		{ "threshold-file", required_argument, NULL, 'f' },
		{ "output", required_argument, NULL, 'o' },
		{ NULL, 0, NULL, 0 }
	};
	
	while( (iota = getopt_long( argc, argv, "vkf:o:", opts, &idx )) != -1 ){
		switch( iota ){
			case 'v' :
				flagger |= VERBOSE;
				break;
			case 'k' :
				flagger |= KLZ_FLAG;
				break;
			case 'o' :
				flagger |= OUT_TO_FILE;
				strncpy( out_fname, optarg, 256 );
				break;
			case 'f' :
				strncpy( thr_name, optarg, 256 );
				break;
			default : exit( 1 );
		}
	}
	
	if( flagger & VERBOSE ) puts( "*** Welcome in the threshold program! ***" );
	
	std::vector<XB::clusterZ> klz, kbuf;
	std::vector<XB::data> data, dbuf;
	
	if( flagger & VERBOSE ) puts( "Reading:" );
	if( flagger & KLZ_FLAG ){
		if( in_fcount ) for( int i=0; i < in_fcount; ++i ){
			if( flagger & VERBOSE ) printf( "\t%s\n", in_fname[i] );
			XB::load( in_fname[i], kbuf );
			klz.insert( klz.end(), kbuf.begin(), kbuf.end() );
			kbuf.clear();
		} else XB::load( stdin, klz );
	} else {
		if( in_fcount ) for( int i=0; i < in_fcount; ++i ){
			if( flagger & VERBOSE ) printf( "\t%s\n", in_fname[i] );
			XB::load( in_fname[i], dbuf );
			data.insert( data.end(), dbuf.begin(), dbuf.end() );
			dbuf.clear();
		} else XB::load( stdin, data );
	}
	
	float cry_thr[162];
	for( int i=0; i < 162; ++i ) cry_thr[i] = INT_MAX;
	if( !strcmp( thr_name, "/dev/null" ) ){
		puts( "No threshold file: nothing to do." );
		exit( 0 );
	}
	if( flagger & VERBOSE ) printf( "Reading thresholds: %s\n", thr_name );
	int nb_thr = get_thresholds( thr_name, cry_thr );
	if( flagger & VERBOSE ) printf( "\t%d thresholds read.\n", nb_thr );
	
	if( flagger & KLZ_FLAG ) apply_thr( klz, cry_thr );
	else apply_thr( data, cry_thr );
	if( flagger & VERBOSE ) puts( "Thresholds applied." );
	
	if( flagger & OUT_TO_FILE ){
		if( flagger & VERBOSE ) printf( "Writing to %s...\n", out_fname );
		if( flagger & KLZ_FLAG ) XB::write( out_fname, klz );
		else XB::write( out_fname, data );
	} else {
		if( flagger & KLZ_FLAG ) XB::write( stdout, klz );
		else XB::write( stdout, data );
	}
	
	//happy thoughts
	if( flagger & VERBOSE ) printf( "*** Done, goodbye. ***\n" );
	return 0;
}

//------------------------------------------------------------------------------------
int get_thresholds( char *thr_name, float *cry_thr ){
	char *pcommand = (char*)calloc( strlen( thr_name ) + 128, 1 );
	strcpy( pcommand, "awk '$1 == \"Crystal\" { print $3 }; $1 == \"Cutoff\" { if( $2 ) print $2; else print 2147483647 }' " );
	strcat( pcommand, thr_name );
	FILE *thr_pipe = popen( pcommand, "r" );
	
	int c, count=0; float t;
	while( !feof( thr_pipe ) ){
		fscanf( thr_pipe, "%d", &c );
		fscanf( thr_pipe, "%f", &t );
		++count;
		cry_thr[c-1] = t;
	}
	
	pclose( thr_pipe );
	return count;
}

//------------------------------------------------------------------------------------
void apply_thr( std::vector<XB::data> &data, float *cry_thr ){
	for( int i=0; i < data.size(); ++i )
		for( int h=0; h < data[i].n; ++h )
			if( data[i].e[h] < cry_thr[data[i].i[h]-1] ) data[i].e[h] = 0;
}

//hideous nidiication of for loops
void apply_thr( std::vector<XB::clusterZ> &klz, float *cry_thr ){
	for( int k=0; k < klz.size(); ++k )
		for( int sk=0; sk < klz[k].clusters.size(); ++sk ){
			for( int c=0; c < klz[k].clusters[sk].crys.size(); ++c )
				if( klz[k].clusters[sk].crys_e[c] > cry_thr[klz[k].clusters[sk].crys[c]-1] )
					klz[k].clusters[sk].crys_e[c] = 0;
			klz[k].clusters[sk].sum_e = 0;
			for( int c=0; c < klz[k].clusters[sk].crys.size(); ++c )
				klz[k].clusters[sk].sum_e += klz[k].clusters[sk].crys[c];
		}
}
