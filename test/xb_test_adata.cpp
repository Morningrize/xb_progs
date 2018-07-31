//the new, expanded version of the test program for the adata structure.

#include <stdio.h>
#include <getopt.h>

#include "xb_arbitrary_data.h"
#include "xb_io.h"

char fields[256][16];
XB::adata_indexer test_indexer( int fcount );

int main( int argc, char **argv ){
    int fcount = 0;
    
	for( int i=1; i < argc; ++i ){
		if( argv[i][0] == '-' ) break;
		strncpy( fields[fcount], argv[i], 16 );
		++fcount;
	}
	
	char iota=0;
	int ns=10; //number of structures to allocate
	while( (iota = getopt( argc, argv, "n:" )) != -1 ){
		switch( iota ){
			case 'n':
				ns = atoi( optarg );
				break;
			default :
				fprintf( stderr, "YOU made a mistake.\n" );
				exit( 1 );
		}
	}
	
	puts( "*** Weclome in XB::adata test program ***" );
    
    //test the indexer.
    puts( ">>>PHASE 1: testing the indexer" );
    XB::adata_indexer idx = test_indexer( fcount );
    
    //TBC
    return 0;
}

//---------------------------------------------------------------------------------
XB::adata_indexer test_indexer( int fcount ){
    XB::adata_indexer idx( fcount );
    printf( "\tConstructed an indexer of size %d\n", fcount );
    for( int i=0; i < idx.names.size(); ++i ){
        strncpy( idx[i].name, fields[i], XB_ADATA_FIELD_NAME_LENGTH );
        idx[i].size = 4*sizeof( float );
        idx.diffs[i] = 4*i*sizeof(float);
    }
    
    XB::adata_indexer idx_again( idx );
    puts( "\tCopy constructed an indexer:" );
    
    for( int i=0; i < fcount; ++i ){
        printf( "\tfield %s:\n", fields[i] );
        printf( "\t\tOriginal: %s:%d:%d; copied %s:%d:%d\n",
                idx[i].name, idx[i].size, idx.diffs[i],
                idx_again[i].name, idx_again[i].size, idx.diffs[i] );
        idx[i].name[0] = 'Z';
        idx.diffs[i] += 100;
    }
    
    idx = idx + idx_again;
    puts( "\tConcatenated indexes" );
    
    for( int i=0; i < idx.size(); ++i )
        printf( "\tField: %s:%d:%d\n", idx[i].name, idx[i].size, idx.diffs[i] );
    
    //reset the structure for sane later usage
    memset( idx_again.diffs, 0, XB_ADATA_NB_FIELDS );

    return idx_again;
}
