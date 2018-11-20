//the new, expanded version of the test program for the adata structure.

#include <stdio.h>
#include <getopt.h>

#include "xb_arbitrary_data.h"
#include "xb_io.h"

char fields[256][16];
XB::adata_indexer test_indexer( int fcount );
XB::adata test_adata( XB::adata_indexer idx );
XB::adata_uniarr test_uniarr( XB::adata_indexer idx );

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
    
    //puts( ">>>PHASE 2: testing the structure" );
    //XB::adata a = test_adata( idx );
    
    puts( ">>>PHASE 3: testing the uniform array" );
    XB::adata_uniarr ua = test_uniarr( idx );
    
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
    
    for( int i=0; i < idx_again.size(); ++i ){
        printf( "\tfield %s:\n", fields[i] );
        printf( "\t\tOriginal: %s:%d:%d; copied %s:%d:%d\n",
                idx[i].name, idx[i].size, idx.diffs[i],
                idx_again[i].name, idx_again[i].size, idx_again.diffs[i] );
        idx[i].name[0] = 'Z';
        idx.diffs[i] = -1;
    }
    
    idx = idx + idx_again;
    puts( "\tConcatenated indexes" );
    
    for( int i=0; i < idx.size(); ++i )
        printf( "\tField: %s:%d:%d\n", idx[i].name, idx[i].size, idx.diffs[i] );
    
    //reset the structure for sane later usage
    for( int i=0; i < 2*fcount; ++i ) idx_again.diffs[i] = 4*i*sizeof( float );
    for( int i=2*fcount; i < XB_ADATA_NB_FIELDS; ++i ) idx_again.diffs[i] = -1;

    return idx_again;
}

//------------------------------------------------------------------------------------
XB::adata test_adata( XB::adata_indexer idx ){
    XB::adata a;
    puts( "\tDefault constructed (a)" );
    XB::adata_put( stdout, a );
    
    XB::adata b( idx.names.data(), idx.size() );
    puts( "\tOld style constructor (b)" );
    XB::adata_put( stdout, b );
    
    XB::adata c( idx );
    puts( "\tIndexer constructor (c)" );
    XB::adata_put( stdout, c );
    
    a = c.copy();
    puts( "\tCopy function (c -> a)" );
    XB::adata_put( stdout, a );
    
    int smt = 42;
    a.dofield( "try", sizeof( int ), &smt );
    puts( "\tPushed a field (a)" );
    XB::adata_put( stdout, a );
    
    a.rmfield( "try" );
    puts( "\tRemoved the field (a)" );
    XB::adata_put( stdout, a );
    
    b.clear();
    puts( "\tClear a structure (b)" );
    XB::adata_put( stdout, b );
    
    std::vector< XB::adata_field > fld = c.lsfields();
    puts( "\tField list (c)" );
    for( int i=0; i < fld.size(); ++i ){
        printf( "\t---%s:%d\n", fld[i].name, fld[i].size );
    }
    
    void *buf;
    XB::adata_getlbuf( &buf, c );
    puts( "\tLinearized c" );
    XB::adata_fromlbuf( a, buf );
    puts( "\tDelinearized in a" );
    XB::adata_put( stdout, a );
    
    a.clear();
    a.dofield( "try", sizeof(int), &smt );
    a.dofield( "this", sizeof(int), &smt );
    XB::adata d = XB::adata_merge( a, c );
    puts( "\tMerged a{ try, this } and c in d" );
    XB::adata_put( stdout, a );
    XB::adata_put( stdout, c );
    XB::adata_put( stdout, d );
    
    return a;
}

//---------------------------------------------------------------------------------
//and now the big one: the uniform array
XB::adata_uniarr test_uniarr( XB::adata_indexer idx ){
    XB::adata_uniarr a;
    puts( "\tUA default constructed" );
    
    XB::adata_uniarr b( idx );
    puts( "\tUA constructed with indexer" );
    XB::adata_indexer idx_again = b.get_indexer();
    for( int i=0; i < idx_again.size(); ++i ){
        printf( "\t\tOriginal: %s:%d:%d; copied %s:%d:%d\n",
                idx[i].name, idx[i].size, idx.diffs[i],
                idx_again[i].name, idx_again[i].size, idx_again.diffs[i] );
    }
    
    XB::adata_uniarr c( 3, idx );
    puts( "\tUA constructed with elements and indexer" );
    for( int i=0; i < c.size(); ++i ) XB::adata_put( stdout, c[i] );
    
    puts( "\tUA pop back" );
    XB::adata popped = c.pop_back();
    printf( "\tsize: %d\n", c.size() );
    puts( "\tpopped element:" );
    XB::adata_put( stdout, popped );
    
    puts( "\tUA push back" );
    //NOTE: this won't work because idx is broken! XB::adata pushed( c.get_indexer() );
    XB::adata pushed( popped );
    XB::adata_put( stdout, pushed );
    c.push_back( pushed );
    printf( "\tsize: %d\n", c.size() );
    for( int i=0; i < c.size(); ++i ) XB::adata_put( stdout, c[i] );
    
    puts( "\tUA push indexer" );
    XB::adata_field f = { "new", 8 };
    c.push_indexer( f );
    XB::adata_put( stdout, c[0] );
    
    puts( "\tUA pop indexer" );
    c.pop_indexer();
    XB::adata_put( stdout, c[0] );
    
    /*NOTE: it works, but idx is broken so here it doesn't
    puts( "\tUA operator+ and =" );
    for( int i=0; i < c.size(); ++i ) b.push_back( XB::adata( idx ) );
    b.push_indexer( f );
    b = b+c;
    printf( "\tsize: %d\n", b.size() );
    for( int i=0; i < b.size(); ++i ) XB::adata_put( stdout, b[i] );
    */
    
    puts( "\tUA clear, isempty" );
    b.clear();
    printf( "\tsize: %d\n", b.size() );
    if( b.empty() ) puts( "\tcorrect" );
    
    return c;
}
    
    
    
    
