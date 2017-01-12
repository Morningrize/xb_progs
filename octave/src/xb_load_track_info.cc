//This GNU Octave function provides an interface to load data
//from files created with the XB toolkit.
//The data must have been generated from the class XB::track_info.
//This means that translated data from experiment and simulation
//can be loaded with this interface.
//NOTE: there's no guarantee whatsoever that this function
//      would be at all usable with MATLAB in a MEX file.
//      It will *not* be tested and it is *not* ment as such.

//stl includes
#include <vector>

//includes from octave
#include <octave/oct.h> //all the gobbins for OCT files
#include <octave/oct-map.h> //data will be reconstructed as a structure (octave_map)
#include <octave/Array.h> //octave arrays arrays
#include <octave/Cell.h> //octave cell arrays
#include <octave/file-stat.h> //file_stat

//includes from the toolkit
#include "xb_io.h" //XB::load
#include "xb_data.h" //XB::data
#include "xb_error.h" //XB::error

//a helper function that does a XB::versor to struct converison.
octave_scalar_map versor2struct( XB::versor &given );

DEFUN_DLD( xb_load_track_info, args, nargout, "XB::load data interface for Octave" ){
	//argument checks: they must be one or more strings
	int nargin = args.length();
	
	//chek on the numerical types
	if( sizeof(octave_uint32) != sizeof(unsigned int) ){
		error( "Quirky types." );
	}
	
	if( nargin == 0 ){
		error( "At least one file name must be provided" );
	}
	
	//loop-load the files
	std::vector<XB::track_info*> data, data_buf;
	char in_fname[256];
	for( int f=0; f < nargin; ++f ){
		if( args(f).is_string() ){ //if the argument is a string
		                           //attempt to load the file
			octave::sys::file_stat fs( args(f).string_value() );
			if( fs.exists() ){
				strcpy( in_fname, args(f).string_value().c_str() );
				XB::load( in_fname, data_buf );
				data.insert( data.end(), data_buf.begin(), data_buf.end() );
			} else {
				octave_stdout << "xb_load_track_info: warning: file \""
				              << args(f).string_value() << "\" doesn't exist.\n";
				continue;
			}
		} else {
			octave_stdout << "xb_load_track_info: warning: argument "
			              << f << " is not a valid filename.\n";
			continue;
		}
	}
	//check that something has been read
	if( !data.size() ){
		octave_stdout << "xb_load_track_info: warning: no data loaded.\n";
		return octave_value_list();
	}
	
	//now, begin the translation into octave structure
	//first, allocate the octave_map that will hold the thing
	dim_vector o_dim_v( data.size(), 1 ), o_dim_null( 0, 0 ); 
	
	//copy the data:
	//prepare the fields:
	Cell o_field_n( o_dim_v );
	Cell o_field_evnt( o_dim_v );
	Cell o_field_in_beta( o_dim_v );
	Cell o_field_beta_0( o_dim_v );
	Cell o_field_in_Z( o_dim_v );
	Cell o_field_in_A_on_Z( o_dim_v );
	Cell o_field_fragment_A( o_dim_v );
	Cell o_field_fragment_Z( o_dim_v );
	Cell o_field_fragment_beta( o_dim_v );
	Cell o_field_incoming( o_dim_v );
	Cell o_field_outgoing( o_dim_v );
	
	//and some buffers
	Array<float> f_buf;
	Cell ov_buf;

	unsigned int current_numel = 0;
	for( int i=0; i < data.size(); ++i ){
		//check for NULLs
		if( data[i] == NULL ) continue;
		current_numel = data[i]->n;
	
		//make the structure:
		//firts, copy the trivially copiable
		o_field_n(i) = current_numel;
		o_field_evnt(i) = data[i]->evnt;
		o_field_in_beta(i) = data[i]->in_beta;
		o_field_beta_0(i) = data[i]->beta_0;
		o_field_in_Z(i) = data[i]->in_Z;
		o_field_in_A_on_Z(i) = data[i]->in_A_on_Z;
				
		//then, copy the arrays
		//sizing
		dim_vector o_dim( current_numel, 1 );
		i_buf.resize( o_dim );	
		f_buf.resize( o_dim );
		ov_buf.resize( o_dim );
		
		//all the rest
		memcpy( f_buf.fortran_vec(), data[i]->fragment_A,
		        current_numel*sizeof(float) );
		o_field_fragment_A(i) = f_buf;
		
		memcpy( f_buf.fortran_vec(), data[i]->fragment_Z,
		        current_numel*sizeof(float) );
		o_field_fragment_Z(i) = f_buf;
		
		memcpy( f_buf.fortran_vec(), data[i]->fragment_beta,
		        current_numel*sizeof(float) );
		o_field_fragment_beta(i) = f_buf;
		
		for( int v=0; v < current_numel; ++v ){
			ov_buf(v) = versor2struct( data[i]->incoming[v] );
		}
		o_field_incoming(i) = ov_buf;
		
		for( int v=0; v < current_numel; ++v ){
			ov_buf(v) = versor2struct( data[i]->outgoing[v] );
		}
		o_field_outgoing(i) = ov_buf;
		
		//finally, deallocate and nullify the copied element
		delete data[i];
		data[i] = NULL;
		f_buf.clear();
		ov_buf.clear();
	}
	
	//make the map
	octave_map o_data_m;
	o_data_m.setfield( "n", o_field_n );
	o_data_m.setfield( "evnt", o_field_evnt );
	o_data_m.setfield( "in_beta", o_field_in_beta );
	o_data_m.setfield( "beta_0", o_field_beta_0 );
	o_data_m.setfield( "in_Z", o_field_in_Z );
	o_data_m.setfield( "in_A_on_Z", o_field_in_A_on_Z );
	o_data_m.setfield( "fragment_A", o_field_fragment_A );
	o_data_m.setfield( "fragment_Z", o_field_fragment_Z );
	o_data_m.setfield( "fragment_beta", o_field_fragment_beta );
	o_data_m.setfield( "incoming", o_field_incoming );
	o_data_m.setfield( "outgoing", o_field_outgoing );
	
	//happy thoughts
	return octave_value_list( octave_value( o_data_m ) );
}

//implementation of the helper function
octave_scalar_map versor2struct( XB::versor &given ){
	octave_scalar_map m;
	m.setfield( "i", given.i );
	m.setfield( "j", given.j );
	m.setfield( "k", given.k );
	
	return m;
}
