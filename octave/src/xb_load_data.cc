//This GNU Octave function provides an interface to load data
//from files created with the XB toolkit.
//The data must have been generated from the class XB::data.
//This means that translated data from experiment and simulation
//can be loaded with this interface.
//NOTE: there's no guarantee whatsoever that this function
//      would be at all usable with MATLAB in a MEX file.
//      It will *not* be tested and it is *not* ment as such.

//octave's documentation string
#define O_DOC_STRING "-*- texinfo -*-\n\
@deftypefn{Function File} {@var{data} =} xb_load_data( @var{filename} )\n\
@deftypefnx{Function File} {@var{data} =} xb_load_data( @var{file_1}, @var{file_2}, ... )\n\
@deftypefnx{Function File} {@var{data} =} xb_load_data( @var{file_1}, ..., @var{nb_events} )\n\
@deftypefnx{Function File} {@var{data} =} xb_load_data( @var{file_1}, ..., [@var{from_event}, @var{to_event}] )\n\
Loads an array of XB::data from a file generated by the program \"xb_data_translator\", among others.\n\
\n\
If one of the @var{file_n} does not exist, a warning message is printed and the file is ignored.\n\
\n\
The format of @var{data} is the same-ish used in the xb_progs toolkit:\n\
@example\n\
@group\n\
data = xb_load_data( 'some_file.xb' );\n\
@result{} structure array data:\n\
     n\n\
     evnt\n\
     i  (array)\n\
     t  (array)\n\
     e  (array)\n\
     he (array)\n\
     sum_e\n\
     in_beta\n\
@end group\n\
@end example\n\
\n\
For more information about the content of the fields, use the documentation of the toolkit.\n\
@end deftypefn"

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

DEFUN_DLD( xb_load_data, args, nargout, O_DOC_STRING ){
	//argument checks: they must be one or more strings
	int nargin = args.length();
	unsigned int load_nb_events[2] = { 0, 0 };
	
	//chek on the numerical types
	if( sizeof(octave_uint32) != sizeof(unsigned int) ){
		error( "Quirky types." );
	}
	
	if( nargin == 0 ){
		error( "At least one file name must be provided" );
	}
	
	//loop-load the files
	std::vector<XB::data> data, data_buf;
	char in_fname[256];
	bool compression_flag = true;
	for( int f=0; f < nargin; ++f ){
		if( args(f).is_string() ){ //if the argument is a string
		                           //attempt to load the file
			octave::sys::file_stat fs( args(f).string_value() );
			if( fs.exists() ){
				strcpy( in_fname, args(f).string_value().c_str() );
				try{
					if( compression_flag ) XB::load( in_fname, data_buf );
					else {
						FILE *input_source = fopen( args(f).string_value().c_str(), "r" );
						XB::load( input_source, data_buf );
						fclose( input_source );
					}
				} catch( XB::error e ){
					error( e.what() );
				}
				data.insert( data.end(), data_buf.begin(), data_buf.end() );
				data_buf.clear();
			} else if( args(f).string_value() == "no-compression" ) compression_flag = false;
			else if( args(f).string_value() == "compression" ) compression_flag = true;
			else {
				octave_stdout << "xb_data_load: warning: file \""
				              << args(f).string_value() << "\" doesn't exist.\n";
				continue;
			}
		} else if( args(f).is_scalar_type() && args(f).is_numeric_type() ){ //just load some events
			load_nb_events[1] = args(f).int_value();
		} else if( !args(f).is_scalar_type() && args(f).is_numeric_type() ){ //load a range of events
			load_nb_events[0] = args(f).int32_array_value()(0);
			load_nb_events[1] = args(f).int32_array_value()(1);
		} else {
			octave_stdout << "xb_load_clusterZ: warning: argument "
			              << f << " is not vaild.\n";
			continue;
		}
	}

	//consistency check on the range: if the range is zero
	//or if it's backward, load everything, silently.
	if( load_nb_events[1] <= load_nb_events[0] ){
		load_nb_events[0] = 0;
		load_nb_events[1] = 0;
		--load_nb_events[1];
	}

	//check that something has been read
	if( !data.size() ){
		octave_stdout << "xb_load_data: warning: no data loaded.\n";
		return octave_value_list();
	}
	
	unsigned int data_size = (( data.size() < load_nb_events[1]-load_nb_events[0] )?
	                          data.size() : load_nb_events[1]-load_nb_events[0] );

	//now, begin the translation into octave structure
	//first, allocate the octave_map that will hold the thing
	dim_vector o_dim_v( data_size, 1 ), o_dim_null( 0, 0 ); 
	
	//copy the data:
	//prepare the fields:
	Cell o_field_n( o_dim_v );
	Cell o_field_evnt( o_dim_v );
	Cell o_field_tpat( o_dim_v );
	Cell o_field_sum_e( o_dim_v );
	Cell o_field_in_beta( o_dim_v );
	Cell o_field_inz( o_dim_v );
	Cell o_field_inaonz( o_dim_v );
	Cell o_field_i( o_dim_v );
	Cell o_field_t( o_dim_v );
	Cell o_field_pt( o_dim_v );
	Cell o_field_e( o_dim_v );
	Cell o_field_he( o_dim_v );
	//and some buffers
	Array<octave_uint32> i_buf;
	Array<float> f_buf;

	unsigned int current_numel = 0;
	for( int i=load_nb_events[0], off_i; i < data.size() && i < load_nb_events[1]; ++i ){
		off_i = i - load_nb_events[0];
		
		current_numel = data[i].n;
	
		//make the structure:
		//firts, copy the trivially copiable
		o_field_n(off_i) = current_numel;
		o_field_evnt(off_i) = data[i].evnt;
		o_field_tpat(off_i) = data[i].tpat;
		o_field_sum_e(off_i) = data[i].sum_e;
		o_field_in_beta(off_i) = data[i].in_beta;
		o_field_inz(off_i) = data[i].in_Z;
		o_field_inaonz(off_i) = data[i].in_A_on_Z;
				
		//then, copy the arrays
		//sizing
		dim_vector o_dim( current_numel, 1 );
		i_buf.resize( o_dim );	
		f_buf.resize( o_dim );
		
		//crystal indexes
		if( current_numel ){
			memcpy( i_buf.fortran_vec(), data[i].i,
			        current_numel*sizeof(unsigned int) );
			o_field_i(off_i) = i_buf;
		} else {
			o_field_i(off_i) = Array<octave_uint32>( o_dim_null );
			continue;
		}
		
		//all the rest
		//TODO: check the i --> off_i
		if( !data[i].empty_t ){
			memcpy( f_buf.fortran_vec(), data[i].t,
			        current_numel*sizeof(float) );
			o_field_t(i) = f_buf;
		} else o_field_t(off_i) = Array<float>( o_dim_null );
		
		if( !data[i].empty_pt ){
			memcpy( f_buf.fortran_vec(), data[i].pt,
			        current_numel*sizeof(float) );
			o_field_pt(i) = f_buf;
		} else o_field_t(off_i) = Array<float>( o_dim_null );
		
		if( !data[i].empty_e ){
			memcpy( f_buf.fortran_vec(), data[i].e,
			        current_numel*sizeof(float) );
			o_field_e(i) = f_buf;
		} else o_field_t(off_i) = Array<float>( o_dim_null );
		
		if( !data[i].empty_he ){
			memcpy( f_buf.fortran_vec(), data[i].he,
			        current_numel*sizeof(float) );
			o_field_he(i) = f_buf;
		} else o_field_t(off_i) = Array<float>( o_dim_null );
		
		//finally, cleanup
		i_buf.clear();
		f_buf.clear();
	}
	
	//make the map
	octave_map o_data_m;
	o_data_m.setfield( "n", o_field_n );
	o_data_m.setfield( "evnt", o_field_evnt );
	o_data_m.setfield( "tpat", o_field_tpat );
	o_data_m.setfield( "i", o_field_i );
	o_data_m.setfield( "t", o_field_t );
	o_data_m.setfield( "pt", o_field_pt );
	o_data_m.setfield( "e", o_field_e );
	o_data_m.setfield( "he", o_field_he );
	o_data_m.setfield( "sum_e", o_field_sum_e );
	o_data_m.setfield( "in_beta", o_field_in_beta );
	o_data_m.setfield( "in_Z", o_field_inz );
	o_data_m.setfield( "in_A_on_Z", o_field_inaonz );
	
	//final cleanup
	data.clear();

	//happy thoughts
	return octave_value_list( octave_value( o_data_m ) );
}		
