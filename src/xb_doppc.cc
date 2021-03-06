#include "xb_doppc.h"

//====================================================================================
//implementation of the function that drives the doppler correction
//in parallel.
//NOTE: for now, the fedault beam out is always used as the direction
//      in the near future, the directional information from the tracker
//      will be used when available
//NOTE2: I do apologize to whoever will read this function, it is a bit
//       more complicated than usual -especially because it's a bit of a fest
//       of double-pointers and whatever.
void apply_doppler_correction( std::vector<XB::data> &xb_book,
                               std::vector<XB::track_info> &xb_track_book,
                               unsigned int default_beam_out, correct_mode mode,
                               int flagger ){
	//ok, now do some actual work:
	//prepare the thing by ordering the vectors according to the event number
	std::sort( xb_track_book.begin(), xb_track_book.end(), evnt_id_comparison );
	std::sort( xb_book.begin(), xb_book.end(), evnt_id_comparison );

	//prepare the beta interpolator
	XB::b_interp interp_beta_0( xb_track_book );
	
	//some useful iterators
	XB::track_info *xbtb_begin = &xb_track_book.front();
	XB::track_info *xbtb_end = &xb_track_book.back()+1;
	XB::track_info *track_iter = xbtb_begin, *track_iter_bak;
	
	//the comparison functional
	is_event_id is_evnt;
	
	//this is a perfect candidate to do in parallel
	#pragma omp parallel
	{
	
	//some useful iterators
	XB::track_info *xbtb_begin = &xb_track_book.front();
	XB::track_info *xbtb_end = &xb_track_book.back()+1;
	XB::track_info *track_iter = xbtb_begin, *track_iter_bak;
	
	//the comparison functional
	is_event_id is_evnt;
	
	//some verbosty cosmetics
	unsigned int thread_num = 0;
	#if !(defined(__APPLE__) && defined(__clang__))
		thread_num = omp_get_thread_num();
	#endif
	if( flagger & VERBOSE && !thread_num ) printf( "Event: 0000000000" );
	
	switch( mode ){
		case RELAXED : case FASTIDIOUS : case VERY_FASTIDIOUS :
			//loop on all of them (cleverly and in parallel)
			#pragma omp for schedule( static, 10 ) 
			for( int i=0; i < xb_book.size(); ++i ){
				if( flagger & VERBOSE && !thread_num ) printf( "\b\b\b\b\b\b\b\b\b\b%010d", i );
				
				//find the track info
				is_evnt = xb_book.at(i); //set up the functional
				track_iter_bak = track_iter; //save the last valid track iter
				//find it. Note that both arrays are sorted, so we can start looking
				//from the last track_iter and reduce the field search as we progress.
				if( xb_book.at(i).evnt < track_iter->evnt ||
				    xb_book.at(i).evnt > (xbtb_end-1)->evnt ) track_iter = xbtb_end;
				else track_iter = std::find_if( track_iter, xbtb_end, is_evnt );
				
				//if the event is in the track bunch (and we didn't reach
				//the end of that vector) use the track data to correct
				if( track_iter != xbtb_end ){
					//if we are very fastidious, check if we need to jump
					if( mode == VERY_FASTIDIOUS && track_iter->n != 1 ){
						xb_book.at(i).evnt = 0;
						continue;
					}
					
					//and correct it
					XB::doppler_correct( xb_book.at(i),
					                     track_iter->beta_0,
					                     default_beam_out );
					
				//else, use the default beam out direction (given as a crystal)
				//and interpolate the beta
				} else {
					if( mode == RELAXED ){
						XB::doppler_correct( xb_book.at(i),
						                     interp_beta_0( xb_book.at(i).in_beta ),
						                     default_beam_out );
					} else {
						xb_book.at(i).evnt = 0;
					}
					track_iter = track_iter_bak;
				}
			}
			break;
		case SIMULATION : //process them all if we aren't fastidious
			//init a random sequence
			srand( time( NULL )+thread_num );
			
			//loop on all of them (cleverly and in parallel)
			#pragma omp for schedule( dynamic ) 
			for( int i=0; i < xb_book.size(); ++i ){
				
				if( flagger & VERBOSE && !thread_num ) printf( "\b\b\b\b\b\b\b\b\b\b%010d", i );
				
				//The simulation has NO in_beta information available, so we need it
				//from the track info of the run we are tailoring the simulation to.
				track_iter = &xb_track_book.at( rand()%xb_track_book.size() ); //set up the functional
				
				//and correct it
				//added some error handling (for now, only for simulations
				//but it will be extended).
				try{
					XB::doppler_correct( xb_book.at(i),
					                     track_iter->beta_0,
					                     default_beam_out );
				} catch( XB::error e ){
					//TODO: figure out how to get rid of the corrupted events
					continue;
				}
			}
			break;
	}
	
	if( flagger & VERBOSE && !thread_num ) printf( "\n" );
	
	} //parallel pragma ends here
	
	//if we were fastidious, clean up the events. In a clever way.
	int n_zeroed = 0;
	if( mode != RELAXED && mode != SIMULATION ){
		if( flagger & VERBOSE ) printf( "Pruning data...\n" );
	
		//first: sort them again, all the 0-ed events will bubble up at the beginning.
		std::sort( xb_book.begin(), xb_book.end(), evnt_id_comparison );
		
		//then, find the first non-zero event. To do so, count the 0-ed elements,
		//which are now at the beginning of the vector.
		is_evnt = 0; //get a functional for event id 0
		n_zeroed = std::count_if( xb_book.begin(), xb_book.end(), is_evnt );
		
		//finally, chop off those elements.
		xb_book.erase( xb_book.begin(), xb_book.begin()+n_zeroed );
	}
}

//------------------------------------------------------------------------------------
//implementation of the function that drives the doppler correction
//in parallel. Cluster edition.
//NOTE: for now, the fedault beam out is always used as the direction
//      in the near future, the directional information from the tracker
//      will be used when available
//NOTE2: I do apologize to whoever will read this function, it is a bit
//       more complicated than usual -especially because it's a bit of a fest
//       of double-pointers and whatever.
void apply_doppler_correction( std::vector<XB::clusterZ> &xb_book,
                               std::vector<XB::track_info> &xb_track_book,
                               unsigned int default_beam_out, correct_mode mode,
                               int flagger ){
	//ok, now do some actual work:
	//prepare the thing by ordering the vectors according to the event number
	std::sort( xb_track_book.begin(), xb_track_book.end(), evnt_id_comparison );
	std::sort( xb_book.begin(), xb_book.end(), evnt_id_comparison );

	//prepare the beta interpolator
	XB::b_interp interp_beta_0( xb_track_book );
	
	//some useful iterators
	XB::track_info *xbtb_begin = &xb_track_book.front();
	XB::track_info *xbtb_end = &xb_track_book.back()+1;
	XB::track_info *track_iter = xbtb_begin, *track_iter_bak;
	
	//the comparison functional
	is_event_id is_evnt;
	
	//this is a perfect candidate to do in parallel
	#pragma omp parallel
	{
	
	//some useful iterators
	XB::track_info *xbtb_begin = &xb_track_book.front();
	XB::track_info *xbtb_end = &xb_track_book.back()+1;
	XB::track_info *track_iter = xbtb_begin, *track_iter_bak;
	
	//the comparison functional
	is_event_id is_evnt;
	
	//some verbosty cosmetics
	unsigned int thread_num = 0;
	#if !(defined(__APPLE__) && defined(__clang__))
		thread_num = omp_get_thread_num();
	#endif
	if( flagger & VERBOSE && !thread_num ) printf( "Event: 0000000000" );
	
	switch( mode ){
		case RELAXED : case FASTIDIOUS : case VERY_FASTIDIOUS :
			//loop on all of them (cleverly and in parallel)
			#pragma omp for schedule( static, 10 ) 
			for( int i=0; i < xb_book.size(); ++i ){
				if( flagger & VERBOSE && !thread_num ) printf( "\b\b\b\b\b\b\b\b\b\b%010d", i );
				
				//find the track info
				is_evnt = xb_book.at(i); //set up the functional
				track_iter_bak = track_iter; //save the last valid track iter
				//find it. Note that both arrays are sorted, so we can start looking
				//from the last track_iter and reduce the field search as we progress.
				if( xb_book.at(i).evnt < track_iter->evnt ||
				    xb_book.at(i).evnt > (xbtb_end-1)->evnt ) track_iter = xbtb_end;
				else track_iter = std::find_if( track_iter, xbtb_end, is_evnt );
				
				//if the event is in the track bunch (and we didn't reach
				//the end of that vector) use the track data to correct
				if( track_iter != xbtb_end ){
					//if we are very fastidious, check if we need to jump
					if( mode == VERY_FASTIDIOUS && track_iter->n != 1 ){
						xb_book.at(i).evnt = 0;
						continue;
					}
					
					//and correct it
					XB::doppler_correct( xb_book.at(i),
					                     track_iter->beta_0,
					                     default_beam_out );
					
				//else, use the default beam out direction (given as a crystal)
				//and interpolate the beta
				} else {
					if( mode == RELAXED ){
						XB::doppler_correct( xb_book.at(i),
						                     interp_beta_0( xb_book.at(i).in_beta ),
						                     default_beam_out );
					} else {
						xb_book.at(i).evnt = 0;
					}
					track_iter = track_iter_bak;
				}
			}
			break;
		case SIMULATION : //process them all if we aren't fastidious
			//init a random sequence
			srand( time( NULL )+thread_num );
			
			//loop on all of them (cleverly and in parallel)
			#pragma omp for schedule( dynamic ) 
			for( int i=0; i < xb_book.size(); ++i ){
				
				if( flagger & VERBOSE && !thread_num ) printf( "\b\b\b\b\b\b\b\b\b\b%010d", i );
				
				//The simulation has NO in_beta information available, so we need it
				//from the track info of the run we are tailoring the simulation to.
				track_iter = &xb_track_book.at( rand()%xb_track_book.size() ); //set up the functional
				
				//and correct it
				//added some error handling (for now, only for simulations
				//but it will be extended).
				try{
					XB::doppler_correct( xb_book.at(i),
					                     track_iter->beta_0,
					                     default_beam_out );
				} catch( XB::error e ){
					//TODO: figure out how to get rid of the corrupted events
					continue;
				}
			}
			break;
	}
	
	if( flagger & VERBOSE && !thread_num ) printf( "\n" );
	
	} //parallel pragma ends here
	
	//if we were fastidious, clean up the events. In a clever way.
	int n_zeroed = 0;
	if( mode != RELAXED && mode != SIMULATION ){
		if( flagger & VERBOSE ) printf( "Pruning data...\n" );
	
		//first: sort them again, all the 0-ed events will bubble up at the beginning.
		std::sort( xb_book.begin(), xb_book.end(), evnt_id_comparison );
		
		//then, find the first non-zero event. To do so, count the 0-ed elements,
		//which are now at the beginning of the vector.
		is_evnt = 0; //get a functional for event id 0
		n_zeroed = std::count_if( xb_book.begin(), xb_book.end(), is_evnt );
		
		//finally, chop off those elements.
		xb_book.erase( xb_book.begin(), xb_book.begin()+n_zeroed );
	}
}

//====================================================================================
//the comparison utility
//for the event holder structure
bool evnt_id_comparison( const XB::event_holder &one, const XB::event_holder &two ){
	return one.evnt < two.evnt;
}
