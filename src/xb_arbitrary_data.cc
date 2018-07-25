//implementation of xb_arbitrary_data.h

#include "xb_arbitrary_data.h"

namespace XB{

	//============================================================================
	//adata_indexer operators

	//----------------------------------------------------------------------------
	//comparison. See the NOTE on field ordering in the header!
	bool adata_indexer::operator!=( adata_indexer &right ){
		if( names.size() != right.names.size() ) return true;
		if( memcmp( diffs, right.diffs, XB_ADATA_NB_FIELDS*sizeof( unsigned short ) ) )
			return true;
		for( int i=0; i < names.size(); ++i )
			if( strcmp( names[i].name, right.names[i].name ) ) return true;
		return false;
	}

	//----------------------------------------------------------------------------
	//concatenation. Remember that order matters!
	adata_indexer &adata_indexer::operator+( adata_indexer &right ){
		//first merge the diff tables. If there are conflicts, then throw
		short unsigned candidate = 0;
		for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ){
			candidate = diffs[i] + right.diffs[i];
			if( candidate != diffs[i] && candidate != right.diffs[i] )
				throw( error( "Colliding fields or integer overflow!", "XB::adata_indexer" ) );
			diffs[i] = candidate;
		}
		names = names + right.names;
		
		return *this;
	}

	//============================================================================
	//_xb_arbitrary_data methods

	//----------------------------------------------------------------------------
	//constructors:
	adata::_xb_arbitrary_data():
		_buf( NULL ),
		_buf_sz( 0 ),
		_fields( new adata_indexer ),
        _is_fields_owned( 1 )
	{
		//just banally init the event_holder members
		n = 0;
		evnt = 0;
		tpat = 0;
		in_Z = 0;
		in_A_on_Z = 0;
		
		for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ) _fields->diffs[i] = 0;
	}
	
	adata::_xb_arbitrary_data( const adata_field *fld_array, size_t n_fld ):
		_buf( NULL ),
		_buf_sz( 0 ),
		_fields( new adata_indexer( n_fld ) ),
		_is_fields_owned( 1 )
	{
		n = 0;
		evnt = 0;
		tpat = 0;
		in_Z = 0;
		in_A_on_Z = 0;
		
		for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ) _fields->diffs[i] = 0;
		for( int i=0; i < n_fld; ++i ) dofield( fld_array[i], NULL );
	}
	
	adata::_xb_arbitrary_data( const adata &given ):
		_buf( NULL ),
		_buf_sz( given._buf_sz ),
		_is_fields_owned( given._is_fields_owned ) //ownership is also copied
	{
		n = given.n;
		evnt = given.evnt;
		tpat = given.tpat;
		in_Z = given.in_Z;
		in_A_on_Z = given.in_A_on_Z;
		
		//copy the indexer
		if( given._is_fields_owned ) //then we also need it owned, copy it
			_fields = new adata_indexer( *given._fields );
		else //if it's not owned, then we just copy the pointer to the external one
			_fields = given._fields;
	
		//copy the buffer
		_buf = malloc( given._buf_sz );
		if( !_buf ) throw error( "Memory error!", "XB::adata::assign" );
		memcpy( _buf, given._buf, given._buf_sz );
	}
	
	//TODO: there's a contructor missing.
	
	//dtor:
	adata::~_xb_arbitrary_data(){
		if( _buf ) free( _buf );
        if( _is_fields_owned ) delete _fields;
	}
	
	//----------------------------------------------------------------------------
	//utils:
	//----------------------------------------------------------------------------
	//now the fun begins: the hash function
	unsigned char adata::phash8( const char *name ) const {
		int len = strlen( name );
		unsigned short h = 0;
		
		for( int i=0; i < len; ++i ) h = adata_pT[ h ^ name[i] ];
		
		return h;
	}
	
	//----------------------------------------------------------------------------
	//operators:
	//----------------------------------------------------------------------------
	//assignment operator
	adata &adata::operator=( const adata &given ){
		n = given.n;
		evnt = given.evnt;
		tpat = given.tpat;
		in_Z = given.in_Z;
		in_A_on_Z = given.in_A_on_Z;
		
		_buf_sz = given._buf_sz;

		//copy the indexer
		if( given._is_fields_owned ) //then we also need it owned, copy it
			_fields = new adata_indexer( *given._fields );
		else //if it's not owned, then we just copy the pointer to the external one
			_fields = given._fields;
		_is_fields_owned = given._is_fields_owned;

		//copy the buffer
		if( _buf ) free( _buf );
		_buf = malloc( given._buf_sz );
		if( !_buf ) throw error( "Memory error!", "XB::adata::assign" );
		memcpy( _buf, given._buf, given._buf_sz );
		
		return *this;
	}
	
	//----------------------------------------------------------------------------
	//comparison ops
	bool adata::operator==( const adata &right ) const {
		if( n != right.n) return false;
		if( evnt != right.evnt) return false;
		if( tpat != right.tpat) return false;
		if( in_Z != right.in_Z) return false;
		if( in_A_on_Z != right.in_A_on_Z ) return false;

		if( _buf_sz != right._buf_sz ) return false;
		if( _fields->size() != right._fields->size() ) return false;
		
		for( int i=0; i < _fields->size(); ++i ){
			if( strcmp( _fields->names[i].name, right._fields->names[i].name ) ||
			    _fields->names[i].size != right._fields->names[i].size ) return false;
		}
		
		if( memcmp( _buf, right._buf, _buf_sz ) ) return false;
		
		return true;
	}
	
	bool adata::operator!=( const adata &right ) const {
		return !( *this == right );
	}
	
	//----------------------------------------------------------------------------
	//accessors:
	
	//----------------------------------------------------------------------------
	//parenthesis operator
	void *adata::operator()( const char *name ) const {
		void *head = (char*)_buf + _fields->diffs[phash8( name )];
		if( !head ) throw error( "Field is empty!", "XB::adata::()" );
		
		int fsize = *(int*)head; //in bytes!
		head = (int*)head + 1;
		
		void *payload = malloc( fsize );
		if( !payload ) throw error( "Memory error!", "XB::adata::()" );
		
		memcpy( payload, head, fsize );
		return payload;
	}
	
	//----------------------------------------------------------------------------
	//access in write a field denoted by adata_field
	void adata::dofield( const adata_field &fld, void *buf ){
		unsigned char i_fld = phash8( fld.name );
		void *head = (char*)_buf + _fields->diffs[i_fld];
		
		if( _fields->diffs[i_fld] || _fields->diffs[i_fld] >= _buf_sz ){
			_buf = realloc( _buf, _buf_sz+fld.size );
			
			//save the new field poiner
			_fields->diffs[i_fld] = _buf_sz;
			head = (char*)_buf + _fields->diffs[i_fld]; //and put it in the head
			*(int*)head = fld.size; //write the size
			head = (int*)head + 1; //move the head
			
			//if there's something to copy, do it
			//else zero the memory (safer)
			if( buf ) memcpy( head, buf, fld.size );
			else memset( head, 0, fld.size );
			
			//finally, update the buffer size
			//and push the new field in the field list
			_buf_sz += fld.size + sizeof(int);
			if( _is_fields_owned ) _fields->names.push_back( fld );
		} else { //the field is populated
			if( !buf ) return; //do nothing
			if( fld.size != *(int*)head ) //freak out
				throw error( "Wrong field size!", "XB::adata::dofield" );
			
			head = (int*)head + 1; //move the head past the size
			memcpy( head, buf, fld.size ); //copy
		}
	}
	
	void adata::dofield( const char *name, short size, void *buf ){
		adata_field fld = { "", size };
		strncpy( fld.name, name, 16 );
		dofield( fld, buf );
	}
	
	//----------------------------------------------------------------------------
	//get the size of a field (in bytes)
	int adata::fsize( const char *name ) const {
		void *head = (char*)_buf + _fields->diffs[phash8( name )];
		if( !head ) return 0;
		return *(int*)head;
		return 0;
	}
	
	//----------------------------------------------------------------------------
	//remove a field (another interesting method)
	void adata::rmfield( const char *name ){
		i_fld = _fields->diffs[phash8( name )];
		if( _fields->diffs[i_fld] < 0 || _fields->diffs[i_fld] >= _buf_sz ) return;
		
		//work out where head is in the buffer
		int fsize = *(int*)head + sizeof(int);
		int from_front = (char*)head - (char*)_buf; //how far from the front
		int to_back = _buf_sz - fsize - from_front; //how far from the back
		
		//find the beginning from the head
		void *rest = (char*)head + fsize;
		memmove( head, rest, to_back );
		
		//move the field pointers
		_fields->diffs[phash8( name )] = 0; //remove the one
		for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ){
			if( !_fields->diffs[i] || _fields->diffs[i] < rest ) continue;
			
			//move them back by fsize (the removed bit)
			_fields->diffs[i] = _fields->diffs[i] - fsize;
		}
		
		//resize the buffer
		_buf_sz -= fsize;
		_buf = realloc( _buf, _buf_sz );
		
		//drum out the field from the field list
		//only if you own the thing
		if( _is_fields_owned ){
			from_front = 0; //recycle
			while( strcmp( (_fields->names.begin()+from_front)->name, name ) ) ++from_front;
			_fields->names.erase( _fields->names.begin() + from_front );
		}
	}
	
	//----------------------------------------------------------------------------
	//clear the structure
	//NOTE: with the shared indexer, this is DANGER-ous
	void adata::clear(){
		n = 0;
		tpat = 0;
		in_Z = 0;
		in_A_on_Z = 0;
		if( _buf ){ free( _buf ); _buf = NULL; }
		_buf_sz = 0;
		//fat fingered safety measure...
		if( _is_fields_owned ){
			for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ) _fields->diffs[i] = 0;
			_fields->names.clear();
		}
	}

	//----------------------------------------------------------------------------
	//un- and subscribe methods
	void adata::subscribe_uniarr( adata_uniarr *ua ){
		if( !_fields ){
			_fields = ua->_indexer;
			_is_fields_owned = 0;
		} else if( _fields && *_fields == ua->_indexer ){
			if( _is_fields_owned ) delete _fields;
			_fields = _indexer;
			_is_fields_owned = 0;
		} else throw( error( "Not uniform to uniform array!", "XB::adata" ) );
	}

	void adata::unsubscribe_uniarr(){
		adata_indexer *old = _fields;
		_fields = new adata_indexer( *old );
		_is_fields_owned = 1;
	}
	
	//============================================================================
	//implementation of the adata_uniarr class
	
	//----------------------------------------------------------------------------
	//constructors:
	
	adata_uniarr::adata_uniarr( const unsigned nb, const adata_indexer &idx ):
		_indexer( idx ),
		_ua( adata( &idx ) )
	{}
	
	adata_uniarr::adata_uniarr( const adata_uniarr &given ):
		_indexer( given._indexer ),
		_ua( given._ua ) //this will copy everything but the pointer to the indexer will be wrong
	{
		//set them right without resetting the data structs.
		for( int i=0; i < _ua.size(); ++i ) _ua[i]._fields = &_indexer;
	}

	//----------------------------------------------------------------------------
	//set the indexer to an existing uniform array
	//note that this will destroy the content if it's different from the
	//previous one
	void adata_uniarr::set_indexer( const adata_indexer &given ){
		if( _indexer != indexer ) clear();
		_indexer = indexer;
    }
	
	//----------------------------------------------------------------------------
	//push a field onto the whole array
	void adata_uniarr::push_indexer( adata_field &given ){
		//indexer shared --> control to this class!
		if( _indexer.diffs[phash8( given.name ) < 0 ) _indexer.names.push_back( given );
		
		//do we touch all the events, so that we don't need to discover at the moment's
		//notice that we need to allocate memory? Yeeah
		//NOTE: this will also set the correct diff.
		for( int i=0; i < size(); ++i ) this->at(i).dofield( given, NULL );
	}
	
	//----------------------------------------------------------------------------
	//remove a field from all entries
	adata_field adata_uniarr::pop_indexer(){
		if( _indexer.diffs[phash8( given.name ) < 0 ) return { "", 0 }; //empty field
		adata_field fld = _indexer.names.back();
		for( int i=0; i < size(); ++i ) this->at(i).rmfield( fld );
		_indexer.names.pop_back();
		_indexer.diffs[phash8( fld.name )] = -1;
		return fld;
	}
	
	//----------------------------------------------------------------------------
	//push and SUBSCRIBE an arbitrary data
    void adata_uniarr::push_back( const adata &given ){
		this->push_back( given );
		back().subscribe_uniarr( this );
	}
	
	//----------------------------------------------------------------------------
	//pop and unsubscribe (pop_back will destroy the element)
    adata adata_uniarr::pop_back(){
		adata ad = this->back();
		this->pop_back();
		ad->unsubscribe_uniarr();
		return ad;
	}
	
	//----------------------------------------------------------------------------
	_ua_iter adata_uniarr::insert( _ua_iter pos, adata &val ){
		val.subscribe_uniarr( this );
		return _ua.insert( pos, val );
	}
	
	//----------------------------------------------------------------------------
	_ua_iter adata_uniarr::insert( _ua_iter pos, _ua_iter first, _ua_iter last ){
		for( _ua_iter i=firt; i != last; ++i ) i->subscribe_uniarr( this );
		return _ua.insert( pos, first, last );
	};
	
	//----------------------------------------------------------------------------
	_ua_iter adata_uniarr::erase( _ua_iter pos ){
		pos->unsubscribe_uniarr();
		return _ua.erase( pos );
	}
	
	//----------------------------------------------------------------------------
	_ua_iter adata_uniarr::erase( _ua_iter first, _ua_iter last ){
		for( _ua_iter i=first; i != last; ++i ) i->unsubscribe_uniarr();
		return _ua.erase( first, last );
	}
	
	//----------------------------------------------------------------------------
	//operators:
	
	//----------------------------------------------------------------------------
	//merge horizontally two uniform arrays (operator+)
	adata_uniarr &adata_uniarr::operator+( const adata_uniarr &right ){
		if( right.empty() || _ua.empty() ) return *this;
		if( right.size() != _ua.size() )
			throw( error e( "Mismatching lengths!", "XB::adata_uniarr" ) );
		
		
		adata_indexer iright = *right._indexer;
		for( int i=0; i < iright.size(); ++i )
			if( iright.diffs[i] >= 0 ) iright.diffs[i] += right._buf_sz;
		_indexer = _indexer + iright;
		
		_ua.resize( ( _ua.size() >= right.size() )? _ua.size() : right.size() );
		unsigned new_buf_sz = _ua[0]._buf_sz + right._ua[0]._buf_sz;
		unsigned right_buf_sz = right._ua[0]._buf_sz;
		
		for( int i=0; i < _ua.size(); ++i ){
			if( !memcmp( &_ua[i].n, &right[i].n, sizeof( event_holder ) ) )
				throw( error( "Structures relate to different events!", "XB::adata_uniarr" ) );
			_ua[i]._buf = realloc( _ua[i]._buf, new_buf_sz );
			memcpy( (char*)_ua[i]._buf + _ua[i]._buf_sz, right[i]._buf, right_buf_sz );
		}
		
		return *this;
	}
	
	//----------------------------------------------------------------------------
	//assignment
	adata_uniarr &adata_uniarr::operator=( const adata_uniarr &right ){
		_indexer = right._indexer;
		_ua = given._ua; //same as the copy ctor, contents will be correct but pointers wrong
		//and set them
		for( int i=0; i < _ua.size(); ++i ) _ua[i]._fields = &_indexer;
		
		return *this;
	}
	
	//============================================================================
	//the three friend functions.

	//----------------------------------------------------------------------------
	//make the linearized buffer:
	//[event_holder| # fields|field list|field pointer deltas|data size|data buffer]
	int adata_getlbuf( void **linbuf, const adata &given ){
		int nf = given._fields->size();
		
		//reorder the deltas (just handy)
		int *deltas = (int*)calloc( nf, sizeof(int) );
        for( int i=0; i < nf; ++i )
            deltas[i] = given._fields->diffs[given.phash8( given._fields->names[i].name )];
		
		//allocate the linear buffer
		int bsize = sizeof(event_holder) + (nf+2)*sizeof(int) +
		            nf*sizeof(adata_field) + given._buf_sz;
		*linbuf = malloc( bsize );
		
		//do the copying
		void *head = *linbuf;
		memcpy( head, &given.n, sizeof(event_holder) ); //the event holder
		head = (event_holder*)head + 1;
		*(int*)head = nf; //# fields
		head = (int*)head + 1;
		memcpy( head, &given._fields->names[0]), nf*sizeof(adata_field) ); //field list
		head = (adata_field*)head + nf;
		memcpy( head, deltas, nf*sizeof(int) ); //deltas
		head = (int*)head + nf;
		*(int*)head = given._buf_sz; //data size
		head = (int*)head + 1;
		memcpy( head, given._buf, given._buf_sz ); //data
		
		free( deltas );
		
		return bsize;
	}
	
	//----------------------------------------------------------------------------
	//now from the linear buffer to the structure
	int adata_fromlbuf( adata &given, const void *buffer, adata_indexer *indexer ){
		void *hdr = (char*)buffer + sizeof(event_holder);
		int nf = *(int*)hdr;
		hdr = (int*)hdr + 1;
		int hdr_sz = nf2hdr_size( nf ) -1*sizeof(int);
		
		void *data = (char*)buffer + hdr_sz;
		int data_sz = *(int*)data;
		data = (int*)data + 1;
		
		//clear the structure
		given.clear();
		
		//copy the event holder
		memcpy( &given.n, buffer, sizeof(event_holder) );
		
		//copy the data
		given._buf_sz = data_sz;
		given._buf = malloc( data_sz );
		memcpy( given._buf, data, data_sz );
		
		//copy the field list
        adata_indexer fields( nf );
        memcpy( &fields.names[0], hdr, nf*sizeof(adata_field) );
        if( !indexer ){
            given._fields = new adata_indexer( fields );
            given._is_fields_owned = 1;
        } else {
            *indexer = fields;
            given._fields = indexer;
            given._is_fields_owned = 0;
        }
		
		//reconstruct the pointer offset map
		hdr = (adata_field*)hdr + nf;
		for( int i=0; i < fields.size(); ++i ){
			given._fields->diffs[given.phash8( fields[i].name )] = *((int*)hdr+i);
		}
		
		return nf; //useless...
	}
	
	//----------------------------------------------------------------------------
	//the structure merger!
	//NOTE: the merger is ORDERED. merges two into one
	adata adata_merge( const adata &one, const adata &two ){
		if( !memcmp( &one.n, &two.n, sizeof( event_holder ) ) )
			throw( error( "structures relate to  events!", "XB::adata_merge" ) );
		
		//make the new indexer
		adata merged( one );
		merged._fields = new adata_indexer( one._fields );
		merged._is_fields_owned = 1;
		adata_indexer itwo = *two._fields;
		
		for( int i=0; i < itwo.names.size(); ++i )
			if( itwo.diffs[i] >= 0 ) itwo.diffs[i] += one._buf_sz;
		
		*merged.fields = *merged._fields + itwo; //NOTE: this will throw if not compatible
		merged._buf = realloc( merged._buf, merged._buf_sz + two._buf_sz );
		memcpy( (char*)merged._buf + merged._buf_sz, two._buf, two._buf_sz );
		
		return merged;
	}
	
} //end of namespace
