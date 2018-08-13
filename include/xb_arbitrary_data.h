//this file defines a sort-of-associatve array, which will be used to store
//data that aren't really crystal ball data but support.
//The structure will inherit from the event_holder.

#ifndef XB_ARBITRARY_DATA__H
#define XB_ARBITRARY_DATA__H

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <algorithm>

#include "xb_error.h"
#include "xb_data.h"

//some defines
#define XB_ADATA_NB_FIELDS 256 //maximum number of fields supported
                               //NOTE: keep it small, but that means
                               //      hash collisions will be a problem
#define XB_ADATA_FIELD_NAME_LENGTH 16

//the hash table
//NOTE: apparently, there's a way to generate one that won't produce collisions
//      on a given set of words. This is NOT it yet, it's just a random one.
//TODO: get a proper one for the 153 fields in the land02 tree.
#define XB_PAERSON_HASH_TABLE { 171, 104, 228, 208, 188, 42, 152, 244, 137, 117, 173, 255, 201, 215, 204, 41, 74, 45, 246, 249, 91, 184, 227, 59, 64, 133, 114, 220, 122, 155, 192, 212, 43, 105, 46, 16, 77, 156, 98, 126, 191, 12, 190, 75, 55, 169, 13, 106, 84, 36, 33, 170, 61, 194, 144, 136, 178, 164, 29, 161, 108, 206, 121, 123, 129, 135, 47, 6, 146, 10, 250, 185, 239, 51, 89, 15, 49, 30, 94, 128, 193, 32, 181, 183, 142, 210, 163, 34, 67, 3, 143, 100, 230, 216, 23, 19, 97, 93, 159, 22, 124, 76, 237, 226, 162, 139, 200, 221, 252, 145, 125, 199, 70, 229, 20, 182, 154, 219, 18, 96, 119, 179, 5, 217, 160, 35, 231, 37, 56, 132, 82, 81, 115, 116, 223, 92, 65, 112, 158, 180, 2, 153, 168, 113, 209, 165, 134, 8, 167, 120, 225, 101, 148, 253, 176, 243, 205, 207, 21, 54, 247, 44, 235, 202, 172, 4, 147, 73, 40, 83, 213, 102, 196, 85, 189, 80, 53, 242, 66, 157, 109, 195, 31, 78, 48, 69, 57, 86, 186, 38, 218, 110, 79, 211, 187, 203, 0, 127, 1, 233, 9, 224, 14, 87, 7, 254, 95, 150, 88, 248, 25, 27, 107, 63, 149, 138, 50, 238, 240, 198, 24, 251, 17, 118, 140, 166, 68, 90, 99, 197, 111, 130, 234, 174, 245, 232, 62, 28, 60, 214, 151, 175, 222, 52, 72, 11, 58, 39, 177, 131, 71, 236, 26, 103, 241, 141 }

//a define to get the header size of the linear buffer
#define nf2hdr_size( nf ) sizeof(event_holder) + ((nf)+2)*sizeof(int) + (nf)*sizeof(adata_field)

/*TODO LIST:
0)   adapt for pointer ofssets instead of pointers [x]
1)   adapt for pointed indexer                     [x]
1.5) implement convenience into indexer            [x]
2)   implement un- and subscribe methods of adata  [x]
3)   implement adata_uniarr                        [x]
4)   implement the merge function                  [x]
*/

namespace XB{
    //----------------------------------------------------------------------------
    //an ugly constant to hold the hash table
    const unsigned char adata_pT[256] = XB_PAERSON_HASH_TABLE;
    
	//----------------------------------------------------------------------------
	//declare the uniform array
	class _xb_arbitrary_data_uniform_array;
	
	//----------------------------------------------------------------------------
	//a data structure representing the field,
	typedef struct _xb_arb_data_field {
		char name[XB_ADATA_FIELD_NAME_LENGTH];
		short size;
	} adata_field;
    
    //a data structure to do indexing (also the offset table can be shared!)
    //NOTE: now an offset of -1 means empty!
    typedef class _xb_arb_data_indexer {
        public:
            _xb_arb_data_indexer(): names( 0 ) {
                for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ) diffs[i] = -1; };
            _xb_arb_data_indexer( const unsigned n ): names( n ) {
                for( int i=0; i < XB_ADATA_NB_FIELDS; ++i ) diffs[i] = -1; };
            _xb_arb_data_indexer( const _xb_arb_data_indexer &given ): names( given.names ) {
                memcpy( diffs, given.diffs, XB_ADATA_NB_FIELDS*sizeof(int) ); };
            _xb_arb_data_indexer &operator=( const _xb_arb_data_indexer &right ){
                names = right.names;
                memcpy( diffs, right.diffs, XB_ADATA_NB_FIELDS*sizeof(int) );
                return *this;
            };
            //NOTE on comparison: to KISS, a difference in the ordering of the
            //     fields is flagged as a difference in the indexers
            //     this also means that the concatenation order _matters_
            bool operator!=( const _xb_arb_data_indexer &right );
            bool operator==( const _xb_arb_data_indexer &right ){ return !( *this != right ); };
            _xb_arb_data_indexer &operator+( const _xb_arb_data_indexer &right );
            
            unsigned size() const { return names.size(); };
            adata_field &operator[]( unsigned i ){ return names[i]; };
            adata_field &at( unsigned i ){ return names.at(i); };
            std::vector< adata_field > names;
            int diffs[XB_ADATA_NB_FIELDS];
    } adata_indexer;

	//----------------------------------------------------------------------------
	//the main thing, the class
	typedef class _xb_arbitrary_data : public event_holder {
		public:
			friend class _xb_arbitrary_data_uniform_array;
			
			//ctors, dtor
			_xb_arbitrary_data();
			_xb_arbitrary_data( const adata_field *fld_array, size_t n_fld );
			_xb_arbitrary_data( adata_indexer *indexer );
			_xb_arbitrary_data( const adata_indexer &indexer );
			_xb_arbitrary_data( const _xb_arbitrary_data &given );
			~_xb_arbitrary_data();
			
			//important operators
			_xb_arbitrary_data &operator=( const _xb_arbitrary_data &right );
			_xb_arbitrary_data &copy() const; //returns a copy of this
			                                  //but with its own indexer & no ua.
			//get data from field, by name.
			//use fsize( char *name ) to ge the returned buffer size
			void *operator()( const char *name ) const;
			bool operator==( const _xb_arbitrary_data &right ) const;
			bool operator!=( const _xb_arbitrary_data &right ) const;
			
			//accessing methods:
			//create/write field
			//if *buf is NULL, the field is just created
			void dofield( const char *name, short size, void *buf );
			void dofield( const adata_field &fld, void *buf );
			//get the size of a field
			//remove a field
			void rmfield( const char *name );
			void clear(); //remove everything.
			//for size-1 fields of a specified type
			//you can use this themplate mehtod, too
			template< class T >
			T tip( const char *name ) const {
                		unsigned char i_fld = phash8( name );
				void *head = (char*)_buf + _fields->diffs[i_fld];
				if( head < _buf || _fields->diffs[i_fld] >= _buf_sz )
					return NULL;
				head = (int*)head + 1;
				return *(T*)head;
			};
			
			//this should be a fast access method
			//since we laugh in the face of danger,
			//no out of bound check is provided.
			//you'll just loop in the array
			template< class T >
			inline T &at( const char *name, int i ){ //note that this is NOT const
				void *head = (char*)_buf + _fields->diffs[phash8( name )];
				if( head < _buf ) throw error( "Not a field!", "XB::adata::getfield" );
				int len = *(int*)head;
				head = (int*)head + 1;
				return ((T*)head)[i%len];
			}
			
			//Indexer accessing and operations
			std::vector< adata_field > lsfields() const { return _fields->names; };
			adata_indexer *get_indexer() { return _fields; };
			int nf() const { return _fields->names.size(); }; //count them
			int fsize( const char *name ) const; //get the size of one
			
			//a couple of friends, for I/O ops
			friend int adata_getlbuf( void **linbuf, const _xb_arbitrary_data &given );
			friend int adata_fromlbuf( _xb_arbitrary_data &here,
			                           const void *buf,
			                           adata_indexer *indexer );
			friend _xb_arbitrary_data adata_merge( const _xb_arbitrary_data &one,
			                                       const _xb_arbitrary_data &two );
			friend void adata_put( FILE *stream, const _xb_arbitrary_data &given );
		private:
			//the data buffer
			//data is stored [int size|data]
			//the field pointer point to the size!!!
			//can contain any type.
			int _buf_sz;
			void *_buf;
            		adata_indexer *_fields; //it's the indexer!
			_xb_arbitrary_data_uniform_array *_parent_array;
            		char _is_fields_owned;
			
			//an utility to has a field name
			//pearson's has, 8 bits.
			static unsigned char phash8( const char *name );
			
			//un- and subscription to an uniform array
			void subscribe_uniarr( _xb_arbitrary_data_uniform_array *ua );
			void unsubscribe_uniarr();
	} adata;

	//----------------------------------------------------------------------------
    	//an array of UNIFORM adata entries --this should be the way to use this
    	//struct since the indexer has been rendered external
    	typedef std::vector< adata > _ua_type;
	typedef _ua_type::iterator _ua_iter;
    	typedef class _xb_arbitrary_data_uniform_array {
		public:
            		friend class _xb_arbitrary_data;
            
			_xb_arbitrary_data_uniform_array() {};
			_xb_arbitrary_data_uniform_array( const adata_indexer &indexer ):
                	_indexer( indexer ) {};
			_xb_arbitrary_data_uniform_array( const unsigned nb, const adata_indexer &idx );
			_xb_arbitrary_data_uniform_array( const _xb_arbitrary_data_uniform_array &given );
			_xb_arbitrary_data_uniform_array &operator=( const _xb_arbitrary_data_uniform_array& );
			_xb_arbitrary_data_uniform_array &operator+( const _xb_arbitrary_data_uniform_array &right );
			
			//NOTE: invoking set_indexer will reset ALL members!
			void set_indexer( const adata_indexer &given );
			void push_indexer( adata_field &given ); //will add a field to all elements
			adata_field pop_indexer(); //pops the last element of indexer
			adata_indexer get_indexer() { return _indexer; };
			//NOTE: where are the signle field operations?
			//      DANGER all members of an array SHARE an indexer and can modify it
			
			//goodies with un- and subscription
			void push_back( const adata &given );
			adata pop_back();
			_ua_iter insert( _ua_iter pos, adata &val );
			_ua_iter insert( _ua_iter pos, _ua_iter first, _ua_iter last );
			_ua_iter erase( _ua_iter pos );
			_ua_iter erase( _ua_iter first, _ua_iter last );
			
			//some forwarding of the underlying vector
			adata &operator[]( unsigned i ){ return _ua[i]; };
			adata operator[]( unsigned i  ) const { return _ua[i]; };
			adata &at( unsigned i ){ return _ua.at(i); };
			adata at( unsigned i ) const { return _ua.at(i); };
			unsigned size() const { return _ua.size(); };
			_ua_iter begin() { return _ua.begin(); };
			_ua_iter end() { return _ua.end(); };
			adata &front(){ return _ua.front(); };
			adata &back(){ return _ua.back(); };
			bool empty() const { return _ua.empty(); };
			void clear(){ _ua.clear(); _indexer = adata_indexer(); };
			void resize( unsigned n ){ _ua.resize( n ); };
		private:
			adata_indexer _indexer;
			_ua_type _ua;
	} adata_uniarr;
	
	//----------------------------------------------------------------------------
	//Two functions to be used in read/write operations (friended by the class)
	//they convert it in and from a linear buffer, ready to be written to file.
	//adata_getlbuf: get the linear buffer in void *buffer and return the size.
	//               buffer will be allocated. An error is thrown if it's already
	//               allocated.
	int adata_getlbuf( void **linbuf, const _xb_arbitrary_data &given );
	int adata_fromlbuf( _xb_arbitrary_data &here, const void *buf,
                        adata_indexer *indexer = NULL );
	_xb_arbitrary_data adata_merge( const _xb_arbitrary_data &one,
	                                const _xb_arbitrary_data &two );
	void adata_put( FILE *stream, const adata &given );

}
#endif
