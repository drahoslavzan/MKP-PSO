/*
 * File:    get_event_component.c
 * Author:  Vince Weaver
 *	        vweaver1@eecs.utk.edu
 */

/*
  This test makes sure PAPI_get_event_component() works
*/

#include "papi_test.h"

int
main( int argc, char **argv )
{
	
    int i;
    int retval;
    PAPI_event_info_t info;
    int numcmp, cid, our_cid;

    /* Set TESTS_QUIET variable */
    tests_quiet( argc, argv );

    /* Init PAPI library */
    retval = PAPI_library_init( PAPI_VER_CURRENT );
    if ( retval != PAPI_VER_CURRENT ) {
       test_fail( __FILE__, __LINE__, "PAPI_library_init", retval );
    }

    numcmp = PAPI_num_components(  );


    /* Loop through all components */
    for( cid = 0; cid < numcmp; cid++ ) {

       i = 0 | PAPI_NATIVE_MASK;
       retval = PAPI_enum_cmp_event( &i, PAPI_ENUM_FIRST, cid );

       do {
          retval = PAPI_get_event_info( i, &info );
	  our_cid=PAPI_get_event_component(i);

	  if (our_cid!=cid) {
             test_fail( __FILE__, __LINE__, "component mismatch", 1 );
	  }

	  if (!TESTS_QUIET) {
	    printf("%d %d %s\n",cid,our_cid,info.symbol);
	  }

	  
       } while ( PAPI_enum_cmp_event( &i, PAPI_ENUM_EVENTS, cid ) == PAPI_OK );

    }

    test_pass( __FILE__, NULL, 0 );
   
    return 0;
}
