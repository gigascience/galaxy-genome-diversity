#!/usr/bin/env python

import sys

def die( message ):
    print >> sys.stderr, message
    sys.exit(1)

def open_or_die( filename, mode='r', message=None ):
    if message is None:
        message = 'Error opening {0}'.format( filename )
    try:
        fh = open( filename, mode )
    except IOError, err:
        die( '{0}: {1}'.format( message, err.strerror ) )
    return fh

class LocationFile( object ):
    def __init__( self, filename, comment_chars=None, delimiter='\t', key_column=0 ):
        self.filename = filename
        if comment_chars is None:
            self.comment_chars = ( '#' )
        else:
            self.comment_chars = tuple( comment_chars )
        self.delimiter = delimiter
        self.key_column = key_column
        self._map = {}
        self._populate_map()

    def _populate_map( self ):
        try:
            with open( self.filename ) as fh:
                line_number = 0
                for line in fh:
                    line_number += 1
                    line = line.rstrip( '\r\n' )
                    if not line.startswith( self.comment_chars ):
                        elems = line.split( self.delimiter )
                        if len( elems ) <= self.key_column:
                            die( 'Location file {0} line {1}: less than {2} columns'.format( self.filename, line_number, self.key_column + 1 ) )
                        else:
                            key = elems.pop( self.key_column )
                            if key in self._map:
                                if self._map[key] != elems:
                                    die( 'Location file {0} line {1}: duplicate key "{2}"'.format( self.filename, line_number, key ) )
                            else:
                                self._map[key] = elems
        except IOError, err:
            die( 'Error opening location file {0}: {1}'.format( self.filename, err.strerror ) )

    def get_values( self, key ):
        if key in self._map:
            rval = self._map[key]
            if len( rval ) == 1:
                return rval[0]
            else:
                return rval
        else:
            die( 'key "{0}" not found in location file {1}'.format( key, self.filename ) )

    def get_values_if_exists( self, key ):
        if key in self._map:
            rval = self._map[key]
            if len( rval ) == 1:
                return rval[0]
            else:
                return rval
        else:
            return None
