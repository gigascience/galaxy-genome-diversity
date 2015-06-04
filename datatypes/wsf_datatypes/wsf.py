"""
SnpFile datatype
"""

import galaxy.datatypes.data
import tempfile
import os
import json
from galaxy import util
from galaxy.datatypes.sniff import *
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.images import Html
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement

class Wped( Html ):
    allow_datatype_change = False
    composite_type = 'basic'
    file_ext = 'gd_ped'

    MetadataElement( name="base_name", desc="base name for all transformed versions of this genetic dataset", default='WpedData', readonly=True, set_in_upload=True )

    def __init__( self, **kwd ):
        Html.__init__( self, **kwd )
        self.add_composite_file( '%s.ped', description = 'Pedigree File', substitute_name_with_metadata = 'base_name', is_binary = False )
        self.add_composite_file( '%s.map', description = 'Map File', substitute_name_with_metadata = 'base_name', is_binary = False )

class Individuals( Tabular ):
    file_ext = 'gd_indivs'
    def __init__(self, **kwd):
        Tabular.__init__( self, **kwd )
        self.column_names = [ 'Column', 'Name', 'Alias' ]

    def display_peek( self, dataset ):
        return Tabular.make_html_table( self, dataset, column_names=self.column_names )

class DatasetComments( object ):
    def __init__( self, dataset, comment_string='#' ):
        self.dataset = dataset
        self.comment_string = comment_string
        self.comment_string_len = len(comment_string)
        self._comments = []
        self._read_comments()

    def _read_comments( self ):
        if self.dataset.has_data():
            try:
                for line in open(self.dataset.file_name, 'rU'):
                    if line.startswith(self.comment_string):
                        comment = line[self.comment_string_len:]
                        self._comments.append(comment)
                    else:
                        break
            except:
                pass

    def __str__( self ):
        return "".join(self._comments)

    @property
    def comments( self ):
        return self._comments

class DatasetCommentMetadata( object ):
    def __init__( self, dataset, comment_string='#' ):
        self.dataset_comments = DatasetComments( dataset, comment_string )
        self._comment_metadata = {}
        self._decode_dataset_comments()

    def _decode_dataset_comments( self ):
        dataset_comment_string = str( self.dataset_comments )
        try:
            self._comment_metadata = simplejson.loads( dataset_comment_string )
        except simplejson.JSONDecodeError as e:
            pass

    @property
    def comment_metadata( self ):
        return self._comment_metadata

class AnnotatedTabular( Tabular ):
    """ Tabular file with optional comment block containing JSON to be imported into metadata """
    MetadataElement( name="comment_metadata", desc="comment metadata", param=metadata.DictParameter, visible=False, readonly=True )

    def set_meta( self, dataset, overwrite = True, **kwd ):
        Tabular.set_meta( self, dataset, overwrite=overwrite, max_data_lines=None, max_guess_type_data_lines=1000, **kwd )
        if dataset.metadata.comment_metadata is None:
            dataset_comment_metadata = DatasetCommentMetadata( dataset )
            dataset.metadata.comment_metadata = dataset_comment_metadata.comment_metadata.copy()
            self.set_dataset_metadata_from_comments( dataset )

    def set_dataset_metadata_from_comments( self, dataset ):
        pass

    def set_peek( self, dataset, line_count=None, is_multi_byte=False ):
        super(Tabular, self).set_peek( dataset, line_count=line_count, is_multi_byte=is_multi_byte, WIDTH='unlimited', skipchars=['#'] )

    def display_peek( self, dataset ):
        """Returns formated html of peek"""
        return Tabular.make_html_table( self, dataset, skipchars=['#'] )

class Fake( AnnotatedTabular ):
    MetadataElement( name="scaffold", desc="scaffold column", param=metadata.ColumnParameter, default=0 )
    MetadataElement( name="pos", desc="pos column", param=metadata.ColumnParameter, default=0 )
    MetadataElement( name="ref", desc="ref column", param=metadata.ColumnParameter, default=0 )
    MetadataElement( name="rPos", desc="rPos column", param=metadata.ColumnParameter, default=0 )
    MetadataElement( name="species", desc="species", default='', no_value='', visible=False, readonly=True )

    def set_dataset_metadata_from_comments( self, dataset ):
        self.set_dataset_column_names_metadata( dataset )
        self.set_dataset_columnParameter_metadata( dataset )
        self.set_dataset_species_metadata( dataset )
        self.set_dataset_dbkey_metadata( dataset )

    def set_dataset_column_names_metadata( self, dataset ):
        value_from_comment_metadata = dataset.metadata.comment_metadata.get( 'column_names', None )
        if isinstance( value_from_comment_metadata, list ):
            dataset.metadata.column_names = value_from_comment_metadata[:]

    def set_dataset_columnParameter_metadata( self, dataset ):
        for name, spec in dataset.metadata.spec.items():
            if isinstance( spec.param, metadata.ColumnParameter ):
                value_from_comment_metadata = dataset.metadata.comment_metadata.get( name, None )
                if value_from_comment_metadata is not None:
                    try:
                        i = int( value_from_comment_metadata )
                    except:
                        i = 0
                    if 0 <= i <= dataset.metadata.columns:
                        setattr( dataset.metadata, name, i )

    def set_dataset_species_metadata( self, dataset ):
        value_from_comment_metadata = dataset.metadata.comment_metadata.get( 'species', None )
        if isinstance( value_from_comment_metadata, basestring ):
            dataset.metadata.species = value_from_comment_metadata

    def set_dataset_dbkey_metadata( self, dataset ):
        value_from_comment_metadata = dataset.metadata.comment_metadata.get( 'dbkey', '?' )
        if isinstance( value_from_comment_metadata, basestring ):
            dataset.metadata.dbkey = value_from_comment_metadata

class GDSnp( Fake ):
    """ Webb's SNP file format """
    file_ext = 'gd_snp'

    MetadataElement( name="individual_names", desc="individual names", visible=False, readonly=True )
    MetadataElement( name="individual_columns", desc="individual columns", visible=False, readonly=True )

    def set_dataset_metadata_from_comments( self, dataset ):
        Fake.set_dataset_metadata_from_comments( self, dataset )
        self.set_dataset_individual_metadata( dataset )

    def set_dataset_individual_metadata( self, dataset ):
        individual_list = dataset.metadata.comment_metadata.get( 'individuals', None )
        if not isinstance( individual_list, list ):
            individual_list = []

        individual_names = []
        individual_columns = []

        for individual in individual_list:
            if not isinstance( individual, list ) or len( individual ) != 2:
                continue
            name, col = individual
            if not isinstance( name, basestring ):
                name = ''
            try:
                c = int( col )
            except:
                c = 0
            if 0 < c <= dataset.metadata.columns:
                individual_names.append( name )
                individual_columns.append( c )

        if individual_names:
            dataset.metadata.individual_names = individual_names[:]
            dataset.metadata.individual_columns = individual_columns[:]

class GDGenotype( GDSnp ):
    """ Webb's genotype file format """
    file_ext = 'gd_genotype'

class GDSap( Fake ):
    """ Webb's SAP file format """
    file_ext = 'gd_sap'

    MetadataElement( name="kegg_gene", desc="KEGG gene code column", param=metadata.ColumnParameter, default=0 )
    MetadataElement( name="kegg_path", desc="KEGG pathway code/name column", param=metadata.ColumnParameter, default=0 )

