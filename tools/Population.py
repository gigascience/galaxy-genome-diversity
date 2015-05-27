#!/usr/bin/env python

import OrderedDict
import base64
import json
import zlib

import sys

class Individual(object):
    __slots__ = ['_column', '_name', '_alias']

    def __init__(self, column, name, alias=None):
        self._column = int(column)
        self._name = name
        self._alias = alias

    @property
    def column(self):
        return self._column

    @property
    def name(self):
        return self._name if self._alias is None else self._alias

    @property
    def alias(self):
        return self._alias

    @alias.setter
    def alias(self, alias):
        self._alias = alias

    @property
    def real_name(self):
        return self._name

    def __eq__(self, other):
        return self._column == other._column and self._name == other._name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return 'Individual: column={0} name={1} alias={2}'.format(self._column, self._name, self._alias)


class Population(object):
    def __init__(self, name=None):
        self._columns = OrderedDict.OrderedDict()
        self._name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    def add_individual(self, individual, alias=None):
        if individual.column not in self._columns:
            self._columns[individual.column] = individual
        elif self._columns[individual.column] == individual:
            # should should this be an error?
            # should we replace the alias using this entry?
            pass
        else:
            raise 'Duplicate column: {0}'.format(individual)

    def is_superset(self, other):
        for column, other_individual in other._columns.items():
            our_individual = self._columns.get(column)
            if our_individual is None or our_individual != other_individual:
                return False
        return True

    def is_disjoint(self, other):
        for column, our_individual in self._columns.items():
            other_individual = other._columns.get(column)
            if other_individual is not None and other_individual == our_individual:
                return False
        return True

    def column_list(self):
        return self._columns.keys()

    def individual_with_column(self, column):
        if column in self._columns:
            return self._columns[column]
        return None

    def tag_list(self, delimiter=':'):
        entries = []
        for column, individual in self._columns.iteritems():
            first_token = individual.name.split()[0]
            entry = '{0}{1}{2}'.format(column, delimiter, first_token)
            entries.append(entry)
        return entries

    def to_string(self, delimiter=':', separator=' ', replace_names_with=None):
        entries = []
        for column, individual in self._columns.items():
            value = individual.name
            if replace_names_with is not None:
                value = replace_names_with
            entry = '{0}{1}{2}'.format(column, delimiter, value)
            entries.append(entry)
        return separator.join(entries)

    def __str__(self):
        return self.to_string()

    def from_population_file(self, filename):
        with open(filename) as fh:
            for line in fh:
                line = line.rstrip('\r\n')
                column, name, alias = line.split('\t')
                alias = alias.strip()
                individual = Individual(column, name)
                if alias:
                    individual.alias = alias
                self.add_individual(individual)

    def from_tag_list(self, tag_list):
        for tag in tag_list:
            column, name = tag.split(':')
            individual = Individual(column, name)
            self.add_individual(individual)

    def from_wrapped_dict(self, wrapped_dict):
        unwraped_dict = self.unwrap_dict(wrapped_dict)
        for name, column in unwraped_dict.iteritems():
            individual = Individual(column, name)
            self.add_individual(individual)

    def unwrap_dict(self, wrapped_dict):
        decoded_value = self.decode_value(wrapped_dict)
        decompressed_value = self.decompress_value(decoded_value)
        def _decode_list(data):
            rv = []
            for item in data:
                if isinstance(item, unicode):
                    item = item.encode('utf-8')
                elif isinstance(item, list):
                    item = _decode_list(item)
                elif isinstance(item, dict):
                    item = _decode_dict(item)
                rv.append(item)
            return rv
        def _decode_dict(data):
            rv = {}
            for key, value in data.iteritems():
                if isinstance(key, unicode):
                    key = key.encode('utf-8')
                if isinstance(value, unicode):
                    value = value.encode('utf-8')
                elif isinstance(value, list):
                    value = _decode_list(value)
                elif isinstance(value, dict):
                    value = _decode_dict(value)
                rv[key] = value
            return rv
        unwrapped_dict = json.loads(decompressed_value, object_hook=_decode_dict)
        return unwrapped_dict

    def decode_value(self, value):
        try:
            return base64.b64decode(value)
        except TypeError, message:
            print >> sys.stderr, 'base64.b64decode: {0}: {1}'.format(message, value)
            sys.exit(1)

    def decompress_value(self, value):
        try:
            return zlib.decompress(value)
        except zlib.error, message:
            print >> sys.stderr, 'zlib.decompress: {0}'.format(message)
            sys.exit(1)

    def individual_names(self):
        for column, individual in self._columns.items():
            first_token = individual.name.split()[0]
            yield first_token

