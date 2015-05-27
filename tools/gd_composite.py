#!/usr/bin/env python

from galaxy import eggs
import pkg_resources
pkg_resources.require( "Cheetah" )
from Cheetah.Template import Template

import errno
import os
from datetime import datetime

################################################################################

def die(message):
    print >> sys.stderr, message
    sys.exit(1)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno <> errno.EEXIST:
            raise

################################################################################

class Display(object):
    def display(self, parameter):
        print parameter

class DisplayFile(Display):
    def display(self, parameter):
        return '<a href="{0}">{1}</a>'.format(parameter.value, parameter.name)

class DisplayValue(Display):
    def display(self, parameter):
        if parameter.value is not None:
            return '{0}: {1}'.format(parameter.description, parameter.value)
        else:
            return '{0}'.format(parameter.description)

class DisplayTagList(Display):
    def display(self, parameter):
        rv = []
        if parameter.name:
            rv.append(parameter.name)
        rv.append('<ol>')
        for tag in parameter.value:
            col, individual_name = tag.split(':')
            rv.append('<li>{0}</li>'.format(individual_name))
        rv.append('</ol>')
        return '\n'.join(rv)

class DisplayPopulationList(Display):
    def display(self, parameter):
        rv = []
        rv.append('Populations')
        rv.append('<ul>')
        for population in parameter.value:
            rv.append('<li>')
            if population.name is not None:
                rv.append(population.name)
            rv.append('<ol>')
            for name in population.individual_names():
                rv.append('<li>{0}</li>'.format(name))
            rv.append('</ol>')
            rv.append('</li>')
        rv.append('</ul>')
        return '\n'.join(rv)

#    def display(self, parameter, name=''):
#        print '<ul> {0}'.format(name)
#        for individual_name in parameter.individual_names():
#            print '<li>{0}>/li>'.format(individual_name)
#        print '</ul>'
        
        
class Parameter(object):
    def __init__(self, name=None, value=None, description=None, display_type=None):
        self.name = name
        self.value = value
        self.description = description
        if display_type is None:
            self.display_type = Display()
        else:
            self.display_type = display_type

    def display(self):
        return self.display_type.display(self)

class InfoPage(object):
    _realpath = os.path.realpath(__file__)
    _script_dir = os.path.dirname(_realpath)
    template_file = os.path.join(_script_dir, 'gd_composite_template.html')
    def __init__(self):
        self.timestamp = datetime.now().strftime('%Y-%m-%d %I:%M:%S %p')
        self.title = 'Genome Diversity Composite Dataset'
        self.inputs = []
        self.outputs = []
        self.misc = ''
        self.template = self.load_template()

    def load_template(self):
        with open(self.template_file) as f:
            return f.read().rstrip('\r\n')

    def set_title(self, title):
        self.title = title

    def add_input_parameter(self, parameter):
        self.inputs.append(parameter)

    def add_output_parameter(self, parameter):
        self.outputs.append(parameter)

    def add_misc(self, misc):
        self.misc = misc

    def render(self):
        return Template(self.template, searchList=[{'tool': self}])
        

            





