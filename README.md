# galaxy-genome-diversity

This repository houses the Galaxy wrappers for the Genome Diversity set of
tools developed by [Webb Miller's group](http://www.bx.psu.edu/miller_lab) at
the Department of Biology at Penn State.

Each Galaxy tool is dependent on one or more associated Tool Shed packages which
provides access to one or more specific binary executables.

## Folder structure

There is one directory for each tool which has been wrapped for use in Galaxy 
and these are contained within the `tools` folder.  The `dependencies` folder 
contains the packages for the Galaxy tool dependency definitions.

The `test-data` directory contains data files which are used for functional
tests by the Galaxy wrappers.

## Installation

These Galaxy tools are to be used from within a Galaxy server. They can be
automatically installed via the
[GigaScience Tool Shed](http://gigatoolshed.net) which handles their required
dependencies. All tools can be manually installed too. Documentation for
automatically and manually installing each Galaxy SOAP tool can be found within
the tools' README file.

## Testing

Functional tests for most of the tools can be found in the tool's XML
configuration file.