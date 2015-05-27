#!/usr/bin/env python

import base64
import errno
import os
import subprocess
import sys
import zlib

################################################################################

def die(message):
    print >> sys.stderr, message
    sys.exit(1)

################################################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno <> errno.EEXIST:
            raise

################################################################################

def run_program(prog, args, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    kwargs = {
        "bufsize": -1,
        "close_fds": False,
        "creationflags": 0,
        "cwd": None,
        "env": None,
        "executable": prog,
        "preexec_fn": None,
        "shell": False,
        "startupinfo": None,
        "stderr": stderr,
        "stdin": None,
        "stdout": stdout,
        "universal_newlines": False
    }

    str_args = [str(x) for x in args]

    p = subprocess.Popen(str_args, **kwargs)
    (stdoutdata, stderrdata) = p.communicate()
    rc = p.returncode

    if rc != 0:
        die('FAILED:\n{0}\nCOMMAND:\n{1}'.format(stderrdata, ' '.join(str_args)))

    return stdoutdata, stderrdata

################################################################################

def unwrap_string(string):
    try:
        decoded_string = base64.b64decode(string)
    except TypeError, message:
        die('base64.b64decode: {0}: {1}'.format(message, string))

    try:
        return zlib.decompress(decoded_string)
    except zlib.error, message:
        die('zlib.decompress: {0}'.format(message))

################################################################################

def wrap_string(string, level=9):
    try:
        compressed_string = zlib.compress(string, level)
    except zlib.error, message:
        die('zlib.compress: {0}'.format(message))
    return base64.b64encode(compressed_string)
    
################################################################################
