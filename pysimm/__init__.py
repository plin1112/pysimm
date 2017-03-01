# ******************************************************************************
# pysimm.__init__ module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from __future__ import print_function
import urllib3

__version__ = '0.1.0'

error = True
warning = True
verbose = True
debug = True

error_print = lambda *a, **k: print('(error) PySIMM:', *a) if error else lambda *a, **k: None
warning_print = lambda *a, **k: print('(warning) PySIMM:', *a) if warning else lambda *a, **k: None
verbose_print = lambda *a, **k: print('PySIMM:', *a) if verbose else lambda *a, **k: None
debug_print = lambda *a, **k: print('(debug) PySIMM:', *a) if debug else lambda *a, **k: None

class PysimmError(Exception):
    pass


def check_version():
    try:
        http = urllib3.PoolManager()
        r = http.request('GET', 'http://pysimm.org/version')
        if r.status == 200:
            remote_version = r.data.split('.')
            local_version = __version__.split('.')
            if remote_version > local_version:
                print('remote version {} is newer than local version {}'.format('.'.join(remote_version), '.'.join(local_version)))
                print('consider checking release notes in case of bug patches or updates http://pysimm.org/release-notes')
        else:
            print('tried to check remote version on pysimm.org server but failed')
    except:
        print('tried to check remote version on pysimm.org server but failed')
        

check_version()
