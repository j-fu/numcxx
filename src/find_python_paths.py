import sys
import os
from distutils import sysconfig
import numpy
import argparse
parser = argparse.ArgumentParser(description=
        """
        Python script to find libs + includes for python along with numpy.
        Works with python2.7 + python3
        """
        )
parser.add_argument('--libs',  action='store_true', default=False,    dest='libs',help='Print out python libs')
parser.add_argument('--includes',  action='store_true', default=False,    dest='includes',help='Print out python includes')
args=parser.parse_args()


python_version = sysconfig.get_config_var('VERSION')
numpy_version=numpy.__version__

for directory in sys.path:
    trydir=directory+"/numpy/"
    if os.path.exists(trydir):
        numpy_root=trydir
        break

if not numpy_root:
    raise Exception("missing numpy path")


python_includes= ['-I' + sysconfig.get_python_inc(),
                  '-I' + sysconfig.get_python_inc(plat_specific=True),
                  '-I' + numpy_root+'core/include']
#
# Need this for python3
#
abiflags=sysconfig.get_config_var('ABIFLAGS')
if not abiflags:
    abiflags=' '

python_libs=(['-lpython' + python_version+abiflags]
             +sysconfig.get_config_var('LIBS').split()
             +sysconfig.get_config_var('SYSLIBS').split())

ldpath=sysconfig.get_config_var('LIBPL')

if args.includes:
    print(  ' '.join(python_includes))
if args.libs:
    print( '-L'+ldpath+' '+' '.join(python_libs))

