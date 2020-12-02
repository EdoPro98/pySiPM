import platform, importlib
version = platform.python_version()
implementation = platform.python_implementation()
compiler = platform.python_compiler()
cur_system = platform.uname()

nversion = 'Not Installed'
cpversion = 'Not Installed'
rootversion = 'Not Installed'
if importlib.util.find_spec('numpy'):
  import numpy
  nversion = numpy.__version__
if importlib.util.find_spec('cupy'):
  import cupy
  cpversion = cupy.__version__
if importlib.util.find_spec('ROOT'):
    try:
        from ROOT import gROOT
        rootversion = gROOT.GetVersion()
    except:
        pass

assert importlib.util.find_spec('numpy'), 'Numpy is needed, please install using: pip install numpy'
assert importlib.util.find_spec('matplotlib'), 'Matplotlib is needed, please install using: pip install matplotlib'
assert importlib.util.find_spec('h5py'), 'H5py is needed, please install using: pip install h5py'
assert importlib.util.find_spec('uproot'), 'Uproot is needed, please install using: pip install uproot'

s1=f'''.__________________________________________________.
|                    ._____ _ .____  .__  __.      |
|     .____ _.  _.  /      (_)|  _ \ |  \/  |      |
|     |  _ | | | | |  (----'_ | |_) || \  / |      |
|     | |_)| |_| |  \  \   | ||  __/ | |\/| |      |
|     | .__/\__  .---)  |  | || |    | |  | |      |
|     | |    __/ \_____/   |_||_|    |_|  |_|      |
|     |_|   |___/                                  |'''
s2=f'''|
|Version:  0.1
|Running on: {cur_system.system} - {cur_system.release}
|Python version: {implementation}  {version}
|C compiler: {compiler}
|Cern ROOT version: {rootversion}
|Numpy version:  {nversion}
|Cupy version: {cpversion}'''
s3='''|
|The IDEA collaboration:
|Edoardo Proserpio
|Massimiliano Antonello
|Romualdo Santoro'''

print(s1)
for l in s2.splitlines():
  if l[0]=='|':
    text = '|'+l[1:].strip().center(50)+'|'
    print(text)
print('|'+'_'*50+'|')
for l in s3.splitlines():
  if l[0]=='|':
    text = '|'+l[1:].strip().center(50)+'|'
    print(text)
print('|'+'_'*50+'|')
