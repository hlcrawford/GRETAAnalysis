import sys
import os

colors = {}
colors['cyan']   = '\033[36m'
colors['purple'] = '\033[95m'
colors['blue']   = '\033[94m'
colors['green']  = '\033[92m'
colors['yellow'] = '\033[93m'
colors['red']    = '\033[91m'
colors['end']    = '\033[0m'

## If the output is not a terminal, remove the colors
if not sys.stdout.isatty():
   for key, value in colors.items():
      colors[key] = ''

compile_source_message = '%sCompiling %s===================> %s$SOURCE%s' % \
   (colors['blue'], colors['purple'], colors['cyan'], colors['end'])

compile_shared_source_message = '%sCompiling shared %s============> %s$SOURCE%s' % \
   (colors['blue'], colors['purple'], colors['cyan'], colors['end'])

link_program_message = '%sLinking Program %s=============> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

link_library_message = '%sLinking Static Library %s=====> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

ranlib_library_message = '%sRanlib Library %s===========> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

link_shared_library_message = '%sLinking Shared Library %s======> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

java_library_message = '%sCreating Java Archive %s======> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

root_dictionary_message = '%sGenerating ROOT dictionary %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])
 
def rootcint(target, source, env):
    """Executes the ROOT dictionary generator over a list of headers. """
    dictname = target[0]
    headers = ""
    for f in source:
    	headers += str(f) + " "

    command = "rootcint -f %s -c -p %s" %(dictname, headers)
    ok = os.system(command)
    return ok

def remember_pcm(target, source, env):
    new_target = os.path.splitext(str(target[0]))[0]+'_rdict.pcm'
    target.append(new_target)
    return target, source

## Create construction environment propagating the external environment
env = Environment(ENV=os.environ, 
      		  CXXCOMSTR = compile_source_message,
  		  CCCOMSTR = compile_source_message,
  		  SHCCCOMSTR = compile_shared_source_message,
  		  SHCXXCOMSTR = compile_shared_source_message,
  		  ARCOMSTR = link_library_message,
  		  RANLIBCOMSTR = ranlib_library_message,
  		  SHLINKCOMSTR = link_shared_library_message,
  		  LINKCOMSTR = link_program_message,
  		  JARCOMSTR = java_library_message,
  		  JAVACCOMSTR = compile_source_message) 

## Create a rootcint builder and attach it to the environment
bld = Builder(action=Action(rootcint,root_dictionary_message), emitter = remember_pcm)
env.Append(BUILDERS = {'RootCint':bld})

## Optimization flags ##################################################
env.Append(CCFLAGS = ['-O2', '-D_FILE_OFFSET_64', '-g'], LINKFLAGS=[])

## Finding dependencies (ROOT)
try:
    env.ParseConfig('root-config --cflags')
    env.ParseConfig('root-config --glibs')
except OSError:
    print ("scons: ROOT not found!")
    exit(1)

env.Append(CPPPATH=['.', './src', './include', './tracking', './src/hfc'])
env.Append(LIBPATH='./lib')

TrackDebug = ARGUMENTS.get('TRACK_DEBUG', '0')
if TrackDebug == '1':
    print('%sBuilding with GRETA_TRACK_DEBUG (tracking trace to stderr)%s' %
          (colors['yellow'], colors['end']))
    env.Append(CXXFLAGS=['-DGRETA_TRACK_DEBUG'])

envUnpack = env.Clone()

## Building GRETINADict and libGRETINA #################################
gretaDictTarget = 'lib/GRETADict.cpp'
gretaDictHeaders = ['include/GRETA.h', 'include/GretaTracked.h', 'include/dptc.h', 'include/SortingStructures.h', 'include/LinkDefGRETA.h']
env.RootCint(gretaDictTarget, gretaDictHeaders)

gretaLibTarget = 'lib/GRETA'
track_cppsources = [
    'tracking/greta_tracking_globals.cpp',
    'tracking/setupTrack.cpp',
    'tracking/track.cpp',
    'tracking/ctkTrackOpt.cpp',
    'tracking/ctkPrTrkPar.cpp',
    'tracking/ctkinit.cpp',
    'tracking/findVector.cpp',
    'tracking/findAngle.cpp',
    'tracking/findCAngle.cpp',
    'tracking/ctksort.cpp',
    'tracking/splitCluster.cpp',
    'tracking/combineCluster.cpp',
    'tracking/matchMaker.cpp',
    'tracking/reCluster.cpp',
    'tracking/str_decomp.cpp',
    'tracking/rotations.cpp',
    'tracking/pairprod.cpp',
    'tracking/ctktk.cpp',
    'tracking/printEvent.cpp',
]

gretaLibSources = ['lib/GRETADict.cpp', 'src/GRETA.cpp', 'src/SortingStructures.cpp',
                   'src/GretaTracked.cpp', 'src/GretaTracking.cpp', 'src/GretaTrackingWorkspace.cpp',
                   'src/GretaTrackingDebug.cpp'] + track_cppsources
envLib = env.Clone()
if TrackDebug == '1':
    envLib.Append(CXXFLAGS=['-DGRETA_TRACK_DEBUG'])
# Legacy trackMain sources as C++; common headers for all tracking TUs
envLib.Append(CXXFLAGS=['-include', 'tracking/greta_tracking_std.hpp',
                        '-DGRETA_TRACKING_USE_WORKSPACE'])
if sys.platform == 'darwin':
    envLib.Append(SHLINKFLAGS=['-install_name', '@loader_path/lib/libGRETA.dylib'])
elif sys.platform.startswith('linux'):
    envLib.Append(SHLINKFLAGS=['-install_name', '$ORIGIN/lib/libGRETA.so'])

envLib.SharedLibrary(target = gretaLibTarget, source = gretaLibSources, 
                  SHLIBPREFIX='lib')

if sys.platform == 'darwin':
    envUnpack.Append(LINKFLAGS=['-Wl,-rpath,@loader_path/lib'])
elif sys.platform.startswith('linux'):
    envUnpack.Append(LINKFLAGS=['-Wl,-rpath,$ORIGIN/lib'])

## Building Unpack executable ###########################################
unpackTarget = 'readGreta'
unpackSources = ['src/readGretaMultiCrystal.cpp', 'src/Globals.cpp']
if TrackDebug == '1':
    envUnpack.Append(CXXFLAGS=['-DGRETA_TRACK_DEBUG'])
envUnpack.Append(LIBS=['GRETA'], LIBPATH=['lib'])
envUnpack.Program(target = unpackTarget, source = unpackSources)

## Pass-2 tracking from mode-2 ROOT tree (OpenMP via -j) #################
def add_openmp(env):
    """Enable -fopenmp for trackGreta (-j). macOS needs libomp (brew install libomp)."""
    if sys.platform == 'darwin':
        env.Append(CXXFLAGS=['-Xpreprocessor', '-fopenmp'])
        env.Append(LIBS=['omp'])
        for prefix in ['/opt/homebrew/opt/libomp', '/usr/local/opt/libomp']:
            if os.path.isdir(prefix):
                env.Append(CPPPATH=[prefix + '/include'])
                env.Append(LIBPATH=[prefix + '/lib'])
                break
    else:
        env.Append(CXXFLAGS=['-fopenmp'])
        env.Append(LINKFLAGS=['-fopenmp'])

trackGretaTarget = 'trackGreta'
trackGretaSources = ['src/trackGretaFromRoot.cpp']
envTrackGreta = envUnpack.Clone()
if TrackDebug == '1':
    envTrackGreta.Append(CXXFLAGS=['-DGRETA_TRACK_DEBUG'])
add_openmp(envTrackGreta)
envTrackGreta.Program(target = trackGretaTarget, source = trackGretaSources)

## Building GRETA HFC executable ########################################
gebTarget = 'GEB_HFC'
gebSources = ['src/hfc/GEB_HFC.cpp', 'src/hfc/HFC.cpp']
env.Program(target = gebTarget, source = gebSources)