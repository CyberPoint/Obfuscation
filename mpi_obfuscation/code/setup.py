#! /usr/bin/sage


from setuptools import setup, Extension


library_dirs = ["./","/usr/lib", '/usr/lib/openmpi/lib','/usr//lib','/usr/local/lib/python2.7/dist-packages/mpi4py/lib-pmpi/']

libraries = [
    'gmp',
    'gomp',
    'mpi_cxx',
    'dl',
    'hwloc'
]
compile_args = [
    '-fopenmp',
    '-O3',
    '-Wall',
    '-pthread'
]

incdirs = [
   # '/usr/src/linux-headers-3.16.0-41/include/linux/',
    '/usr/lib/openmpi/include',
    '/usr/lib/openmpi/include/openmpi',
    '/usr/local/lib/python2.7/dist-packages/mpi4py/include/'
]

zobfuscator = Extension(
    'obf._zobfuscator',
    library_dirs = library_dirs,
    include_dirs = incdirs,
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/circuit.cpp',
        'src/clt_mlm.cpp',
        'src/_zobfuscator.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
        'src/utils.cpp'
    ]
)

obfuscator = Extension(
    'obf._obfuscator',
   library_dirs = library_dirs,
    include_dirs=incdirs,
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/clt_mlm.cpp',
        'src/_obfuscator.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
        'src/utils.cpp',
    ]
)

setup(name='obfuscator',
      author='Alex J. Malozemoff',
      author_email='amaloz@cs.umd.edu',
      version='1.0alpha',
      description='Implementation of cryptographic program obfuscation',
      license='GPLv2',
      url='https://amaloz.github.io/obfuscation',
      packages=['obf'],
      ext_modules=[obfuscator, zobfuscator],
      scripts=['obfuscator'],
      test_suite='t',
      classifiers=[
          'Topic :: Security :: Cryptography',
          'Environment :: Console',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Operating System :: Unix',
          'License :: OSI Approved :: Free For Educational Use',
          'Programming Language :: C',
          'Programming Language :: Sage',
      ])
