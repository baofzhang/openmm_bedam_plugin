# OpenMM BEDAM Plugin

A plugin to add the BEDAM model to OpenMM

Baofeng Zhang <BZhang@brooklyn.cuny.edu>

Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>

## License

This software is released under the LGPL license. See LICENSE.

## Credits

This software is written and maintained by Baofeng Zhang <BZhang@brooklyn.cuny.edu> and Emilio Gallicchio <egallicchio@brooklyn.cuny.edu> with support from a grant from the National Science Foundation (ACI 1440665).

The plugin interface is based on the [openmmexampleplugin](https://github.com/peastman/openmmexampleplugin) by Peter Eastman.

## Installation

Locate the OpenMM installation directory, otherwise it will default to `/usr/local/openmm`.

Download the package from github:

```
git clone https://github.com/baofzhang/openmm_bedam_plugin.git
```


Build and install the plugin with cmake. Assuming a unix system:

```
mkdir build_openmm_bedam_plugin
cd build_openmm_bedam_plugin
ccmake -i ../openmm_bedam_plugin
```

Hit `c` (configure) until all variables are correctly set, then `g` to generate the makefiles. `OPENMM_DIR` should point to an existing OpenMM installation. `CMAKE_INSTALL_PREFIX` normally is the same as `OPENMM_DIR`. The BEDAM plugin requires the python API. You need `python` and `swig` to install it.

Once the configuration is done do:

```
make
make install
make PythonInstall
```

The last two steps may need superuser access depending on the installation target. It is recommended to to build the plugin under a `virtualenv` environment to install the python modules without superuser access.

## Test


```

cd example
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<openmm_dir>/lib:<openmm_dir>/lib/plugins
python test.py

Notes:

1 the dmsreader.py in the example/ is the required file, this file can read two dms files into the openmm system.
2 the example calculates the binding energy of a beta-cyclodextrin/benzene complex in implicit solvent. The dms files in the example/ have already been added the HCT parameters. 

```

`<openmm_dir>` is the OpenMM installation directory. Again, the last step is best accomplished under the same `virtualenv` environment used to build the python modules.

## async-bedam-test-example

In this example, the GPU-accelerated version of BEDAM absolute binding free energy calculations has been implemented in the OpemMM MD engine and the Asynchronous Replica Exchange job coordination system. The absolute binding free energies of one inhibitor of MCL1 has been computed. Total 18 intermediate steps at lambda=0.0, 0.002, 0.004, 0.008, 0.01, 0.02, 0.04, 0.07, 0.1, 0.17, 0.25, 0.35, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0, were used to conduct the lambda-biased replica exchanges.

cd async-bedam-example
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<openmm_dir>/lib
python <async_re-openmm_dir>/bedamtempt_async_re.py mcl1-23.cntl

Notes:

1 <async_re-openmm_dir> is the async_re-openmm installation directory. Refer to [the async_re-openmm_dir github page](https://github.com/baofzhang/async_re-openmm) for download and installation instructions.
2 the dmsreaderwscwrestraint.py is a new version of dmsreader.py. In this script, softcore, c-alpha restraint functions are added. In order to apply c-alpha restraint, one has to prepare the receptor dms file in advance so that there are x0, y0 ,z0 values, which are exactly the copies of original x,y,z values.
3 One can control the soft-core value and other parameters in mcl1-23.py. 
4 in this example, we use four GPUs to run the jobs, each GPU running one replica. The nodefile sets up the GPU arrangement.

# openmm_bedam_plugin
