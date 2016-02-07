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

Important Notes for Installation

Right now, the OpenMM 6.2.0 shows several errors when one installs the BEDAM plugin.

1 After installing OpenMM, when one install BEDAM plugin, it shows the sfmt/sfmt.h is missing. One has to manually copy the sfmt folder to openmm/include/openmm/reference/

2 After installing OpenMM, when one install BEDAM plugin, it shows the lepton/lepton.h is missing. One has to manually copy the lepton folder to openmm/include/openmm/opencl/

3 When one runs 'make PythonInstall', it shows RPMD has errors. One has to manually disable several statements in the file BEDAMPluginWrapper.cpp in the folder build_openmm_bedam_plugin/python/ .

The required disabled statements are as follows
1 

//static void *_p_OpenMM__RPMDIntegratorTo_p_OpenMM__Integrator(void *x, int *SWIGUNUSEDPARM(newmemory)) {
  //return (void *)((OpenMM::Integrator *)  ((OpenMM::RPMDIntegrator *) x));
  //}


2

/*  {&_swigt__p_OpenMM__RPMDIntegrator, _p_OpenMM__RPMDIntegratorTo_p_OpenMM__Integrator, 0, 0},*/

Otherwise, one can not successfully install BEDAM plugin


## Test


```

cd example
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<openmm_dir>/lib
python test.py

Notes:

1 the dmsreader.py in the example/ is the required file, this file can read two dms files into the openmm system.
2 the example is calculating the binding energy for bcd and benzene with lambda equals to 1.0 and with HCT implicit solvent model. The dms files in the example/ have already been added the HCT parameters. 

```

`<openmm_dir>` is the OpenMM installation directory. Again, the last step is best accomplished under the same `virtualenv` environment used to build the python modules.

# openmm_bedam_plugin
