# Package Management with Spack

  [Spack](https://spack.readthedocs.io/en/latest/) is a python-based package manager
   that provides access to both precompiled binaries and locally-compilable source code.
   Spack is non-destructive: installing a new version does not break existing installations,
   so many configurations can coexist on the same system.

## Getting Started

   Spack can be installed by directly cloning from the git repository:  

   ```bash
   git clone https://github.com/spack/spack && cd spack
   git checkout releases/v0.15
   ```

   You may also need some additional python dependencies, for example:

   ```bash
   pip3 install botocore boto3  
   ```

   Next add Spack to your path by sourcing a shell script in your ~/.bashrc:  

   ```bash
   source <PATH-TO-SPACK>/share/spack/setup-env.sh
   ```

   This provides additional command-line features. Please refer to the
   [website](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html)
   for the most current information.

## Installing Dependencies with the BASH Scripts

   Two BASH scripts have been written to aid installation. The `linux-GCC` script
   installs GCC 8 and all external libraries, and should be fine to execute on any
   linux machine.  

   The `sol-intel` script has been written with exciting developer machines in mind,
   where the Intel compiler suite is available as a module. If you wish to adapt this
   script, remove `module load intel/2019` and ensure you have an Intel compiler
   installed. Intel is available via spack or via
   [apt-get](https://software.intel.com/content/www/us/en/develop/articles/installing-intel-oneapi-toolkits-via-apt.html)
   on linux systems.  

## Managing Build Stacks with Virtual Environments

|  Description                        |   Commands                     |
|-------------------------------------|--------------------------------|
|  List virtual envs                  |  `spack env list`              |
|  Show current virtual env (if any)  |  `spack env status`            |
|  Create virtual env 'project'       |  `spack env create project`    |
|  Activate virtual env 'project'     |  `spack env activate project`  |
|  Deactivate current virtual env     |  `spack env deactivate`        |
