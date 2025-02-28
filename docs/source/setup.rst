Setup
=====

This document describes how to set up a development environment for the course. 

Using COSMOS at LUNARC
----------------------

If you have applied for an account at LUNARC you can login to our remote desktop environment or connect using a SSH terminal. More information on this can be found here:

`Login into remote desktop environment <https://lunarc-documentation.readthedocs.io/en/latest/getting_started/using_hpc_desktop/>`_

`Login using SSH <https://lunarc-documentation.readthedocs.io/en/latest/getting_started/login_howto/>`_

A suitable environment can be loaded using the following commands:

.. code-block:: bash

    module load foss/2024a
    module load Eigen
    module load CMake

To validate the environment you can run the following commands:

.. code-block:: bash

    g++ --version
    cmake --version

This shoudl give the following output:

.. code-block:: bash

    $ g++ --version
    g++ (GCC) 13.3.0
    Copyright (C) 2023 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    $ cmake --version
    cmake version 3.29.3

CMake suite maintained and supported by Kitware (kitware.com/cmake).
or similar.

Linux
-----

macOS
-----

Windows
-------