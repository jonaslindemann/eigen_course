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

To validate that Eigen is working properly, create a file called ``ex0.cpp`` with the following content:

.. code-block:: cpp

    #include <iostream>
    #include <Eigen/Dense>

    int main()
    {
        Eigen::Matrix3d m = Eigen::Matrix3d::Random();
        std::cout << "Here is the matrix m:" << std::endl;
        std::cout << m << std::endl;
        return 0;
    }

Compile the file using the following command:

.. code-block:: bash

    g++ ex0.cpp -o ex0

Run the compiled file using the following command:

.. code-block:: bash

    ./ex0

This should give the following output:

.. code-block:: bash

    Here is the matrix m:
    0.680375   0.59688 -0.329554
    -0.211234  0.823295  0.536459
    0.566198 -0.604897 -0.444451


Linux
-----

The easiest way to set up a development environment on Linux is to use a package manager. The following commands can be used to install the necessary tools on Ubuntu:

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install g++
    sudo apt-get install cmake
    sudo apt-get install libeigen3-dev

To validate the environment you can run the following commands:

.. code-block:: bash

    g++ --version
    cmake --version

This shoudl give the following output (or similar):

.. code-block:: bash

    g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0
    Copyright (C) 2023 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    (base) lindemann@GRAYSON:~$ cmake --version
    cmake version 3.28.3

    CMake suite maintained and supported by Kitware (kitware.com/cmake).

To validate that Eigen is working properly, create a file called ``ex0.cpp`` with the following content:

.. code-block:: cpp

    #include <iostream>
    #include <Eigen/Dense>

    int main()
    {
        Eigen::Matrix3d m = Eigen::Matrix3d::Random();
        std::cout << "Here is the matrix m:" << std::endl;
        std::cout << m << std::endl;
        return 0;
    }

Compile the file using the following command:

.. code-block:: bash

    g++ ex0.cpp -I/usr/include/eigen3 -o ex0

Run the compiled file using the following command:

.. code-block:: bash

    ./ex0

This should give the following output:

.. code-block:: bash

    Here is the matrix m:
    0.680375   0.59688 -0.329554
    -0.211234  0.823295  0.536459
    0.566198 -0.604897 -0.444451


macOS
-----

The easiest way to set up a development environment on macOS is to use a package manager. The following commands can be used to install the necessary tools using Homebrew:

.. code-block:: bash

    brew install gcc
    brew install cmake
    brew install eigen

To validate that Eigen is working properly, create a file called ``ex0.cpp`` with the following content:

.. code-block:: cpp

    #include <iostream>
    #include <Eigen/Dense>

    int main()
    {
        Eigen::Matrix3d m = Eigen::Matrix3d::Random();
        std::cout << "Here is the matrix m:" << std::endl;
        std::cout << m << std::endl;
        return 0;
    }

Compile the file using the following command:

.. code-block:: bash

    g++ -std=c++11 -I/usr/local/include/eigen3 ex0.cpp -o ex0

Run the compiled file using the following command:

.. code-block:: bash

    ./ex0

This should give the following output:

.. code-block:: bash

    Here is the matrix m:
    0.680375   0.59688 -0.329554
    -0.211234  0.823295  0.536459
    0.566198 -0.604897 -0.444451


Windows
-------