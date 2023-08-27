===========
EMAN plugin
===========

This plugin provides wrappers around several programs of `EMAN <https://blake.bcm.edu/emanwiki/EMAN2>`_ software suite.

.. image:: https://img.shields.io/pypi/v/scipion-em-eman2.svg
        :target: https://pypi.python.org/pypi/scipion-em-eman2
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-eman2.svg
        :target: https://pypi.python.org/pypi/scipion-em-eman2
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-eman2.svg
        :target: https://pypi.python.org/pypi/scipion-em-eman2
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-eman2?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-eman2
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-eman2
        :target: https://pypi.python.org/pypi/scipion-em-eman2
        :alt: Downloads


Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-eman2

b) Developer's version

    * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-eman2.git

    * install

    .. code-block::

        scipion installp -p /path/to/scipion-em-eman2 --devel

EMAN software will be installed automatically with the plugin but you can also use an existing installation by providing *EMAN_ENV_ACTIVATION* (see below).

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*EMAN_ENV_ACTIVATION* (default = conda activate eman-2.99.52):
Command to activate the EMAN environment.

The default scratch directory is assumed */tmp/*. You can change it by setting *EMAN2SCRATCHDIR* in ``scipion.conf`` or your shell environment.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

   scipion test eman2.tests.test_protocols_eman.TestEmanTiltValidate
   scipion test eman2.tests.test_protocols_eman.TestEmanRefineEasy
   scipion test eman2.tests.test_protocols_eman.TestEmanRefine2DBispec
   scipion test eman2.tests.test_protocols_eman.TestEmanRefine2D
   scipion test eman2.tests.test_protocols_eman.TestEmanReconstruct
   scipion test eman2.tests.test_protocols_eman.TestEmanInitialModelMda
   scipion test eman2.tests.test_protocols_eman.TestEmanInitialModelGroel
   scipion test eman2.tests.test_protocols_eman.TestEmanInitialModelSGD
   scipion test eman2.tests.test_protocols_eman.TestEmanCtfAuto
   scipion test eman2.tests.test_protocols_eman.TestEmanAutopick

A complete list of tests can also be seen by executing ``scipion test --show --grep eman``

Supported versions
------------------

2.99.47, 2.99.52

Protocols
---------

* boxer
* boxer auto
* ctf auto
* initial model
* initial model SGD
* reconstruct
* refine 2d
* refine 2d bispec
* refine easy
* tilt validate

References
----------

1. \G. Tang, L. Peng, P.R. Baldwin, D.S. Mann, W. Jiang, I. Rees & S.J. Ludtke. (2007) EMAN2: an extensible image processing suite for electron microscopy. J Struct Biol. 157, 38-46. PMID: 16859925
