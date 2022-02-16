============
EMAN2 plugin
============

This plugin provides wrappers around several programs of `EMAN2 <https://blake.bcm.edu/emanwiki/EMAN2>`_ software suite.

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

        git clone https://github.com/scipion-em/scipion-em-eman2.git

    * install

    .. code-block::

        scipion installp -p /path/to/scipion-em-eman2 --devel

**Important: starting from the plugin v3.2 all tomo protocols have been removed, they are now in https://github.com/scipion-em/scipion-em-emantomo**

EMAN2 binaries will be installed automatically with the plugin, but you can also link an existing installation.

    * Default installation path assumed is ``software/em/eman-2.91``, if you want to change it, set *EMAN2_HOME* in ``scipion.conf`` file pointing to the folder where the EMAN2 is installed.
    * The default scratch directory is assumed */tmp/*. You can change it by setting *EMAN2SCRATCHDIR* in ``scipion.conf`` or your shell environment.

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

2.31, 2.9, 2.91

Protocols
---------

* boxer (new interactive e2boxer.py)
* boxer auto (fully automated new boxer in EMAN >= 2.21)
* ctf auto
* `initial model <https://github.com/scipion-em/scipion-em-eman2/wiki/EmanProtInitModel>`_
* initial model SGD
* reconstruct
* refine 2d
* refine 2d bispec
* refine easy
* tilt validate

References
----------

1. \G. Tang, L. Peng, P.R. Baldwin, D.S. Mann, W. Jiang, I. Rees & S.J. Ludtke. (2007) EMAN2: an extensible image processing suite for electron microscopy. J Struct Biol. 157, 38-46. PMID: 16859925
