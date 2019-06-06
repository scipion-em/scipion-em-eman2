============
EMAN2 plugin
============

This plugin provide wrappers around several programs of `EMAN2 <https://blake.bcm.edu/emanwiki/EMAN2>`_ software suite.

.. figure:: http://scipion-test.cnb.csic.es:9980/badges/eman2_devel.svg
   :align: left
   :alt: build status

Installation
------------

You will need to use `2.0 <https://github.com/I2PC/scipion/releases/tag/V2.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-eman

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-eman.git

    * install

    .. code-block::

        scipion installp -p path_to_scipion-em-eman --devel

EMAN2 binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is ``software/em/eman-2.3``, if you want to change it, set *EMAN2DIR* in ``scipion.conf`` file to the folder where the EMAN2 is installed. To check the installation, simply run one of the following Scipion tests:

.. code-block::

   scipion test eman2.tests.test_protocols_eman.TestEmanTiltValidate
   scipion test eman2.tests.test_protocols_eman.TestEmanRefineEasy
   scipion test eman2.tests.test_protocols_eman.TestEmanRefine2DBispec
   scipion test eman2.tests.test_protocols_eman.TestEmanRefine2D
   scipion test eman2.tests.test_protocols_eman.TestEmanReconstruct
   scipion test eman2.tests.test_protocols_eman.TestEmanInitialModelMda
   scipion test eman2.tests.test_protocols_eman.TestEmanInitialModelGroel
   scipion test eman2.tests.test_protocols_eman.TestEmanCtfAuto
   scipion test eman2.tests.test_protocols_eman.TestEmanAutopick

A complete list of tests can also be seen by executing ``scipion test --show --grep eman``

Supported versions
------------------

2.12, 2.21, 2.3

In 2018 the plugin was updated to support the latest (at that moment) EMAN2: 2.21. This required a lot of code refactoring and the support of old EMAN2 version 2.11 had to be discontinued. Several new protocols were added: 2D refinements, tilt validation, ctf-auto and new e2boxer. The full changelog since Scipion-1.x is available `here <https://github.com/scipion-em/scipion-em-eman2/issues/1`_. In the following releases we plan to deprecate 2.12 support as well.

Protocols
---------

* [boxer (old and new e2boxer.py)](EmanProtBoxing)
* boxer auto (fully automated new boxer in 2.21)
* ctf auto
* [initial model](EmanProtInitModel)
* reconstruct
* refine 2d
* refine 2d bispec
* refine easy
* [sparx picker](SparxGaussianProtPicking)
* tilt validate

References
----------

1. \G. Tang, L. Peng, P.R. Baldwin, D.S. Mann, W. Jiang, I. Rees & S.J. Ludtke. (2007) EMAN2: an extensible image processing suite for electron microscopy. J Struct Biol. 157, 38-46. PMID: 16859925