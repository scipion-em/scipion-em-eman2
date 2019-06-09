# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import re
from os.path import exists, basename, join
from glob import glob

import pyworkflow.em as em
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam, BooleanParam,
                                        LabelParam)
from pyworkflow.utils import makePath, createLink, cleanPath


import eman2
from eman2.constants import *


class EmanProtRefine2DBispec(em.ProtClassify2D):
    """
    This protocol wraps *e2refine2d_bispec.py* EMAN2 program.

    This program is used to produce reference-free class averages
    from a population of mixed, unaligned particle images. These averages
    can be used to generate initial models or assess the structural
    variability of the data. They are not normally themselves used as part
    of the single particle reconstruction refinement process, which
    uses the raw particles in a reference-based classification
    approach. However, with a good structure, projections of the
    final 3-D model should be consistent with the results of
    this reference-free analysis.

    This variant of the program uses rotational/translational
    invariants derived from the bispectrum of each particle.
"""
    _label = 'refine 2D bispec'

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {
            'partSetFlipFullRes': self._getExtraPath('sets/all__ctf_flip_fullres.lst'),
            'partSetFlipLp5': self._getExtraPath('sets/all__ctf_flip_lp5.lst'),
            'partSetFlipLp7': self._getExtraPath('sets/all__ctf_flip_lp7.lst'),
            'partSetFlipLp12': self._getExtraPath('sets/all__ctf_flip_lp12.lst'),
            'partSetFlipLp20': self._getExtraPath('sets/all__ctf_flip_lp20.lst'),
            'partBispecSet': self._getExtraPath('sets/all__ctf_flip_bispec.lst'),
            'partInvarSet': self._getExtraPath('sets/all__ctf_flip_invar.lst'),
            'classes_scipion': self._getExtraPath('classes_scipion_it%(iter)02d.sqlite'),
            'classes': 'r2db_%(run)02d/classes_%(iter)02d.hdf',
            'cls': 'r2db_%(run)02d/classmx_%(iter)02d.hdf',
            'results': self._getExtraPath('results_it%(iter)02d.txt'),
            'basis': self._getExtraPath('r2db_%(run)02d/basis_%(iter)02d.hdf')
        }
        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self, currRun):
        """ Setup the regex on how to find iterations. """
        clsFn = self._getExtraPath(self._getFileName('classes', run=currRun, iter=1))
        self._iterTemplate = clsFn.replace('classes_01', 'classes_??')
        # Iterations will be identify by classes_XX_ where XX is the iteration
        #  number and is restricted to only 2 digits.
        self._iterRegex = re.compile('classes_(\d{2})')

    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputBispec', PointerParam,
                      label='Choose e2ctf auto protocol',
                      pointerClass='EmanProtCTFAuto',
                      help='Select EMAN CTF auto protocol that has '
                           'generated bispectra.')
        form.addParam('numberOfClassAvg', IntParam, default=32,
                      label='Number of class-averages',
                      help='Number of class-averages to generate. Normally you '
                           'would want a minimum of ~10-20 particles per class on '
                           'average, but it is fine to have 100-200 for a large data '
                           'set. If you plan on making a large number (>100) of '
                           'classes, you should use the *Fast seed* option. Note '
                           'that these averages are not used for final 3-D '
                           'refinement, so generating a very large number is not '
                           'useful in most situations.')
        form.addParam('numberOfIterations', IntParam, default=3,
                      label='Number of iterations',
                      help='Number of iterations of the overall 2-D refinement '
                           'process to run. For high contrast data, 4-5 iterations '
                           'may be more than enough, but for low contrast data '
                           'it could take 10-12 iterations to converge well.')
        form.addParam('nbasisfp', IntParam, default=8,
                      label='Number of MSA vectors to use',
                      help='Number of MSa basis vectors to use when '
                           'classifying particles.')
        if self._isVersion23():
            form.addParam('alignSort', BooleanParam, default=True,
                          label='Align and sort?',
                          help='This will align and sort the final class-averages '
                               'based on mutual similarity.')

        line = form.addLine('Centering: ',
                            help="If the default centering algorithm "
                                 "(xform.center) doesn't work well, "
                                 "you can specify one of the others "
                                 "here (e2help.py processor center)")
        line.addParam('centerType', EnumParam,
                      choices=['nocenter', 'xform.center',
                               'xform.centeracf', 'xform.centerofmass', 'None'],
                      label="", default=XFORM_CENTER,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('centerParams', StringParam, default='',
                      label='params')

        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "e2refine2d_bispec.py -h" will print a list '
                           'of possible options.')

        form.addSection(label='Class averaging')
        form.addParam('paramsMsg2', LabelParam, default=True,
                      label='These parameters are for advanced users only!\n',
                      help='For help please address to EMAN2 %s or run:\n'
                           '*scipion run e2help.py cmp -v 2* or\n'
                           '*scipion run e2help.py averagers -v 2*' % WIKI_URL)
        form.addParam('classIter', IntParam, default=4,
                      label='Number of iterations for class-averages',
                      help='Number of iterations to use when making '
                           'class-averages (default=5)')
        form.addParam('classKeep', FloatParam, default=0.8,
                      label='Fraction of particles to keep',
                      help='The fraction of particles to keep in each class, '
                      'based on the similarity score generated by cmp '
                      '(default=0.8)')
        form.addParam('classKeepSig', BooleanParam, default=False,
                      label='Keep particles based on sigma?',
                      help='Change the *keep* criterion from fraction-based '
                      'to sigma-based')
        form.addParam('classAveragerType', EnumParam,
                      choices=['absmaxmin', 'ctf.auto', 'ctf.weight',
                               'ctf.weight.autofilt', 'ctfw.auto', 'iteration',
                               'localweight', 'mean', 'mean.tomo',
                               'minmax', 'sigma', 'weightedfourier'],
                      label='Class averager: ',
                      default=AVG_CTF_WEIGHT_AUTOFILT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='The averager used to generated class-averages')

        line = form.addLine('classnormproc: ',
                            help='Normalization applied during class-averaging')
        line.addParam('classnormprocType', EnumParam,
                      choices=['normalize', 'normalize.bymass',
                               'normalize.circlemean', 'normalize.edgemean',
                               'normalize.local', 'normalize.lredge',
                               'normalize.mask', 'normalize.maxmin',
                               'normalize.ramp.normvar', 'normalize.rows',
                               'normalize.toimage', 'normalize.unitlen',
                               'normalize.unitsum', 'None'],
                      label='',
                      default=PROC_NORMALIZE_EDGEMEAN,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classnormprocParams', StringParam,
                      default='', label='params')

        line = form.addLine('classcmp: ')
        line.addParam('classcmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classcmpParams', StringParam,
                      default='', label='params',
                      help='The name of a cmp to be used in class averaging '
                           '(default=ccc)')

        group = form.addGroup('First stage aligner (clsavg)')
        line = group.addLine('classalign: ')
        line.addParam('classalignType', EnumParam,
                      choices=['frm2d', 'rotate_flip',
                               'rotate_flip_iterative', 'rotate_precenter',
                               'rotate_trans_flip_scale',
                               'rotate_trans_flip_scale_iter',
                               'rotate_trans_scale_iter',
                               'rotate_translate', 'rotate_translate_bispec',
                               'rotate_translate_flip',
                               'rotate_translate_flip_iterative',
                               'rotate_translate_flip_resample',
                               'rotate_translate_iterative',
                               'rotate_translate_resample',
                               'rotate_translate_scale', 'rotate_translate_tree',
                               'rotational', 'rotational_bispec',
                               'rotational_iterative', 'rtf_exhaustive',
                               'rtf_slow_exhaustive', 'scale', 'symalign',
                               'symalignquat', 'translational', 'None'],
                      label='', default=ALN_ROTATE_TRANSLATE_TREE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classalignParams', StringParam,
                      default='flip=0', label='params')
        line = group.addLine('classaligncmp: ')
        line.addParam('classaligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classaligncmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('Second stage aligner (clsavg)')
        line = group.addLine('classralign: ')
        line.addParam('classralignType', EnumParam,
                      choices=['None', 'refine'],
                      label='', default=RALN_NONE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classralignParams', StringParam,
                      default='', label='params')
        line = group.addLine('classraligncmp: ')
        line.addParam('classraligncmpType', EnumParam,
                      choices=['ccc', 'dot', 'frc', 'lod', 'optsub',
                               'optvariance', 'phase', 'quadmindot',
                               'sqeuclidean', 'vertical', 'None'],
                      label='', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classraligncmpParams', StringParam,
                      default='', label='params')

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._createIterTemplates(currRun=1)
        self._insertFunctionStep('createLinksStep')
        args = self._prepareParams()
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions ------------------------------
    def createLinksStep(self):
        prot = self._inputProt()
        prevPartDir = prot._getExtraPath("particles")
        currPartDir = self._getExtraPath("particles")
        prevSetsDir = prot._getExtraPath("sets")
        currSetsDir = self._getExtraPath("sets")

        createLink(prevPartDir, currPartDir)
        createLink(prevSetsDir, currSetsDir)

    def refineStep(self, args):
        """ Run the EMAN program to refine 2d. """
        program = eman2.Plugin.getProgram('e2refine2d_bispec.py')
        # mpi and threads are handled by EMAN itself
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        partSet = self._getInputParticles()
        classes2D = self._createSetOfClasses2D(partSet)
        self._fillClassesFromIter(classes2D, self._lastIter())

        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(partSet, classes2D)

    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input CTF protocol: %s" % self.getObjectTag('inputBispec'))
            summary.append("Classified into *%d* classes." % self.numberOfClassAvg)
            summary.append("Output set: %s" % self.getObjectTag('outputClasses'))

        summary.append('\n\n*Note:* output particles are not '
                       'aligned when using this classification method.')
        return summary

    def _methods(self):
        methods = "We classified input particles from %s" % (
            self.getObjectTag('inputBispec'))
        methods += "into %d classes using e2refine2d_bispec.py " %\
                   self.numberOfClassAvg
        return [methods]

    #--------------------------- UTILS functions ------------------------------
    def _prepareParams(self):
        args1 = " --input=%s" % self._getParticlesStack()
        args2 = self._commonParams()
        args = args1 + args2

        return args

    def _commonParams(self):
        args = " --ncls=%(ncls)d --iter=%(numberOfIterations)d --nbasisfp=%(nbasisfp)d"
        args += " --classkeep=%(classKeep)f --classiter=%(classiter)d "
        args += " --classaverager=%s" % self.getEnumText('classAveragerType')

        if self._isVersion23() and self.alignSort:
            args += " --alignsort"

        if self.classKeepSig:
            args += " --classkeepsig"

        for param in ['classnormproc', 'classcmp', 'classalign', 'center',
                      'classaligncmp', 'classralign', 'classraligncmp']:
            args += self._getOptsString(param)

        if self.numberOfMpi > 1:
            args += " --parallel=mpi:%(mpis)d:%(scratch)s --threads=%(threads)d"
        else:
            args += " --parallel=thread:%(threads)d --threads=%(threads)d"

        params = {'ncls': self.numberOfClassAvg.get(),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'nbasisfp': self.nbasisfp.get(),
                  'classKeep': self.classKeep.get(),
                  'classiter': self.classIter.get(),
                  'threads': self.numberOfThreads.get(),
                  'mpis': self.numberOfMpi.get(),
                  'scratch': eman2.SCRATCHDIR
                  }
        args = args % params

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        return args

    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key, **args))

    def _getParticlesStack(self):
        protType = self._inputProt().type.get()
        if protType == HIRES:
            return "sets/" + basename(self._getFileName("partSetFlipLp5"))
        elif protType == MIDRES:
            return "sets/" + basename(self._getFileName("partSetFlipLp7"))
        else:
            return "sets/" + basename(self._getFileName("partSetFlipLp12"))

    def _iterTextFile(self, iterN):
        f = open(self._getFileName('results', iter=iterN))

        for line in f:
            if '#' not in line:
                yield map(float, line.split())

        f.close()

    def _getRun(self):
        return 1

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 2 digits iteration number

        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)

        if clean:
            cleanPath(data_classes)

        if not exists(data_classes):
            clsSet = em.SetOfClasses2D(filename=data_classes)
            clsSet.setImages(self._getInputParticles())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _getInputParticles(self):
        protType = self._inputProt().type.get()
        if protType == HIRES:
            return self._inputProt().outputParticles_flip_lp5
        elif protType == MIDRES:
            return self._inputProt().outputParticles_flip_lp7
        else:
            return self._inputProt().outputParticles_flip_lp12

    def _fillClassesFromIter(self, clsSet, iterN):
        self._execEmanProcess(iterN)
        params = {'orderBy': ['_micId', 'id'],
                  'direction': 'ASC'}
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=self._iterTextFile(iterN),
                             iterParams=params)

    def _execEmanProcess(self, iterN):
        clsFn = self._getFileName("cls", run=1, iter=iterN)
        classesFn = self._getFileName("classes", run=1, iter=iterN)

        proc = eman2.Plugin.createEmanProcess(args='read %s %s %s %s 2d'
                                                   % (self._getParticlesStack(), clsFn, classesFn,
                                                      self._getBaseName('results', iter=iterN)),
                                              direc=self._getExtraPath())
        proc.wait()

        self._classesInfo = {}  # store classes info, indexed by class id
        for classId in range(self.numberOfClassAvg.get()):
            self._classesInfo[classId + 1] = (classId + 1,
                                              self._getExtraPath(classesFn))

    def _getOptsString(self, option):
        optionType = "optionType = self.getEnumText('" + option + "Type')"
        optionParams = 'optionParams = self.' + option + 'Params.get()'
        exec(optionType)
        exec(optionParams)

        if optionType == 'None':
            return ''
        if optionParams != '':
            argStr = ' --%s=%s:%s' % (option, optionType, optionParams)
        else:
            argStr = ' --%s=%s' % (option, optionType)

        return argStr

    def _updateParticle(self, item, row):
        if row[1] == 1:  # enabled
            item.setClassId(row[2] + 1)
        else:
            setattr(item, "_appendItem", False)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn = self._classesInfo[classId]
            item.getRepresentative().setLocation(classId, fn)

    def _inputProt(self):
        return self.inputBispec.get()

    def _isVersion23(self):
        return eman2.Plugin.isVersion('2.3')
