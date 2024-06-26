# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from glob import glob
from enum import Enum

from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtRefine3D
from pwem.objects.data import Volume, SetOfParticles
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam, BooleanParam)
from pyworkflow.utils.path import cleanPattern, makePath, createLink

from .. import Plugin
from ..convert import rowToAlignment, writeSetOfParticles
from ..constants import TOPHAT_NONE, SPEED_5, AMP_AUTO, EMAN2SCRATCHDIR


class outputs(Enum):
    outputVolume = Volume
    outputParticles = SetOfParticles


class EmanProtRefine(ProtRefine3D):
    """
    This protocol wraps *e2refine_easy.py* EMAN2 program.

This is the primary single particle refinement program in EMAN2.1+.
It replaces earlier programs such as e2refine.py and e2refine_evenodd.py.

Major features of this program:

 * While a range of command-line options still exist. You should not
 normally specify more than a few basic requirements. The rest will
 be auto-selected for you.
 * This program will split your data in half and automatically
 refine the halves independently to produce a gold standard resolution
 curve for every step in the refinement.
 * An HTML report file will be generated as this program runs,
 telling you exactly what it decided to do and why, as well as giving
 information about runtime, etc while the job is still running.
 * The gold standard FSC also permits us to automatically filter the
 structure at each refinement step. The resolution you specify is
 a target, NOT the filter resolution.
    """
    _label = 'refine easy'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {
            'partSet': 'sets/inputSet.lst',
            'partFlipSet': 'sets/inputSet__ctf_flip.lst',
            'data_scipion': self._getExtraPath('data_scipion_it%(iter)02d.sqlite'),
            'projections': self._getExtraPath('projections_it%(iter)02d_%(half)s.sqlite'),
            'classes': 'refine_%(run)02d/classes_%(iter)02d',
            'classesEven': self._getExtraPath('refine_%(run)02d/classes_%(iter)02d_even.hdf'),
            'classesOdd': self._getExtraPath('refine_%(run)02d/classes_%(iter)02d_odd.hdf'),
            'cls': 'refine_%(run)02d/cls_result_%(iter)02d',
            'clsEven': self._getExtraPath('refine_%(run)02d/cls_result_%(iter)02d_even.hdf'),
            'clsOdd': self._getExtraPath('refine_%(run)02d/cls_result_%(iter)02d_odd.hdf'),
            'angles': self._getExtraPath('projectionAngles_it%(iter)02d.txt'),
            'mapEven': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d_even.hdf'),
            'mapOdd': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d_odd.hdf'),
            'mapFull': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d.hdf'),
            'mapEvenUnmasked': self._getExtraPath('refine_%(run)02d/threed_even_unmasked.hdf'),
            'mapOddUnmasked': self._getExtraPath('refine_%(run)02d/threed_odd_unmasked.hdf'),
            'fscUnmasked': self._getExtraPath('refine_%(run)02d/fsc_unmasked_%(iter)02d.txt'),
            'fscMasked': self._getExtraPath('refine_%(run)02d/fsc_masked_%(iter)02d.txt'),
            'fscMaskedTight': self._getExtraPath('refine_%(run)02d/fsc_maskedtight_%(iter)02d.txt'),
            'reportHtml': self._getExtraPath('refine_%(run)02d/report/index.html')
        }
        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self, currRun):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('mapFull', run=currRun,
                                               iter=1).replace('threed_01',
                                                               'threed_??')
        # Iterations will be identify by threed_XX_ where XX is the iteration
        # number and is restricted to only 2 digits.
        self._iterRegex = re.compile(r'threed_(\d{2})')

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous '
                           'run of type *%s* class and most of the input parameters '
                           'will be taken from it.' % self.getClassName())
        form.addParam('continueRun', PointerParam,
                      pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('inputParticles', PointerParam,
                      label="Input particles",
                      important=True, pointerClass='SetOfParticles',
                      condition='not doContinue', allowsNull=True,
                      help='Select the input particles.')
        form.addParam('input3DReference', PointerParam,
                      important=True,
                      pointerClass='Volume', allowsNull=True,
                      label='Initial 3D reference volume',
                      condition='not doContinue',
                      help='Input 3D reference reconstruction.')
        form.addParam('skipctf', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Skip ctf estimation?',
                      help='Use this if you want to skip running e2ctf.py. '
                           'It is not recommended to skip this step unless CTF '
                           'estimation was already done with EMAN2.')
        form.addParam('numberOfIterations', IntParam, default=6,
                      label='Number of iterations',
                      help='The total number of refinement iterations to '
                           'perform.')
        form.addParam('tophat', EnumParam,
                      choices=['none', 'local', 'localwiener', 'global'],
                      label="Tophat filter?", default=TOPHAT_NONE,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Instead of imposing a final '
                           'Wiener filter (tophat = none)), use a tophat '
                           'filter (global similar to Relion). local '
                           'determines local resolution and '
                           'filters. Danger of feature exaggeration.')
        form.addParam('symmetry', StringParam, default='c1',
                      condition='not doContinue',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then the '
                           'model is assumed to have no symmetry. \n'
                           'Choices are: c(n), d(n), tet, icos, or oct.\n'
                           'See https://blake.bcm.edu/emanwiki/EMAN2/Symmetry '
                           'for a detailed description of symmetry in Eman.')
        form.addParam('doBreaksym', BooleanParam, default=False,
                      label='Break symmetry?',
                      help='If set True, reconstruction will be asymmetric '
                           'with *Symmetry group* parameter specifying a '
                           'known pseudosymmetry, not an imposed symmetry.')
        form.addParam('resol', FloatParam, default=25.0,
                      label='Target resolution (A)',
                      help='Target resolution in A of this refinement run. '
                           'Usually works best in at least two steps '
                           '(low/medium) resolution, then final resolution) '
                           'when starting with a poor starting model. '
                           'Usually 3-4 iterations is sufficient.')
        form.addParam('molMass', FloatParam, default=500.0,
                      label='Molecular mass (kDa)',
                      help='Approximate molecular mass of the particle, in kDa. '
                           'This is used to run normalize.bymass. Due to '
                           'resolution effects, not always the true mass.')
        form.addParam('useE2make3d', BooleanParam, default=False,
                      label='Use old e2make3d?',
                      help='Use the traditional e2make3d program instead of '
                           'the new e2make3dpar program.')
        form.addParam('maskExpand', IntParam, default=-1,
                      label='Expand mask by (px)',
                      help='Default=boxsize/20. Specify number of voxels to '
                           'expand mask before soft edge. Use this if low '
                           'density peripheral features are cut off by the mask.')
        form.addParam('noRandPhase', BooleanParam, default=False,
                      label='Supress phase randomization',
                      help='Suppress independent phase randomization '
                           'of input map. Only appropriate if input '
                           'map has been preprocessed in some suitable '
                           'fashion.')

        form.addSection(label='Advanced')
        form.addParam('speed', EnumParam,
                      choices=['1', '2', '3', '4', '5', '6', '7'],
                      label="Speed", default=SPEED_5,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Balances speed vs precision. '
                           'Larger values sacrifice a bit of potential '
                           'resolution for significant speed increases. Set '
                           'to 1 when really pushing resolution. Set to 7 for '
                           'initial refinements.')
        form.addParam('classKeep', FloatParam, default=0.9,
                      label='Fraction of particles to use in final average',
                      help='The fraction of particles to keep in each class,'
                           'based on the similarity score.')
        form.addParam('m3dKeep', FloatParam, default=0.8,
                      label='Fraction of class-averages to use in 3-D map',
                      help='The fraction of slices to keep in reconstruction.')
        form.addParam('useBispec', BooleanParam, default=False,
                      label='Use bispectra? (experimental)',
                      help='Will use bispectra for orientation '
                           'determination (EXPERIMENTAL).')
        form.addParam('useSetsfref', BooleanParam, default=False,
                      label='Use the setsfref option in class averaging?',
                      help='This matches the filtration of the class-averages '
                           'to the projections for easier comparison. May '
                           'also improve convergence. '
                           'Disabled when ampcorrect=flatten is used.')
        form.addParam('doAutomask', BooleanParam, default=False,
                      label='Do automask to the class-average?',
                      help='This will apply an automask to the class-average '
                           'during iterative alignment for better accuracy. '
                           'The final class averages are unmasked.')
        form.addParam('doThreshold', BooleanParam, default=False,
                      label='Apply threshold before project the volume?',
                      help='Applies a threshold to the volume just before '
                           'generating projections. A sort of aggressive '
                           'solvent flattening for the reference.')
        form.addParam('m3dPostProcess', StringParam,
                      default='none',
                      label='Postprocess parameters',
                      help="<name>:<parm>=<value>:...  An arbitrary processor "
                           "(e2help.py processors -v2) to apply to the 3-D map "
                           "after each iteration. Default=none")
        form.addParam('ampCorrect', EnumParam,
                      choices=['auto', 'strucfac', 'flatten', 'none'],
                      label="Amplitude correction:", default=AMP_AUTO,
                      display=EnumParam.DISPLAY_COMBO,
                      help="Will perform amplitude correction via the specified "
                           "method. 'flatten' requires a target resolution better "
                           "than 8 angstroms (experimental). 'none' will disable "
                           "amplitude correction (experimental).")
        form.addParam('extraParams', StringParam,
                      default='',
                      label='Additional parameters',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. "
                           "See e2refine_easy.py -h.")
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._createIterTemplates(self._getRun())
        if self.doContinue:
            self.input3DReference.set(None)
            self.inputParticles.set(None)
            self._insertFunctionStep('createLinkSteps')
            args = self._prepareContinueParams()
        else:
            self._insertFunctionStep('convertImagesStep')
            args = self._prepareParams()
        self._insertFunctionStep('refineStep', args)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def createLinkSteps(self):
        continueRun = self.continueRun.get()
        prevPartDir = continueRun._getExtraPath("particles")
        currPartDir = self._getExtraPath("particles")
        runN = self._getRun()
        prevRefDir = continueRun._getExtraPath("refine_%02d" % runN)
        currRefDir = self._getExtraPath("refine_%02d" % runN)
        prevSetsDir = continueRun._getExtraPath("sets")
        currSetsDir = self._getExtraPath("sets")

        createLink(prevPartDir, currPartDir)
        createLink(prevRefDir, currRefDir)
        createLink(prevSetsDir, currSetsDir)

    def convertImagesStep(self):
        partSet = self._getInputParticles()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)
        if not self.skipctf:
            program = Plugin.getProgram('e2ctf.py')
            acq = partSet.getAcquisition()

            args = " --voltage %d" % acq.getVoltage()
            args += " --cs %f" % acq.getSphericalAberration()
            args += " --ac %f" % (100 * acq.getAmplitudeContrast())
            args += " --threads=%d" % self.numberOfThreads.get()
            if not partSet.isPhaseFlipped():
                args += " --phaseflip"
            args += " --computesf --apix %f" % partSet.getSamplingRate()
            args += " --allparticles --autofit --curdefocusfix --storeparm -v 8"
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

        program = Plugin.getProgram('e2buildsets.py')
        args = " --setname=inputSet --allparticles"
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def refineStep(self, args):
        """ Run the EMAN program to refine a volume. """
        if not self.doContinue:
            cleanPattern(self._getExtraPath('refine_01'))
        program = Plugin.getProgram('e2refine_easy.py')
        # mpi and threads are handled by EMAN itself
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        iterN = self.numberOfIterations.get()
        partSet = self._getInputParticles(pointer=True)
        numRun = self._getRun()

        vol = Volume()
        vol.setFileName(self._getFileName("mapFull", run=numRun, iter=iterN))
        halfMap1 = self._getFileName("mapEvenUnmasked", run=numRun)
        halfMap2 = self._getFileName("mapOddUnmasked", run=numRun)
        vol.setHalfMaps([halfMap1, halfMap2])
        vol.copyInfo(partSet.get())

        newPartSet = self._createSetOfParticles()
        newPartSet.copyInfo(partSet.get())
        self._fillDataFromIter(newPartSet, iterN)

        self._defineOutputs(**{outputs.outputVolume.name: vol,
                               outputs.outputParticles.name: newPartSet})
        self._defineSourceRelation(partSet, vol)
        self._defineTransformRelation(partSet, newPartSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        particles = self._getInputParticles()
        samplingRate = particles.getSamplingRate()

        if self.resol < 2 * samplingRate:
            errors.append("\nTarget resolution is smaller than Nyquist limit.")

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume is not ready yet.")
        else:
            inputSize = self._getInputParticles().getSize()
            outputSize = self.outputParticles.getSize()
            diff = inputSize - outputSize
            if diff > 0:
                summary.append("Warning!!! %d particles "
                               "were discarded during refinement." % diff)

        summary.append("To see progress report, click "
                       "*Analyze Results*  and choose *Show "
                       "HTML report*.")
        return summary

    def _citations(self):
        return ['Bell2016']

    # --------------------------- UTILS functions -----------------------------
    def _prepareParams(self):
        args1 = " --input=%(imgsFn)s --model=%(volume)s"
        args2 = self._commonParams()

        refVolFn = "ref_vol.hdf"
        origVol = os.path.relpath(self.input3DReference.get().getFileName(),
                                  self._getExtraPath()).replace(":mrc", "")
        args = "%s %s --apix=%0.3f" % (origVol, refVolFn,
                                       self.input3DReference.get().getSamplingRate())
        self.runJob(Plugin.getProgram('e2proc3d.py'), args,
                    cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

        params = {'imgsFn': self._getParticlesStack(),
                  'volume': refVolFn}

        args = args1 % params + args2
        return args

    def _prepareContinueParams(self):
        runN = self._getRun()
        args1 = "--startfrom=refine_%02d" % runN
        args2 = self._commonParams()
        args = args1 + args2
        return args

    def _commonParams(self):
        args = " --targetres=%(resol)f --speed=%(speed)d --sym=%(sym)s "
        args += " --iter=%(numberOfIterations)d --mass=%(molMass)f "
        args += " --apix=%(samplingRate)f --classkeep=%(classKeep)f"
        if self.numberOfMpi > 1:
            args += " --m3dkeep=%(m3dKeep)f --parallel=mpi:%(mpis)d:%(scratch)s"
        else:
            args += " --m3dkeep=%(m3dKeep)f --parallel=thread:%(threads)d"
        args += " --threads=%(threads)d"

        samplingRate = self._getInputParticles().getSamplingRate()
        params = {'resol': self.resol.get(),
                  'speed': int(self.getEnumText('speed')),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'sym': self.symmetry.get(),
                  'molMass': self.molMass.get(),
                  'samplingRate': samplingRate,
                  'classKeep': self.classKeep.get(),
                  'm3dKeep': self.m3dKeep.get(),
                  'threads': self.numberOfThreads.get(),
                  'mpis': self.numberOfMpi.get(),
                  'scratch': Plugin.getVar(EMAN2SCRATCHDIR)
                  }
        args %= params

        if self.doBreaksym:
            args += " --breaksym"
        if self.useE2make3d:
            args += " --m3dold"
        if self.maskExpand.get() != -1:
            args += " --automaskexpand=%d" % self.maskExpand.get()
        if self.useSetsfref:
            args += " --classrefsf"
        if self.doAutomask:
            args += " --classautomask"
        if self.doThreshold:
            args += " --prethreshold"
        if self.m3dPostProcess.get() != 'none':
            args += " --m3dpostprocess=%s" % self.m3dPostProcess.get()

        args += " --ampcorrect=%s" % self.getEnumText('ampCorrect')

        if self.tophat != TOPHAT_NONE:
            args += " --tophat=%s" % self.getEnumText('tophat')

        if self.useBispec:
            args += " --invar"

        if self.noRandPhase:
            args += " --norandomaphase"

        if self.extraParams.hasValue():
            args += ' ' + self.extraParams.get()
        return args

    def _getRun(self):
        if not self.doContinue:
            return 0
        else:
            files = sorted(glob(self.continueRun.get()._getExtraPath("refine*")))
            if files:
                f = files[-1]
                refineNumber = int(f.split("_")[-1]) + 1
                return refineNumber

    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key, **args))

    def _getParticlesStack(self):
        if not self.inputParticles.get().isPhaseFlipped() and not self.skipctf:
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")

    def _iterTextFile(self, iterN):
        with open(self._getFileName('angles', iter=iterN)) as f:
            for line in f:
                if '#' not in line and line.strip():
                    yield [float(x) for x in line.split()]

    def _createItemMatrix(self, item, rowList):
        if rowList[1] == 1:
            item.setTransform(rowToAlignment(rowList[2:],
                                             alignType=ALIGN_PROJ))
        else:
            setattr(item, "_appendItem", False)

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 3 digits iteration number

        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterData(self, it):
        data_sqlite = self._getFileName('data_scipion', iter=it)
        if not os.path.exists(data_sqlite):
            iterImgSet = SetOfParticles(filename=data_sqlite)
            iterImgSet.copyInfo(self._getInputParticles())
            self._fillDataFromIter(iterImgSet, it)
            iterImgSet.write()
            iterImgSet.close()

        return data_sqlite

    def _getInputParticles(self, pointer=False):
        if self.doContinue:
            self.inputParticles.set(self.continueRun.get().inputParticles.get())
        if pointer:
            return self.inputParticles
        else:
            return self.inputParticles.get()

    def _fillDataFromIter(self, imgSet, iterN):
        numRun = self._getRun()
        self._execEmanProcess(numRun, iterN)
        initPartSet = self._getInputParticles()
        imgSet.setAlignmentProj()
        partIter = iter(initPartSet.iterItems(orderBy=['_micId', 'id'],
                                              direction='ASC'))

        imgSet.copyItems(partIter,
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=self._iterTextFile(iterN))

    def _execEmanProcess(self, numRun, iterN):
        clsFn = self._getFileName("cls", run=numRun, iter=iterN)
        classesFn = self._getFileName("classes", run=numRun, iter=iterN)
        angles = self._getFileName('angles', iter=iterN)

        if not os.path.exists(angles) and os.path.exists(self._getFileName('clsEven',
                                                                           run=numRun, iter=iterN)):
            proc = Plugin.createEmanProcess(args='read %s %s %s %s 3d'
                                                 % (self._getParticlesStack(), clsFn, classesFn,
                                                    self._getBaseName('angles', iter=iterN)),
                                            direc=self._getExtraPath())
            proc.wait()
