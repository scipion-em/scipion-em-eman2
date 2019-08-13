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

from glob import glob

from pyworkflow.utils.path import cleanPattern
from pyworkflow.protocol.params import (PointerParam, TextParam, IntParam,
                                        BooleanParam, StringParam,
                                        EnumParam, FloatParam)
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.em.data import SetOfClasses2D, SetOfAverages, Volume

import eman2
from eman2.constants import *


class EmanProtInitModelSGD(ProtInitialVolume):
    """
    This protocol wraps *e2initialmodel_sgd.py* EMAN2 program.

    This program makes initial models using a (kind of) stochastic gradient
    descent approach. It is recommended that the box size of
    particles is around 100.
    """

    _label = 'initial model SGD'

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', EnumParam,
                      choices=['Averages', 'Particles'],
                      default=SGD_INPUT_AVG,
                      label='Select input type',
                      help='You can choose either class averages '
                      'or particles as input.')
        form.addParam('inputAvg', PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      condition='inputType==%d' % SGD_INPUT_AVG,
                      # pointerCondition='hasRepresentatives',
                      label="Input averages", important=True,
                      help='Select the class averages to build your '
                           '3D model.\nYou can select SetOfAverages or '
                           'SetOfClasses2D as input.')
        form.addParam('inputPart', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='inputType==%d' % SGD_INPUT_PTCLS,
                      label="Input particles", important=True,
                      help='Select the particles to build your '
                           '3D model.')
        form.addParam('symmetry', TextParam, default='c1',
                      label='Symmetry group',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('batchSize', IntParam, default=10,
                      label='Batch size',
                      help='Batch size of stochastic gradient descent. '
                           'N particles are randomly selected to '
                           'generate an initial model at each step.')
        form.addParam('numberOfIterations', IntParam, default=20,
                      label='Number of iterations to perform',
                      help='The total number of refinement to perform.')
        form.addParam('numberOfModels', IntParam, default=10,
                      label='Number of different initial models',
                      help='The number of different initial models to '
                           'generate in search of a good one.')
        form.addParam('targetRes', FloatParam, default='20.0',
                      label='Target resolution (A)',
                      help='Target resolution in A of the model.')
        form.addParam('shrink', IntParam, default=1,
                      label='Shrink factor',
                      help='Using a box-size >64 is not optimal for making '
                           'initial models. Suggest using this option to '
                           'shrink the input particles by an integer amount '
                           'prior to reconstruction. Default = 1, no shrinking')

        form.addSection('Advanced')
        form.addParam('learnRate', FloatParam, default='0.3',
                      label='Learning rate',
                      help='Learning rate is how much the initial model changes '
                           'toward the gradient direction in each iteration. '
                           'Ranges from 0.0 to 1.0. Default is 0.3')
        form.addParam('lrDecay', FloatParam, default='1.0',
                      label='Learning decay',
                      help='Learning rate multiplier after each iteration.')
        form.addParam('addNoise', FloatParam, default='3.0',
                      label='Add noise',
                      help='Add noise on particles at each iteration. '
                           'Stablize convergence for some reason.')
        form.addParam('fullCov', BooleanParam, default=False,
                      label='Full coverage',
                      help='Assume the input particles covers most of the '
                           'orientation of the model. This gives better '
                           'performance when the model is relatively featureless, '
                           'but is more likely to fail when there are incorrect '
                           'particles in the input.')
        form.addParam('writeTmp', BooleanParam, default=False,
                      label='Write tmp output?',
                      help='Write output for each iteration.')
        form.addParam('extraParams', StringParam, default='',
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "e2initialmodel_sgd.py -h" will print a list '
                           'of possible options.')

        form.addParallelSection(threads=10, mpi=0)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._prepareDefinition()
        self._insertFunctionStep('createStackImgsStep')
        self._insertInitialModelStep()
        self._insertFunctionStep('createOutputStep')

    def _insertInitialModelStep(self):
        args = '--ptcls=input_set.spi'
        if self.shrink > 1:
            args += ' --shrink=%(shrink)d'

        args += ' --ntry=%(numberOfModels)d --niter=%(numberOfIterations)d'
        args += ' --batchsize=%(batchSize)d --targetres=%(targetRes)f'
        args += ' --learnrate=%(learnRate)f --lrdecay=%(lrDecay)f'
        args += ' --addnoise %(addNoise)f --sym=%(symmetry)s'

        if self.writeTmp:
            args += ' --writetmp'
        if self.fullCov:
            args += ' --fullcov'

        args += ' --threads=%(threads)d'

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        self._insertFunctionStep('createInitialModelStep', args % self._params)

    # --------------------------- STEPS functions -----------------------------
    def createStackImgsStep(self):
        imgsFn = self._params['imgsFn']
        inputSet = self._getInputSet()
        if isinstance(inputSet, SetOfClasses2D):
            imgSet = self._createSetOfParticles("_averages")
            for i, cls in enumerate(self.inputAvg.get()):
                img = cls.getRepresentative()
                img.setSamplingRate(cls.getSamplingRate())
                img.setObjId(i + 1)
                imgSet.append(img)
        elif isinstance(inputSet, SetOfAverages):
            imgSet = self.inputAvg.get()
        else:
            imgSet = self.inputPart.get()

        imgSet.writeStack(imgsFn)

    def createInitialModelStep(self, args):
        """ Run the EMAN program to create the initial model. """
        cleanPattern(self._getExtraPath('initmodel_??'))
        program = eman2.Plugin.getProgram('e2initialmodel_sgd.py')
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        volumes = self._createSetOfVolumes()
        shrink = self.shrink.get()
        inputSet = self._getInputSet()
        if isinstance(inputSet, SetOfClasses2D):
            volumes.setSamplingRate(inputSet.getImages().getSamplingRate() * shrink)
        elif isinstance(inputSet, SetOfAverages):
                volumes.setSamplingRate(inputSet.getSamplingRate() * shrink)
        else:
            volumes.setSamplingRate(inputSet.getSamplingRate() * shrink)

        outputVols = self._getVolumes()
        for k, volFn in enumerate(outputVols):
            vol = Volume()
            vol.setFileName(volFn)
            vol.setObjComment('eman initial model %02d' % (k + 1))
            volumes.append(vol)

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(inputSet, volumes)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output initial volumes: %s" % self.outputVolumes.getSize())
        return summary

    # --------------------------- UTILS functions -----------------------------

    def _prepareDefinition(self):
        self._params = {'imgsFn': self._getExtraPath('input_set.spi'),
                        'numberOfIterations': self.numberOfIterations.get(),
                        'numberOfModels': self.numberOfModels.get(),
                        'shrink': self.shrink.get(),
                        'symmetry': self.symmetry.get(),
                        'threads': self.numberOfThreads.get(),
                        'batchSize': self.batchSize.get(),
                        'targetRes': self.targetRes.get(),
                        'learnRate': self.learnRate.get(),
                        'lrDecay': self.lrDecay.get(),
                        'addNoise': self.addNoise.get()}

    def _getVolumes(self):
        outputVols = glob(self._getExtraPath('initmodel_??/model_??.hdf'))
        outputVols.sort()

        return outputVols

    def _getInputSet(self):
        if self.inputType.get() == SGD_INPUT_AVG:
            return self.inputAvg.get()
        else:
            return self.inputPart.get()
