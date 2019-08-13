# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
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

from pyworkflow.utils.properties import Message
from pyworkflow.protocol import params
import pyworkflow.em as pwem
from pyworkflow.em.data import Volume
from eman2.convert import writeSetOfParticles
from pyworkflow.utils.path import makePath

from tomo.protocols import ProtTomoBase

from glob import glob
import eman2

class EmanProtTomoInitialModel(pwem.EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_sgd.py.py* EMAN2 program.

    It will take a set of class-averages/projections and build a set
    of 3-D models suitable for use as initial models in single
    particle reconstruction. The output set is theoretically sorted
    in order of quality (best one is numbered 1), though it's best
    to look at the other answers as well.
    """
    _label = 'tomo initial model'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('particles', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      # pointerCondition='hasRepresentatives',
                      label="Particles", important=True,
                      help='Select the set of subtomograms to build an initial model')

        form.addParam('reference', params.PointerParam,
                      pointerClass='SubTomogram',
                      # pointerCondition='hasRepresentatives',
                      label="Reference", important=True,
                      help='Select the subtomogram to use as reference')

        form.addParam('symmetry', params.TextParam, default='c1',
                      label='Symmetry',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('mask', params.PointerParam,
                      label='Mask',
                      allowsNull=True,
                      pointerClass='VolumeMask',
                      help='Select the subtomogram to use as reference')
        form.addParam('gaussFilter', params.FloatParam, default=-1,
                      label='Gauss',
                      help='The Gaussian filter level')
        form.addParam('filterto', params.FloatParam, default=0.02,
                      label='Filterto',
                      help='Filter map to frequency after each iteration. Default is 0.02')
        form.addParam('fourier', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Fourier',
                      help='Gradient descent in fourier space')
        form.addParam('batchSize', params.IntParam, default=12,
                      label='Batch Size',
                      help='SGD batch size')
        form.addParam('learningRate', params.FloatParam, default=0.1,
                      label='Learn Rate',
                      help='Learning Rate. Default is 0.1')
        form.addParam('numberOfIterations', params.IntParam, default=5,
                      label='Number of iterations to perform',
                      help='The total number of refinement to perform.')
        form.addParam('numberOfBatches', params.IntParam, default=10,
                      label='Number of batches',
                      help='Number of batches per iteration')
        form.addParam('shrink', params.IntParam, default=1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Shrink factor',
                      help='Using a box-size >64 is not optimal for making '
                           'initial models. Suggest using this option to '
                           'shrink the input particles by an integer amount '
                           'prior to reconstruction. Default = 1, no shrinking')
        form.addParam('applySim', params.BooleanParam, default=False,
                      label='Apply Symmetry',
                      help='Apply Symmetry')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # self.convertImagesStep()
        # self.prepareDefinition()
        # self.launchInitialModelStep()
        self._getOutput()

    # Eman2 output is saved on ScipionUserData/projects/{PROJECT_NAME}
    # must fetch from there, read and create scipion output in extra file.
    def _getOutput(self):
        toBaseDir = self.getWorkingDir() + '/../..'

        print({'workingDic': self._enterDir(toBaseDir), 'baseDir': toBaseDir })

    # --------------------------- STEPS functions -----------------------------
    # Get Scipion references to subtomograms and write hdf files for eman2 to process.
    def convertImagesStep(self):
        partSet = self.particles.get()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)

    def launchInitialModelStep(self):
        command_params = self._params
        args = '%s/*.hdf' % self._getExtraPath("particles")
        if command_params['reference']:
            args += ' --reference=%(reference)s'

        args += (' --sym=%(symmetry)s --gaussz=%(gaussFilter)f --filterto=%(filterto)f'
                 ' --batchsize=%(batchSize)d --learnrate=%(learningRate)f --niter=%(numberOfIterations)d'
                 ' --nbatch=%(numberOfBatches)d')

        if command_params['shrink'] > 1:
            args += ' --shrink=%(shrink)d'
        if self.fourier.get():
            args += ' --fourier'
        if self.applySim.get():
            args += ' --applysim'
        if command_params['mask']:
            args += '--mask=%(mask)s'

        program = eman2.Plugin.getProgram("e2spt_sgd.py")
        self._log.info('Launching: ' + program + ' ' + args % command_params)
        self.runJob(program, args % command_params)

        #self._insertFunctionStep('createInitialModelStep', args % self._params)

    ## Missing what to do with output
    def createOutputStep(self):
        setOfSubtomograms = self.particles.get()
        volumes = self._createSetOfVolumes()
        shrink = self.shrink.get()
        samplingRate = self.samplingRate.get()
        volumes.setSamplingRate(samplingRate * shrink)
        outputVols = self._getVolumes()
        for k, volFn in enumerate(outputVols):
            vol = Volume()
            vol.setFileName(volFn)
            vol.setObjComment('eman initial model %02d' % (k + 1))
            volumes.append(vol)

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputSet, volumes)


    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        ProtTomoBase._runSteps(self, startIndex)


    def _getVolumes(self):
        outputVols = glob(self._getExtraPath('initial_models/model_??_??.hdf'))
        outputVols.sort()
        return outputVols

    def prepareDefinition(self):
        self._params = {
            'symmetry': self.symmetry.get(),
            'gaussFilter': self.gaussFilter.get(),
            'filterto': self.filterto.get(),
            'batchSize': self.batchSize.get(),
            'learningRate': self.learningRate.get(),
            'numberOfIterations': self.numberOfIterations.get(),
            'numberOfBatches': self.numberOfBatches.get(),
            'mask': self.mask.get(),
            'shrink': self.shrink.get(),
            'reference': self.reference.get().getFileName()
         }

