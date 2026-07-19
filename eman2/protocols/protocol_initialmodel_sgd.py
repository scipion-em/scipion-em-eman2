# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from glob import glob
from enum import Enum

from pyworkflow.utils.path import cleanPattern
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, IntParam,
                                        BooleanParam, StringParam,
                                        EnumParam, FloatParam)
from pwem.protocols import ProtInitialVolume
from pwem.objects.data import SetOfClasses2D, SetOfAverages, Volume, SetOfVolumes

from .. import Plugin
from ..constants import SGD_INPUT_AVG, SGD_INPUT_PTCLS


class outputs(Enum):
    outputVolumes = SetOfVolumes


class EmanProtInitModelSGD(ProtInitialVolume):
    """
    Generates initial 3D cryo-EM models from particle images or class
    averages using a stochastic gradient descent strategy implemented in
    EMAN2. The protocol is intended for ab initio structure determination
    when no reliable starting reference is available and provides multiple
    candidate volumes for exploratory structural analysis.

    AI Generated:

    Initial Model SGD (EmanProtInitModelSGD) — User Manual
        Overview

        The Initial Model SGD protocol creates three-dimensional initial
        models from cryo-EM particle images or two-dimensional class
        averages using a stochastic gradient descent optimization strategy.
        Its main purpose is to generate plausible low-resolution structural
        references that can serve as starting points for downstream
        refinement and classification workflows.

        In single-particle cryo-EM, obtaining a reliable initial model is
        often one of the most critical and challenging stages of the entire
        reconstruction process. A good initial volume helps guide later
        refinement toward biologically meaningful solutions, while poor or
        biased starting models may lead to incorrect structural
        interpretations. This protocol is designed to reduce reference bias
        by constructing models directly from experimental image data.

        Inputs and General Workflow

        The protocol accepts either particle images or class averages as
        input. Class averages are often preferred because they contain
        improved signal-to-noise ratio and partially reduce image
        heterogeneity. However, particles may also be used directly when
        averaging has not yet been performed or when preserving structural
        diversity is important.

        During execution, the protocol repeatedly samples subsets of the
        input data and incrementally updates candidate three-dimensional
        models. Multiple independent models can be generated in parallel,
        allowing users to compare alternative structural solutions and
        identify reproducible features across reconstructions.

        This strategy is particularly valuable in difficult datasets where
        orientation assignment is uncertain or where the specimen may contain
        flexible or heterogeneous conformations.

        Stochastic Gradient Descent Strategy

        The reconstruction process relies on iterative optimization in which
        the current model is gradually adjusted according to information
        extracted from randomly selected image subsets. This stochastic
        behavior helps avoid strong dependence on any individual subset of
        particles and can improve convergence stability in noisy cryo-EM
        datasets.

        The batch size determines how many particles contribute to each
        optimization step. Smaller batches introduce more stochastic
        variation and may help explore a broader solution space, whereas
        larger batches generally produce more stable but potentially less
        exploratory convergence behavior.

        The number of iterations controls how extensively the optimization is
        refined. More iterations may improve convergence quality, although
        excessive refinement at the initial-model stage is rarely necessary
        because subsequent high-resolution refinement protocols will further
        optimize the reconstruction.

        Multiple Initial Models

        One of the biologically important features of this protocol is the
        ability to generate several independent initial models during the
        same execution. This approach is highly recommended in exploratory
        cryo-EM projects because it helps assess reconstruction robustness
        and detect potential convergence artifacts.

        When multiple independently generated models converge toward similar
        structural features, confidence in the biological validity of the
        reconstruction increases substantially. Conversely, strongly
        different outcomes may indicate insufficient data quality,
        orientation bias, structural heterogeneity, or optimization
        instability.

        In practical cryo-EM workflows, users often compare these candidate
        models visually and select the most biologically plausible structure
        for downstream refinement.

        Symmetry Considerations

        The protocol supports symmetry imposition during initial model
        generation. Correct symmetry specification can significantly improve
        reconstruction quality by reinforcing equivalent structural views and
        increasing effective signal strength.

        Common biological examples include cyclic symmetry in ring-like
        assemblies, dihedral symmetry in multimeric complexes, or
        icosahedral symmetry in viral particles. Applying the correct
        symmetry often accelerates convergence and improves structural
        consistency.

        However, users should apply symmetry carefully. Incorrect symmetry
        assumptions may artificially distort the reconstruction, obscure
        asymmetric features, or create misleading structural artifacts. When
        structural asymmetry or compositional variability is suspected,
        reconstructing without symmetry is usually the safest initial
        strategy.

        Resolution and Shrinking

        The target resolution parameter defines the approximate resolution
        goal for the generated models. Since initial models are intended only
        as low-resolution starting references, aggressive high-resolution
        reconstruction is generally unnecessary at this stage.

        The protocol also allows shrinking the particle images before
        reconstruction. Reducing image size substantially decreases
        computational cost and often improves optimization stability,
        especially for large particle boxes. In many practical workflows,
        moderate shrinking is beneficial for rapid exploratory model
        generation.

        Biologically, shrinking sacrifices fine structural detail but
        preserves the large-scale architecture required for reliable initial
        orientation estimation and downstream refinement initialization.

        Learning Rate and Optimization Stability

        The learning rate determines how strongly the model changes during
        each optimization step. Larger values accelerate convergence but may
        increase instability, while smaller values provide smoother but
        slower optimization.

        Learning decay gradually reduces the update magnitude over time,
        helping stabilize convergence as the reconstruction progresses. This
        behavior is often beneficial when approaching a consistent structural
        solution.

        The protocol also allows the addition of controlled noise during
        optimization. Although counterintuitive from a biological
        perspective, introducing noise can sometimes improve convergence
        robustness by preventing optimization from becoming trapped in poor
        local solutions.

        Orientation Coverage

        The protocol includes an option to assume broad orientation coverage
        within the dataset. This assumption may improve performance for
        relatively featureless particles when the input images sample most
        orientations evenly.

        However, this strategy becomes less reliable when the dataset
        contains substantial contamination, incorrect particles, or severe
        preferred orientation effects. In such cases, more conservative
        reconstruction assumptions are generally safer.

        Outputs and Interpretation

        The protocol produces one or more initial three-dimensional volumes
        representing candidate structural solutions derived from the input
        images. These volumes are typically low to intermediate resolution
        and are intended primarily for subsequent refinement rather than
        direct biological interpretation.

        Biologically meaningful features at this stage usually include global
        shape, major domains, overall symmetry, and large conformational
        organization. Fine structural details should not be overinterpreted
        because the models remain highly dependent on limited orientation
        accuracy and low-resolution optimization.

        Practical Recommendations

        For most cryo-EM datasets, class averages provide more stable inputs
        than raw particles and are often the preferred starting point for
        initial model generation. Using moderate shrinking and generating
        multiple candidate models are generally recommended practices.

        Users should carefully compare independent models and verify that
        major structural features remain reproducible across solutions before
        proceeding to high-resolution refinement. Visual inspection and
        biological plausibility remain essential evaluation criteria at this
        stage.

        In difficult datasets with preferred orientations or substantial
        heterogeneity, reducing optimization aggressiveness and increasing
        the diversity of candidate models may improve robustness.

        Final Perspective

        Initial model generation is one of the most biologically sensitive
        stages in cryo-EM processing because it establishes the structural
        framework for all subsequent refinement. Reliable initial models
        emerge not only from computational optimization but also from careful
        experimental data selection, thoughtful symmetry assumptions, and
        critical biological interpretation. Producing multiple reproducible
        candidate structures is often the best strategy for ensuring robust
        and trustworthy cryo-EM reconstructions.
    """

    _label = 'initial model SGD'
    _devStatus = PROD
    _possibleOutputs = outputs

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
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See https://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
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
        form.addParam('targetRes', FloatParam, default=20.0,
                      label='Target resolution (A)',
                      help='Target resolution in A of the model.')
        form.addParam('shrink', IntParam, default=1,
                      label='Shrink factor',
                      help='Using a box-size >64 is not optimal for making '
                           'initial models. Suggest using this option to '
                           'shrink the input particles by an integer amount '
                           'prior to reconstruction. Default = 1, no shrinking')

        form.addSection('Advanced')
        form.addParam('learnRate', FloatParam, default=0.3,
                      label='Learning rate',
                      help='Learning rate is how much the initial model changes '
                           'toward the gradient direction in each iteration. '
                           'Ranges from 0.0 to 1.0. Default is 0.3')
        form.addParam('lrDecay', FloatParam, default=1.0,
                      label='Learning decay',
                      help='Learning rate multiplier after each iteration.')
        form.addParam('addNoise', FloatParam, default=3.0,
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
        self._insertFunctionStep('createStackImgsStep', needsGPU=False)
        self._insertInitialModelStep()
        self._insertFunctionStep('createOutputStep', needsGPU=False)

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

        self._insertFunctionStep('createInitialModelStep', args % self._params,
                                 needsGPU=False)

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
        program = Plugin.getProgram('e2initialmodel_sgd.py')
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        volumes = self._createSetOfVolumes()
        shrink = self.shrink.get()
        inputSet = self._getInputSet()
        if isinstance(inputSet, SetOfClasses2D):
            volumes.setSamplingRate(inputSet.getImages().getSamplingRate() * shrink)
        else:
            volumes.setSamplingRate(inputSet.getSamplingRate() * shrink)

        outputVols = self._getVolumes()
        for k, volFn in enumerate(outputVols):
            vol = Volume()
            vol.setFileName(volFn)
            vol.setObjComment('eman initial model %02d' % (k + 1))
            volumes.append(vol)

        self._defineOutputs(**{outputs.outputVolumes.name: volumes})
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
