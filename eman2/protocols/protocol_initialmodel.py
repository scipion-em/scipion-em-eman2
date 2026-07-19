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
from glob import glob
from enum import Enum

from pyworkflow.utils.path import cleanPattern
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, IntParam,
                                        BooleanParam, LEVEL_ADVANCED,
                                        StringParam)
from pwem.protocols import ProtInitialVolume
from pwem.objects.data import SetOfClasses2D, Volume, SetOfVolumes

from .. import Plugin
from ..constants import EMAN2SCRATCHDIR


class outputs(Enum):
    outputVolumes = SetOfVolumes


class EmanProtInitModel(ProtInitialVolume):
    """
    Generates initial 3D models for single-particle cryo-EM analysis using
    EMAN2 ab initio reconstruction strategies. The protocol is designed to
    transform a collection of 2D class averages or averaged particle images
    into one or more candidate 3D volumes suitable for downstream refinement
    and structural interpretation. Multiple independent models can be created
    in parallel to increase the probability of obtaining a biologically
    meaningful starting structure. More info:
    https://blake.bcm.edu/emanwiki/EMAN2/Programs/e2initialmodel

    AI Generated:

    Initial Model Generation (EmanProtInitModel) - User Manual
        Overview

        The Initial Model protocol provides an entry point for ab initio
        three-dimensional reconstruction in cryo-EM workflows. Its purpose is
        to estimate one or several plausible 3D maps directly from 2D class
        averages or averaged particle projections, without requiring a prior
        structural reference. This stage is one of the most biologically
        important moments in single-particle analysis because the quality and
        reliability of the initial model strongly influence all subsequent
        refinement and interpretation steps.

        The protocol is especially useful when studying macromolecular
        assemblies whose structure is unknown or when validating whether a
        reconstruction can emerge consistently from the experimental data
        itself. By generating multiple candidate maps, the procedure allows
        users to compare alternative structural solutions and identify the
        most stable or biologically realistic reconstruction.

        Inputs and Biological Context

        The protocol accepts either 2D class averages or sets of averaged
        particle projections. In practice, these inputs should already
        represent meaningful structural views with reduced noise and improved
        signal quality. Well-defined classes with broad angular coverage are
        critical for successful model generation.

        From a biological perspective, the quality of the input data often
        determines the interpretability of the resulting volumes more than any
        advanced parameter adjustment. Poorly aligned classes, preferred
        orientations, or strong structural heterogeneity may produce unstable
        or ambiguous initial models. For this reason, users are encouraged to
        inspect class averages carefully before attempting ab initio
        reconstruction.

        Symmetry and Structural Assumptions

        The protocol supports several symmetry groups commonly encountered in
        cryo-EM studies, including cyclic, dihedral, tetrahedral, octahedral,
        and icosahedral symmetries. Correct symmetry assignment is biologically
        essential because it constrains the reconstruction and strongly
        influences the resulting map quality.

        Applying the correct symmetry can dramatically improve signal and
        structural consistency. However, imposing an incorrect symmetry may
        introduce artificial features or distort biologically meaningful
        asymmetries. When uncertainty exists, it is generally safer to begin
        with low symmetry assumptions and increase constraints only after
        validation.

        Highly symmetric particles such as viral capsids or symmetric protein
        cages benefit particularly from dedicated reconstruction strategies.
        These systems often converge more rapidly and produce cleaner maps due
        to the large amount of redundant structural information provided by
        symmetry.

        Iterative Model Exploration

        The protocol can generate multiple candidate models through repeated
        reconstruction attempts. This strategy is valuable because ab initio
        reconstruction may converge toward different structural solutions,
        particularly when the experimental data are noisy or structurally
        heterogeneous.

        Producing several models allows the user to compare consistency across
        independent reconstructions. If multiple runs converge toward similar
        volumes, confidence in the biological validity of the structure
        increases substantially. Conversely, strong variability among models
        may indicate insufficient angular coverage, excessive heterogeneity,
        or poor data quality.

        The number of iterations controls how extensively the protocol refines
        candidate structures. Increasing iterations may improve convergence,
        but excessive refinement at this early stage can sometimes reinforce
        noise or overfit unstable features.

        Shrink Factor and Computational Efficiency

        The protocol supports optional downsampling of the input images before
        reconstruction. This approach is commonly used when working with large
        particle boxes because low-resolution structural information is often
        sufficient for generating an initial model.

        Biologically, shrinking is usually acceptable during early exploratory
        reconstruction because the goal is to establish the overall molecular
        architecture rather than recover fine structural detail. Reduced image
        sizes also decrease computational cost and improve convergence speed.

        However, excessive downsampling may remove weak but biologically
        relevant features, particularly in flexible complexes or assemblies
        with small domains. Users should therefore balance computational
        efficiency against the structural complexity of the target system.

        Randomization and Model Diversity

        The protocol can introduce randomized orientation strategies during
        model initialization. This helps avoid reconstruction bias and
        encourages broader exploration of possible structural solutions.

        In practical cryo-EM analysis, this option is especially useful when
        the orientation distribution of the particles is incomplete or when
        the data contain substantial uncertainty. Randomization may increase
        robustness against local convergence artifacts, although it can also
        produce greater variability between candidate maps.

        Biological users should interpret model diversity carefully. Consistent
        convergence across randomized attempts often indicates that the data
        contain strong structural information. Highly divergent solutions may
        suggest heterogeneity, preferred orientations, or insufficient data
        quality.

        Mask Expansion and Peripheral Features

        Advanced masking controls help preserve peripheral densities and avoid
        truncation of extended structural regions. This becomes particularly
        important for elongated assemblies, membrane-associated complexes, or
        flexible macromolecular systems where biologically relevant densities
        may extend toward the edges of the reconstruction.

        Appropriate mask expansion improves continuity of peripheral regions
        while reducing the risk of artificially clipping flexible domains.
        Nevertheless, excessively permissive masking may retain noise and
        destabilize convergence.

        Outputs and Interpretation

        The protocol produces one or more candidate initial volumes that can
        be inspected, validated, and used as starting references for later
        refinement procedures. The generated maps are typically ranked by
        overall reconstruction quality, although users are encouraged to
        examine all candidate solutions rather than relying exclusively on the
        first result.

        From a biological perspective, the initial model should be evaluated
        based on structural plausibility, internal consistency, and agreement
        with known biochemical information. Features such as unrealistic
        disconnected densities, strong directional artifacts, or unstable
        conformations may indicate unreliable reconstruction.

        Initial models are not final biological structures. Instead, they
        serve as starting hypotheses that guide subsequent high-resolution
        refinement and validation workflows.

        Practical Recommendations

        In routine cryo-EM workflows, users should begin with carefully
        selected class averages that cover the widest possible range of
        orientations. Generating several candidate models is generally
        recommended, particularly for novel targets or heterogeneous samples.

        Moderate downsampling is often beneficial for large particles, while
        highly symmetric systems should always use the appropriate symmetry
        definition to maximize reconstruction quality. When results appear
        unstable, increasing the diversity of reconstruction attempts may help
        identify reproducible structural features.

        Visual inspection remains critical throughout the process. The most
        biologically meaningful model is not always the numerically highest
        ranked one, especially in difficult datasets containing flexibility or
        compositional variability.

        Final Perspective

        Initial model generation represents the transition from two-dimensional
        experimental observations to a three-dimensional structural hypothesis.
        The quality of this transition strongly affects all downstream cryo-EM
        analysis steps. Careful selection of input averages, appropriate
        symmetry assumptions, and critical interpretation of candidate models
        are essential for obtaining reliable and biologically meaningful
        reconstructions.
    """

    _label = 'initial model'
    _devStatus = PROD
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label="Input averages", important=True,
                      help='Select the your class averages to build your '
                           '3D model.\nYou can select SetOfAverages or '
                           'SetOfClasses2D as input.')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See https://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('numberOfIterations', IntParam, default=8,
                      label='Number of iterations to perform',
                      help='The total number of refinement to perform.')
        form.addParam('numberOfModels', IntParam, default=10,
                      label='Number of different initial models',
                      help='The number of different initial models to '
                           'generate in search of a good one.')
        form.addParam('shrink', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Shrink factor',
                      help='Using a box-size >64 is not optimal for making '
                           'initial models. Suggest using this option to '
                           'shrink the input particles by an integer amount '
                           'prior to reconstruction. Default = 1, no shrinking')
        form.addParam('randOrient', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use random orientations?',
                      help='Instead of seeding with a random volume, '
                           'seeds by randomizing input orientations')
        form.addParam('autoMaskExp', IntParam, default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Automask expand (px)',
                      help='Number of voxels of post-threshold expansion '
                           'in the mask, for use when peripheral '
                           'features are truncated '
                           '(default=shrunk boxsize/20)')
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "e2initialmodel.py -h" will print a list '
                           'of possible options.')

        form.addParallelSection(threads=8, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._prepareDefinition()
        self._insertFunctionStep('createStackImgsStep', needsGPU=False)
        self._insertInitialModelStep()
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    def _insertInitialModelStep(self):
        args = '--input=%(relImgsFn)s --sym=%(symmetry)s'
        if self.shrink > 1:
            args += ' --shrink=%(shrink)d'
        if not self._isHighSym():
            args += ' --tries=%(numberOfModels)d --iter=%(numberOfIterations)d'
            if self.randOrient:
                args += ' --randorient'
            if self.autoMaskExp.get() != -1:
                args += '--automaskexpand %d'
            if self.numberOfMpi > 1:
                args += ' --parallel=mpi:%(mpis)d:%(scratch)s'
            else:
                args += ' --parallel=thread:%(threads)d'
        else:
            args += ' --threads=%(threads)d'
        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        self._insertFunctionStep('createInitialModelStep', args % self._params,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def createStackImgsStep(self):
        if isinstance(self.inputSet.get(), SetOfClasses2D):
            pixSize = self.inputSet.get().getImages().getSamplingRate()
            imgSet = self._createSetOfParticles("_averages")
            for i, cls in enumerate(self.inputSet.get()):
                img = cls.getRepresentative()
                img.setSamplingRate(pixSize)
                img.setObjId(i + 1)
                imgSet.append(img)
        else:
            imgSet = self.inputSet.get()
            pixSize = imgSet.getSamplingRate()

        tmpStack = self._getTmpPath("averages.spi")
        imgSet.writeStack(tmpStack)
        orig = os.path.relpath(tmpStack,
                               self._getExtraPath())
        args = "%s %s --apix=%0.3f" % (orig, self._params['relImgsFn'], pixSize)
        self.runJob(Plugin.getProgram('e2proc2d.py'), args,
                    cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createInitialModelStep(self, args):
        """ Run the EMAN program to create the initial model. """
        cleanPattern(self._getExtraPath('initial_models'))
        if self._isHighSym():
            program = Plugin.getProgram('e2initialmodel_hisym.py')
        else:
            program = Plugin.getProgram('e2initialmodel.py')

        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        classes2DSet = self.inputSet.get()
        volumes = self._createSetOfVolumes()
        shrink = self.shrink.get()
        if isinstance(self.inputSet.get(), SetOfClasses2D):
            volumes.setSamplingRate(classes2DSet.getImages().getSamplingRate() * shrink)
        else:
            volumes.setSamplingRate(self.inputSet.get().getSamplingRate() * shrink)
        outputVols = self._getVolumes()
        for k, volFn in enumerate(outputVols):
            vol = Volume()
            vol.setFileName(volFn)
            vol.setObjComment('eman initial model %02d' % (k + 1))
            volumes.append(vol)

        self._defineOutputs(**{outputs.outputVolumes.name: volumes})
        self._defineSourceRelation(self.inputSet, volumes)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Input images: %s" % self.getObjectTag('inputSet'))
            summary.append("Output initial volumes: %s" % self.outputVolumes.getSize())
            if self._isHighSym():
                summary.append("Used e2initialmodel_hisym.py for high symmetry reconstruction.")
        return summary

    # --------------------------- UTILS functions -----------------------------

    def _prepareDefinition(self):
        self._params = {'imgsFn': self._getExtraPath('representatives.hdf'),
                        'relImgsFn': 'representatives.hdf',
                        'numberOfIterations': self.numberOfIterations.get(),
                        'numberOfModels': self.numberOfModels.get(),
                        'shrink': self.shrink.get(),
                        'symmetry': self.symmetry.get(),
                        'threads': self.numberOfThreads.get(),
                        'mpis': self.numberOfMpi.get(),
                        'scratch': Plugin.getVar(EMAN2SCRATCHDIR)}

    def _isHighSym(self):
        return self.symmetry.get() in ["oct", "tet", "icos"]

    def _getVolumes(self):
        if self._isHighSym():
            outputVols = [self._getExtraPath('final.hdf')]
        else:
            outputVols = glob(self._getExtraPath('initial_models/model_??_??.hdf'))
            outputVols.sort()
        return outputVols
