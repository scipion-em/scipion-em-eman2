# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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

import os

from pyworkflow.protocol.params import (IntParam, FloatParam,
                                        EnumParam, PointerParam,
                                        StringParam, USE_GPU,
                                        GPU_LIST, BooleanParam)
from pyworkflow.utils import makePath, createLink
from pyworkflow.constants import PROD
from pwem.protocols import ProtParticlePickingAuto

from .. import Plugin
from ..convert import readSetOfCoordinates, convertReferences
from ..constants import AUTO_CONVNET, AUTO_GAUSS


class EmanProtAutopick(ProtParticlePickingAuto):
    """
    Provides automated particle picking for single particle analysis workflows
    using the EMAN2 Boxer framework. The protocol supports several automated
    detection strategies, including reference-based searching, neural-network
    classification, local alignment approaches, and simple reference-free
    Gaussian picking for rapid exploratory analyses.

    AI Generated:

    EMAN Automatic Particle Picking (EmanProtAutopick) - User Manual
        Overview

        The EMAN automatic particle picking protocol is designed to identify
        particle coordinates directly from cryo-EM micrographs with minimal
        manual intervention. In single particle analysis workflows, automated
        picking is essential for processing large datasets efficiently while
        maintaining consistency across thousands of micrographs.

        The protocol offers several picking strategies adapted to different
        biological situations and data qualities. These methods range from
        simple reference-free approaches suitable for highly contrasted samples
        to advanced neural-network-based detection systems capable of handling
        heterogeneous and noisy datasets.

        For most cryo-EM projects, automated picking represents the transition
        from exploratory manual inspection to large-scale data production. The
        quality of the selected particles strongly influences downstream
        classification, refinement, and final reconstruction quality.

        Inputs and General Workflow

        The protocol requires a set of input micrographs together with particle
        size and box size estimates. Depending on the selected strategy, the
        workflow may also require reference averages or a previously trained
        particle picking model obtained from an earlier interactive boxing
        session.

        During execution, the protocol processes each micrograph individually
        and generates particle coordinates automatically. The resulting
        coordinates can later be inspected, cleaned, and used for particle
        extraction and downstream reconstruction workflows.

        In practical biological applications, users commonly perform an initial
        manual or semi-automated picking session to establish representative
        particles and then apply this automated protocol to process the complete
        dataset consistently.

        Choice of Picking Strategy

        The protocol provides multiple particle detection approaches, each with
        distinct biological advantages and limitations.

        The local search strategy performs reference-guided matching using
        strongly downsampled particles and references. This method is often
        useful when particles exhibit moderate heterogeneity but still maintain
        recognizable structural features. It generally produces fewer false
        positives than simpler correlation-based approaches, although it may be
        computationally slower.

        The reference-based strategy relies on high-quality reference images
        representing the expected particle projections. This approach performs
        well when reliable references are available and the particles display
        limited structural variability. It is especially useful for symmetric
        or structurally stable complexes. However, biologically important rare
        conformations may be missed if they are poorly represented in the
        references.

        The neural-network strategy is the most advanced and generally the most
        accurate option when properly trained. It is particularly effective for
        difficult datasets containing contamination, low contrast, preferred
        orientations, or substantial heterogeneity. By combining multiple neural
        discrimination stages, the protocol attempts to distinguish true
        particles from both background noise and non-particle contaminants.

        The Gaussian strategy provides a rapid reference-free alternative for
        simple datasets with strong contrast and highly homogeneous particles.
        It may perform adequately for large symmetric particles such as viruses
        or ribosomes, but it is usually less reliable for small, flexible, or
        heterogeneous molecular assemblies.

        Particle Size and Box Size Considerations

        Accurate particle sizing is one of the most biologically important
        aspects of automated picking. The particle size parameter should reflect
        the approximate largest visible dimension of the molecular complex in
        the micrographs. If the value is too small, peripheral structural
        information may be excluded. If it is too large, contamination and
        background noise may interfere with particle discrimination.

        The box size determines the extraction region surrounding each detected
        particle. A biologically appropriate box should contain the complete
        particle together with sufficient surrounding solvent to support later
        alignment and classification procedures.

        Users should visually verify the resulting picks carefully, particularly
        when processing datasets with flexible conformations, aggregation, or
        crowded particle distributions.

        Neural-Network-Based Picking

        The neural-network mode is especially valuable for challenging cryo-EM
        datasets where classical template matching methods may fail. In these
        situations, training data generated from a previous supervised boxing
        session can significantly improve particle detection accuracy.

        The protocol allows reuse of previously trained particle selection
        models together with curated examples of valid particles, contaminants,
        and background regions. This enables biologically consistent particle
        selection across multiple datasets or acquisition sessions.

        GPU acceleration can optionally be used in neural-network mode to
        accelerate particle detection. This becomes increasingly important for
        large cryo-EM projects involving thousands of micrographs.

        Threshold Selection and Particle Quality

        Detection thresholds strongly influence the balance between sensitivity
        and specificity. Lower thresholds generally increase particle recovery
        but may introduce more false positives. Higher thresholds improve purity
        at the risk of discarding rare or low-contrast particles.

        Biological interpretation should guide threshold optimization. In many
        projects, retaining a broader particle population during initial picking
        is preferable because downstream two-dimensional classification can
        remove contaminants later. However, excessively permissive thresholds
        may overwhelm classification procedures with noise and artifacts.

        Outputs and Their Interpretation

        The protocol produces particle coordinate sets associated with the input
        micrographs. These coordinates define the particle centers that will be
        used during extraction and subsequent cryo-EM analysis stages.

        The biological quality of the output should always be validated through
        visual inspection. Automated methods can mistakenly identify ice,
        carbon edges, aggregation artifacts, or contamination as particles,
        especially in low-contrast datasets.

        Careful evaluation of coordinate distributions and extracted particle
        images is therefore essential before proceeding to classification and
        refinement.

        Practical Recommendations

        In routine cryo-EM workflows, it is often advisable to begin with the
        neural-network strategy if representative training data are available.
        This approach generally provides the best balance between particle
        recovery and contamination rejection.

        For simple, high-contrast datasets, the Gaussian or reference-based
        methods may provide sufficiently accurate results with lower
        computational cost. Users working with heterogeneous or flexible
        complexes should inspect particle diversity carefully to avoid losing
        biologically meaningful conformational states.

        Thresholds should be adjusted iteratively while visually inspecting
        representative micrographs. It is generally safer to tolerate moderate
        levels of false positives during initial picking than to exclude rare
        structural states irreversibly.

        Final Perspective

        Automated particle picking is a biologically critical selection step
        that determines which molecular observations contribute to the final
        cryo-EM reconstruction. Appropriate choice of picking strategy,
        realistic particle sizing, careful threshold tuning, and continuous
        visual validation are essential for generating reliable datasets suitable
        for high-quality structural interpretation.
    """
    _label = 'boxer auto'
    _devStatus = PROD

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {'goodRefsFn': self._getExtraPath('info/boxrefs.hdf'),
                  'badRefsFn': self._getExtraPath('info/boxrefsbad.hdf'),
                  'bgRefsFn': self._getExtraPath('info/boxrefsbg.hdf'),
                  'nnetFn': self._getExtraPath('nnet_pickptcls.hdf'),
                  'nnetClFn': self._getExtraPath('nnet_classify.hdf'),
                  'trainoutFn': self._getExtraPath('trainout_pickptcl.hdf'),
                  'trainoutClFn': self._getExtraPath('trainout_classify.hdf')
                  }
        self._updateFilenamesDict(myDict)

    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addHidden(USE_GPU, BooleanParam, default=False,
                       label="Use GPU?",
                       help="Set to Yes if you want to run Neural Net "
                            "boxer on GPU. Default is CPU.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       label="Choose GPU ID",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on.\n"
                            "Eman boxer can use only one GPU.")
        form.addParam('boxSize', IntParam, default=128,
                      allowsPointers=True,
                      label='Box size (px)',
                      help="Box size in pixels. See https://eman2.org/BoxSize")
        form.addParam('particleSize', IntParam, default=100,
                      label='Particle size (px)',
                      help="Longest axis of particle in pixels (diameter, "
                           "not radius).")
        form.addParam('boxerMode', EnumParam,
                      choices=['local search', 'by ref', 'neural net', 'gauss'],
                      label="Autopicker mode:", default=AUTO_CONVNET,
                      display=EnumParam.DISPLAY_COMBO,
                      help="Choose autopicker mode:\n\n"
                           " _local search_ - Heavily downsamples the "
                           "particles and references, and actually performs "
                           "a 2-D alignment of each putative particle to each "
                           "reference to identify the best particles in the "
                           "image. In theory this should produce fewer "
                           "false positives. Reference requirements are "
                           "similar to By Ref, though a smaller number of "
                           "references may be fine.\n"
                           " _by ref_ - This is a classic reference based "
                           "particle picker. To use it, you need to have "
                           "high quality good references in all possible "
                           "3-D orientations. The algorithm will do in-plane "
                           "rotation of the references and cross-correlate "
                           "to look for peaks. It is recommended that you "
                           "use projections of a 3-D map (low resolution) "
                           "as references.\n"
                           " _neural net_ - This is the most accurate boxer "
                           "by far, both in terms of false positives and "
                           "false negatives, when trained properly. It is "
                           "based on a pair of neural networks, one to "
                           "discriminate between putative particles and "
                           "the background, and a second to discriminate "
                           "between real particles and contamination or "
                           "other high contrast non-particles, which is "
                           "why it has two thresholds.\n"
                           " _gauss_ - This is a simple and fast "
                           "reference-free picker, which provides a simple "
                           "solution for easy particle picking cases, "
                           "where the particles have good contrast and are "
                           "monodisperse. This may work well for things "
                           "like viruses or ribosomes. It is very fast "
                           "and since it requires no references, it's "
                           "easy to try. It likely won't work well for "
                           "most projects. It was ported from the old "
                           "boxer program by a volunteer (Vadim Kotov).")
        form.addParam('threshold', FloatParam, default=5.0,
                      label='Threshold',
                      condition='boxerMode!=%d' % AUTO_GAUSS)
        form.addParam('threshold2', FloatParam, default=-5.0,
                      condition='boxerMode==%d' % AUTO_CONVNET,
                      label='Threshold2')
        form.addParam('gaussLow', FloatParam, default=1.,
                      condition='boxerMode==%d' % AUTO_GAUSS,
                      label='Threshold low')
        form.addParam('gaussHigh', FloatParam, default=2.,
                      condition='boxerMode==%d' % AUTO_GAUSS,
                      label='Threshold high')
        form.addParam('gaussWidth', FloatParam, default=1.,
                      condition='boxerMode==%d' % AUTO_GAUSS,
                      label='Gaussian width')

        form.addSection('References')
        form.addParam('boxerProt', PointerParam,
                      pointerClass='EmanProtBoxing',
                      condition='boxerMode==%d' % AUTO_CONVNET,
                      label='Previous e2boxer protocol',
                      help='Provide previously executed e2boxer protocol '
                           'that has all 3 types of references and '
                           'pre-trained neural network.')
        form.addParam('goodRefs', PointerParam,
                      pointerClass='SetOfAverages',
                      condition='boxerMode<%d' % AUTO_CONVNET,
                      label="Good references",
                      help="Good particle references.")

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertInitialSteps(self):
        self._createFilenameTemplates()
        initId = self._insertFunctionStep('convertInputStep', needsGPU=False)
        return [initId]

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        goodRefs = self.goodRefs.get() if self.goodRefs.hasValue() else None
        boxerProt = self.boxerProt.get() if self.boxerProt.hasValue() else None
        storePath = self._getExtraPath("info")
        makePath(storePath)

        if goodRefs is not None:
            convertReferences(goodRefs, self._getFileName('goodRefsFn'))

        if boxerProt is not None:
            boxerProt._createFilenameTemplates()
            keys = ['goodRefsFn', 'badRefsFn', 'bgRefsFn',
                    'nnetFn', 'nnetClFn',
                    'trainoutFn', 'trainoutClFn']

            for fn in keys:
                if os.path.exists(boxerProt._getFileName(fn)):
                    createLink(boxerProt._getFileName(fn),
                               self._getFileName(fn))

    def _pickMicrograph(self, mic, *args):
        micFile = os.path.relpath(mic.getFileName(), self.getCoordsDir())
        params = " --apix=%f --no_ctf" % self.inputMicrographs.get().getSamplingRate()
        params += " --boxsize=%d" % self.boxSize.get()
        params += " --ptclsize=%d" % self.particleSize.get()
        params += " --threads=%d" % self.numberOfThreads.get()

        modes = ['auto_local', 'auto_ref', 'auto_convnet', 'auto_gauss']
        params += " --autopick=%s" % modes[self.boxerMode.get()]

        if self.boxerMode.get() == AUTO_GAUSS:
            params += ":gauss_width=%0.3f:thr_low=%0.3f:thr_high=%0.3f:boxsize=%d" % (
                self.gaussWidth.get(), self.gaussLow.get(),
                self.gaussHigh.get(), self.boxSize.get())
        else:
            params += ":threshold=%0.2f" % self.threshold.get()

        if self.boxerMode.get() == AUTO_CONVNET:
            params += ":threshold2=%0.2f" % self.threshold2.get()

            if self.useGpu:
                params += " --device=gpu%s" % self.gpuList.get().strip()
            else:
                params += " --device=cpu"

        params += ' %s' % micFile
        program = Plugin.getProgram('e2boxer.py')

        self.runJob(program, params, cwd=self.getCoordsDir())

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.useGpu and (self.boxerMode.get() != AUTO_CONVNET):
            errors.append("You can use GPU only for neural net picker!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def getCoordsDir(self):
        return self._getExtraPath()

    def getFiles(self):
        return (self.inputMicrographs.get().getFiles() |
                ProtParticlePickingAuto.getFiles(self))

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        coordSet.setBoxSize(self.boxSize.get())
        readSetOfCoordinates(workingDir, micList, coordSet)
