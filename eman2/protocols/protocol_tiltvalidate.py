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

import os

from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, EnumParam, StringParam,
                                        BooleanParam, IntParam, LEVEL_ADVANCED)
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pwem.protocols import ProtAnalysis3D

from .. import Plugin
from ..constants import (WIKI_URL, SIMCMP_CHOICES, CMP_CCC, SIMALIGN_CHOICES,
                         ALN_ROTATE_TRANSLATE, RALN_NONE, CMP_DOT)
from ..convert import writeSetOfParticles


class EmanProtTiltValidate(ProtAnalysis3D):
    """
    Validates cryo-EM tilt pair data against a reconstructed 3D volume
    using EMAN2 tilt geometry analysis methods derived from the approach
    described by Rosenthal and Henderson. The protocol evaluates whether
    experimentally observed tilted and untilted particle projections are
    geometrically consistent with the provided reconstruction, helping
    assess the reliability and angular correctness of the map.

    AI Generated:

    Tilt Validate (EmanProtTiltValidate) - User Manual
        Overview

        The Tilt Validate protocol performs angular validation of a 3D
        reconstruction by comparing experimental tilt pair particles
        against projections derived from an input volume. In cryo-EM,
        tilt pair analysis is an important strategy for independently
        confirming that a reconstruction possesses the correct handedness,
        angular assignment consistency, and overall orientation accuracy.

        The protocol is especially valuable in workflows where the
        reliability of orientation determination must be carefully
        verified, such as during initial model validation, publication
        preparation, or the analysis of challenging datasets with low
        signal-to-noise ratios. By comparing tilted and untilted particle
        images against projections from the reconstructed map, the method
        provides an experimental consistency check that is independent of
        many assumptions used during refinement.

        Biological and Experimental Context

        In single-particle cryo-EM, tilt pair experiments involve imaging
        the same particles at two known specimen tilts. Because the
        relative orientation between the tilted and untilted images is
        experimentally constrained, these data provide a powerful way to
        verify whether the reconstructed volume and its assigned particle
        orientations are physically meaningful.

        From a biological perspective, tilt validation increases
        confidence that structural features observed in the reconstruction
        correspond to real molecular organization rather than alignment
        artifacts or reconstruction bias. This is particularly important
        for asymmetric complexes, flexible assemblies, or datasets
        reconstructed near the limits of achievable resolution.

        Inputs and General Workflow

        The protocol requires two main inputs: a reconstructed 3D volume
        and a set of experimentally paired tilted and untilted particles.
        The volume acts as the structural reference against which the
        particle orientations are evaluated.

        During processing, the protocol generates projections from the
        volume and compares them with the experimental particle images
        across a defined angular search space. The agreement between
        predicted and observed tilt relationships is then analyzed to
        determine whether the reconstruction is geometrically consistent
        with the experimental tilt data.

        For best results, the input particles should correspond closely
        to the particles used during reconstruction. Significant
        heterogeneity, poor particle quality, or inaccurate particle
        pairing may reduce the reliability of the validation.

        Symmetry Considerations

        The protocol allows the user to define the symmetry of the input
        reconstruction. Correct symmetry assignment is biologically
        critical because symmetry directly influences the interpretation
        of angular relationships between projections.

        Symmetric complexes such as viral capsids or highly ordered
        oligomers benefit from the use of the appropriate symmetry group,
        since it improves the consistency of orientation comparisons.
        However, assigning incorrect symmetry can artificially improve or
        distort validation statistics. When uncertainty exists regarding
        the true symmetry of the specimen, conservative choices are often
        preferable.

        Tilt Geometry and Angular Constraints

        The maximum tilt angle parameter determines the allowable angular
        range during validation. This should reflect the experimental
        acquisition conditions as closely as possible. Excessively broad
        ranges may introduce ambiguous matches, whereas overly narrow
        limits can exclude valid solutions.

        The angular projection step controls the sampling precision used
        during orientation comparison. Smaller angular increments improve
        precision but increase computational cost. For exploratory
        analyses, moderate angular sampling is often sufficient, while
        publication-quality validation may benefit from finer searches.

        Quaternion-Based Orientation Analysis

        The protocol optionally supports quaternion-based angular
        calculations. Quaternion representations can improve numerical
        stability and provide smoother orientation handling in complex
        rotational searches. This option is particularly useful in cases
        where conventional angular parameterizations may introduce
        ambiguities or discontinuities.

        In most routine workflows, standard orientation handling is
        adequate, but advanced users working with difficult angular
        distributions or highly heterogeneous datasets may benefit from
        quaternion analysis.

        Particle Shrinking and Computational Efficiency

        For large datasets or preliminary analyses, particles may be
        computationally reduced in size before similarity evaluation.
        This can significantly accelerate validation while preserving the
        overall angular trends necessary for interpretation.

        Biological users should recognize that aggressive shrinking may
        reduce sensitivity to fine structural details. Therefore, coarse
        analyses may be appropriate during parameter exploration, whereas
        final validation is generally more reliable when performed with
        minimally reduced data.

        Similarity Metrics and Alignment Strategies

        The protocol supports multiple comparison metrics and alignment
        strategies for evaluating particle similarity. These options
        influence how projections are matched against experimental
        particles and can affect robustness under different imaging
        conditions.

        Cross-correlation based approaches are commonly suitable for most
        cryo-EM datasets and provide a good balance between sensitivity
        and computational efficiency. More advanced refinement or
        alignment methods may improve performance in difficult datasets
        containing noise, partial occupancy, or substantial orientation
        uncertainty.

        In practice, default comparison settings are often sufficient for
        routine validation. Advanced optimization is generally reserved
        for specialized analyses or problematic datasets.

        Contour Plot Analysis

        The protocol optionally produces contour plots similar to those
        traditionally used in tilt pair validation studies. These plots
        provide a visual representation of angular consistency across the
        dataset and can help identify systematic orientation deviations
        or reconstruction problems.

        From a biological interpretation standpoint, well-defined contour
        distributions generally indicate good agreement between the
        reconstruction and experimental tilt geometry. Broad or irregular
        distributions may suggest alignment instability, preferred
        orientation artifacts, incorrect symmetry assignment, or particle
        heterogeneity.

        Outputs and Interpretation

        The protocol produces validation measurements describing the
        agreement between the input volume and the experimental tilt pair
        particles. These results can be used to assess angular accuracy,
        reconstruction consistency, and overall map reliability.

        When contour visualization is enabled, additional graphical
        outputs are generated to facilitate interpretation of angular
        distributions and tilt agreement quality.

        Validation results should always be interpreted together with
        complementary cryo-EM quality indicators such as FSC curves,
        angular distribution analysis, local resolution estimation, and
        visual inspection of reconstructed features.

        Practical Recommendations

        For most biological applications, it is advisable to begin with
        conservative angular sampling and standard similarity metrics.
        Once stable behavior is confirmed, finer angular searches and
        advanced alignment refinements may be introduced if necessary.

        Correct pairing of tilted and untilted particles is essential for
        meaningful validation. Datasets with substantial heterogeneity,
        particle damage, or strong preferred orientation effects may
        produce weaker validation statistics even when the reconstruction
        itself is accurate.

        In difficult cases, careful inspection of contour distributions,
        symmetry assignments, and tilt angle constraints often provides
        the most useful diagnostic information.

        Final Perspective

        Tilt pair validation provides an experimentally grounded
        assessment of cryo-EM reconstruction reliability. Rather than
        relying solely on internal refinement statistics, this approach
        tests whether the reconstructed volume is consistent with known
        experimental geometry. For many cryo-EM studies, especially those
        involving novel structures or challenging datasets, tilt
        validation represents an important step toward ensuring
        biologically trustworthy structural interpretation.
    """

    _label = 'tilt validate'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtAnalysis3D.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        myDict = {
            'untiltPartSet': 'sets/untilted_ptcls.lst',
            'tiltPartSet': 'sets/tilted_ptcls.lst',
            'outputAngles': self._getExtraPath('TiltValidate_01/perparticletilts.json'),
            'outputContourPlot': self._getExtraPath('TiltValidate_01/contour.hdf')
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",
                      help='Select the input volume that will be validated.')
        form.addParam('inputTiltPair', PointerParam,
                      label="Input tilt pair particles",
                      pointerClass='ParticlesTiltPair',
                      help='Select the input set of tilt pair particles.')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then '
                           'the model is assumed to have no symmetry. \n'
                           'Choices are: *i, c, d, tet, icos, or oct* \n'
                           'See https://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('maxtilt', FloatParam, default=180.0,
                      label='Max tilt angle',
                      help='Maximum tilt angle permitted when finding tilt '
                           'distances.')
        form.addParam('quaternion', BooleanParam, default=False,
                      label='Use quaternions', expertLevel=LEVEL_ADVANCED,
                      help='Use quaternions for tilt distance computation')
        form.addParam('delta', FloatParam, default=5.0,
                      label='Projection step (deg.)',
                      help='Angular step size for alignment')
        form.addParam('shrink', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Shrink particles',
                      help='Optionally shrink the input particles by an integer '
                           'amount prior to computing similarity scores. '
                           'For speed purposes.')
        form.addParam('doContourPlot', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do contour plot?',
                      help='Also make a contour plot similar to fig. 6 '
                           'in Henderson paper')
        form.addParam('tiltRange', IntParam, default=15,
                      expertLevel=LEVEL_ADVANCED,
                      condition='doContourPlot',
                      label='Tilt range',
                      help='The angular tilt range to search')
        form.addParam('verbose', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Verbose level',
                      help='Verbose level from 0 to 9. ')

        form.addSection(label='Similarity matrix')
        form.addParam('paramsMsg', LabelParam, default=True,
                      label='These parameters are for advanced users only!\n',
                      help='For help please address to EMAN2 %s or run:\n'
                           '*e2help.py cmp -v 2* or\n'
                           '*e2help.py aligners -v 2*' % WIKI_URL)
        line = form.addLine('simcmp: ',
                            help='The name of a cmp to be used in comparing '
                                 'the aligned images (default=ccc)')
        line.addParam('simcmpType', EnumParam,
                      choices=list(SIMCMP_CHOICES.values()),
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simcmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('First stage aligner')
        line = group.addLine('simalign: ')
        line.addParam('simalignType', EnumParam,
                      choices=list(SIMALIGN_CHOICES.values()),
                      label='type', default=ALN_ROTATE_TRANSLATE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simalignParams', StringParam,
                      default='', label='params')
        line = group.addLine('simaligncmp: ')
        line.addParam('simaligncmpType', EnumParam,
                      choices=list(SIMCMP_CHOICES.values()),
                      label='type', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simaligncmpParams', StringParam,
                      default='', label='params')

        group = form.addGroup('Second stage aligner')
        line = group.addLine('simralign: ')
        line.addParam('simralignType', EnumParam,
                      choices=['None', 'refine', 'refine_3d',
                               'refine_3d_grid', 'refinecg'],
                      label='type', default=RALN_NONE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simralignParams', StringParam,
                      default='', label='params')
        line = group.addLine('simraligncmp: ')
        line.addParam('simraligncmpType', EnumParam,
                      choices=list(SIMCMP_CHOICES.values()),
                      label='type', default=CMP_DOT,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('simraligncmpParams', StringParam,
                      default='', label='params')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertImagesStep', needsGPU=False)
        args = self._prepareParams()
        self._insertFunctionStep('runValidateStep', args, needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def convertImagesStep(self):
        part = self.inputTiltPair.get()
        partUnt = part.getUntilted()
        partTilt = part.getTilted()
        storePath = self._getExtraPath("particles")
        pwutils.makePath(storePath)
        self.info("Converting input particle set..")

        for partSet, suffix in zip([partUnt, partTilt],
                                   ['_untilted_ptcls', '_tilted_ptcls']):
            partAlign = partSet.getAlignment()
            writeSetOfParticles(partSet, storePath,
                                alignType=partAlign, suffix=suffix)

            setName = suffix.split('_')[1]
            program = Plugin.getProgram('e2buildsets.py')
            args = " particles/*%s.hdf --setname=%s" % (
                suffix, setName)
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

    def runValidateStep(self, args):
        program = Plugin.getProgram('e2tiltvalidate.py')
        self.runJob(program, args, cwd=self._getExtraPath(), numberOfThreads=1)

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        self._validateDim(self.inputTiltPair.get().getUntilted(),
                          self.inputVolume.get(), errors,
                          'Input tilt pair particles', 'Input volume')

        return errors

    def _summary(self):
        summary = list()
        summary.append("Max. tilt angle: *%0.2f*" % self.maxtilt.get())
        summary.append("Projection step: *%d deg.*" % self.delta.get())
        summary.append("Symmetry: *%s*" % self.symmetry.get())

        return summary

    # --------------------------- UTILS functions -----------------------------

    def _prepareParams(self):
        args = " --untiltdata=%(untilt)s --tiltdata=%(tilt)s --volume=%(volume)s"
        args += " --maxtiltangle=%(maxtilt)f --sym=%(sym)s --delta=%(delta)f"
        args += " --verbose=%(verb)d"

        if self.shrink.get() != 1:
            args += " --shrink=%d" % self.shrink.get()
        if self.quaternion:
            args += " --quaternion"
        if self.doContourPlot:
            args += " --docontourplot --tiltrange %d" % self.tiltRange.get()
        if self.numberOfThreads.get() > 1:
            args += " --parallel=thread:%d" % self.numberOfThreads.get()

        params = {'untilt': self._getFileName("untiltPartSet"),
                  'tilt': self._getFileName("tiltPartSet"),
                  'volume': os.path.relpath(self.inputVolume.get().getFileName(),
                                            self._getExtraPath()).replace(":mrc", ""),
                  'maxtilt': self.maxtilt.get(),
                  'sym': self.symmetry.get(),
                  'delta': self.delta.get(),
                  'verb': self.verbose.get()
                  }

        args %= params

        for param in ['simcmp', 'simalign', 'simaligncmp',
                      'simralign', 'simraligncmp']:
            args += self._getSimmxOpts(param)

        return args

    def _getSimmxOpts(self, option):
        optionType = self.getEnumText(option + 'Type')
        optionParams = getattr(self, option + 'Params').get()

        if optionType == 'None':
            return ''
        if optionParams != '':
            argStr = ' --%s=%s:%s' % (option, optionType, optionParams)
        else:
            argStr = ' --%s=%s' % (option, optionType)

        return argStr
