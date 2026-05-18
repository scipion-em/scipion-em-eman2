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
from enum import Enum

from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam, BooleanParam,
                                        LEVEL_ADVANCED)
from pyworkflow.constants import PROD
from pyworkflow.utils.path import cleanPattern, makePath
from pwem.objects.data import Volume
from pwem.protocols import ProtReconstruct3D

from .. import Plugin
from ..convert import writeSetOfParticles
from ..constants import (RECON_FOURIER, FOURIER_GAUSS2, KEEP_PERCENTAGE,
                         KEEP_STDDEV, KEEP_ABSQUAL)


class outputs(Enum):
    outputVolume = Volume


class EmanProtReconstruct(ProtReconstruct3D):
    """
    Reconstructs 3D cryo-EM volumes from aligned 2D particle images using
    EMAN2 reconstruction methods. The protocol generates a three-dimensional
    density map by combining particle projections with their associated
    orientation information while optionally applying symmetry and CTF-aware
    corrections to improve structural consistency.

    AI Generated:

    Reconstruct 3D Volume (EmanProtReconstruct) — User Manual
        Overview

        The Reconstruct 3D Volume protocol generates a three-dimensional
        cryo-EM map from a set of aligned two-dimensional particle images.
        Each particle contributes information about the structure from a
        specific viewing direction, and the protocol combines these views
        into a coherent volumetric reconstruction. This step is one of the
        central stages in single-particle cryo-EM workflows because it
        transforms image alignment information into an interpretable 3D map.

        The protocol is designed for biological studies ranging from
        exploratory structural analysis to high-resolution reconstruction
        projects. It supports several EMAN2 reconstruction approaches,
        including Fourier-based methods and alternative interpolation
        strategies, allowing users to adapt the reconstruction process to
        different datasets and image qualities.

        Inputs and Reconstruction Workflow

        The protocol requires a set of particles containing projection
        alignment information. These orientations are used to determine how
        each particle contributes to the final three-dimensional structure.
        Accurate angular assignments are critical because reconstruction
        quality depends directly on the consistency and precision of the
        input alignments.

        During execution, particles may undergo additional preprocessing,
        including CTF estimation and phase correction. These operations help
        compensate for microscope-induced contrast distortions and generally
        improve the interpretability of the reconstructed map. In most
        biological workflows, applying proper CTF handling is strongly
        recommended unless the dataset has already been fully processed in a
        compatible EMAN2 environment.

        The protocol then combines all particle projections into a single
        volume while enforcing the selected symmetry and reconstruction
        strategy. The resulting map can subsequently be used for refinement,
        classification, atomic modeling, or structural interpretation.

        Symmetry Considerations

        Symmetry plays an important biological and computational role in
        three-dimensional reconstruction. When a macromolecular complex is
        known to possess rotational or point-group symmetry, imposing the
        correct symmetry improves signal-to-noise ratio and often enhances
        achievable resolution.

        Common examples include cyclic symmetries in ring-shaped assemblies,
        dihedral symmetries in multimeric complexes, or icosahedral symmetry
        in viral capsids. Correct symmetry assignment can substantially
        strengthen weak structural features by averaging equivalent views.

        However, biological caution is required. Applying incorrect symmetry
        may introduce severe artifacts, obscure asymmetric features, or
        produce misleading structural interpretations. If structural
        asymmetry or conformational heterogeneity is suspected, reconstruction
        without symmetry is usually the safest initial approach.

        Reconstruction Methods

        The protocol supports several reconstruction strategies with
        different numerical characteristics and computational behaviors.
        For most biological applications, the Fourier-based reconstructor is
        the recommended starting point because it offers a strong balance
        between reconstruction quality, robustness, and computational
        efficiency.

        Different Fourier insertion modes are available to control how
        particle information is interpolated into Fourier space. Simpler
        interpolation methods are generally faster but may reduce accuracy,
        while smoother Gaussian-based approaches usually provide improved
        reconstruction quality for demanding datasets.

        Alternative reconstruction modes are also available for specialized
        applications, including methods optimized for CTF-aware insertion,
        weighted reconstruction, or experimental interpolation schemes.
        Advanced users may select these approaches when adapting the protocol
        to unusual datasets or testing alternative reconstruction behaviors.

        Iterative Reconstruction and Particle Quality

        The protocol can perform iterative reconstruction procedures in which
        particles are repeatedly incorporated into the volume while image
        normalization and quality assessment are refined. Iterative
        reconstruction can improve consistency across heterogeneous particle
        populations and reduce the impact of poor-quality images.

        Several strategies are available for retaining or excluding particles
        according to reconstruction quality. Users may define particle
        retention using percentages, statistical thresholds, or absolute
        quality criteria. These controls are especially valuable for
        heterogeneous or noisy datasets where low-quality projections may
        degrade the final map.

        In biological practice, retaining all particles is often acceptable
        for highly homogeneous datasets, while more selective filtering may
        improve results for challenging samples containing contaminants,
        preferred orientations, or damaged particles.

        Volume Dimensions and Padding

        The protocol allows control over reconstruction dimensions and image
        padding. Padding increases the effective working area around particles
        during reconstruction and can reduce edge-related artifacts. This is
        particularly useful when particles occupy a large fraction of the
        image box or when interpolation stability becomes important.

        Users may also explicitly define the dimensions of the reconstructed
        volume and the final written map. These controls are useful when
        matching reconstruction sizes across different software packages or
        preparing maps for downstream refinement and visualization workflows.

        In most routine cases, default dimensions derived from the particle
        box size are sufficient. However, advanced users working with large
        complexes or specialized processing pipelines may benefit from
        customized volume dimensions.

        Weighting and Particle Contributions

        By default, the protocol automatically weights particles according to
        their contribution quality. This behavior is particularly beneficial
        when particles originate from class averages or when some projections
        contain stronger signal than others.

        Automatic weighting generally improves reconstruction robustness by
        reducing the influence of weaker or noisier images. Nevertheless,
        users may disable weighting when equal particle contributions are
        desired for methodological testing or highly controlled datasets.

        Outputs and Interpretation

        The protocol produces a reconstructed three-dimensional density map
        representing the consensus structure derived from the input particles.
        The resulting volume preserves the sampling information and metadata
        associated with the original dataset, facilitating integration with
        downstream cryo-EM workflows.

        Biologically, the reconstructed map should always be interpreted in
        the context of dataset quality, alignment accuracy, imposed symmetry,
        and particle heterogeneity. High-resolution features generally
        indicate strong structural consistency, while blurred regions may
        reflect flexibility, conformational variability, or reconstruction
        limitations.

        Practical Recommendations

        For most cryo-EM projects, the Fourier reconstruction method with
        standard settings provides an effective starting point. Applying the
        correct symmetry and maintaining accurate particle alignments are the
        most critical factors influencing reconstruction quality.

        Users working with noisy or heterogeneous datasets should carefully
        evaluate particle retention parameters and consider iterative
        reconstruction approaches to improve map stability. Conversely,
        highly homogeneous datasets often reconstruct successfully with
        relatively conservative filtering and default settings.

        It is generally advisable to inspect reconstructed volumes visually
        after processing to verify that expected structural features are
        preserved and that no obvious symmetry artifacts or reconstruction
        distortions are present.

        Final Perspective

        Three-dimensional reconstruction is the stage at which aligned
        particle images become an interpretable structural model of the
        biological specimen. The quality of the final map depends not only on
        reconstruction parameters but also on the consistency of the input
        particles, the correctness of symmetry assumptions, and the biological
        homogeneity of the sample. Careful selection of reconstruction
        strategies and thoughtful interpretation of the resulting volume are
        essential for obtaining biologically meaningful cryo-EM structures.
    """

    _label = 'reconstruct'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {
            'partSet': 'sets/inputSet.lst',
            'partFlipSet': 'sets/inputSet__ctf_flip.lst',
            'volume': self._getExtraPath('volume.hdf'),
        }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles",
                      pointerCondition='hasAlignmentProj',
                      help='Select the input images from the project.')
        form.addParam('skipctf', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Skip ctf estimation?',
                      help='Use this if you want to skip running e2ctf.py. '
                           'It is not recommended to skip this step unless CTF '
                           'estimation was already done with EMAN2.')
        form.addParam('useE2make3d', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use old e2make3d?',
                      help='Use the traditional e2make3d program instead of '
                           'the new e2make3dpar program.')
        form.addParam('numberOfIterations', IntParam, default=2,
                      condition='useE2make3d',
                      label='Number of iterations:',
                      help='Set the number of iterations. Iterative '
                           'reconstruction improves the overall normalization '
                           'of the 2D images as they are inserted into the '
                           'reconstructed volume, and allows for the '
                           'exclusion of the poorer quality images.')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then the '
                           'model is assumed to have no symmetry. \n'
                           'Choices are: *i, c, d, tet, icos, or oct* \n'
                           'See https://blake.bcm.edu/emanwiki/EMAN2/Symmetry \n'
                           'for a detailed description of symmetry in Eman.')
        line = form.addLine('Padding to Reconstruct: ',
                            expertLevel=LEVEL_ADVANCED,
                            help='Will zero-pad images to the specifed size '
                                 '(x,y) or (x,x) prior to reconstruction. '
                                 'If not specified no padding occurs.')
        line.addParam('padX', IntParam, default=0, label='X ')
        line.addParam('padY', IntParam, default=0, label='Y ')

        line = form.addLine('Dimensions Volume: ',
                            expertLevel=LEVEL_ADVANCED,
                            help='Defines the dimensions (x,y,z) or (x,x,x) '
                                 'of the reconstructed volume. If omitted, '
                                 'implied value based on padded 2D images '
                                 'is used. ')
        line.addParam('dimVolX', IntParam, default=0, label='X')
        line.addParam('dimVolY', IntParam, default=0, label='Y')
        line.addParam('dimVolZ', IntParam, default=0, label='Z')

        line = form.addLine('Dimensions to Write Volume: ',
                            expertLevel=LEVEL_ADVANCED,
                            help='Defines the dimensions (x,y,z) or (x,x,x) '
                                 'of the final volume written to disk, if '
                                 'omitted, size will be based on unpadded '
                                 'input size. ')
        line.addParam('dimWriteVolX', IntParam, default=0, label='X')
        line.addParam('dimWriteVolY', IntParam, default=0, label='Y')
        line.addParam('dimWriteVolZ', IntParam, default=0, label='Z')
        form.addParam('reconstructionMethod', EnumParam,
                      choices=['back_projection', 'fourier', 'fourier_iter',
                               'fouriersimple2D', 'nn4', 'nn4_ctf',
                               'nn4_ctf_rect', 'nn4_ctfw', 'nn4_ctfws', 'nn4_rect',
                               'nnSSNR', 'nnSSNR_ctf', 'real_median', 'wiener_fourier'],
                      label="Reconstruction Method:", default=RECON_FOURIER,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Reconstructor to use. See e2help.py reconstructors '
                           '-v 9. Default is fourier:mode=gauss_2.')
        form.addParam('fourierMode', EnumParam,
                      condition="reconstructionMethod==1 or reconstructionMethod==2 or reconstructionMethod==12",
                      choices=['nearest_neighbor', 'gauss_2', 'gauss_3',
                               'gauss_5', 'gauss_5_slow', 'gypergeom_5',
                               'experimental'],
                      label="Mode to Fourier method:", default=FOURIER_GAUSS2,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Fourier pixel insertion mode. See e2help.py reconstructors '
                           'fourier -v 9. Default mode is gauss_2.')
        form.addParam('keepSense', EnumParam, expertLevel=LEVEL_ADVANCED,
                      choices=['percentage', 'standard deviation',
                               'absolute quality'],
                      label="Sense of keep:", default=KEEP_PERCENTAGE,
                      display=EnumParam.DISPLAY_COMBO,
                      help="If *percentage* is selected, *keep* parameter "
                           "will be interpreted as a percentage. It is the "
                           "default option.\nIf *standard deviation* is "
                           "selected, *keep* parameter will be interpreted "
                           "as a standard deviation coefficient instead of "
                           "as a percentage.\n"
                           "If *absolute quality* is selected, *keep* "
                           "parameter will refer to the absolute quality "
                           "of the class-average, not a local quality "
                           "relative to other similar sized classes.")
        form.addParam('keep', FloatParam, default=1.0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Fraction of slices to keep",
                      help='The fraction of slices to keep, in fraction,'
                           ' based on quality scores (1.0 = use all slices).')
        form.addParam('doNotAutoWt', BooleanParam, default=False,
                      label='Disable automatic weighting?',
                      help='This argument turns automatic weighting off '
                           'causing all images to be weighted by 1. If this '
                           'argument is False images inserted into the '
                           'reconstructed volume are weighted by the number '
                           'of particles that contributed to them (i.e. as '
                           'in class averages), which is extracted from the '
                           'image header.')
        form.addParam('extraParams', StringParam,
                      expertLevel=LEVEL_ADVANCED,
                      default='',
                      label='Additional parameters',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. "
                           "See e2make3dpar.py -h.")

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertImagesStep', needsGPU=False)
        self._insertFunctionStep('reconstructVolumeStep',
                                 self._prepareParams(), needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def convertImagesStep(self):
        partSet = self.inputParticles.get()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)
        if not self.skipctf:
            program = Plugin.getProgram('e2ctf.py')
            acq = partSet.getAcquisition()

            args = " --voltage %3d" % acq.getVoltage()
            args += " --cs %f" % acq.getSphericalAberration()
            args += " --ac %f" % (100 * acq.getAmplitudeContrast())
            if not partSet.isPhaseFlipped():
                args += " --phaseflip"
            args += " --computesf --apix %f " % partSet.getSamplingRate()
            args += " --allparticles --autofit --curdefocusfix --storeparm -v 8"
            args += " --threads=%d" % self.numberOfThreads.get()
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfThreads=1)

        program = Plugin.getProgram('e2buildsets.py')
        args = " --setname=inputSet --allparticles"
        self.runJob(program, args, cwd=self._getExtraPath(), numberOfThreads=1)

    def reconstructVolumeStep(self, args):
        """ Run the EMAN program to reconstruct a volume. """
        cleanPattern(self._getFileName("volume"))
        if self.useE2make3d:
            program = Plugin.getProgram('e2make3d.py')
        else:
            program = Plugin.getProgram('e2make3dpar.py')
        self.runJob(program, args, cwd=self._getExtraPath(), numberOfThreads=1)

    def createOutputStep(self):
        partSet = self.inputParticles.get()
        vol = Volume()
        vol.setFileName(self._getFileName("volume"))
        vol.copyInfo(partSet)
        self._defineOutputs(**{outputs.outputVolume.name: vol})
        self._defineSourceRelation(self.inputParticles, vol)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if not self.useE2make3d and self.reconstructionMethod != RECON_FOURIER:
            errors.append('e2make3dpar.py program can use only Fourier method '
                          'for reconstruction!')

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Input images: %s" % self.getObjectTag('inputParticles'))
            summary.append("Output volume: %s" % self.getObjectTag('outputVolume'))
        return summary

    # --------------------------- UTILS functions -----------------------------

    def _prepareParams(self):
        args = "--input %(imgsFn)s --output %(outputVol)s --sym %(sym)s"

        if self.useE2make3d:
            args += " --iter %(numberOfIterations)d "
            args += " --recon %(reconsMethod)s"
        else:
            args += " --mode %s" % self.getEnumText('fourierMode')
            args += " --threads=%d" % self.numberOfThreads.get()

        if self.extraParams.hasValue():
            args += ' ' + self.extraParams.get()

        reconsMethod = self.getEnumText('reconstructionMethod')
        if (reconsMethod == 'fourier' or reconsMethod == 'fourier_plane' or
                reconsMethod == 'fouriersimple2D' or
                reconsMethod == 'wiener_fourier'):
            reconsMethod = reconsMethod + ':mode=' + self.getEnumText('fourierMode')

        params = {'imgsFn': self._getParticlesStack(),
                  'outputVol': self._getBaseName("volume"),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'sym': self.symmetry.get(),
                  'reconsMethod': reconsMethod
                  }

        args %= params

        if self.padX.get() > 0:
            if self.padY.get() <= 0 or self.padX.get() == self.padY.get():
                args += " --pad %d" % self.padX.get()
            else:
                args += " --pad %d,%d" % (self.padX.get(), self.padY.get())

        if self.dimVolX.get() > 0:
            if ((self.dimVolY.get() <= 0 and self.dimVolZ.get() <= 0) or
                    (self.dimVolY.get() == self.dimVolX.get() and
                     self.dimVolZ.get() == self.dimVolX.get())):
                args += " --padvol %d" % self.dimVolX.get()
            else:
                args += " --padvol %d,%d,%d" % (self.dimVolX.get(),
                                                self.dimVolY.get(),
                                                self.dimVolZ.get())

        if self.dimWriteVolX.get() > 0:
            if ((self.dimWriteVolY.get() <= 0 and self.dimWriteVolZ.get() <= 0) or
                    (self.dimWriteVolY.get() == self.dimWriteVolX.get() and
                     self.dimWriteVolZ.get() == self.dimWriteVolX.get())):
                args += " --outsize %d" % self.dimWriteVolX.get()
            else:
                args += " --outsize %d,%d,%d" % (self.dimWriteVolX.get(),
                                                 self.dimWriteVolY.get(),
                                                 self.dimWriteVolZ.get())

        if self.keepSense == KEEP_STDDEV:
            args += "  --keep %f --keepsig" % self.keep.get()
        elif self.keepSense == KEEP_ABSQUAL:
            args += "  --keep %f --keepabs" % self.keep.get()

        if self.keep.get() != 1.0 and self.keepSense == KEEP_PERCENTAGE:
            args += " --keep %f" % self.keep.get()

        if self.doNotAutoWt:
            args += " --no_wt"

        return args

    def _getBaseName(self, key):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key))

    def _getParticlesStack(self):
        if not self.inputParticles.get().isPhaseFlipped() and not self.skipctf:
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")
