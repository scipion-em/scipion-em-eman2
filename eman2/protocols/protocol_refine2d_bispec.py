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
import re
from glob import glob
from enum import Enum

from pwem.objects import SetOfClasses2D
from pwem.protocols import ProtClassify2D
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam, BooleanParam,
                                        LabelParam)
from pyworkflow.utils import createLink, cleanPath

from .. import Plugin
from ..constants import *


class outputs(Enum):
    outputClasses = SetOfClasses2D


class EmanProtRefine2DBispec(ProtClassify2D):
    """
    Performs reference-free 2D classification of cryo-EM particle images
    using bispectrum-based invariants within the EMAN2 framework. The
    protocol generates representative class averages from heterogeneous
    particle populations, allowing researchers to evaluate data quality,
    identify structural variability, and obtain biologically meaningful
    views of macromolecular complexes without requiring an initial
    structural reference. These class averages are commonly used during
    exploratory stages of single-particle analysis and can support the
    generation of initial 3D models or the interpretation of conformational
    diversity. More info:
    https://blake.bcm.edu/emanwiki/EMAN2

    AI Generated:

    Refine 2D Bispectrum Classification (EmanProtRefine2DBispec) —
        User Manual

        Overview

        The Refine 2D Bispectrum Classification protocol performs
        reference-free classification of particle images using the
        EMAN2 bispectrum-based refinement strategy. Its primary goal
        is to organize heterogeneous particle datasets into groups of
        structurally similar views while avoiding dependence on an
        external reference model. This makes the protocol especially
        valuable during the early stages of cryo-EM processing, when
        the quality, homogeneity, and orientation distribution of the
        dataset are still being evaluated.

        Unlike reference-based refinement approaches, this protocol
        relies on rotational and translational invariants derived from
        particle bispectra. These mathematical descriptors help compare
        particles independently of their initial orientation, allowing
        robust classification even when the dataset contains strongly
        misaligned images or multiple structural states. From a
        biological perspective, this enables the identification of
        dominant particle views, contaminants, damaged particles, or
        flexible conformations without introducing model bias.

        Inputs and Biological Context

        The protocol requires particles that have already undergone
        CTF estimation and bispectrum preprocessing through a compatible
        EMAN2 workflow. These particles are typically extracted from
        cryo-EM micrographs and may represent multiple orientations,
        conformational states, or compositional assemblies.

        In practical biological workflows, the resulting class averages
        serve several important purposes. They help determine whether
        the dataset contains recognizable structural features, whether
        particles are sufficiently homogeneous for high-resolution
        refinement, and whether additional cleaning or classification
        steps are required. Researchers frequently inspect the classes
        visually to identify preferred orientations, flexibility,
        aggregation, or the presence of contaminants.

        Since the classification is reference free, the protocol is
        particularly suitable for novel macromolecular complexes,
        poorly characterized assemblies, or datasets where avoiding
        reference bias is essential.

        Number of Classes and Dataset Balance

        One of the most biologically important parameters is the number
        of class averages to generate. This choice determines how the
        structural variability of the dataset will be represented.

        Using too few classes may merge distinct particle views or
        conformations into overly broad averages, potentially masking
        important biological heterogeneity. Conversely, using too many
        classes can fragment the dataset excessively, producing noisy
        averages with insufficient particle support.

        In routine cryo-EM practice, the optimal number of classes
        depends on dataset size and structural complexity. Large and
        heterogeneous datasets often benefit from a larger number of
        classes, while smaller or highly homogeneous datasets usually
        require fewer classes to maintain strong signal quality.

        Bispectrum and Invariant-Based Classification

        The defining characteristic of this protocol is its use of
        bispectrum-derived invariants. These descriptors reduce the
        sensitivity of classification to rotational and translational
        differences between particles. Biologically, this is useful
        because particles in cryo-EM datasets are rarely perfectly
        aligned during the initial stages of processing.

        By emphasizing invariant image features, the protocol can more
        effectively identify structurally related particles even in
        noisy datasets. This improves the robustness of exploratory
        classification and can reveal underlying structural organization
        before precise alignment procedures are introduced.

        However, users should remember that invariant-based approaches
        prioritize classification robustness rather than final alignment
        precision. The resulting classes are intended primarily for
        structural interpretation and dataset assessment rather than
        direct high-resolution refinement.

        Multivariate Analysis and Basis Selection

        The protocol uses multivariate statistical analysis to project
        particle images into a lower-dimensional feature space. The
        number of basis vectors controls how much structural information
        is retained during this representation.

        From a biological perspective, increasing the number of basis
        vectors may improve sensitivity to subtle conformational
        differences or fine structural details. However, excessively
        large values can also increase sensitivity to noise, especially
        in low-quality datasets.

        For most routine cryo-EM datasets, moderate values provide a
        good balance between discrimination power and robustness.
        Highly heterogeneous or exceptionally large datasets may benefit
        from more detailed representations.

        Alignment and Sorting of Class Averages

        The protocol optionally aligns and sorts the resulting class
        averages according to their mutual similarity. This operation
        improves visual interpretability and facilitates downstream
        analysis by organizing related views together.

        Biologically, sorted classes can reveal smooth angular coverage
        of the particle or expose missing orientations and preferred
        views. This information is often critical when evaluating
        whether the dataset is suitable for isotropic 3D reconstruction.

        Alignment of class averages also improves the interpretability
        of flexible assemblies by emphasizing shared structural regions.
        Nevertheless, users should interpret highly variable or weakly
        populated classes cautiously, since apparent structural features
        may sometimes reflect noise or alignment instability.

        Centering and Averaging Strategies

        The protocol provides multiple centering and averaging options
        that influence the appearance and stability of the final class
        averages. Proper centering is especially important for elongated,
        asymmetric, or flexible particles where automatic centering may
        otherwise drift toward solvent regions or peripheral density.

        Different averaging strategies can also affect how structural
        variability is represented. Conservative averaging approaches
        tend to preserve only strongly reproducible features, whereas
        more permissive approaches may retain weaker or more flexible
        densities.

        In biological practice, users often compare several averaging
        configurations when working with difficult datasets to determine
        which settings best preserve meaningful structural information.

        Outputs and Interpretation

        The protocol produces a set of 2D class averages together with
        particle assignments for each class. These outputs provide a
        compact visual summary of the structural content of the dataset.

        Strong class averages with recognizable secondary-structure
        features generally indicate good particle quality and alignment
        potential. Conversely, diffuse or poorly resolved classes may
        suggest excessive heterogeneity, inaccurate particle picking,
        ice contamination, or insufficient signal-to-noise ratio.

        The classes generated by this protocol are not intended as the
        final step of refinement. Instead, they serve as an intermediate
        biological interpretation tool that guides subsequent processing
        decisions such as particle cleaning, initial model generation,
        or selection of homogeneous subsets for high-resolution
        reconstruction.

        Practical Recommendations

        In routine workflows, it is often advisable to begin with a
        moderate number of classes and inspect the results visually.
        If the classes appear overly broad or contain mixed views,
        increasing the class count may reveal additional heterogeneity.
        If classes become excessively noisy, reducing the number of
        classes or increasing particle counts may improve stability.

        Datasets containing strong preferred orientations may produce
        many highly similar classes, whereas highly flexible complexes
        may generate diffuse or structurally inconsistent averages.
        In these situations, additional particle cleaning or focused
        classification strategies may be required.

        When evaluating the results, users should prioritize biological
        interpretability rather than numerical metrics alone. The most
        useful class averages are those that reveal reproducible
        structural features and meaningful conformational organization.

        Final Perspective

        Reference-free 2D classification is one of the most important
        exploratory stages in cryo-EM single-particle analysis because
        it provides the first direct visualization of structural
        reproducibility within the dataset. The bispectrum-based
        strategy implemented in this protocol offers a robust approach
        for organizing heterogeneous particle populations while reducing
        sensitivity to initial particle orientation.

        For most cryo-EM researchers, the quality and interpretability
        of the resulting class averages strongly influence downstream
        decisions regarding dataset quality, structural heterogeneity,
        and refinement strategy. Careful selection of classification
        parameters and thoughtful biological interpretation are therefore
        essential for obtaining reliable and meaningful results.
    """
    _label = 'refine 2D bispec'
    _devStatus = PROD
    _possibleOutputs = outputs

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
        self._iterRegex = re.compile(r'classes_(\d{2})')

    # --------------------------- DEFINE param functions ----------------------
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
        form.addParam('nbasisfp', IntParam, default=8,
                      label='Number of MSA vectors to use',
                      help='Number of MSa basis vectors to use when '
                           'classifying particles.')
        form.addParam('doNormProj', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Normalize projected vectors?',
                      help='Normalizes each projected vector into the MSA '
                           'subspace. Note that this is different from normalizing '
                           'the input images since the subspace is not expected to '
                           'fully span the image')
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
                      choices=list(CENTER_CHOICES.values()),
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
                      choices=list(AVG_CHOICES.values()),
                      label='Class averager: ',
                      default=AVG_CTF_WEIGHT_AUTOFILT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='The averager used to generated class-averages')

        line = form.addLine('classnormproc: ',
                            help='Normalization applied during class-averaging')
        line.addParam('classnormprocType', EnumParam,
                      choices=list(NORM_CHOICES.values()),
                      label='',
                      default=PROC_NORMALIZE_EDGEMEAN,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classnormprocParams', StringParam,
                      default='', label='params')

        line = form.addLine('classcmp: ')
        line.addParam('classcmpType', EnumParam,
                      choices=list(SIMCMP_CHOICES.values()),
                      label='', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classcmpParams', StringParam,
                      default='', label='params',
                      help='The name of a cmp to be used in class averaging '
                           '(default=ccc)')

        group = form.addGroup('First stage aligner (clsavg)')
        line = group.addLine('classalign: ')
        line.addParam('classalignType', EnumParam,
                      choices=list(SIMALIGN_CHOICES.values()),
                      label='', default=ALN_ROTATE_TRANSLATE_TREE,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classalignParams', StringParam,
                      default='flip=1', label='params')
        line = group.addLine('classaligncmp: ')
        line.addParam('classaligncmpType', EnumParam,
                      choices=list(SIMCMP_CHOICES.values()),
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
                      choices=list(SIMCMP_CHOICES.values()),
                      label='', default=CMP_CCC,
                      display=EnumParam.DISPLAY_COMBO)
        line.addParam('classraligncmpParams', StringParam,
                      default='', label='params')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._createIterTemplates(currRun=self._getRun())
        self._insertFunctionStep('createLinksStep', needsGPU=False)
        args = self._prepareParams()
        self._insertFunctionStep('refineStep', args, needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
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
        program = Plugin.getProgram('e2refine2d_bispec.py')
        # mpi and threads are handled by EMAN itself
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        partSet = self._getInputParticles()
        classes2D = self._createSetOfClasses2D(partSet)
        self._fillClassesFromIter(classes2D, self._lastIter())

        self._defineOutputs(**{outputs.outputClasses.name: classes2D})
        self._defineSourceRelation(partSet, classes2D)

    # --------------------------- INFO functions ------------------------------
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
        methods += "into %d classes using e2refine2d_bispec.py " % \
                   self.numberOfClassAvg
        return [methods]

    # --------------------------- UTILS functions -----------------------------
    def _prepareParams(self):
        args1 = " --input=%s" % self._getParticlesStack()
        args2 = self._commonParams()
        args = args1 + args2

        return args

    def _commonParams(self):
        args = " --ncls=%(ncls)d --nbasisfp=%(nbasisfp)d"
        args += " --classkeep=%(classKeep)f --classiter=%(classiter)d "
        args += " --classaverager=%s" % self.getEnumText('classAveragerType')

        if self.doNormProj:
            args += " --normproj"

        if self.alignSort:
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
                  'nbasisfp': self.nbasisfp.get(),
                  'classKeep': self.classKeep.get(),
                  'classiter': self.classIter.get(),
                  'threads': self.numberOfThreads.get(),
                  'mpis': self.numberOfMpi.get(),
                  'scratch': Plugin.getVar(EMAN2SCRATCHDIR)
                  }
        args %= params

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        return args

    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key, **args))

    def _getParticlesStack(self):
        protType = self._inputProt().type.get()
        if protType == HIRES:
            return "sets/" + os.path.basename(self._getFileName("partSetFlipLp5"))
        elif protType == MIDRES:
            return "sets/" + os.path.basename(self._getFileName("partSetFlipLp7"))
        else:
            return "sets/" + os.path.basename(self._getFileName("partSetFlipLp12"))

    def _iterTextFile(self, iterN):
        with open(self._getFileName('results', iter=iterN)) as f:
            for line in f:
                if '#' not in line and line.strip():
                    yield [float(x) for x in line.split()]

    def _getRun(self):
        return 0

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
        If the file doesn't exist, it will be created by
        converting from this iteration data.star file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)

        if clean:
            cleanPath(data_classes)

        if not os.path.exists(data_classes):
            clsSet = SetOfClasses2D(filename=data_classes)
            clsSet.setImages(self._getInputParticles())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _getInputParticles(self):
        prot = self._inputProt()
        protType = prot.type.get()
        if protType == HIRES:
            output = getattr(prot, "outputParticles_flip_lp5")
        elif protType == MIDRES:
            output = getattr(prot, "outputParticles_flip_lp7")
        else:
            output = getattr(prot, "outputParticles_flip_lp12")

        return output

    def _fillClassesFromIter(self, clsSet, iterN):
        self._execEmanProcess(iterN)
        params = {'orderBy': ['_micId', 'id'],
                  'direction': 'ASC'}
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=self._iterTextFile(iterN),
                             iterParams=params)

    def _execEmanProcess(self, iterN):
        runN = self._getRun()
        clsFn = self._getFileName("cls", run=runN, iter=iterN)
        classesFn = self._getFileName("classes", run=runN, iter=iterN)

        proc = Plugin.createEmanProcess(args='read %s %s %s %s 2d'
                                             % (self._getParticlesStack(), clsFn, classesFn,
                                                self._getBaseName('results', iter=iterN)),
                                        direc=self._getExtraPath())
        proc.wait()

        self._classesInfo = {}  # store classes info, indexed by class id
        for classId in range(self.numberOfClassAvg.get()):
            self._classesInfo[classId + 1] = (classId + 1,
                                              self._getExtraPath(classesFn))

    def _getOptsString(self, option):
        optionType = self.getEnumText(option + 'Type')
        optionParams = getattr(self, option + 'Params').get()

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
            _, fn = self._classesInfo[classId]
            item.getRepresentative().setLocation(classId, fn)

    def _inputProt(self):
        return self.inputBispec.get()
