# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from glob import glob

from pyworkflow.utils.path import cleanPattern
from pyworkflow.protocol.params import (PointerParam, TextParam, IntParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.em.data import SetOfClasses2D, Volume

import eman2


class EmanProtInitModel(ProtInitialVolume):
    """
    This protocol wraps *e2initialmodel.py* Eman2 program.

    It will take a set of class-averages/projections and build a set
    of 3-D models suitable for use as initial models in single
    particle reconstruction. The output set is theoretically sorted
    in order of quality (best one is numbered 1), though it's best
    to look at the other answers as well.

    See more details in:
    http://blake.bcm.edu/emanwiki/EMAN2/Programs/e2initialmodel
    """

    _label = 'initial model'

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      # pointerCondition='hasRepresentatives',
                      label="Input averages", important=True,
                      help='Select the your class averages to build your '
                           '3D model.\nYou can select SetOfAverages or '
                           'SetOfClasses2D as input.')
        form.addParam('symmetry', TextParam, default='c1',
                      label='Symmetry group',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
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
        if eman2.Plugin.isNewVersion():
            form.addParam('randOrient', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label='Use random orientations?',
                          help='Instead of seeding with a random volume, '
                               'seeds by randomizing input orientations')
        form.addParallelSection(threads=8, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._prepareDefinition()
        self._insertFunctionStep('createStackImgsStep')
        self._insertInitialModelStep()
        self._insertFunctionStep('createOutputStep')

    def _insertInitialModelStep(self):
        args = '--input=%(relImgsFn)s --sym=%(symmetry)s'
        if self.shrink > 1:
            args += ' --shrink=%(shrink)d'
        if not self._isHighSym():
            args += ' --tries=%(numberOfModels)d --iter=%(numberOfIterations)d'
            if eman2.Plugin.isNewVersion() and self.randOrient:
                args += ' --randorient'
            if self.numberOfMpi > 1:
                args += ' --parallel=mpi:%(mpis)d:%(scratch)s'
            else:
                args += ' --parallel=thread:%(threads)d'
        else:
            args += ' --threads=%(threads)d'

        self._insertFunctionStep('createInitialModelStep', args % self._params)

    # --------------------------- STEPS functions -----------------------------
    def createStackImgsStep(self):

        imgsFn = self._params['imgsFn']
        if isinstance(self.inputSet.get(), SetOfClasses2D):
            imgSet = self._createSetOfParticles("_averages")
            for i, cls in enumerate(self.inputSet.get()):
                img = cls.getRepresentative()
                img.setSamplingRate(cls.getSamplingRate())
                img.setObjId(i + 1)
                imgSet.append(img)
        else:
            imgSet = self.inputSet.get()
        imgSet.writeStack(imgsFn)

    def createInitialModelStep(self, args):
        """ Run the EMAN program to create the initial model. """
        cleanPattern(self._getExtraPath('initial_models'))
        if self._isHighSym():
            program = eman2.Plugin.getProgram('e2initialmodel_hisym.py')
        else:
            program = eman2.Plugin.getProgram('e2initialmodel.py')

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

        self._defineOutputs(outputVolumes=volumes)
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
        imgsFn = self._getPath('representatives.stk')

        self._params = {'imgsFn': imgsFn,
                        'relImgsFn': os.path.relpath(imgsFn, self._getExtraPath()),
                        'numberOfIterations': self.numberOfIterations.get(),
                        'numberOfModels': self.numberOfModels.get(),
                        'shrink': self.shrink.get(),
                        'symmetry': self.symmetry.get(),
                        'threads': self.numberOfThreads.get(),
                        'mpis': self.numberOfMpi.get(),
                        'scratch': eman2.SCRATCHDIR
                        }

    def _isHighSym(self):
        return (eman2.Plugin.getActiveVersion() != '2.11' and
                self.symmetry.get() in ["oct", "tet", "icos"])

    def _getVolumes(self):
        if self._isHighSym():
            outputVols = [self._getExtraPath('final.hdf')]
        else:
            outputVols = glob(self._getExtraPath('initial_models/model_??_??.hdf'))
            outputVols.sort()
        return outputVols
