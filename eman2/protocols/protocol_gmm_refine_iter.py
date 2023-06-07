# **************************************************************************
# *
# * Authors:
# *
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

from enum import Enum
from os.path import basename, abspath, join

from pwem import ALIGN_PROJ, ALIGN_NONE
from pyworkflow.object import String
from pyworkflow.utils import Message
from pyworkflow.utils.path import createLink, makePath
from pyworkflow.constants import BETA
from pyworkflow.protocol.params import (PointerParam, StringParam, FloatParam)
from pwem.protocols import EMProtocol
from pwem.objects.data import Volume, AtomStruct
from .. import Plugin
import logging
logger = logging.getLogger(__name__)


class outputs(Enum):
    volume = Volume


class EmanProtGmmRefineIter(EMProtocol):
    """
    This protocol wraps *e2gmm_refine_iter.py* EMAN2 program.

    It will tstart from the initial orientation assignment and run five iterations of
    refinement using GMMs as references.

    See more details in:
    https://blake.bcm.edu/emanwiki/EMAN2/e2gmm_refine
    """

    _label = 'GMM-based global refinement'
    _devStatus = BETA
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------
    # DEV TUTORIALS
    # https://scipion-em.github.io/docs/release-3.0.0/docs/developer/tutorials/dev-tutorials.html

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.refFile = None
        self.resultFile = None

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label='Particles',
                      important=True)
        form.addParam('initRef', PointerParam,
                      pointerClass='Volume',
                      label="Reference volume",
                      important=True)
        form.addParam('pdbObj', PointerParam,
                      pointerClass=AtomStruct,
                      allowsNull=True,
                      label='Initial points (opt.)')
        form.addParam('startres', FloatParam,
                      default=4,
                      label='Starting resolution (Å)')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then '
                           'the model is assumed to have no symmetry. \n'
                           'Choices are: *i, c, d, tet, icos, or oct* \n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        # TODO 1 (Muyuan): finish the form with the required params. See https://docs.google.com/presentation/d/1sACaNZFgH0qWeXE6BLUWEDW3cjYTS4kbojrKvvRp78s/present?slide=id.g9f7805ad80_0_538

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.runGmmRefineIterStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.refFile = self._getExtraPath(basename(self.initRef.get().getFileName()))

    def convertInputStep(self):
        # Convert the reference to HDF if necessary
        program = Plugin.getProgram('e2proc3d.py')
        inRef = self.initRef.get()
        inFile = inRef.getFileName()
        sRate = inRef.getSamplingRate()
        if inFile.endswith('.hdf'):
            createLink(inFile, self.refFile)
        else:
            args = '%s %s --apix %.2f ' % (abspath(inFile), self.refFile, sRate)
            self.runJob(program, args, cwd=self._getExtraPath())

        # Generate a star file using the input particles
        starFile = self._getExtraPath('inParticles.star')
        inParticles = self.inParticles.get()
        self.particles2Star(inParticles, starFile)

        # Use EMAN to generate a .lst file with the location of all particles and their initial orientation assignment
        program = Plugin.getProgram('e2convertrelion.py')
        acquisition = inParticles.getAcquisition()
        particlesDir = 'particles'
        makePath(self._getExtraPath(particlesDir))
        args = [
            f'{starFile}',
            f'--voltage {acquisition.getVoltage()}',
            f'--cs {acquisition.getSphericalAberration():.2f}',
            f'--apix {inParticles.getSamplingRate():.3f}',
            f'--amp {acquisition.getAmplitudeContrast():-2f}',
            f'--onestack {join(particlesDir, "particles_all.hdf")}',
            f'--sym {self.symmetry.get()}',
            '--make3d'
        ]
        self.runJob(program, args, cwd=self._getExtraPath())

        # TODO 3 (Jorge) [DONE]: convert Scipion set to star file. Same as here --> https://raw.githubusercontent.com/scipion-em/scipion-em-cryosparc2/devel/cryosparc2/convert/convert.py#:~:text=logger.info(%27Trying%20to%20generate%20the%20star%20file%20with%20Relion,fileName%2C%20**args)%0A%20%20%20%20%20%20%20%20logger.info(%27The%20star%20file%20was%20generate%20successfully%20...%27)
        # TODO 4 (Muyuan / Jorge) [DONE]: call e2convertrelion.py to get the Eman format (check relion SPA plugin)

    def runGmmRefineIterStep(self):
        # TODO 2 (Muyuan): finish the arguments generation for the calling command to EMAN (method _getGmmRefineArgs)
        program = Plugin.getProgram('e2gmm_refine_iter.py')
        args = self._getGmmRefineArgs()
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        # TODO 5 (Jorge): Convert the resulting gmm_00/threed_xx.hdf into mrc

        # Create output object
        vol = Volume()
        vol.copyInfo(self.initRef.get())
        vol.setFileName(self.resultFile)
        # Set FSC data for the viewer
        self.fscUnmasked = String()
        self.fscMasked = String()
        self.fscMaskedTight = String()

        # Define outputs and relations
        self._defineOutputs(**{outputs.volume.name: vol})
        self._defineSourceRelation(self.initRef, vol)
        if self.pdbObj:
            self._defineSourceRelation(self.pdbObj, vol)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        pass

    def _summary(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _getGmmRefineArgs(self):
        atomStruct = self.pdbObj.get()
        args = [
            f'{self.refFile}',
            f'--startres {self.startres.get():.2f}',
            f'--sym {self.symmetry.get()}'
        ]
        if atomStruct:
            args.append(f'--initpts {atomStruct.getFileName()}')
        return ' '.join(args)

    def particles2Star(self, imgSet, fileName):
        logger.info('Trying to generate the star file with Relion convert...')
        args = {'outputDir': self._getExtraPath(),
                'fillMagnification': True,
                'fillRandomSubset': True}
        from relion import convert
        alignType = ALIGN_PROJ if imgSet.hasAlignmentProj() else ALIGN_NONE
        args['alignType'] = alignType
        args['incompatibleExtensions'] = ['hdf', 'stk']
        convert.writeSetOfParticles(imgSet, fileName, **args)
        logger.info('The star file was generate successfully ...')
