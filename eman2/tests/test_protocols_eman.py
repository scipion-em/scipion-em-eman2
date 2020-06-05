# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

from pyworkflow.tests import (BaseTest, setupTestProject, DataSet)
from pwem.protocols import (ProtImportMicrographs, ProtImportParticles, ProtImportVolumes,
                            ProtImportAverages, ProtImportCoordinatesPairs,
                            ProtImportMicrographsTiltPairs)
from pwem import Domain
from pyworkflowtests.protocols import ProtOutputTest
from pwem.objects.data import Pointer
from pyworkflow.utils import magentaStr

from ..protocols import *
import tomo.protocols


class TestEmanBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')

    @classmethod
    def setData(cls, projectData='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(projectData)
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        if samplingRate is not None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=0, filesPath=pattern,
                                             samplingRate=samplingRate, magnification=magnification,
                                             voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=1, filesPath=pattern,
                                             scannedPixelSize=scannedPixelSize,
                                             voltage=voltage, magnification=magnification,
                                             sphericalAberration=sphericalAberration)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        cls.assertIsNotNone(cls.protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")

        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern, samplingRate=samplingRate,
                                         checkStack=checkStack)
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return cls.protImport

    @classmethod
    def runImportParticlesSqlite(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         importFrom=4,
                                         sqliteFile=pattern, samplingRate=samplingRate)
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes,
                                     filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def runImportAverages(cls, pattern, samplingRate):
        """ Run an Import averages protocol. """
        cls.protImportAvg = cls.newProtocol(ProtImportAverages,
                                            filesPath=pattern,
                                            samplingRate=samplingRate,
                                            checkStack=True)
        cls.launchProtocol(cls.protImportAvg)
        return cls.protImportAvg


class TestEmanInitialModelMda(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.averages = cls.dataset.getFile('averages')
        cls.symmetry = 'd6'
        cls.numberOfIterations = 5
        cls.numberOfModels = 2
        print(magentaStr("\n==> Importing data - class averages:"))
        cls.protImportAvg = cls.runImportAverages(cls.averages, 3.5)

    def test_initialmodel(self):
        print(magentaStr("\n==> Testing eman2 - initial model:"))
        protIniModel = self.newProtocol(EmanProtInitModel,
                                        symmetry=self.symmetry,
                                        numberOfIterations=self.numberOfIterations,
                                        numberOfModels=self.numberOfModels,
                                        numberOfThreads=4)
        protIniModel.inputSet.set(self.protImportAvg.outputAverages)
        self.launchProtocol(protIniModel)
        self.assertIsNotNone(protIniModel.outputVolumes,
                             "There was a problem with eman initial model protocol")


class TestEmanInitialModelGroel(TestEmanInitialModelMda):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('groel')
        cls.averages = cls.dataset.getFile('averages')
        cls.symmetry = 'd7'
        cls.numberOfIterations = 10
        cls.numberOfModels = 10
        print(magentaStr("\n==> Importing data - class averages:"))
        cls.protImportAvg = cls.runImportAverages(cls.averages, 2.1)


class TestEmanInitialModelSGD(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('groel')
        cls.averages = cls.dataset.getFile('averages')
        cls.symmetry = 'd7'
        cls.numberOfIterations = 20
        cls.numberOfModels = 2
        print(magentaStr("\n==> Importing data - class averages:"))
        cls.protImportAvg = cls.runImportAverages(cls.averages, 2.1)

    def test_initialmodel(self):
        print(magentaStr("\n==> Testing eman2 - initial model SGD:"))
        protIniModel = self.newProtocol(EmanProtInitModelSGD,
                                        symmetry=self.symmetry,
                                        numberOfIterations=self.numberOfIterations,
                                        numberOfModels=self.numberOfModels,
                                        numberOfThreads=4)
        protIniModel.inputType.set(0)  # averages
        protIniModel.inputAvg.set(self.protImportAvg.outputAverages)
        self.launchProtocol(protIniModel)
        self.assertIsNotNone(protIniModel.outputVolumes,
                             "There was a problem with eman initial model SGD protocol")


class TestEmanReconstruct(TestEmanBase):
    def test_ReconstructEman(self):
        print(magentaStr("\n==> Importing data - particles:"))
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from scipion (to-reconstruct)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dsRelion.getFile('import/case2/particles.sqlite'),
                                 magnification=10000,
                                 samplingRate=7.08
                                 )
        self.launchProtocol(prot1)

        print(magentaStr("\n==> Testing eman2 - recostruct:"))
        protReconstruct = self.newProtocol(EmanProtReconstruct)
        protReconstruct.inputParticles.set(prot1.outputParticles)
        self.launchProtocol(protReconstruct)
        self.assertIsNotNone(protReconstruct.outputVolume,
                             "There was a problem with eman reconstruction protocol")


class TestEmanRefineEasy(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData('mda')
        print(magentaStr("\n==> Importing data - particles:"))
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        print(magentaStr("\n==> Importing data - volume:"))
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)

    def test_RefineEasyEman(self):
        print(magentaStr("\n==> Testing eman2 - refine easy:"))
        protRefine = self.newProtocol(EmanProtRefine,
                                      symmetry="d6",
                                      speed=6,
                                      numberOfIterations=1)
        protRefine.inputParticles.set(self.protImport.outputParticles)
        protRefine.input3DReference.set(self.protImportVol.outputVolume)
        self.launchProtocol(protRefine)
        self.assertIsNotNone(protRefine.outputVolume,
                             "There was a problem with eman refine easy protocol")


class TestEmanRefine2D(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData('mda')
        print(magentaStr("\n==> Importing data - particles:"))
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)

    def test_Refine2DEman(self):
        print(magentaStr("\n==> Testing eman2 - refine 2d:"))
        protRefine = self.newProtocol(EmanProtRefine2D,
                                      numberOfIterations=2, numberOfClassAvg=5,
                                      classIter=2, nbasisfp=3)
        protRefine.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protRefine)
        self.assertIsNotNone(protRefine.outputClasses,
                             "There was a problem with eman refine2d protocol")


class TestEmanRefine2DBispec(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData('relion_tutorial')
        cls.partsFn = cls.dataset.getFile('import/case2/particles.sqlite')
        print(magentaStr("\n==> Importing data - particles:"))
        cls.protImport = cls.runImportParticlesSqlite(cls.partsFn, 3.5)

    def test_Refine2DBispecEman(self):
        print(magentaStr("\n==> Running eman2 - ctf auto:"))
        protCtf = self.newProtocol(EmanProtCTFAuto,
                                   numberOfThreads=3)
        protCtf.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protCtf)
        self.assertIsNotNone(protCtf.outputParticles_flip_fullRes,
                             "There was a problem with eman ctf auto protocol")

        print(magentaStr("\n==> Testing eman2 - refine 2d bispec:"))
        protRefine = self.newProtocol(EmanProtRefine2DBispec,
                                      inputBispec=protCtf,
                                      numberOfIterations=2, numberOfClassAvg=5,
                                      classIter=2, nbasisfp=5)
        self.launchProtocol(protRefine)
        self.assertIsNotNone(protRefine.outputClasses,
                             "There was a problem with eman refine2d bispec protocol")


class TestEmanTiltValidate(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('eman')
        cls.vol = cls.dataset.getFile('volume')
        cls.micsUFn = cls.dataset.getFile('micU')
        cls.micsTFn = cls.dataset.getFile('micT')
        cls.patternU = cls.dataset.getFile("coords/ip3r10252011-0005_0-2_info.json")
        cls.patternT = cls.dataset.getFile("coords/ip3r10252011-0005_10_info.json")
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.6)

    def test_RefineEman(self):
        print(magentaStr("\n==> Importing data - micrograph pairs:"))
        protImportMicsPairs = self.newProtocol(ProtImportMicrographsTiltPairs,
                                               patternUntilted=self.micsUFn,
                                               patternTilted=self.micsTFn,
                                               samplingRate=1.88, voltage=200,
                                               sphericalAberration=2.0)
        self.launchProtocol(protImportMicsPairs)
        self.assertIsNotNone(protImportMicsPairs.outputMicrographsTiltPair,
                             "There was a problem with the import of mic pairs")

        print(magentaStr("\n==> Importing data - coordinate pairs:"))
        protImportCoords = self.newProtocol(ProtImportCoordinatesPairs,
                                            importFrom=1,  # from eman
                                            patternUntilted=self.patternU,
                                            patternTilted=self.patternT,
                                            boxSize=256)
        protImportCoords.inputMicrographsTiltedPair.set(protImportMicsPairs.outputMicrographsTiltPair)
        self.launchProtocol(protImportCoords)
        self.assertIsNotNone(protImportCoords.outputCoordinatesTiltPair,
                             "There was a problem with the import of coord pairs")

        print(magentaStr("\n==> Running xmipp3 - extract particles:"))
        XmippProtExtractParticlesPairs = Domain.importFromPlugin('xmipp3.protocols',
                                                                 'XmippProtExtractParticlesPairs')
        protExtractPairs = self.newProtocol(XmippProtExtractParticlesPairs,
                                            downFactor=2.0,
                                            boxSize=128,
                                            doInvert=True)

        protExtractPairs.inputCoordinatesTiltedPairs.set(protImportCoords.outputCoordinatesTiltPair)
        self.launchProtocol(protExtractPairs)
        self.assertIsNotNone(protExtractPairs.outputParticlesTiltPair,
                             "There was a problem with particle pair extraction")

        print(magentaStr("\n==> Testing eman2 - tilt validate:"))
        protValidate = self.newProtocol(EmanProtTiltValidate, symmetry="c4",
                                        maxtilt=60.0, delta=2.0, shrink=2,
                                        quaternion=True,
                                        simcmpType=2,  # frc
                                        simcmpParams='maxres=60',
                                        simalignType=7,  # rotate_translate
                                        simralignType=1,  # refine
                                        numberOfThreads=4)
        protValidate.inputTiltPair.set(protExtractPairs.outputParticlesTiltPair)
        protValidate.inputVolume.set(self.protImportVol.outputVolume)
        protValidate._createFilenameTemplates()
        outputAngles = protValidate._getFileName('outputAngles')
        self.launchProtocol(protValidate)
        self.assertIsNotNone(outputAngles, "Missing some output files!")


class TestEmanCtfAuto(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData('relion_tutorial')
        cls.partsFn = cls.dataset.getFile('import/case2/particles.sqlite')
        print(magentaStr("\n==> Importing data - particles:"))
        cls.protImport = cls.runImportParticlesSqlite(cls.partsFn, 3.5)

    def test_CtfAutoEman(self):
        print(magentaStr("\n==> Testing eman2 - ctf auto:"))
        protCtf = self.newProtocol(EmanProtCTFAuto,
                                   numberOfThreads=3)
        protCtf.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protCtf)
        self.assertIsNotNone(protCtf.outputParticles_flip_fullRes,
                             "There was a problem with eman ctf auto protocol")


class TestEmanAutopick(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData('igbmc_gempicker')
        cls.micsFn = cls.dataset.getFile('micrographs/*.mrc')
        cls.avgFn = cls.dataset.getFile('templates/templates_white.stk')
        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.runImportMicrograph(cls.micsFn,
                                                     samplingRate=4.4,
                                                     voltage=120,
                                                     sphericalAberration=2.0,
                                                     scannedPixelSize=None,
                                                     magnification=60000)
        print(magentaStr("\n==> Importing data - class averages:"))
        cls.protImportAvg = cls.runImportAverages(cls.avgFn, 4.4)

    def test_AutopickEman(self):
        print(magentaStr("\n==> Testing eman2 - e2boxer auto:"))
        protPick = self.newProtocol(EmanProtAutopick,
                                    boxerMode=1,  # by_ref
                                    goodRefs=self.protImportAvg.outputAverages,
                                    threshold=12.0,
                                    numberOfThreads=2)
        protPick.inputMicrographs.set(self.protImportMics.outputMicrographs)
        self.launchProtocol(protPick)
        self.assertIsNotNone(protPick.outputCoordinates,
                             "There was a problem with e2boxer auto protocol")

    def test_AutopickSparx(self):
        print(magentaStr("\n==> Testing eman2 - e2boxer gauss:"))
        protPick2 = self.newProtocol(SparxGaussianProtPicking,
                                     boxSize=128,
                                     lowerThreshold=0.004,
                                     higherThreshold=0.1,
                                     gaussWidth=0.525,
                                     useVarImg=False,
                                     doInvert=True)
        protPick2.inputMicrographs.set(self.protImportMics.outputMicrographs)
        self.launchProtocol(protPick2)
        self.assertIsNotNone(protPick2.outputCoordinates,
                             "There was a problem with e2boxer gauss auto protocol")

    def test_AutopickSparxPointer(self):
        print(magentaStr("\n==> Running auto boxsize simulator:"))
        protAutoBoxSize = self.newProtocol(ProtOutputTest,
                                           iBoxSize=64,  # output is twice
                                           objLabel='auto boxsize simulator')
        self.launchProtocol(protAutoBoxSize)

        print(magentaStr("\n==> Testing eman2 - e2boxer gauss:"))
        protPick2 = self.newProtocol(SparxGaussianProtPicking,
                                     lowerThreshold=0.004,
                                     higherThreshold=0.1,
                                     gaussWidth=0.525,
                                     useVarImg=False,
                                     doInvert=True)
        protPick2.inputMicrographs.set(self.protImportMics.outputMicrographs)
        protPick2.boxSize.setPointer(Pointer(protAutoBoxSize,
                                             extended="oBoxSize"))
        self.launchProtocol(protPick2)
        self.assertIsNotNone(protPick2.outputCoordinates,
                             "There was a problem with e2boxer gauss auto protocol")


class TestEmanTomoBase(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def setData(cls, projectData='tomo-em'):
        from tomo.tests import DataSet
        cls.dataset = DataSet.getDataSet(projectData)
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
        cls.inputSetOfSubTomogram = cls.dataset.getFile('subtomo')
        cls.smallTomogram = cls.dataset.getFile('coremask_normcorona.mrc')


class TestEmanTomoExtraction(TestEmanTomoBase):
    """This class check if the protocol to extract subtomograms
    in Eman works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runTomoExtraction(self, tomoSource=0, doInvert=False, doNormalize=False, boxSize=32, downFactor=1):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")

        if tomoSource == 1:
            protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                                  inputTomograms=protImportTomogram.outputTomograms,
                                                  inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                                  tomoSource=tomoSource,
                                                  doInvert=doInvert,
                                                  doNormalize=doNormalize,
                                                  boxSize=boxSize,
                                                  downFactor=downFactor)
        else:
            protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                                  inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                                  tomoSource=tomoSource,
                                                  doInvert=doInvert,
                                                  doNormalize=doNormalize,
                                                  boxSize=boxSize,
                                                  downFactor=downFactor)
        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(protTomoExtraction.outputSetOfSubtomogram,
                             "There was a problem with SetOfSubtomogram output")
        return protTomoExtraction

    def test_extractParticlesSameAsPicking(self):
        protTomoExtraction = self._runTomoExtraction()
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction

    def test_extractParticlesOther(self):
        protTomoExtraction = self._runTomoExtraction(tomoSource=1)
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction

    def test_extractParticlesWithDoInvert(self):
        protTomoExtraction = self._runTomoExtraction(doInvert=True)
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction

    def test_extractParticlesWithDoNormalize(self):
        protTomoExtraction = self._runTomoExtraction(doNormalize=True)
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction

    def test_extractParticlesModifiedDownFactor(self):
        protTomoExtraction = self._runTomoExtraction(downFactor=2)
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction

    def test_extractParticlesModifiedBoxSize(self):
        protTomoExtraction = self._runTomoExtraction(boxSize=64)
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction

    def test_extractParticlesWithAllOptions(self):
        protTomoExtraction = self._runTomoExtraction(boxSize=64, downFactor=2, doNormalize=True, doInvert=True)
        output = getattr(protTomoExtraction, 'outputSetOfSubtomogram', None)
        self.assertTrue(output)
        self.assertTrue(output.hasCoordinates3D())
        self.assertTrue(output.getCoordinates3D().getObjValue())
        return protTomoExtraction


class TestEmanTomoInitialModel(TestEmanTomoBase):
    """This class check if the protocol to extract particles
    in Relion works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                              inputTomograms=protImportTomogram.outputTomograms,
                                              inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                              tomoSource=0,
                                              doInvert=False,
                                              doNormalize=False,
                                              boxSize=32)

        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(protTomoExtraction.outputSetOfSubtomogram,
                             "There was a problem with SetOfSubtomogram output")
        return protTomoExtraction

    def _performFinalValidation(self, protInitialModel):
        averageSubTomogram = protInitialModel.averageSubTomogram
        self.assertEqual(os.path.basename(averageSubTomogram.getFirstItem().getFileName()), "output.hdf")
        self.assertEqual(averageSubTomogram.getFirstItem().getSamplingRate(), 20.0)
        self.assertEqual(averageSubTomogram.getSamplingRate(), 20.0)

        setOfSubTomograms = protInitialModel.outputParticles
        self.assertEqual(setOfSubTomograms.getSize(), 5)
        self.assertEqual(setOfSubTomograms.getCoordinates3D().getObjValue().getSize(), 5)
        self.assertEqual(setOfSubTomograms.getSamplingRate(), 20.0)

        for subTomogram in setOfSubTomograms:
            self.assertEqual(subTomogram.getSamplingRate(), 20)
            self.assertTrue(hasattr(subTomogram, "coverage"))
            self.assertTrue(hasattr(subTomogram, "score"))
            matrix = subTomogram.getTransform().getMatrix()
            self.assertEqual(matrix.shape, (4, 4))

        summary = protInitialModel.summary()
        self.assertIsNotNone(summary)
        self.assertEqual(len(summary), 2)
        self.assertEqual(summary[0], "Particles: 5")
        self.assertTrue(summary[1].startswith("Reference file used:"))

    def _runTomoSubtomogramInitialModelWithSubtomo(self):
        protTomoExtraction = self._runPreviousProtocols()

        particles = protTomoExtraction.outputSetOfSubtomogram

        protInitialModel = self.newProtocol(EmanProtTomoInitialModel,
                                            particles=particles,
                                            reference=protTomoExtraction,
                                            symmetry="c1",
                                            gaussFilter=-1.5,
                                            filterto=0.03,
                                            fourier=False,
                                            batchSize=20,
                                            learningRate=1,
                                            numberOfIterations=2,
                                            numberOfBatches=1,
                                            shrink=4,
                                            applySim=False)
        protInitialModel.reference.setExtended("outputSetOfSubtomogram.1")

        self.launchProtocol(protInitialModel)

        self.assertIsNotNone(protInitialModel.averageSubTomogram,
                             "There was a problem with subTomograms output")
        self.assertIsNotNone(protInitialModel.outputParticles,
                             "There was a problem with particles output")

        return protInitialModel

    def _runTomoSubtomogramInitialModelWithVolume(self):
        protTomoExtraction = self._runPreviousProtocols()

        particles = protTomoExtraction.outputSetOfSubtomogram

        self.dataset = DataSet.getDataSet('eman')
        self.vol = self.dataset.getFile('volume')
        self.protImportVol = self.runImportVolumes(self.vol, 3.5)

        self.assertIsNotNone(self.protImportVol.outputVolume,
                             "There was a problem with SetOfSubtomogram output")

        protInitialModel = self.newProtocol(EmanProtTomoInitialModel,
                                            particles=particles,
                                            reference=self.protImportVol.outputVolume,
                                            symmetry="c1",
                                            gaussFilter=-1.5,
                                            filterto=0.03,
                                            fourier=False,
                                            batchSize=20,
                                            learningRate=1,
                                            numberOfIterations=2,
                                            numberOfBatches=1,
                                            shrink=4,
                                            applySim=False)

        self.launchProtocol(protInitialModel)

        self.assertIsNotNone(protInitialModel.averageSubTomogram,
                             "There was a problem with subTomograms output")
        self.assertIsNotNone(protInitialModel.outputParticles,
                             "There was a problem with particles output")

        return protInitialModel

    def test_initialModelOutputWithSubtomo(self):
        protInitialModel = self._runTomoSubtomogramInitialModelWithSubtomo()

        self._performFinalValidation(protInitialModel)

        return protInitialModel

    def test_initialModelOutputWithVolume(self):
        protInitialModel = self._runTomoSubtomogramInitialModelWithVolume()

        self._performFinalValidation(protInitialModel)

        return protInitialModel


class TestEmanTomoSubtomogramRefinement(TestEmanTomoBase):
    """This class check if the protocol Subtomogram refinement works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        doInvert = False
        doNormalize = False
        boxSize = 32
        protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                              inputTomograms=protImportTomogram.outputTomograms,
                                              inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                              tomoSource=0,
                                              doInvert=doInvert,
                                              doNormalize=doNormalize,
                                              boxSize=boxSize)

        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(protTomoExtraction.outputSetOfSubtomogram,
                             "There was a problem with SetOfSubtomogram output")

        return protTomoExtraction

    def _performFinalValidation(self, protTomoSubtomogramRefinement):
        outputSetOfSubTomograms = protTomoSubtomogramRefinement.outputParticles
        averageSubTomogram = protTomoSubtomogramRefinement.averageSubTomogram

        self.assertTrue("threed" in averageSubTomogram.getFirstItem().getFileName())
        self.assertEqual(averageSubTomogram.getFirstItem().getSamplingRate(), 5.0)
        self.assertEqual(averageSubTomogram.getSamplingRate(), 5.0)

        self.assertEqual(outputSetOfSubTomograms.getDimensions(), (32, 32, 32))
        self.assertEqual(outputSetOfSubTomograms.getSize(), 5)
        self.assertEqual(outputSetOfSubTomograms.getCoordinates3D().getObjValue().getSize(), 5)

        for subTomogram in outputSetOfSubTomograms:
            self.assertEqual(subTomogram.getSamplingRate(), 5)
            self.assertTrue(hasattr(subTomogram, "coverage"))
            self.assertTrue(hasattr(subTomogram, "score"))
            matrix = subTomogram.getTransform().getMatrix()
            self.assertEqual(matrix.shape, (4, 4))

    def _runTomoSubtomogramRefinementWithSubtomo(self, niter=2, mass=500.0, threads=1, pkeep=0.8, goldstandard=-1,
                                                 goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0):
        protTomoExtraction = self._runPreviousProtocols()
        protTomoRefinement = self.newProtocol(EmanProtTomoRefinement,
                                              inputSetOfSubTomogram=protTomoExtraction.outputSetOfSubtomogram,
                                              inputRef=protTomoExtraction,
                                              niter=niter,
                                              mass=mass,
                                              threads=threads,
                                              pkeep=pkeep,
                                              goldstandard=goldstandard,
                                              goldcontinue=goldcontinue,
                                              sym=sym,
                                              localfilter=localfilter,
                                              maxtilt=maxtilt)
        protTomoRefinement.inputRef.setExtended("outputSetOfSubtomogram.1")

        self.launchProtocol(protTomoRefinement)

        self.assertIsNotNone(protTomoRefinement.averageSubTomogram,
                             "There was a problem with subTomograms output")
        self.assertIsNotNone(protTomoRefinement.outputParticles,
                             "There was a problem with particles output")

        return protTomoRefinement

    def _runTomoSubtomogramRefinementWithVolume(self, niter=2, mass=500.0, threads=1, pkeep=0.8, goldstandard=-1,
                                                goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0):
        protTomoExtraction = self._runPreviousProtocols()

        particles = protTomoExtraction.outputSetOfSubtomogram

        self.dataset = DataSet.getDataSet('eman')
        self.vol = self.dataset.getFile('volume')
        self.protImportVol = self.runImportVolumes(self.vol, 3.5)

        self.assertIsNotNone(self.protImportVol.outputVolume,
                             "There was a problem with SetOfSubtomogram output")

        protTomoRefinement = self.newProtocol(EmanProtTomoRefinement,
                                              inputSetOfSubTomogram=particles,
                                              inputRef=self.protImportVol.outputVolume,
                                              niter=niter,
                                              mass=mass,
                                              threads=threads,
                                              pkeep=pkeep,
                                              goldstandard=goldstandard,
                                              goldcontinue=goldcontinue,
                                              sym=sym,
                                              localfilter=localfilter,
                                              maxtilt=maxtilt)

        self.launchProtocol(protTomoRefinement)

        self.assertIsNotNone(protTomoRefinement.averageSubTomogram,
                             "There was a problem with subTomograms output")
        self.assertIsNotNone(protTomoRefinement.outputParticles,
                             "There was a problem with particles output")

        return protTomoRefinement

    def test_defaultSubTomogramRefinementWithSubTomo(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinementWithSubtomo()
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement

    def test_defaultSubTomogramRefinementWithVolume(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinementWithVolume()
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement


class TestEmanTomoTempMatch(TestEmanTomoBase):
    """This class check if the program Template Matching
    from Eman works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runTomoTempMatch(self):
        protImportTomogramBig = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                                 filesPath=self.tomogram,
                                                 samplingRate=5)

        protImportTomogramSmall = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                                   filesPath=self.smallTomogram,
                                                   samplingRate=5)

        self.launchProtocol(protImportTomogramBig)

        self.launchProtocol(protImportTomogramSmall)

        self.assertSetSize(protImportTomogramBig.outputTomograms, size=1,
                           msg="There was a problem with tomogram output")

        self.assertSetSize(protImportTomogramSmall.outputTomograms, size=1, msg="There was a problem with tomogram "
                                                                                "output")

        self.dataset = DataSet.getDataSet('eman')
        self.vol = self.dataset.getFile('volume')
        self.protImportVol = self.runImportVolumes(self.vol, 3.5)

        self.assertIsNotNone(self.protImportVol.outputVolume,
                             "There was a problem with SetOfSubtomogram output")

        protTomoTempMatchBig = self.newProtocol(EmanProtTomoTempMatch,
                                                inputSet=protImportTomogramBig.outputTomograms,
                                                ref=self.protImportVol.outputVolume,
                                                boxSize=128,
                                                sym="c1")

        protTomoTempMatchSmall = self.newProtocol(EmanProtTomoTempMatch,
                                                  inputSet=protImportTomogramSmall.outputTomograms,
                                                  ref=self.protImportVol.outputVolume,
                                                  boxSize=128,
                                                  sym="d1")

        self.launchProtocol(protTomoTempMatchBig)
        self.launchProtocol(protTomoTempMatchSmall)
        return protTomoTempMatchBig, protTomoTempMatchSmall

    def test_TempMatch(self):
        protTomoTempMatch = self._runTomoTempMatch()

        outputCoordsBig = protTomoTempMatch[0].output3DCoordinates
        self.assertEqual(outputCoordsBig.getSize(), 19)
        self.assertEqual(outputCoordsBig.getBoxSize(), 128)
        self.assertEqual(outputCoordsBig.getSamplingRate(), 5)

        outputCoordsSmall = protTomoTempMatch[1].output3DCoordinates
        self.assertEqual(outputCoordsSmall.getSize(), 2)
        self.assertEqual(outputCoordsSmall.getBoxSize(), 128)
        self.assertEqual(outputCoordsSmall.getSamplingRate(), 5)

        return protTomoTempMatch


class TestEmanTomoReconstruction(TestEmanTomoBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _getProtImportTs(self):
        return self.newProtocol(
            tomo.protocols.ProtImportTs,
            filesPath=TestEmanTomoReconstruction.dataset.getPath(),
            filesPattern='tomo{TS}.hdf',
            minAngle=-55,
            maxAngle=60,
            stepAngle=2,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=1.35,
            doseInitial=0,
            dosePerFrame=0.3)

    def _getProtReconstruct(self, protImportTs):
        return self.newProtocol(
            EmanProtTomoReconstruction,
            tiltSeries=protImportTs.outputTiltSeries,
            tiltStep=2.0,
            niter='1,1,1,1',
            bxsz=64)

    def _runPreviousProtocols(self):
        protImportTs = self._getProtImportTs()
        self.launchProtocol(protImportTs)
        self.assertIsNotNone(protImportTs.outputTiltSeries, "Output tilt series not found")

        protReconstruct = self._getProtReconstruct(protImportTs)
        self.launchProtocol(protReconstruct)

        return protReconstruct

    def _validateOutput(self, protReconstruct):
        tomograms = list(protReconstruct.tomograms)
        # 1 tomogram per input file
        self.assertEqual(len(tomograms), 1)
        for tomogram in tomograms:
            self.assertEqual(tomogram.getSamplingRate(), 1.35)
            self.assertEqual(tomogram.getDimensions(), (1024, 1024, 256))

    def test_protocol(self):
        protTomoExtraction = self._runPreviousProtocols()
        self._validateOutput(protTomoExtraction)
        self.assertTrue(protTomoExtraction.summary())
        self.assertTrue(protTomoExtraction.methods())
