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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import (ProtImportMicrographs, ProtImportParticles, ProtImportVolumes,
                            ProtImportAverages, ProtImportCoordinatesPairs,
                            ProtImportMicrographsTiltPairs)
from pwem import Domain
from pyworkflowtests.protocols import ProtOutputTest
from pwem.objects.data import Pointer
from pyworkflow.utils import magentaStr

from ..protocols import *


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
        # We have two options: pass the SamplingRate or
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
