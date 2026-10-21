# **************************************************************************
# *
# * EMAN2 auto-picking in streaming: a micrograph e2boxer fails on is left
# * out and picked again on Continue, a micrograph is published whole, and
# * the run fails only when every micrograph failed.
# *
# **************************************************************************

import os
import tempfile
import unittest
from collections import OrderedDict
from unittest import mock

from eman2.protocols.protocol_autopick_boxer import EmanProtAutopick
import eman2.protocols.protocol_autopick_boxer as module


class _Mic:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return "mic_%03d" % self._objId

    def getFileName(self):
        return "/data/mic_%03d.mrc" % self._objId


class _Coords:
    def __init__(self):
        self.items = []
        self.boxSize = None

    def append(self, item):
        self.items.append(item)

    def setBoxSize(self, boxSize):
        self.boxSize = boxSize

    def getSize(self):
        return len(self.items)


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _Picking(EmanProtAutopick):
    def __init__(self, root):
        self._root = root
        os.makedirs(os.path.join(root, 'extra', 'DONE'), exist_ok=True)
        self.micDict = OrderedDict()
        self.errors = []
        self.boxSize = _Value(64)

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    def error(self, message):
        self.errors.append(message)


class _Base(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmp.cleanup()


class TestAFailedMicrographIsLeftOut(_Base):
    def test_AnE2boxerFailureDoesNotFailTheStep(self):
        protocol = _Picking(self._tmp.name)
        protocol.runJob = mock.Mock(side_effect=RuntimeError("e2boxer crashed"))
        protocol.inputMicrographs = mock.Mock()
        protocol.inputMicrographs.get.return_value.getSamplingRate.return_value = 1.0
        protocol.particleSize = _Value(50)
        protocol.numberOfThreads = _Value(1)
        protocol.boxerMode = _Value(0)
        protocol.threshold = _Value(5.0)

        with mock.patch.object(module.Plugin, 'getProgram', return_value='e2boxer.py'):
            protocol._pickMicrograph(_Mic(3))

        self.assertEqual({3}, protocol._failedPickingMicIds())
        self.assertTrue(protocol.errors)

    def test_AFailedMicrographIsNotRead(self):
        protocol = _Picking(self._tmp.name)
        protocol._recordFailedPicking([2])
        read = []

        def reader(workDir, mics, coordSet, invertY=False):
            read.extend(mic.getObjId() for mic in mics)

        with mock.patch.object(module, 'readSetOfCoordinates', reader):
            protocol.readCoordsFromMics(self._tmp.name, [_Mic(1), _Mic(2)],
                                        _Coords())

        self.assertEqual([1], read)

    def test_AMicrographThatFailsHalfWayPublishesNothing(self):
        protocol = _Picking(self._tmp.name)

        def reader(workDir, mics, coordSet, invertY=False):
            coordSet.append(('coord', mics[0].getObjId()))
            coordSet.setBoxSize(70)
            if mics[0].getObjId() == 2:
                raise ValueError("truncated json")

        coords = _Coords()

        with mock.patch.object(module, 'readSetOfCoordinates', reader):
            protocol.readCoordsFromMics(self._tmp.name, [_Mic(1), _Mic(2)],
                                        coords)

        self.assertEqual([('coord', 1)], coords.items)
        self.assertEqual({2}, protocol._failedPickingMicIds())
        self.assertEqual(70, coords.boxSize, "The box size e2boxer recorded wins.")


class TestFailedMicrographsArePickedAgain(_Base):
    def test_ContinueForgetsWhatMadeThemLookDone(self):
        protocol = _Picking(self._tmp.name)
        protocol.micDict = OrderedDict([("mic_001", _Mic(1)), ("mic_002", _Mic(2))])
        results = protocol._pickingStepResults(["mic_001", "mic_002"])
        for micId in (1, 2):
            open(protocol._getExtraPath('DONE', 'mic_%06d.TXT' % micId), 'w').close()
        with open(protocol._getExtraPath('DONE', 'all.TXT'), 'w') as handle:
            handle.write("1\n2\n")
        protocol._recordFailedPicking([2])

        protocol._prepareFailedForRetry()

        self.assertEqual([True, False], [os.path.exists(r) for r in results])
        self.assertTrue(os.path.exists(protocol._getExtraPath('DONE', 'mic_000001.TXT')))
        self.assertFalse(os.path.exists(protocol._getExtraPath('DONE', 'mic_000002.TXT')))
        self.assertEqual([1], protocol._readDoneList())

    def test_AMicrographPickedWellIsNoLongerFailed(self):
        protocol = _Picking(self._tmp.name)
        protocol._recordFailedPicking([1])

        protocol._forgetFailedPicking([1])

        self.assertEqual(set(), _Picking(self._tmp.name)._failedPickingMicIds())


class TestPickingFailsOnlyWhenEveryMicrographFailed(_Base):
    def test_EveryMicrographFailedFailsTheRun(self):
        protocol = _Picking(self._tmp.name)
        protocol._recordFailedPicking([1, 2])

        with self.assertRaisesRegex(Exception, "failed for every one of the 2"):
            protocol._failIfEveryMicrographFailed()

    def test_SomePickedMicrographsAreNoFailure(self):
        protocol = _Picking(self._tmp.name)
        protocol._recordFailedPicking([2])
        coords = _Coords()
        coords.append('coord')
        protocol.outputCoordinates = coords

        protocol._failIfEveryMicrographFailed()


if __name__ == '__main__':
    unittest.main()
