# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
import subprocess

import pwem
import pyworkflow.utils as pwutils

from .constants import EMAN2_HOME, V2_3, V2_31, V3_0_0


_logo = "eman2_logo.png"
_references = ['Tang2007']


SCRATCHDIR = pwutils.getEnvVariable('EMAN2SCRATCHDIR', default='/tmp/')


class Plugin(pwem.Plugin):
    _homeVar = EMAN2_HOME
    _pathVars = [EMAN2_HOME]
    _supportedVersions = [V2_3, V2_31, V3_0_0]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMAN2_HOME, 'eman-' + cls.getActiveVersion())

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = pwutils.Environ(os.environ)
        pathList = [cls.getHome(d) for d in ['lib', 'bin']]

        # This environment variable is used to setup OpenGL (Mesa)
        # library in remote desktops
        if 'REMOTE_MESA_LIB' in os.environ:
            pathList.append(os.environ['REMOTE_MESA_LIB'])

        environ.update({'PATH': cls.getHome('bin')},
                       position=pwutils.Environ.BEGIN)

        environ.update({
            'LD_LIBRARY_PATH': os.pathsep.join(pathList),
            'PYTHONPATH': os.pathsep.join(pathList),
            'SCIPION_MPI_FLAGS': os.environ.get('EMANMPIOPTS', '')
        }, position=pwutils.Environ.REPLACE)

        return environ

    @classmethod
    def isVersion(cls, version='2.31'):
        return cls.getActiveVersion().startswith(version)

    @classmethod
    def getEmanActivation(cls):
        return "conda activate eman" + V3_0_0

    @classmethod
    def getProgram(cls, program, python=False):
        """ Return the program binary that will be used. """
        if cls.isVersion(V3_0_0):
            program = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getEmanActivation(), program)
        else:
            program = os.path.join(cls.getHome('bin'), program)

        if python:
            python = cls.getHome('bin/python')
            return '%(python)s %(program)s ' % locals()
        else:
            return '%(program)s ' % locals()

    @classmethod
    def getEmanCommand(cls, program, args, python=False):
        return cls.getProgram(program, python) + args

    @classmethod
    def getBoxerCommand(cls, boxerVersion='new'):
        cmd = 'e2boxer.py' if boxerVersion == 'new' else 'e2boxer_old.py'

        return os.path.join(cls.getHome('bin'), cmd)

    @classmethod
    def createEmanProcess(cls, script='e2converter.py', args=None, direc="."):
        """ Open a new Process with all EMAN environment (python...etc)
        that will server as an adaptor to use EMAN library
        """
        program = os.path.join(__path__[0], script)
        cmd = cls.getEmanCommand(program, args, python=True)

        print("** Running: '%s'" % cmd)
        cmd = cmd.split()
        proc = subprocess.Popen(cmd, env=cls.getEnviron(),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                cwd=direc,
                                universal_newlines=True)

        # Python 2 to 3 conversion: iterating over lines in subprocess stdout -> If universal_newlines is False the file
        # objects stdin, stdout and stderr will be opened as binary streams, and no line ending conversion is done.
        # If universal_newlines is True, these file objects will be opened as text streams in universal newlines mode
        # using the encoding returned by locale.getpreferredencoding(False). For stdin, line ending characters '\n' in
        # the input will be converted to the default line separator os.linesep. For stdout and stderr, all line endings
        # in the output will be converted to '\n'. For more information see the documentation of the io.TextIOWrapper
        # class when the newline argument to its constructor is None.

        return proc

    @classmethod
    def defineBinaries(cls, env):
        SW_EM = env.getEmFolder()

        shell = os.environ.get("SHELL", "bash")
        eman23_commands = [
            (shell + ' ./eman2.3.linux64.sh -b -p "%s/eman-2.3"' %
             SW_EM, '%s/eman-2.3/bin/python' % SW_EM)]
        eman231_commands = [
            (shell + ' ./eman2.31_sphire1.3.linux64.sh -b -p "%s/eman-2.31"' %
             SW_EM, '%s/eman-2.31/bin/python' % SW_EM)]

        # For Eman3.0.0-alpha
        installationCmd = cls.getCondaActivationCmd()
        installationCmd += 'conda create -n eman' + V3_0_0 + ' eman-deps-dev=22.1 -c cryoem -c defaults -c conda-forge && '
        installationCmd += 'cd .. && mv eman2-* eman-source && '
        installationCmd += 'mkdir eman-build && '
        installationCmd += 'conda activate eman' + V3_0_0 + ' && '
        installationCmd += 'cd eman-build && '
        installationCmd += 'cmake ../eman-source/ -DENABLE_OPTIMIZE_MACHINE=ON && '
        installationCmd += 'make -j 4 && make install'
        eman3_commands = [(installationCmd, "")]

        env.addPackage('eman', version='2.3',
                       tar='eman2.3.linux64.tgz',
                       commands=eman23_commands)

        env.addPackage('eman', version='2.31',
                       tar='eman2.31.linux64.tgz',
                       commands=eman231_commands,
                       default=False)

        env.addPackage('eman', version='3.0.0-alpha',
                       # url='https://github.com/cryoem/eman2/tarball/master/',
                       url='https://github.com/cryoem/eman2/archive/8170d34.tar.gz',
                       buildDir='eman2-8170d345255c39a2441109562cccf4cb59e7e014',
                       commands=eman3_commands,
                       targetDir="eman-source",
                       default=True)
