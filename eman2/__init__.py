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

from .constants import EMAN2_HOME, EMAN2SCRATCHDIR, V2_9


__version__ = '3.3'
_logo = "eman2_logo.png"
_references = ['Tang2007']


class Plugin(pwem.Plugin):
    _homeVar = EMAN2_HOME
    _pathVars = [EMAN2_HOME]
    _supportedVersions = [V2_9]
    _url = "https://github.com/scipion-em/scipion-em-eman2"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMAN2_HOME, 'eman-' + V2_9)
        cls._defineVar(EMAN2SCRATCHDIR, '/tmp')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = pwutils.Environ(os.environ)
        environ.update({'PATH': cls.getHome('bin')},
                       position=pwutils.Environ.BEGIN)

        return environ

    @classmethod
    def isVersion(cls, version='2.9'):
        return cls.getActiveVersion().startswith(version)

    @classmethod
    def getProgram(cls, program, python=False):
        """ Return the program binary that will be used. """
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
    def createEmanProcess(cls, script='e2converter.py', args=None, direc="."):
        """ Open a new Process with all EMAN environment (python...etc)
        that will serve as an adaptor to use EMAN library
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

        return proc

    @classmethod
    def defineBinaries(cls, env):
        SW_EM = env.getEmFolder()
        shell = os.environ.get("SHELL", "bash")
        url = 'https://cryoem.bcm.edu/cryoem/static/software/release-2.9/eman2.9_sphire1.4_sparx.linux64.sh'
        install_cmd = 'cd %s && wget %s && ' % (SW_EM, url)
        install_cmd += '%s ./eman2.9_sphire1.4_sparx.linux64.sh -b -f -p "%s/eman-2.9"' % (shell, SW_EM)
        eman29_commands = [(install_cmd, '%s/eman-2.9/bin/python' % SW_EM)]

        env.addPackage('eman', version=V2_9,
                       tar='void.tgz',
                       commands=eman29_commands, default=True)
