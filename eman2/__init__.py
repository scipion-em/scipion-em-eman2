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
import logging
logger = logging.getLogger(__name__)

import pwem
import pyworkflow.utils as pwutils

from .constants import EMAN2_HOME, EMAN2SCRATCHDIR, VERSIONS


__version__ = '3.4.2'
_logo = "eman2_logo.png"
_references = ['Tang2007']


class Plugin(pwem.Plugin):
    _homeVar = EMAN2_HOME
    _pathVars = [EMAN2_HOME]
    _supportedVersions = VERSIONS
    _url = "https://blake.bcm.edu/emanwiki/EMAN2"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMAN2_HOME, 'eman-' + VERSIONS[-1])
        cls._defineVar(EMAN2SCRATCHDIR, '/tmp')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = pwutils.Environ(os.environ)
        environ.update({'PATH': cls.getHome('bin')},
                       position=pwutils.Environ.BEGIN)

        for v in ['PYTHONPATH', 'PYTHONHOME']:
            if v in environ:
                del environ[v]

        return environ

    @classmethod
    def versionGE(cls, version):
        """ Return True if current version of eman is newer
         or equal than the input argument.
         Params:
            version: string version (semantic version, e.g 2.91)
        """
        v1 = cls.getActiveVersion()
        if v1 not in VERSIONS:
            raise Exception("This version of EMAN is not supported: ", v1)

        if VERSIONS.index(v1) < VERSIONS.index(version):
            return False
        return True

    @classmethod
    def getActiveVersion(cls, *args):
        """ Reimplemented here, assumes EMAN2_HOME = eman-xxx """
        ver = os.path.basename(cls.getHome()).split("-")[-1]
        versions = cls.getSupportedVersions()
        for v in versions:
            if v == ver:
                return v

        return ''

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
        logger.info(f"\tRunning: {cmd}")
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
        urls = ['https://cryoem.bcm.edu/cryoem/static/software/release-2.91/eman2.91_sphire1.4_sparx.linux64.sh',
                'https://cryoem.bcm.edu/cryoem/static/software/continuous_build/eman2_sphire_sparx.linux.unstable.sh']

        for ver, url in zip(VERSIONS, urls):
            install_cmd = 'cd %s && wget --no-check-certificate -q --show-progress %s && ' % (SW_EM, url)
            install_cmd += '%s ./%s -b -f -p "%s/eman-%s"' % (shell, url.split('/')[-1], SW_EM, ver)
            eman_commands = [(install_cmd, '%s/eman-%s/bin/python' % (SW_EM, ver))]

            env.addPackage('eman', version=ver,
                           tar='void.tgz',
                           commands=eman_commands, default=ver == VERSIONS[-1])
