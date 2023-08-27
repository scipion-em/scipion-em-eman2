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
from pyworkflow import Config

from .constants import (EMAN2SCRATCHDIR, VERSIONS, EMAN_ENV_ACTIVATION,
                        DEFAULT_ACTIVATION_CMD, EMAN_DEFAULT_VER_NUM)


__version__ = '3.6'
_logo = "eman2_logo.png"
_references = ['Tang2007']


class Plugin(pwem.Plugin):
    _supportedVersions = VERSIONS
    _url = "https://github.com/scipion-em/scipion-em-eman2"

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(EMAN_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)
        cls._defineVar(EMAN2SCRATCHDIR, '/tmp')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = pwutils.Environ(os.environ)
        for v in ['PYTHONPATH', 'PYTHONHOME']:
            if v in environ:
                del environ[v]

        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

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
        """ Return the env name that is currently active. """
        envVar = cls.getVar(EMAN_ENV_ACTIVATION)
        return envVar.split()[-1].split("-")[-1]

    @classmethod
    def getEmanEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(EMAN_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        return '%s %s' % (cls.getCondaActivationCmd(),
                          cls.getEmanEnvActivation())

    @classmethod
    def getProgram(cls, program, python=False):
        """ Create EMAN command line. """
        cmd = "python" if python else ""
        return f"{cls.getActivationCmd()} && {cmd} {program} "

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
        proc = subprocess.Popen(cmd, shell=True,  # required for "eval" to work
                                env=cls.getEnviron(),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                cwd=direc,
                                universal_newlines=True)

        return proc

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addEmanPackage(env, ver,
                               default=ver == EMAN_DEFAULT_VER_NUM)

    @classmethod
    def addEmanPackage(cls, env, version, default=False):
        ENV_NAME = f"eman-{version}"
        FLAG = f"eman_{version}_installed"

        # try to get CONDA activation command
        installCmds = [
            cls.getCondaActivationCmd(),
            f'conda create -y -n {ENV_NAME} eman-dev={version} -c cryoem -c conda-forge &&',
            f'touch {FLAG}'  # Flag installation finished
        ]
        emanCmds = [(" ".join(installCmds), FLAG)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage('eman', version=version,
                       tar='void.tgz',
                       commands=emanCmds,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)
