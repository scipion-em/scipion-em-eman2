# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
import subprocess

import pyworkflow.em
import pyworkflow.utils as pwutils

from .constants import EMAN2DIR, V2_11, V2_12, V2_21


_logo = "eman2_logo.png"
_references = ['Tang2007']


SCRATCHDIR = pwutils.getEnvVariable('EMAN2SCRATCHDIR', default='/tmp/')


class Plugin(pyworkflow.em.Plugin):
    _homeVar = EMAN2DIR
    _pathVars = [EMAN2DIR]
    _supportedVersions = [V2_11, V2_12, V2_21]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMAN2DIR, 'eman-2.21')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = pwutils.Environ(os.environ)
        pathList = [cls.getHome(d) for d in ['lib', 'bin']]

        if not cls.isNewVersion():
            pathList.append(cls.getHome('extlib', 'site-packages'))

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

        if not cls.isNewVersion():
            environ.update({'EMAN_PYTHON': cls.getHome('Python/bin/python')},
                           position=pwutils.Environ.END)
        return environ

    @classmethod
    def isNewVersion(cls):
        return not cls.getActiveVersion().startswith("2.1")

    @classmethod
    def getProgram(cls, program, python=False):
        """ Return the program binary that will be used. """
        program = os.path.join(cls.getHome('bin'), program)

        if python:
            if cls.isNewVersion():
                python = cls.getHome('bin/python')
            else:
                python = cls.getHome('Python/bin/python')

            return '%(python)s %(program)s ' % locals()
        else:
            return '%(program)s ' % locals()

    @classmethod
    def getEmanCommand(cls, program, args, python=False):
        return cls.getProgram(program, python) + args

    @classmethod
    def getBoxerCommand(cls, emanVersion=None, boxerVersion='new'):
        """ Returns the Boxer command depending on Eman2 version.
         If emanVersion is None, the current active version will be used.
        """
        emanVersion = emanVersion or cls.getActiveVersion()
        new = emanVersion in [V2_11, V2_12] or boxerVersion == 'new'
        cmd = 'e2boxer.py' if new else 'e2boxer_old.py'

        return os.path.join(cls.getHome('bin'), cmd)

    @classmethod
    def createEmanProcess(cls, script='e2converter.py', args=None, direc="."):
        """ Open a new Process with all EMAN environment (python...etc)
        that will server as an adaptor to use EMAN library
        """
        program = os.path.join(__path__[0], script)
        cmd = cls.getEmanCommand(program, args, python=True)

        print ("** Running: '%s'" % cmd)
        proc = subprocess.Popen(cmd, shell=True, env=cls.getEnviron(),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, cwd=direc)

        return proc

    @classmethod
    def defineBinaries(cls, env):
        eman2_commands = [('./eman2-installer',
                           'eman2.*rc')]

        env.addPackage('eman', version='2.11',
                       tar='eman2.11.linux64.tgz',
                       commands=eman2_commands)

        env.addPackage('eman', version='2.12',
                       tar='eman2.12.linux64.tgz',
                       commands=eman2_commands)

        SW_EM = env.getEmFolder()

        eman22_commands = [
            ('./eman2.21.linux64.centos7.sh -b -p "%s/eman-2.21"' %
             SW_EM, '%s/eman-2.21/bin/python' % SW_EM)]

        env.addPackage('eman', version='2.21',
                       tar='eman2.21.linux64.centos7.tgz',
                       commands=eman22_commands,
                       default=True)


pyworkflow.em.Domain.registerPlugin(__name__)
