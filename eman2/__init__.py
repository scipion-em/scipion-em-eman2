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
from pyworkflow.utils import Environ, getEnvVariable, join


_logo = "eman2_logo.jpg"
_references = ['Tang2007']

EMAN_DIR_VAR = 'EMAN2DIR'
SCRATCHDIR = getEnvVariable('EMAN2SCRATCHDIR', default='/tmp/')


def hello():
    print ("a different hello from Eman2....")

# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    #__metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = Environ(os.environ)
        EMAN2DIR = os.environ[('%s' % EMAN_DIR_VAR)]
        pathList = [os.path.join(EMAN2DIR, d) for d in ['lib', 'bin']]

        if cls.getActiveVersion() in ['2.11', '2.12']:
            pathList.append(os.path.join(EMAN2DIR, 'extlib/site-packages'))

        # This environment variable is used to setup OpenGL (Mesa)
        # library in remote desktops
        if 'REMOTE_MESA_LIB' in os.environ:
            pathList.append(os.environ['REMOTE_MESA_LIB'])

        environ.update({'PATH': join(EMAN2DIR, 'bin')}, position=Environ.BEGIN)

        environ.update({
            'LD_LIBRARY_PATH': os.pathsep.join(pathList),
            'PYTHONPATH': os.pathsep.join(pathList),
            'SCIPION_MPI_FLAGS': os.environ.get('EMANMPIOPTS', '')
        }, position=Environ.REPLACE)

        if cls.getActiveVersion() in ['2.11', '2.12']:
            environ.update({'EMAN_PYTHON': os.path.join(EMAN2DIR, 'Python/bin/python')
                            }, position=Environ.END)
        return environ

    @classmethod
    def getActiveVersion(cls):
        """ Return the version of the Eman2 binaries that is currently active.
        In the current implementation it will be inferred from the EMAN2DIR
        variable, so it should contain the version number in it. """
        path = os.environ['EMAN2DIR']
        for v in cls.getSupportedVersions():
            if v in path:
                return v
        return ''

    @classmethod
    def isNewVersion(cls):
        return not cls.getActiveVersion().startswith("2.1")

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return ['2.11', '2.12', '2.21']

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if RELION binaries are
        properly installed. """
        environ = cls.getEnviron()
        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [EMAN_DIR_VAR]
                        if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []

    @classmethod
    def validateVersion(protocol, errors):
        """ Validate if eman version is set properly according
         to installed version and the one set in the config file.
         Params:
            protocol: the input protocol calling to validate
            errors: a list that will be used to add the error message.
        """
        protocol.validatePackageVersion('EMAN2DIR', errors)

    @classmethod
    def getEmanProgram(cls, program):
        if not 'EMAN_PYTHON' in os.environ:
            if cls.getActiveVersion() in ['2.11', '2.12']:
                pyPath = 'Python/bin/python'
            else:
                pyPath = 'bin/python'
            os.environ['EMAN_PYTHON'] = os.path.join(os.environ['EMAN2DIR'], pyPath)
        # For EMAN2 python scripts, join the path to bin
        program = os.path.join(os.environ['EMAN2DIR'], 'bin', program)
        python = os.environ['EMAN_PYTHON']
        return '%(python)s %(program)s ' % locals()

    @classmethod
    def getEmanCommand(cls, program, args):
        return cls.getEmanProgram(program) + args

    @classmethod
    def getBoxerCommand(cls, emanVersion=None, boxerVersion='new'):
        """ Returns the Boxer command depending on Eman2 version.
         If emanVersion is None, the current active version will be used.
        """
        emanVersion = emanVersion or cls.getActiveVersion()
        # returns boxer program depending on Eman version
        new = emanVersion in ['2.11', '2.12'] or boxerVersion == 'new'
        return 'e2boxer.py' if new else 'e2boxer_old.py'

    @classmethod
    def createEmanProcess(cls, script='e2converter.py', args=None, direc="."):
        """ Open a new Process with all EMAN environment (python...etc)
        that will server as an adaptor to use EMAN library
        """
        print cls.__path__
        program = join(cls.__path__, script)
        cmd = cls.getEmanCommand(program, args)

        print ("** Running: '%s'" % cmd)
        proc = subprocess.Popen(cmd, shell=True, env=cls.getEnviron(),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, cwd=direc)

        return proc


pyworkflow.em.Domain.registerPlugin(__name__)