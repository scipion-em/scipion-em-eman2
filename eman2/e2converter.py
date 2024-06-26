# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
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
"""
This script should be launched using the EMAN2 python interpreter 
and its own environment.
This scripts will convert any SetOfImages to an EMAN2 .hdf stack
It will read from the stdin the number of images and then
for each image will read: index, filename and a possible transformation.
As parameters will receive the output filename for the hdf stack
"""

import json
import os
import sys

import EMAN2 as eman


MODE_WRITE = 'write'
MODE_READ = 'read'


def writeParticles():
    line = sys.stdin.readline()
    fnHdf = ""
    while line:
        objDict = json.loads(line)
        index = int(objDict['_index'])
        filename = str(objDict['_filename'])
        imageData = eman.EMData()
        imageData.read_image(filename, index)

        if '_ctfModel._defocusU' in objDict:
            ctf = eman.EMAN2Ctf()
            defU = objDict['_ctfModel._defocusU']
            defV = objDict['_ctfModel._defocusV']

            ctf.from_dict({"defocus": (defU + defV) / 20000.0,
                           "dfang": objDict['_ctfModel._defocusAngle'],
                           "dfdiff": (defU - defV) / 10000.0,
                           "voltage": objDict['_acquisition._voltage'],
                           "cs": max(objDict['_acquisition._sphericalAberration'], 0.0001),
                           "ampcont": objDict['_acquisition._amplitudeContrast'] * 100.0,
                           "apix": objDict['_samplingRate']})
            imageData.set_attr('ctf', ctf)

        imageData.set_attr('apix_x', objDict['_samplingRate'])
        imageData.set_attr('apix_y', objDict['_samplingRate'])

        transformation = None
        if '_angles' in objDict:
            # TODO: convert to vector not matrix
            angles = objDict['_angles']
            shifts = objDict['_shifts']
            transformation = eman.Transform({"type": "spider",
                                             "phi": angles[0],
                                             "theta": angles[1],
                                             "psi": angles[2],
                                             "tx": -shifts[0],
                                             "ty": -shifts[1],
                                             "tz": -shifts[2],
                                             "mirror": 0,  # TODO: test flip
                                             "scale": 1.0})

        if transformation is not None:
            imageData.set_attr('xform.projection', transformation)

        outputFile = str(objDict['hdfFn'])
        if outputFile != fnHdf:
            i = 0
            fnHdf = outputFile

        imageData.write_image(outputFile, i,
                              eman.EMUtil.ImageType.IMAGE_HDF, False)
        i += 1
        print("OK")  # it is necessary to add newline
        sys.stdout.flush()
        line = sys.stdin.readline()


def readParticles(inputParts, inputCls, inputClasses, outputTxt, alitype='3d'):
    imgs = eman.EMUtil.get_image_count(inputParts)
    clsClassDict = {}
    shiftXList = {}
    shiftYList = {}
    dAlphaList = {}
    flipList = {}
    with open(outputTxt, 'w') as f:

        if alitype == '2d':
            # reading 2d refinement results
            clsImgs = eman.EMData.read_images(inputCls)
            classes = eman.EMData.read_images(inputClasses)
            clsClassList = clsImgs[0]
            f.write('#index, enable, cls, rot, tilt, psi, shiftX, shiftY\n')

            for i in range(clsClassList.get_attr_dict()['ny']):
                clsClassDict[i] = clsClassList[0, i]
                shiftXList[i] = clsImgs[2][0, i]
                shiftYList[i] = clsImgs[3][0, i]
                dAlphaList[i] = clsImgs[4][0, i]
                flipList[i] = clsImgs[5][0, i]

            # now convert eman orientation to scipion
            for index in range(imgs):
                classNum = clsClassDict[index]
                imgRotation = classes[int(classNum)].get_attr_dict().get('xform.projection', None)

                if imgRotation is not None:
                    enable = 1
                    transform = eman.Transform({"type": "2d",
                                                "alpha": dAlphaList[index],
                                                "mirror": True if flipList[index] else False,
                                                "tx": shiftXList[index],
                                                "ty": shiftYList[index],
                                                })
                    if flipList[index]:
                        tilt = 180 - transform.get_rotation("spider")['theta']
                        psi = transform.get_rotation("spider")['phi'] * -1
                    else:
                        tilt = transform.get_rotation("spider")['theta']
                        psi = transform.get_rotation("spider")['phi']

                    rot = transform.get_rotation("spider")['psi']
                    shifts = transform.get_trans()
                    shiftX, shiftY = shifts[0], shifts[1]

                    f.write(('{} '*8 + '\n').format(index, enable, int(classNum), rot, tilt, psi, shiftX, shiftY))
                else:
                    # disabled image
                    f.write(('{} '*2 + '\n').format(index, 0))

        else:
            # reading 3d refinement results
            clsImgsEven = eman.EMData.read_images(inputCls + "_even.hdf")
            clsImgsOdd = eman.EMData.read_images(inputCls + "_odd.hdf")
            classesEven = eman.EMData.read_images(inputClasses + "_even.hdf")
            classesOdd = eman.EMData.read_images(inputClasses + "_odd.hdf")

            clsClassListEven = clsImgsEven[0]
            clsClassListOdd = clsImgsOdd[0]
            f.write('#index, enable, rot, tilt, psi, shiftX, shiftY\n')

            for i in range(clsClassListEven.get_attr_dict()['ny']):
                clsClassDict[2 * i] = clsClassListEven[0, i]
                shiftXList[2 * i] = clsImgsEven[2][0, i]
                shiftYList[2 * i] = clsImgsEven[3][0, i]
                dAlphaList[2 * i] = clsImgsEven[4][0, i]
                flipList[2 * i] = clsImgsEven[5][0, i]

            for i in range(clsClassListOdd.get_attr_dict()['ny']):
                clsClassDict[2 * i + 1] = clsClassListOdd[0, i]
                shiftXList[2 * i + 1] = clsImgsOdd[2][0, i]
                shiftYList[2 * i + 1] = clsImgsOdd[3][0, i]
                dAlphaList[2 * i + 1] = clsImgsOdd[4][0, i]
                flipList[2 * i + 1] = clsImgsOdd[5][0, i]

            # now convert eman orientation to scipion
            for index in range(imgs):
                classNum = clsClassDict[index]
                if index % 2 == 0:
                    classes = classesEven
                else:
                    classes = classesOdd

                imgRotation = classes[int(classNum)].get_attr_dict().get('xform.projection', None)

                if imgRotation is not None:
                    enable = 1
                    az = imgRotation.get_rotation("eman")['az']
                    alt = imgRotation.get_rotation("eman")['alt']

                    transform = eman.Transform({"type": "eman",
                                                "az": az,
                                                "alt": alt,
                                                "phi": dAlphaList[index],
                                                "mirror": True if flipList[index] else False,
                                                "tx": shiftXList[index],
                                                "ty": shiftYList[index],
                                                })
                    transform = transform.inverse()

                    if flipList[index]:
                        tilt = 180 - transform.get_rotation("spider")['theta']
                        psi = transform.get_rotation("spider")['phi'] * -1
                    else:
                        tilt = transform.get_rotation("spider")['theta']
                        psi = transform.get_rotation("spider")['phi']

                    rot = transform.get_rotation("spider")['psi']
                    shifts = transform.get_trans()
                    shiftX, shiftY = shifts[0], shifts[1]

                    f.write(('{} '*7 + '\n').format(index, enable, rot, tilt, psi, shiftX, shiftY))
                else:
                    # disabled image
                    f.write(('{} '*2 + '\n').format(index, 0))


if __name__ == '__main__':
    if len(sys.argv) > 0:
        mode = sys.argv[1]

        if mode == MODE_WRITE:
            writeParticles()
        elif mode == MODE_READ:
            inputParts = sys.argv[2]
            inputCls = sys.argv[3]
            inputClasses = sys.argv[4]
            outputTxt = sys.argv[5]
            alitype = sys.argv[6]
            readParticles(inputParts, inputCls, inputClasses, outputTxt, alitype)
        else:
            raise ValueError("e2converter: Unknown mode '%s'" % mode)
    else:
        print("usage: %s outputFile" % os.path.basename(sys.argv[0]))
