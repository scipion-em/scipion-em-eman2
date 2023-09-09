# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

EMAN2SCRATCHDIR = 'EMAN2SCRATCHDIR'

# Supported versions
VERSIONS = ['2.99.47', '2.99.52']
EMAN_DEFAULT_VER_NUM = VERSIONS[-1]

DEFAULT_ENV_NAME = f"eman-{EMAN_DEFAULT_VER_NUM}"
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
EMAN_ENV_ACTIVATION = 'EMAN_ENV_ACTIVATION'

# ------------------ Constants values -----------------------------------------

# ctf processing type
HIRES = 0
MIDRES = 1
LORES = 2

# ctf invar type
INVAR_AUTO = 0
INVAR_BISPEC = 1
INVAR_HARMONIC = 2

# centering algorithms
XFORM_NOCENTER = 0
XFORM_CENTER = 1
XFORM_CENTERACF = 2
XFORM_CENTEROFMASS = 3
XFORM_CENTER_NONE = 4

CENTER_CHOICES = {
    XFORM_NOCENTER: 'nocenter',
    XFORM_CENTER: 'xform.center',
    XFORM_CENTERACF: 'xform.centeracf',
    XFORM_CENTEROFMASS: 'xform.centerofmass',
    XFORM_CENTER_NONE: 'None'
}

# comparators
CMP_CCC = 0
CMP_DOT = 1
CMP_FRC = 2
CMP_FRC_FREQ = 3
CMP_LOD = 4
CMP_OPTSUB = 5
CMP_OPTVARIANCE = 6
CMP_PHASE = 7
CMP_QUADMINDOT = 8
CMP_SQEUCLIDEAN = 9
CMP_VERTICAL = 10
CMP_NONE = 11

SIMCMP_CHOICES = {
    CMP_CCC: 'ccc',
    CMP_DOT: 'dot',
    CMP_FRC: 'frc',
    CMP_FRC_FREQ: 'frc.freq',
    CMP_LOD: 'lod',
    CMP_OPTSUB: 'optsub',
    CMP_OPTVARIANCE: 'optvariance',
    CMP_PHASE: 'phase',
    CMP_QUADMINDOT: 'quadmindot',
    CMP_SQEUCLIDEAN: 'sqeuclidean',
    CMP_VERTICAL: 'vertical',
    CMP_NONE: 'None'
}

# aligners
ALN_FRM2D = 0
ALN_ROTATE_FLIP = 1
ALN_ROTATE_FLIP_ITERATIVE = 2
ALN_ROTATE_PRECENTER = 3
ALN_ROTATE_TRANS_FLIP_SCALE = 4
ALN_ROTATE_TRANS_FLIP_SCALE_ITER = 5
ALN_ROTATE_TRANS_SCALE_ITER = 6
ALN_ROTATE_TRANSLATE = 7
ALN_ROTATE_TRANSLATE_BISPEC = 8
ALN_ROTATE_TRANSLATE_FLIP = 9
ALN_ROTATE_TRANSLATE_FLIP_ITERATIVE = 10
ALN_ROTATE_TRANSLATE_FLIP_RESAMPLE = 11
ALN_ROTATE_TRANSLATE_ITERATIVE = 12
ALN_ROTATE_TRANSLATE_RESAMPLE = 13
ALN_ROTATE_TRANSLATE_SCALE = 14
ALN_ROTATE_TRANSLATE_TREE = 15
ALN_ROTATIONAL = 16
ALN_ROTATIONAL_BISPEC = 17
ALN_ROTATIONAL_ITERATIVE = 18
ALN_RTF_EXHAUSTIVE = 19
ALN_RTF_SLOW_EXHAUSTIVE = 20
ALN_SCALE = 21
ALN_SYMALIGN = 22
ALN_SYMALIGNQUAT = 23
ALN_TRANSLATIONAL = 24
ALN_NONE = 25

SIMALIGN_CHOICES = {
    ALN_FRM2D: 'frm2d',
    ALN_ROTATE_FLIP: 'rotate_flip',
    ALN_ROTATE_FLIP_ITERATIVE: 'rotate_flip_iterrative',
    ALN_ROTATE_PRECENTER: 'rotate_precenter',
    ALN_ROTATE_TRANS_FLIP_SCALE: 'rotate_trans_flip_scale',
    ALN_ROTATE_TRANS_FLIP_SCALE_ITER: 'rotate_trans_flip_scale_iter',
    ALN_ROTATE_TRANS_SCALE_ITER: 'rotate_trans_scale_iter',
    ALN_ROTATE_TRANSLATE: 'rotate_translate',
    ALN_ROTATE_TRANSLATE_BISPEC: 'rotate_translate_bispec',
    ALN_ROTATE_TRANSLATE_FLIP: 'rotate_translate_flip',
    ALN_ROTATE_TRANSLATE_FLIP_ITERATIVE: 'rotate_translate_flip_iterative',
    ALN_ROTATE_TRANSLATE_FLIP_RESAMPLE: 'rotate_translate_flip_resample',
    ALN_ROTATE_TRANSLATE_ITERATIVE: 'rotate_translate_iterative',
    ALN_ROTATE_TRANSLATE_RESAMPLE: 'rotate_translate_resample',
    ALN_ROTATE_TRANSLATE_SCALE: 'rotate_translate_scale',
    ALN_ROTATE_TRANSLATE_TREE: 'rotate_translate_tree',
    ALN_ROTATIONAL: 'rotational',
    ALN_ROTATIONAL_BISPEC: 'rotational_bispec',
    ALN_ROTATIONAL_ITERATIVE: 'rotational_iterative',
    ALN_RTF_EXHAUSTIVE: 'rtf_exhaustive',
    ALN_RTF_SLOW_EXHAUSTIVE: 'rtf_slow_exhaustive',
    ALN_SCALE: 'scale',
    ALN_SYMALIGN: 'symalign',
    ALN_SYMALIGNQUAT: 'symalignquat',
    ALN_TRANSLATIONAL: 'translational',
    ALN_NONE: 'None'
}

RALN_NONE = 0
RALN_REFINE = 1
RALN_REFINE_3D = 2
RALN_REFINE_3D_GRID = 3
RALN_REFINECG = 4

# averagers
AVG_CTF_AUTO = 0
AVG_CTF_WEIGHT = 1
AVG_CTF_WEIGHT_AUTOFILT = 2
AVG_CTFW_AUTO = 3
AVG_ITERATIVE = 4
AVG_LOCALWEIGHT = 5
AVG_MEAN = 6
AVG_MEAN_TOMO = 7
AVG_MEDIAN = 8
AVG_MINMAX = 9
AVG_SIGMA = 10
AVG_WEIGHTEDFOURIER = 11

AVG_CHOICES = {
    AVG_CTF_AUTO: 'ctf.auto',
    AVG_CTF_WEIGHT: 'ctf.weight',
    AVG_CTF_WEIGHT_AUTOFILT: 'ctf.weight.autofilt',
    AVG_CTFW_AUTO: 'ctfw.auto',
    AVG_ITERATIVE: 'iterative',
    AVG_LOCALWEIGHT: 'localweight',
    AVG_MEAN: 'mean',
    AVG_MEAN_TOMO: 'mean.tomo',
    AVG_MEDIAN: 'median',
    AVG_MINMAX: 'minmax',
    AVG_SIGMA: 'sigma',
    AVG_WEIGHTEDFOURIER: 'weightedfourier'
}

# processors normalize
PROC_NORMALIZE = 0
PROC_NORMALIZE_BYMASS = 1
PROC_NORMALIZE_CIRCLEMEAN = 2
PROC_NORMALIZE_EDGEMEAN = 3
PROC_NORMALIZE_HISTPEAK = 4
PROC_NORMALIZE_LOCAL = 5
PROC_NORMALIZE_LREDGE = 6
PROC_NORMALIZE_MASK = 7
PROC_NORMALIZE_MAXMIN = 8
PROC_NORMALIZE_RAMP_NORMVAR = 9
PROC_NORMALIZE_ROWS = 10
PROC_NORMALIZE_TOIMAGE = 11
PROC_NORMALIZE_UNITLEN = 12
PROC_NORMALIZE_UNITSUM = 13
PROC_NONE = 14

NORM_CHOICES = {
    PROC_NORMALIZE: 'normalize',
    PROC_NORMALIZE_BYMASS: 'normalize.bymass',
    PROC_NORMALIZE_CIRCLEMEAN: 'normalize.circlemean',
    PROC_NORMALIZE_EDGEMEAN: 'normalize.edgemean',
    PROC_NORMALIZE_HISTPEAK: 'normalize.histpeak',
    PROC_NORMALIZE_LOCAL: 'normalize.local',
    PROC_NORMALIZE_LREDGE: 'normalize.lredge',
    PROC_NORMALIZE_MASK: 'normalize.mask',
    PROC_NORMALIZE_MAXMIN: 'normalize.maxmin',
    PROC_NORMALIZE_RAMP_NORMVAR: 'normalize.ramp.normvar',
    PROC_NORMALIZE_ROWS: 'normalize.rows',
    PROC_NORMALIZE_TOIMAGE: 'normalize.toimage',
    PROC_NORMALIZE_UNITLEN: 'normalize.unitlen',
    PROC_NORMALIZE_UNITSUM: 'normalize.unitsum',
    PROC_NONE: 'None'
}

# Reconstruction methods
RECON_BACKPROJ = 0
RECON_FOURIER = 1
RECON_FOURIER_ITER = 2
RECON_FOURIER_SIMPLE = 3
RECON_NN4 = 4
RECON_NN4_CTF = 5
RECON_NN4_CTF_RECT = 6
RECON_NN4_CTFW = 7
RECON_NN4_CTFWS = 8
RECON_NN4_RECT = 9
RECON_NNSSNR = 10
RECON_NNSSNR_CTF = 11
RECON_REAL_MEDIAN = 12
RECON_WIENER_FOURIER = 13

# modes to reconstruct with fourier method
FOURIER_NEIGHBOR = 0
FOURIER_GAUSS2 = 1
FOURIER_GAUSS3 = 2
FOURIER_GAUSS5 = 3
FOURIER_GAUSS5_SLOW = 4
FOURIER_GYPERGEOM5 = 5
FOURIER_EXPERIMENTAL = 6

# speed
SPEED_1 = 0
SPEED_2 = 1
SPEED_3 = 2
SPEED_4 = 3
SPEED_5 = 4
SPEED_6 = 5
SPEED_7 = 6

# Keep parameter for e2make3d.py
KEEP_PERCENTAGE = 0
KEEP_STDDEV = 1
KEEP_ABSQUAL = 2

# Amplitude correction type for e2refine_easy
AMP_AUTO = 0
AMP_STRUCFAC = 1
AMP_FLATTEN = 2
AMP_NONE = 3

# tophat filter for e2refine_easy
TOPHAT_NONE = 0
TOPHAT_LOCAL = 1
TOPHAT_LOCALWIENER = 2
TOPHAT_GLOBAL = 3

# e2boxer autopick modes
AUTO_LOCAL = 0
AUTO_REF = 1
AUTO_CONVNET = 2
AUTO_GAUSS = 3

WIKI_URL = "[[http://blake.bcm.edu/emanwiki/EMAN2][Wiki]]"

# viewer.py constants
LAST_ITER = 0
ALL_ITERS = 1
SELECTED_ITERS = 2

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

TILT_SCATTER = 0
TILT_CONTOUR = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

FSC_UNMASK = 0
FSC_MASK = 1
FSC_MASKTIGHT = 2
FSC_ALL = 3

HALF_EVEN = 0
HALF_ODD = 1
FULL_MAP = 2
ALL_MAPS = 3

OBJCMD_CLASSAVG_PROJS = 'Show class-averages/projections'
OBJCMD_PROJS = 'Show only projections'
OBJCMD_INITVOL = 'Show initial random volume'

# SGD input types
SGD_INPUT_AVG = 0
SGD_INPUT_PTCLS = 1
