#!/bin/bash
# Bash script to export all toolbox paths

# Change a possible relative path to an absolute one
# This is just a hack
cd ${TOOLBOX_SRC}
TOOLBOX_SRC=`pwd`
cd -

# Framework paths
export DRIVER_SRC=${TOOLBOX_SRC}"/driver"
export IO_SRC=${TOOLBOX_SRC}"/io"
export POSTPR_SRC=${TOOLBOX_SRC}"/postprocessing"
export TOOLS_SRC=${TOOLBOX_SRC}"/tools"
export UTILITY_SRC=${TOOLBOX_SRC}"/utility"

export FRAME_SRC=${DRIVER_SRC}"/frame"
TOOLBOX_INC=" -I"${FRAME_SRC}

export MNTR_SRC=${FRAME_SRC}"/mntr"
TOOLBOX_INC+=" -I"${MNTR_SRC}
export RPRM_SRC=${FRAME_SRC}"/rprm"
TOOLBOX_INC+=" -I"${RPRM_SRC}

export CHKPT_SRC=${IO_SRC}"/chkpt"
TOOLBOX_INC+=" -I"${CHKPT_SRC}
export CHKPTDMM_SRC=${CHKPT_SRC}"/chkptdm"
export CHKPTMS_SRC=${CHKPT_SRC}"/chkptms"
TOOLBOX_INC+=" -I"${CHKPTMS_SRC}

export IO_TOOLS_SRC=${IO_SRC}"/io_tools"
TOOLBOX_INC+=" -I"${IO_TOOLS_SRC}

# Statistics post-processing include file has to be copied to setup directory
export PSTAT2D_SRC=${POSTPR_SRC}"/pstat2d"
TOOLBOX_INC+=" -I"${PSTAT2D_SRC}
export PSTAT3D_SRC=${POSTPR_SRC}"/pstat3d"
TOOLBOX_INC+=" -I"${PSTAT3D_SRC}

export BASEFLOW_SRC=${TOOLS_SRC}"/baseflow"
# Statistics include file has to be copied to setup directory
export STAT_SRC=${TOOLS_SRC}"/stat"
# Time series include file has to be copied to setup directory
export TSRS_SRC=${TOOLS_SRC}"/tsrs"
export TSTPR_SRC=${TOOLS_SRC}"/tstpr"
TOOLBOX_INC+=" -I"${TSTPR_SRC}

export SFD_SRC=${BASEFLOW_SRC}"/sfd"
TOOLBOX_INC+=" -I"${SFD_SRC}

# Arnoldi include file has to be copied to setup directory
export ARNA_SRC=${TSTPR_SRC}"/arna"
export PWIT_SRC=${TSTPR_SRC}"/pwit"
TOOLBOX_INC+=" -I"${PWIT_SRC}

export BCND_SRC=${UTILITY_SRC}"/bcnd"
export COMM_SRC=${UTILITY_SRC}"/comm"
TOOLBOX_INC+=" -I"${COMM_SRC}
export CNHT_SRC=${UTILITY_SRC}"/cnht"
TOOLBOX_INC+=" -I"${CNHT_SRC}
export FORCING_SRC=${UTILITY_SRC}"/forcing"
export GRID_SRC=${UTILITY_SRC}"/grid"
export MATH_SRC=${UTILITY_SRC}"/math"
TOOLBOX_INC+=" -I"${MATH_SRC}

# Generalised synthetic eddy method include file has to be copied to setup directory
export GSYEM_SRC=${BCND_SRC}"/gsyem"

export NSEB_SRC=${FORCING_SRC}"/nseb"
TOOLBOX_INC+=" -I"${NSEB_SRC}
export SPNB_SRC=${FORCING_SRC}"/spnb"
TOOLBOX_INC+=" -I"${SPNB_SRC}
# Tripping line include file has to be copied to setup directory
export TRIPL_SRC=${FORCING_SRC}"/tripl"

export MAP2D_SRC=${GRID_SRC}"/map2D"
TOOLBOX_INC+=" -I"${MAP2D_SRC}

# Include path
export TOOLBOX_INC

