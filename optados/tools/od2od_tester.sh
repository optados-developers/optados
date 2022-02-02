#/bin/bash
# AJM 6/12/19

# Put in the examples directory
# un comment the castep running lines if you haven't run it before

OD2OD_EXE="../od2od" # /Users/ajm143/src/optados/optados/src/od2od"
CASTEP_EXE="castep.serial" #"/Users/ajm143/bin/castep.serial"

OME_SEEDNAME=Si2_OPTICS
DOME_SEEDNAME=Si2_CORE
PDOS_SEEDNAME=Si2
ELNES_SEEDNAME=Si2_CORE

OME_TESTDIR=Si2_OPTICS
DOME_TESTDIR=Si2_CORE
PDOS_TESTDIR=Si2_PDOS
ELNES_TESTDIR=Si2_CORE

WD=`pwd`

# OME TESTS
cd ${OME_TESTDIR}
#${CASTEP_EXE} ${OME_SEEDNAME}
## TEST 1 -- Size of the bin file
echo " TEST 1"
${OD2OD_EXE} ome_bin ome_fmt  ${OME_SEEDNAME} ${OME_SEEDNAME}_test1
cp ${OME_SEEDNAME}.bands ${OME_SEEDNAME}_test1.bands
${OD2OD_EXE} ome_fmt ome_bin ${OME_SEEDNAME}_test1 ${OME_SEEDNAME}_test2
ls -lth ${OME_SEEDNAME}_test2.ome_bin
ls -lth ${OME_SEEDNAME}.ome_bin
echo " Any differences in the BIN files?"
diff ${OME_SEEDNAME}.ome_bin ${OME_SEEDNAME}_test2.ome_bin
## TEST 2 -- Format files match 
echo " TEST 2"
cp ${OME_SEEDNAME}.bands ${OME_SEEDNAME}_test2.bands
${OD2OD_EXE} ome_bin ome_fmt ${OME_SEEDNAME}_test2 ${OME_SEEDNAME}_test3
ls -lth ${OME_SEEDNAME}_test1.ome_fmt
ls -lth ${OME_SEEDNAME}_test3.ome_fmt
echo " Any differences in the FMT files?"
diff ${OME_SEEDNAME}_test1.ome_fmt ${OME_SEEDNAME}_test3.ome_fmt
echo " TEST 3"


cd ${WD}

# DOME TESTS
cd ${DOME_TESTDIR}
#${CASTEP_EXE} ${DOME_SEEDNAME}
## TEST 1 -- Size of the bin file
echo " TEST 1"
${OD2OD_EXE} dome_bin dome_fmt  ${DOME_SEEDNAME} ${DOME_SEEDNAME}_test1
cp ${DOME_SEEDNAME}.bands ${DOME_SEEDNAME}_test1.bands
${OD2OD_EXE} dome_fmt dome_bin ${DOME_SEEDNAME}_test1 ${DOME_SEEDNAME}_test2
ls -lth ${DOME_SEEDNAME}_test2.dome_bin
ls -lth ${DOME_SEEDNAME}.dome_bin
echo " Any differences in the BIN files?"
diff ${DOME_SEEDNAME}.dome_bin ${DOME_SEEDNAME}_test2.dome_bin
## TEST 2 -- Format files match 
echo " TEST 2"
cp ${DOME_SEEDNAME}.bands ${DOME_SEEDNAME}_test2.bands
${OD2OD_EXE} dome_bin dome_fmt ${DOME_SEEDNAME}_test2 ${DOME_SEEDNAME}_test3
ls -lth ${DOME_SEEDNAME}_test1.dome_fmt
ls -lth ${DOME_SEEDNAME}_test3.dome_fmt
echo " Any differences in the FMT files?"
diff ${DOME_SEEDNAME}_test1.dome_fmt ${DOME_SEEDNAME}_test3.dome_fmt
cd ${WD}


# PDOS TESTS
cd ${PDOS_TESTDIR}
#${CASTEP_EXE} ${PDOS_SEEDNAME}
## TEST 1 -- Size of the bin file
echo " TEST 1"
${OD2OD_EXE} pdos_bin pdos_fmt  ${PDOS_SEEDNAME} ${PDOS_SEEDNAME}_test1
cp ${PDOS_SEEDNAME}.bands ${PDOS_SEEDNAME}_test1.bands
${OD2OD_EXE} pdos_fmt pdos_bin ${PDOS_SEEDNAME}_test1 ${PDOS_SEEDNAME}_test2
ls -lth ${PDOS_SEEDNAME}_test2.pdos_bin
ls -lth ${PDOS_SEEDNAME}.pdos_bin
echo " Any differences in the BIN files?"
diff ${PDOS_SEEDNAME}.pdos_bin ${PDOS_SEEDNAME}_test2.pdos_bin
## TEST 2 -- Format files match 
echo " TEST 2"
cp ${PDOS_SEEDNAME}.bands ${PDOS_SEEDNAME}_test2.bands
${OD2OD_EXE} pdos_bin pdos_fmt ${PDOS_SEEDNAME}_test2 ${PDOS_SEEDNAME}_test3
ls -lth ${PDOS_SEEDNAME}_test1.pdos_fmt
ls -lth ${PDOS_SEEDNAME}_test3.pdos_fmt
echo " Any differences in the FMT files?"
diff ${PDOS_SEEDNAME}_test1.pdos_fmt ${PDOS_SEEDNAME}_test3.pdos_fmt
cd ${WD}


# ELNES TESTS
cd ${ELNES_TESTDIR}
#${CASTEP_EXE} ${ELNES_SEEDNAME}
# TEST 1 -- Size of the bin file
echo " TEST 1"
${OD2OD_EXE} elnes_bin elnes_fmt  ${ELNES_SEEDNAME} ${ELNES_SEEDNAME}_test1
cp ${ELNES_SEEDNAME}.bands ${ELNES_SEEDNAME}_test1.bands
${OD2OD_EXE} elnes_fmt elnes_bin ${ELNES_SEEDNAME}_test1 ${ELNES_SEEDNAME}_test2
ls -lth ${ELNES_SEEDNAME}_test2.elnes_bin
ls -lth ${ELNES_SEEDNAME}.elnes_bin
echo " Any differences in the BIN files?"
diff ${ELNES_SEEDNAME}.elnes_bin ${ELNES_SEEDNAME}_test2.elnes_bin
## TEST 2 -- Format files match 
echo " TEST 2"
cp ${ELNES_SEEDNAME}.bands ${ELNES_SEEDNAME}_test2.bands
${OD2OD_EXE} elnes_bin elnes_fmt ${ELNES_SEEDNAME}_test2 ${ELNES_SEEDNAME}_test3
ls -lth ${ELNES_SEEDNAME}_test1.elnes_fmt
ls -lth ${ELNES_SEEDNAME}_test3.elnes_fmt
echo " Any differences in the FMT files?"
diff ${ELNES_SEEDNAME}_test1.elnes_fmt ${ELNES_SEEDNAME}_test3.elnes_fmt
cd ${WD}



# Pad-slice test
