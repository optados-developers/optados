#!/bin/bash
# Runs OptaDOS for the test suite where its important to 
# store gzipped fmt files but give OptaDOS bin files
#
# Andrew J. Morris
# Feb 2022
SEEDNAME=$1
ARG2=$2

DEBUG="F"
if (${ARG2} == "--debug"); then
 DEBUG=T
fi

UNZIPPER="gunzip"
CONVERTER="../od2od"
EXE="../optados.x"

UNZIPPER_ARGS="-k" # keep

which ${EXE}
which ${CONVERTER}
which ${UNZIPPER}

echo " Seedname:  ${SEEDNAME}"

# Extract the formatted data 
for ZIPFILE in ${ARGS}*.gz
do
 echo " Untarring ${ZIPFILE}"
	${UNZIPPER} ${UNZIPPER_ARGS} ${ZIPFILE}
done

# Convert the formatted data to a bin
for FMTFILE in ome dome pdos elnes
do
  #check the file exits
  if test -f ${SEEDNAME}.${FMTFILE}_fmt; then
    echo " Creating ${SEEDNAME}.${FMTFILE}_bin from ${SEEDNAME}.${FMTFILE}_fmt"
    ${CONVERTER} ${FMTFILE}_fmt ${FMTFILE}_bin ${SEEDNAME}
  fi
done

# Run Optados
echo " Running OptaDOS: ${EXE} ${SEEDNAME}"
${EXE} ${SEEDNAME} 

# Clean up
if [ ${DEBUG}=="F" ]; then
  rm *_bin
  rm *_fmt
fi
