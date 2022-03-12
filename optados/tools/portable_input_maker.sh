#!/bin/bash

EXE="../../od2od"

if [[ "$1" == "pack" ]]; then

	for FILE in *.dome_bin *.elnes_bin *.ome_bin *.pdos_bin;
	do
		if [ -f $FILE ];
		then
			EXT=${FILE##*.} 
			# Match after last dot	
			SEEDNAME=`basename $FILE .$EXT`
			FILETYPE=${EXT%_*} 
			echo $EXE ${FILETYPE}_bin ${FILETYPE}_fmt $SEEDNAME
			$EXE ${FILETYPE}_bin ${FILETYPE}_fmt $SEEDNAME
			bzip2 -9 ${SEEDNAME}.${FILETYPE}_fmt	
		fi	
	done
	bzip2 -k -9 *.bands

	if [ -f ${SEEDNAME}-out.cell ]; 
	then
		bzip2 -k -9 ${SEEDNAME}-out.cell
	fi

elif [[ "$1" == "unpack" ]]; then
	for FILE in *.bz;
	do
		bunzip $FILE
	done	

	for FILE in *.dome_fmt *.elnes_fmt *.ome_fmt *.pdos_fmt;
        do
                if [ -f $FILE ];
                then
                        EXT=${FILE##*.}
                        # Match after last dot
                        SEEDNAME=`basename $FILE .$EXT`
                        FILETYPE=${EXT%_*}
                        echo $EXE ${FILETYPE}_fmt ${FILETYPE}_bin $SEEDNAME
                        $EXE ${FILETYPE}_fmt ${FILETYPE}_bin $SEEDNAME
                  
                fi
        done
else
	echo ' Command should be either "pack" or "unpack" '
fi
