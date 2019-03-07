#!/bin/tcsh

#Saves Bichsel functions from the given STAR Libary version to a file.
#If the file exists, then the functions are added to the file.

#Get the Library version from the user argument
set libVersion = $1

#Make sure there is just one argument
if ( "$#" != 1 ) then
    echo "ERROR: This script requires only one argument - the star lib version."
    echo "EXITING."
    exit 1
    endif

#Make sure the argument is not empty
if ( "$libVersion" == "" ) then
    echo "ERROR: No library version given. You must pass a library version!"
    echo "EXITING."
    exit 1
endif

#Make sure the argument starts with SL
set temp = `echo $libVersion | cut -c1-2`
if ( "$temp" != "SL" ) then
    echo "WARNING: Your library version does not begin with the characters SL per STAR standards!"
endif

#Get the Current Libary Version so we can change back to it after we are done
set currentLibVersion = ` echo $STAR | cut -d"/" -f6`

#Change to the requested library version
starver $libVersion

#Get the Bichsel Function
root4star -l -b -q BichselFunctionCollector.C\(\""$libVersion"\"\)

#Chech to make sure there is a value for the currentLibVersion, 
#if there is then change back to it, if not warn the user and do nothing
if ( "$currentLibVersion" != "" ) then
    echo "Returning you to starver $currentLibVersion where you started."
    starver $currentLibVersion
else
    echo "No value of starver was set when you ran this script so that is how you are being returned."
    starver
endif

exit 0
