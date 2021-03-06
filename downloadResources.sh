#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to download
# the SKIRT 9 resource files provided on the public SKIRT server, and
# place them in the 'resources' directory next to the git directory.
#

# select download command: wget (Linux) or curl (Mac OS X)
if which wget >/dev/null
then
    DOWNLOAD="wget --no-check-certificate https://sciences.ugent.be/skirtextdat/SKIRT9/Resources/"
elif which curl >/dev/null
then
    DOWNLOAD="curl --insecure -O https://sciences.ugent.be/skirtextdat/SKIRT9/Resources/"
else
    echo error: no wget or curl available to download files
    exit 1
fi

# loop over the list of expected archives and versions given in the text file in the SKIRT repository
while read -u 3 LINE        # use explicit file descriptor 3 to allow nested read from terminal
do
    read NAME VERSION <<< $LINE
    RESOURCENAME=SKIRT9_Resources_$NAME
    VERSIONPATH=../resources/${RESOURCENAME}/version.txt
    ZIPFILENAME=${RESOURCENAME}_v${VERSION}.zip

    # find out whether the expected version is already installed, and if not, ask the user what to do
    PROCEED=0
    if [ -e $VERSIONPATH ]
    then
        read INSTALLEDVERSION <<< $(<$VERSIONPATH)
        if [ $INSTALLEDVERSION = $VERSION ]
        then
            echo $RESOURCENAME version $VERSION is already installed -- skipping
        else
            echo $RESOURCENAME version $INSTALLEDVERSION is installed while version $VERSION is expected
            read -p "Do you want to download and install $RESOURCENAME version $VERSION? [y/n] " RESPONSE
            if [[ $RESPONSE =~ ^[Yy] ]]
            then
                PROCEED=1
            fi
        fi
    else
        echo $RESOURCENAME is not installed
        read -p "Do you want to download and install $RESOURCENAME version $VERSION? [y/n] " RESPONSE
        if [[ $RESPONSE =~ ^[Yy] ]]
        then
            PROCEED=1
        fi
    fi

    # if confirmed by the user, download the archive
    if [ $PROCEED = 1 ]
    then
        echo Downloading $ZIPFILENAME ...
        echo "------------------------------------------------"
        mkdir -p ../resources
        cd ../resources
        $DOWNLOAD$ZIPFILENAME
        unzip -o $ZIPFILENAME  # overwrite exisiting files
        rm $ZIPFILENAME
        cd ../git
        echo "------------------------------------------------"
    fi

done 3< "SKIRT/resources/ExpectedResources.txt"

echo Done.
