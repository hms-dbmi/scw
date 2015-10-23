# This file should be 'source'd from
# the directory in which it resides.
domain=$( domainname )
if [ "$domain"=="Orchestra" ]; then
    module load stats/R/3.2.1
    export R_LIBS=/groups/pklab/scw2014/lib
    echo "Using pre-installed R libraries."
else 
    echo "You will need to install the following R libraries:"
    cat required_libraries.txt
fi
export SCW_HOME=$( pwd )

