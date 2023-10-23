#!/bin/bash

# Search all files to find first one with data
file_with_data=""

for file in solidification/*.csv
do
    if [ -f "$file" ]
    then
        if [ -s "$file" ]
        then
            file_with_data="$file"
            break
        fi
    fi
done

if [ -z "$file_with_data" ]
then
    echo "No CSV file found with solidification fields"
    exit 1
fi

# Count the number of fields to determine header format
line=$(head -n 1 "$file_with_data")
nfields=$(echo ${line} | tr -cd , | wc -c)

if [ ${nfields} == 8 ]
then
   echo "x,y,z,tm,tl,cr,Gx,Gy,Gz" > solidification_data.csv
elif [ ${nfields} == 5 ]
then
   echo "x,y,z,tm,tl,cr" > solidification_data.csv
else
   echo "Unknown header format"
   exit 1
fi

# Append the contents of the processed file to solidification_data.csv
cat solidification/*csv >> solidification_data.csv
