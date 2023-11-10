#!/bin/bash

while getopts i:o: flag
do
    case "${flag}" in
        i) input_directory=${OPTARG};;
        o) output_filename=${OPTARG};;
    esac
done
echo "Combining solidification data in: $input_directory";

# Search all files to find first one with data
file_with_data=""

for file in ${input_directory}/*.csv
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

echo "Creating the combined data file: $output_filename";
if [ ${nfields} == 8 ]
then
   echo "x,y,z,tm,tl,cr,Gx,Gy,Gz" > ${output_filename}
elif [ ${nfields} == 5 ]
then
   echo "x,y,z,tm,tl,cr" > ${output_filename}
else
   echo "Unknown header format"
   exit 1
fi

# Append the contents of the processed file
cat ${input_directory}/*csv >> ${output_filename}
