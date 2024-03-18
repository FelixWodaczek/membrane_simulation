#!/usr/bin/env bash

if [ $# -eq 0 ];
  then
    file_path="directories_in_directory.dat"
else
    file_path="$1"directories_in_directory.dat
fi

pardir="$(dirname $file_path)"

if [ ! -f $file_path ];
  then
    echo "File $file_path does not exist"
    exit 1
fi

if [ -f "$pardir"/rerun_directories_in_directory.dat ];
  then
    rm "$pardir"/rerun_directories_in_directory.dat
fi
touch "$pardir"/rerun_directories_in_directory.dat

array_index=1
for subdir in $(cat $file_path);
do
    n_steps=$(find "$pardir"/"$subdir"/data/ -name "*.txt" -printf '.' | wc -m) # $(ls "$pardir"/"$subdir"/data/*.txt | wc -l)
    # echo $n_steps
    if [ $n_steps -eq 0 ];
      then
        echo "$subdir" >> "$pardir"/rerun_directories_in_directory.dat
    fi
    array_index+=1
done