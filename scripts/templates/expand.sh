#!/usr/bin/env bash

filename_rules=$1
filename_templ=$2
rm -f tmp
sh ${filename_rules} | while read line; do
    eval "$line"
    cat ${filename_templ} >> tmp
    for key in `echo "$line" | awk -F= '{for (i=1; i<NF; i++) {print $i}}' | awk '{print $NF}'`;
    do
      sed -i -e 's/${'${key}'}/'"${!key}"'/g' tmp
    done
done 
echo "SUCCESS: The generated files is called <tmp>!"
