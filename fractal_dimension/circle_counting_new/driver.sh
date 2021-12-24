#!/bin/bash

cmd="./target/release/circle_counting dimension"
flags="--max=$1 --n=100 --time"

for dir in $(find data -type d)
do
    echo $dir | sed -e "s/data/output$1/" | xargs mkdir
done

count=0
total=$(find data -type f | wc -l)
for file in $(find data -type f)
do
    newfile=$(echo $file | sed -e "s/data/output$1/")
    if test -f $newfile; then
        echo "$newfile already exists"
    else
        $cmd $file $flags | tee $newfile
    fi
    count=$(($count + 1))
    echo "$(bc <<< "(100 * $count) / $total")% done"
done
