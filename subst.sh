#!/bin/bash
# Remove the -pedantic option for the sake of old compilers
cd $(dirname $0)
find . -name configure.ac | while read FILE; do
  sed -i.deleteme -e 's|-pedantic||g' $FILE
done
find . -name '*.deleteme' -delete
