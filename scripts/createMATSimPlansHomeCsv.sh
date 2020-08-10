#!/usr/bin/env bash

dir=$(dirname "$0")

epsg=$(cat $dir/../scenarios/mount-alexander-shire/castlemaine-region/archetypes_matsim_main.xml | grep coordinateSystem | cut -f4 -d'"' | tr '[:upper:]' '[:lower:]')

outfile=$dir/../data/population-archetypes.home.csv
cmd="echo \"$epsg\" > $outfile"
echo $cmd && eval $cmd

#cmd="zless $dir/../scenarios/mount-alexander-shire/castlemaine-region/population-archetypes.xml.gz | grep 'activity type=\"home\"' | cut -d'\"' -f4,6 | sed -e 's/\"/,/g' >> $outfile"
cmd="zless $dir/../scenarios/mount-alexander-shire/castlemaine-region/population-archetypes.xml.gz | grep "Geographical.Coordinate" | cut -f2 -d'[' | cut -f1 -d']' | sed -e 's/ //g' >> $outfile"

echo $cmd && eval $cmd
