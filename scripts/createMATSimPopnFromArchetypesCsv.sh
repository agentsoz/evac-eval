#!/usr/bin/env bash
INCSV=../scenarios/mount-alexander-shire/castlemaine-region/population-archetypes.csv.gz
OUTXML=./population-archetypes.xml
EESDIR=../bin/ees-2.1.1-SNAPSHOT

Rscript --vanilla './create-archetypes-csv.R'

#CMD="unzip -o ../../../target/ees-2.1.1-SNAPSHOT-release.zip -d ."
#echo $CMD; eval $CMD

CMD="java -cp $EESDIR/libs/*:$EESDIR/ees-2.1.1-SNAPSHOT.jar \
  io.github.agentsoz.ees.util.SynthpopToMatsim \
  --incsv $INCSV \
  --outxml $OUTXML \
  --incrs \"EPSG:4326\" \
  --outcrs \"EPSG:28355\" \
  "
echo $CMD; eval $CMD

CMD="gzip -f -9 $OUTXML"
echo $CMD; eval $CMD
