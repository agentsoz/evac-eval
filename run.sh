#!/usr/bin/env bash

SECONDS=0 # builtin BASh var to count seconds

EES_JAR=./bin/ees-2.1.1-SNAPSHOT/ees-2.1.1-SNAPSHOT.jar
EES_LIBS=./bin/ees-2.1.1-SNAPSHOT/libs/*

if [ ! -f $EES_JAR ]; then
    printf "\nJAR $EES_JAR not found.\nPlease download and install the EES distribution first as per ./bin/DOWNLOAD.md\n\n"
    exit 0
fi

if [ "$#" -ne 1 ]; then
    echo "usage: $0 CONFIG"
    echo "  CONFIG  path to scenario config XML, e.g., "
    echo "          ./scenarios/mount-alexander-shire/castlemaine-region/archetypes.xml"
    exit 0
fi
CONFIG=$1

CMD="java -Xms8g -Xmx8g \
  -cp $EES_LIBS:$EES_JAR \
  io.github.agentsoz.ees.Run \
  --config $CONFIG"

echo $CMD && eval $CMD

DURATION=$SECONDS
printf "\nFinished (in $(($DURATION / 60)) mins $(($DURATION % 60)) secs)\n\n"
