#!/bin/bash

inputFile=$1
runid=$2

sed 's, ,_,g' -i $inputFile
sed "s,runid=$runid,_,g" -i $inputFile
sed 's/_ch.*//' -i $inputFile
sed 's,___,_,g' -i $inputFile
