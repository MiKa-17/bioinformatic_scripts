#!/bin/bash

inputFile=$1

fileName=$(basename $inputFile .sam)
outputPath=$(dirname $inputFile)

samtools view $inputFile -bo $outputPath/$fileName.bam
samtools sort $outputPath/$fileName.bam -o $outputPath/$fileName.sorted.bam
samtools index $outputPath/$fileName.sorted.bam

