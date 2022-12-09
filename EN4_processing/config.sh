#!/bin/bash

## Set the machine to be used. Pick one.
export MACHINE="LOTUS"  # resource on JASMIN
#export MACHINE="SPICE"  # resource at MO


## SETTINGS FOR JASMIN LOTUS PROCESSING
if [ $MACHINE = "LOTUS" ]; then
  source lotus_config.sh

## SETTINGS FOR MET OFFICE SPICE PROCESSING
elif [ $MACHINE = "SPICE" ]; then
  source spice_config.sh
