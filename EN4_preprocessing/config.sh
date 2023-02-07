#!/bin/bash

## Set the machine to be used. Pick one.
#export MACHINE="LOTUS"  # resource on JASMIN
export MACHINE="SPICE"  # resource at MO
#export MACHINE="LIVMAZ"  # local resource for testing


## SETTINGS FOR JASMIN LOTUS PROCESSING
if [ $MACHINE = "LOTUS" ]; then
  source lotus_config.sh

## SETTINGS FOR MET OFFICE SPICE PROCESSING
elif [ $MACHINE = "SPICE" ]; then
  source spice_config.sh

elif [ $MACHINE = "LIVMAZ" ]; then
  source livmaz_config.sh
fi
