#!/usr/bin/env bash

for SIGNAL in 4 5; do
	for COVARIATE_TYPE in dmc; do
		    for KNOCKOFF_TYPE in ${COVARIATE_TYPE} gaussian; do
 		 	    Rscript main.R $SIGNAL $COVARIATE_TYPE $KNOCKOFF_TYPE coefdiff configs/config_desktop.yml
 		 	done
	done
done
