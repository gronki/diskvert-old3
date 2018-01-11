#!/usr/bin/env bash

rm -rfv par/* datac/* dataw/*
python par.py
cat jobs.lst | parallel bash job.sh
python plot.py &
