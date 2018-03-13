#!/usr/bin/env bash

parallel python ::: plotpar.py maps.py instabil.py cooling.py
