#!/usr/bin/env bash

parallel python ::: cooling.py maps.py plotpar.py instabil.py multi.py
