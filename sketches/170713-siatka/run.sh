#!/usr/bin/env bash

ls *.par | sed s/.par// | parallel bash job.sh {}
echo DONE
