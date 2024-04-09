#!/bin/bash

seq 0 23 1791 | xargs -n1 -P80 "$@"