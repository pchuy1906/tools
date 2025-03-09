#!/bin/bash
scontrol show job=$1 | grep WorkDir  | sed -e 's/=/ /g'
scontrol show job=$1 | grep JobId
scontrol show job=$1 | grep Command  | sed -e 's/=/ /g'
scontrol show job=$1 | grep SubmitTime
