#!/usr/bin/env bash

REPO=~/dev/poseidon-framework/published_data

# stack run xerxes -- ras --maxSnps 100000 -d $REPO --popConfigFile popConfig.yml -f testOut
stack run xerxes -- csfs --maxSnps 100000 -d $REPO --popConfigFile popConfig.yml -f testOutCSFS