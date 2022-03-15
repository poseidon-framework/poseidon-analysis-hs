#!/usr/bin/env bash

REPO=~/dev/poseidon-framework/published_data

stack run xerxes -- ras --maxSnps 100000 -k 100 -d $REPO --popConfigFile popConfig.yml -f testOut