#!/usr/bin/env bash

IMAGE_NAME=lr-cartographer

########################################
# Prepare environment: 

## Build the LRMA python package and copy it here:
#cd ../../packages/python/lrma
#python3 setup.py sdist
#cd -
#lrma_package=$(ls -rt ../../packages/python/lrma/dist/*.gz | tail -n1)
#cp -va ${lrma_package} .

########################################
# Build the container:
time docker build -t ${IMAGE_NAME} .

########################################
# Clean local files:

########################################
# Instructions for user:

echo "You must now tag these images as the next version and \`latest\` then upload them to us.gcr.io/broad-dsp-lrma/${IMAGE_NAME}:"
echo ""
echo "  docker tag ${IMAGE_NAME}:latest us.gcr.io/broad-dsp-lrma/${IMAGE_NAME}:latest"
echo "  docker tag ${IMAGE_NAME}:latest us.gcr.io/broad-dsp-lrma/${IMAGE_NAME}:<NEXT_VERSION>"
echo ""
echo "  docker push us.gcr.io/broad-dsp-lrma/${IMAGE_NAME}:latest"
echo "  docker push us.gcr.io/broad-dsp-lrma/${IMAGE_NAME}:<NEXT_VERSION>"

