#!/usr/bin/env bash

set -e

IMAGE_NAME=lr-jupyter
INTERACTIVE_SUFFIX=_interactive

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

# Now build a container we can run interactively so we can 
# create our HTML reports programmatically:
mv Dockerfile Dockerfile.bak

cat Dockerfile.bak | sed 's#CMD \["start-notebook.sh"\]#CMD \["/bin/bash"\]#g' > Dockerfile
time docker build -t ${IMAGE_NAME}${INTERACTIVE_SUFFIX} .

mv Dockerfile Dockerfile.bak
cat Dockerfile.bak | sed 's#CMD \["/bin/bash"\]#CMD \["start-notebook.sh"\]#g' > Dockerfile

rm Dockerfile.bak

echo ""
echo "Built two images:"
echo -e "\t${IMAGE_NAME}:latest"
echo -e "\t${IMAGE_NAME}${INTERACTIVE_SUFFIX}:latest"

########################################
# Clean local files:

########################################
# Instructions for user:

echo "You must now tag these images as the next version and \`latest\` then upload them to us.gcr.io/broad-dsp-lrma/lr-jupyter:"
echo ""
echo "  docker tag lr-jupyter:latest us.gcr.io/broad-dsp-lrma/lr-jupyter:latest"
echo "  docker tag lr-jupyter_interactive:latest us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:latest"
echo "  docker tag lr-jupyter:latest us.gcr.io/broad-dsp-lrma/lr-jupyter:<NEXT_VERSION>"
echo "  docker tag lr-jupyter_interactive:latest us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:<NEXT_VERSION>"
echo ""
echo "  docker push us.gcr.io/broad-dsp-lrma/lr-jupyter:latest"
echo "  docker push us.gcr.io/broad-dsp-lrma/lr-jupyter:<NEXT_VERSION>"
echo "  docker push us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:latest"
echo "  docker push us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:<NEXT_VERSION>"

