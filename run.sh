#!/bin/bash

cd wdl/; rm lr_wdls.zip; zip -r lr_wdls.zip * > /dev/null; cd ..

cromshell submit wdl/PBFlowcell.wdl \
    PBFlowcell.json \
    resources/workflow_options/default.json \
    wdl/lr_wdls.zip
