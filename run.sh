#!/bin/bash

cd wdl/; rm -f lr_wdls.zip; zip -r lr_wdls.zip * > /dev/null; cd ..

/Users/ewan/Documents/cromshell/cromshell submit wdl/PBFlowcell.wdl \
    PBFlowcell.json \
    resources/workflow_options/default.json \
    wdl/lr_wdls.zip
