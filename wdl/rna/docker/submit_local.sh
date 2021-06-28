curl -X POST\
 -v "localhost:8000/api/workflows/v1" \
 -F workflowSource=@celescope/wdl/rna.wdl \
 -F workflowInputs=@celescope/wdl/rna/docker/rna_input_local.json\
 -F workflowDependencies=@celescope/wdl/wdl.zip