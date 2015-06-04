PDBProcessing
=============

Scripts to process the PDB mmCIF files and run the redundancy removal.

If running the processing, BLASTing and culling in one go (or locally) is not feasible, then the steps can be split up. For example,

1) Run the controller.py script with a dummy BLAST directory argument. This will run the processing, but prevent any BLASTing or culling (and will exit with an error).

2) Run the BLASTing by directly using the BLAST executables and the output from step 1.

3) Run the culling using the generateculledsubsets.py script and the outputs of the preceding steps.
