This test data is based on the XChem lb27995-1 visit.

To run this, setup the environment as described in the main USER_GUIDE.md.
Then:

mkdir test-data/outputs/upload_1
python -m xchemalign.collator -c test-data/config_1.yaml
python -m xchemalign.aligner -d test-data/outputs/upload_1