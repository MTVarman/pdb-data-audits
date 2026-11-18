# Protein–ligand bond records in the Protein Data Bank follow a distance table in its processing software, the bonds outside it are listed as close contacts, and the clash check gives boron a 0.20 Å radius

Michael T. Varman · [ORCID 0009-0001-9453-0226](https://orcid.org/0009-0001-9453-0226)

The scripts compare protein–ligand distances with the MAXIT processing windows, separate a recorded bond from a close-contact row, and check the boron clash radius. `verify_counts.py` recounts the note from the tables in `data/` and `tables/`.

Note (all versions): [10.5281/zenodo.24295881](https://doi.org/10.5281/zenodo.24295881).

Scripts and uncompressed tables at or under 1 MiB are copied byte-for-byte from `Varman_PDB_bond_records_v1.1_supplementary.zip` on Zenodo record 24295951, the latest version of the note. The gzipped archive scans stay in the zip.
