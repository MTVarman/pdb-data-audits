# Undeclared substitutions in PDB entries do not track source organism; reference curation does

Michael T. Varman · [ORCID 0009-0001-9453-0226](https://orcid.org/0009-0001-9453-0226)

`scan.py` aligns each polymer entity to its SIFTS UniProt reference and, when that record is unreviewed, to a reviewed member of its UniRef cluster, then counts substitutions the entry does not declare. `entries.py` and `entries2.py` build the entry lists the scan reads.

Note (all versions): [10.5281/zenodo.24194208](https://doi.org/10.5281/zenodo.24194208).

Scripts and tables at or under 1 MiB are copied byte-for-byte from `orgscan_2026-10-28_full.zip` on Zenodo record 24194209, the latest version of the note. The larger per-organism tables stay in the zip.
