# Release checklist (preprint / submission)

The repository stays private until the preprint. At release:

1. Confirm no study data or identifiers are tracked: `git ls-files` should list only code, docs, tests and `src/edgeuq/data/isotonic_1mm.json`.
2. Fill the remaining dataset accessions in `docs/reproduce.md`.
3. Bump `version` in `pyproject.toml`, `src/edgeuq/__init__.py`, `CITATION.cff` and `.zenodo.json`.
4. Make the repository public (GitHub → Settings → General → Danger zone → Change visibility).
5. Make the GHCR package public (GitHub → Packages → edgeuq → Package settings).
6. Log in to Zenodo with GitHub, enable this repository under *GitHub* in Zenodo.
7. Tag and publish a GitHub release (`git tag v1.0.0 && git push origin v1.0.0`, then *Create release*). Zenodo archives it and mints a DOI; the container workflow publishes `ghcr.io/phindagijimana/edgeuq:1.0.0`.
8. Add the DOI to `CITATION.cff`, the README and the manuscripts.
