# Release checklist (preprint / submission)

**Publication split (2026-10-21):** `edgeuq` is the shared library for **Paper 1** (shift/jitter/rescan boundary) and **Paper 2** (`instability` / `calibrate-leftover`). **DKT-Connectome** may call these CLIs as an optional post-step and publish as a **software** note citing both papers — do not fold all three into one manuscript. See TrackTBI-Sub `dwi_pipeline/docs/papers/mapping_uhat/PUBLICATION_PLAN.md`.

The repository stays private until the preprint. At release:

1. Confirm no study data or identifiers are tracked: `git ls-files` should list only code, docs, tests and `src/edgeuq/data/isotonic_1mm.json`.
2. Fill the remaining dataset accessions in `docs/reproduce.md`.
3. Bump `version` in `pyproject.toml`, `src/edgeuq/__init__.py`, `CITATION.cff` and `.zenodo.json`.
4. Make the repository public (GitHub → Settings → General → Danger zone → Change visibility).
5. Make the GHCR package public (GitHub → Packages → edgeuq → Package settings).
6. Log in to Zenodo with GitHub, enable this repository under *GitHub* in Zenodo.
7. Tag and publish a GitHub release (`git tag v0.2.0 && git push origin v0.2.0`, then *Create release* from the tag). Zenodo archives it and mints a DOI; the container workflow publishes `ghcr.io/phindagijimana/edgeuq:0.2.0`.
8. Add the DOI to `CITATION.cff`, the README and the manuscripts.

## v0.2.0 ready checklist (lab)

- [x] Version bump `0.2.0` in `pyproject.toml`, `__init__.py`, `CITATION.cff`
- [x] `.zenodo.json` present
- [x] Interim Kirby leftover map in `src/edgeuq/data/`
- [x] `docs/CHANGELOG_0.2.0.md`
- [ ] Repo public + Zenodo GitHub enabled (browser)
- [ ] `git push origin v0.2.0` + GitHub Release
- [ ] Paste Zenodo DOI back into CITATION/README
