# HRTF Explorer — EarSim vs BTE vs P0006

Interactive comparison of KEMAR **EarSim** (blocked ear canal), **BTE** (Meta Quest 3 + Cochlear Nucleus 5 mic ports), and SONICOM **P0006** measured human HRTF (reference).

- EarSim / BTE / P0006: Windowed with ITD, 48 kHz  
- P0006: SONICOM measured human HRTF (benchmark)  
- 793 matched directions; direction overlay, slice spectrograms, ITD/ILD, binaural preview  
- UI adapted from [yoyolicoris/hrtf-explorer](https://yoyolicoris.github.io/hrtf-explorer/)

## Live site

**https://ludovic-pirard.github.io/hrtf-explorer-earsim-bte-p0006/**

## Local

```bash
python3 -m http.server 8766
# open http://127.0.0.1:8766/
```

## Regenerate data

```bash
python3 export_explorer_data.py
```
