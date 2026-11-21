# Fit folder: non-growing-antimicrobial (eos1soi)

This folder holds the data, code and results needed to rebuild the model. It is stored with
[eosvc](https://github.com/ersilia-os/eosvc) in the **public** bucket; only this README is
tracked in git.

## Data processing

`src/00_dataanalysis.py` merges the screened compounds (`data/R4A_input.xlsx`, `data/Prestwick.csv`) with the actives (`data/actives.csv`) and standardises them into `data/train.csv`. `src/01_sim_search.py` then drops inactives with Tanimoto similarity above 0.90 to an active, giving `data/train_filtered.csv`, the training set.

## Model building

`src/01_fit_oof.py` builds one LazyQSAR 3.9.3 classifier per task, in slow mode. For each task
it runs five-fold stratified out-of-fold (OOF) cross-validation (seed 42), so every compound
gets one prediction from a model that never saw it. It writes the OOF predictions and metrics to
`results/`, then trains the final model on all compounds into `model/checkpoints/<task>/`.
`src/02_performance.py` plots the OOF ROC curve and the OOF rank of actives and inactives per
task, and collects all metrics in `results/oof_metrics.json`.

```bash
cd src
OMP_NUM_THREADS=1 python 01_fit_oof.py        # or one task: python 01_fit_oof.py model
python 02_performance.py
```

| Task | Output column | Training file | Compounds | Actives | OOF AUROC | Descriptors in final model |
|---|---|---|---|---|---|---|
| `model` | `nongrowing_antimicrobial_rank` | `data/train_filtered.csv` | 5,584 | 37 | 0.970 | chemeleon, clamp |
