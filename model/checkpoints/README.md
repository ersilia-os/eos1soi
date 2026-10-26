# Model pretrained parameters

These models have been trained using the [LazyQSAR package](https://github.com/ersilia-os/lazy-qsar) v3.6.3, in `slow` mode. The ensemble selected the CheMeleon and CLAMP descriptors.
Upon 5-fold crossvalidation on train test splits (80-20%), we obtain the following performance:


| **Model**       | **Data** | **Frac Actives (%)** | **AUROC ± STDev** | **AUPRC ± STDev** |
|-----------------|----------|------------------|-------------------|-------------------|
| antimicrobial   | 5584   | 0.66   |  0.959 ± 0.050 | 0.599 ± 0.060 |
