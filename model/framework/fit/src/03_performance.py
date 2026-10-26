import json
import os

import numpy as np
import stylia
from sklearn.metrics import auc, roc_curve

# Format: slide | Style: ersilia — change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

root = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.abspath(os.path.join(root, "..", "results"))

TITLE = "Antimicrobial potential in non-growing bacteria"


def read_report(name):
    with open(os.path.join(results_dir, f"report_{name}.json"), "r") as f:
        return json.load(f)


def plot_roc_curve(ax, results):
    colors = stylia.ErsiliaColors()
    mean_fpr = np.linspace(0, 1, 100)
    tprs = []
    aucs = []
    for fold in results.keys():
        y_hat = results[fold]["y_hat"]
        y_true = results[fold]["y_true"]
        fpr, tpr, _ = roc_curve(y_true, y_hat)
        aucs += [auc(fpr, tpr)]
        tpr_interp = np.interp(mean_fpr, fpr, tpr)
        tpr_interp[0] = 0.0
        tprs += [tpr_interp]
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.fill_between(
        mean_fpr,
        np.maximum(mean_tpr - std_tpr, 0),
        np.minimum(mean_tpr + std_tpr, 1),
        color=colors.get("gray", lighten=0.5),
        label=r"$\pm$ 1 std. dev.",
    )
    ax.plot(mean_fpr, mean_tpr, color=colors.plum, label="Mean ROC")
    ax.plot([0, 1], [0, 1], linestyle="--", color=colors.gray)
    ax.legend(loc="lower right")
    stylia.label(
        ax,
        xlabel="False positive rate",
        ylabel="True positive rate",
        title=r"AUROC {0:.2f} $\pm$ {1:.2f}".format(mean_auc, std_auc),
    )


def plot_yhat(ax, results, max_points=100):
    colors = stylia.ErsiliaColors()
    rng = np.random.default_rng(42)
    # Pool the five out-of-fold predictions: one fold holds only 7 or 8 actives.
    y_hat = np.array([v for fold in results for v in results[fold]["y_hat"]])
    y_true = np.array([v for fold in results for v in results[fold]["y_true"]])
    for label, color in [(0, colors.gray), (1, colors.pink)]:
        values = y_hat[y_true == label]
        rng.shuffle(values)
        values = values[:max_points]
        jitter = rng.uniform(-0.1, 0.1, size=len(values))
        ax.scatter(label + jitter, values, color=color, alpha=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Neg", "Pos"])
    ax.set_xlim(-0.5, 1.5)
    # Scores are calibrated to the 0.66% prevalence of the training set, so they span
    # orders of magnitude near zero; a linear axis would flatten both classes onto it.
    ax.set_yscale("log")
    stylia.label(ax, xlabel="", ylabel="Score", title=TITLE)


results = read_report("crossval")

fig, axs = stylia.create_figure(1, 2)
plot_roc_curve(axs.next(), results)
plot_yhat(axs.next(), results)
stylia.save_figure(os.path.join(results_dir, "performance_5fold.png"))

roc_aucs = [v["roc_auc"] for v in results.values()]
print("MEAN AUC: ", np.mean(roc_aucs), np.std(roc_aucs))
