"""
Plotting utilities for pySCENIC AUCell results.

Author: Jung Hyun Park
"""

import matplotlib.pyplot as plt
import pandas as pd


def plot_regulon_activity(merged_df, regulon, groupby=None):
    """
    Plot AUCell activity for a single regulon.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Metadata + AUCell matrix
    regulon : str
        Regulon name
    groupby : str, optional
        Metadata column to group by
    """
    if regulon not in merged_df.columns:
        raise ValueError(f"{regulon} not found in DataFrame")

    if groupby is None:
        plt.hist(merged_df[regulon], bins=50)
        plt.xlabel("AUCell score")
        plt.ylabel("Cell count")
        plt.title(regulon)
    else:
        merged_df.boxplot(column=regulon, by=groupby, rot=90)
        plt.title(regulon)
        plt.suptitle("")

    plt.tight_layout()
    plt.show()


def plot_top_regulons(auc_df, n_top=10):
    """
    Plot mean activity of top regulons.

    Parameters
    ----------
    auc_df : pd.DataFrame
        AUCell matrix
    n_top : int
    """
    mean_activity = auc_df.mean(axis=0).sort_values(ascending=False)
    top = mean_activity.head(n_top)

    top.plot(kind="bar")
    plt.ylabel("Mean AUCell score")
    plt.title("Top regulons")
    plt.tight_layout()
    plt.show()

