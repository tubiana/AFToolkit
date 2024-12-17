import pandas as pd
import numpy as np
import json
import glob
import os
from tqdm.notebook import tqdm
from af2_analysis import analysis
from af2_analysis import docking

import matplotlib.pyplot as plt

class MurailWrapper:
    def __init__(self, workdir):
        self.workdir = workdir
        self.data = None

    def plot_single_PAE(self, index, ax=None, cmap='bwr'):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(4, 4))

        json_file = self.data.df["json"][index]
        with open(json_file) as f:
            json_data = json.load(f)

        query = self.data.df.iloc[index]["query"]

        borders = self.data.chain_length[query]
        res_max = sum(borders)

        PAE = json_data["pae"]

        ax.imshow(PAE, cmap=cmap,
                  vmin=0.0,
                  interpolation='nearest',
                  vmax=30.0,)

        ax.hlines(
            np.cumsum(borders[:-1]) - 0.5,
            xmin=-0.5,
            xmax=res_max,
            colors="black",
        )

        ax.vlines(
            np.cumsum(borders[:-1]) - 0.5,
            ymin=-0.5,
            ymax=res_max,
            colors="black",
        )

        ax.set_xlim(-0.5, res_max - 0.5)
        ax.set_ylim(res_max - 0.5, -0.5)

        modelNumber = self.data.df["model"][index]
        ax.set_title(f"Rank {modelNumber}")
        return ax

    def save_all_PAE(self, save=True):
        fig, axes = plt.subplots(1, 5, figsize=(20, 4))
        for i in range(len(self.data.df)):
            self.plot_single_PAE(i, cmap='bwr', ax=axes[i])

        plt.tight_layout()

        if save:
            plt.savefig(f"{self.data.dir}/PAE.png", dpi=300)

    def save_plddt(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        for i in range(len(self.data.df)):
            query = self.data.df.iloc[i]["query"]

            borders = self.data.chain_length[query]
            res_max = sum(borders)

            plddt = self.data.get_plddt(i)
            ax.plot(plddt, label=f"Model {i}", linewidth=0.5)

            ax.vlines(
                np.cumsum(borders[:-1]) - 0.5,
                ymin=-0.5,
                ymax=100,
                colors="black",
            )

            ax.set_ylim(0, 100)

            # add legend for every plot, on the side
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            plt.tight_layout()
            plt.savefig(f"{self.data.dir}/plddt.png", dpi=300)

    def get_all_data(self):
        os.chdir(self.workdir)

        folders = [f for f in glob.glob("*") if os.path.isdir(f)]

        print(folders)

        errors = []
        list_of_scores = []
        models_list = {}

        for model in tqdm(folders):

            try:
                self.data = af2.Data(model + "/", verbose=False)
            except:
                errors.append(model)
                continue

            print(model)

            try:
                analysis.pdockq(self.data, verbose=False)
            except:
                print("pdockQ failed")
            try:
                analysis.pdockq2(self.data, verbose=False)
            except:
                print("pdockq2 failed")

            try:
                analysis.mpdockq(self.data, verbose=False)
            except:
                print("mpdockq failed")

            try:
                analysis.inter_chain_pae(self.data, verbose=False)
            except:
                print("inter_chain_pae failed")

            try:
                analysis.LIS_matrix(self.data, verbose=False)
            except:
                print("LIS_matrix failed")

            try:
                docking.LIS_pep(self.data, verbose=False)
            except:
                print("LIS PEP failed")

            for i in range(5):

                json_file = self.data.df["json"][i]
                with open(json_file) as f:
                    json_data = json.load(f)
                pae_array = json_data["pae"]

                pae_mean = np.mean(pae_array)
                plddt = self.data.get_plddt(i)
                plddt_mean = np.mean(plddt)

                self.data.df.loc[i, "pae_mean"] = pae_mean
                self.data.df.loc[i, "plddt_mean"] = plddt_mean

            name = str(self.data.df.iloc[0]["query"])
            models_list[name] = self.data

            list_of_scores.append(self.data.df)

        results = pd.concat(list_of_scores)

        return results
