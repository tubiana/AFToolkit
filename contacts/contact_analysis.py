from glob import glob
import pandas as pd
import mdtraj as md
from collections import defaultdict
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a mega plot of all PAE matrices from AlphaFold2 predictions")
    parser.add_argument("--Folder", "-f", help="Path to the directory containing the predictions", required=True)
    parser.add_argument("--selection1", "-s1", help="Selection 1",  default="chainid 0")
    parser.add_argument("--selection2", "-s2", help="Selection 2", default="chainid 1")
    parser.add_argument("--add1", "-a1", help="Add to the residue number of selection 1",  default=0)
    parser.add_argument("--add2", "-a2", help="Add to the residue number of selection 2", default=0)
    parser.add_argument("--cutoff", "-c", help="Cutoff for the contacts", default=0.5)
    parser.add_argument("--minimumContacts", "-m", help="Minimum number of contacts to display", default=5)
    parser.add_argument("--outputname", "-o", help="Output name for the plot", default="contacts.png")
    parser.add_argument("--xlabel", "-xl", help="X axis label", default="chain B")
    parser.add_argument("--ylabel", "-yl", help="Y axis label", default="chain A")
    return parser.parse_args()

def get_residue_label(traj, atom_index, add=0):
    resname = traj.topology.atom(atom_index).residue.name
    resid = traj.topology.atom(atom_index).residue.resSeq + add
    return f"{resname} {resid}"

def compute_contacts(pdbs, selection1, selection2, add1, add2, cutoff=0.4):
    traj = md.load(pdbs, top=pdbs[0])

    chainA = traj.topology.select(selection1)
    chainB = traj.topology.select(selection2)

    chainAB = np.concatenate((chainA,chainB))

    contactsB = md.compute_neighbors(traj, cutoff, query_indices=chainA, haystack_indices=chainB)
    contactsA = md.compute_neighbors(traj, cutoff, query_indices=chainB, haystack_indices=chainA)


    contacts_linearB = [item for sublist in contactsB for item in sublist]
    contacts_linearA = [item for sublist in contactsA for item in sublist]
    contacts_atoms = np.unique(np.concatenate((contacts_linearB, contacts_linearA)))

    subtraj = traj.atom_slice(contacts_atoms)

    subtraj.save_pdb("subtraj.pdb")

    #get the number of chains 
    nchains = subtraj.top.n_chains

    chain1selsubtraj = " or ".join([f"chainid {i}" for i in range(nchains-1)])
    chain2selsubtraj = f"chainid {nchains-1}"

    chain1_atoms = subtraj.topology.select(chain1selsubtraj)
    chain2_atoms = subtraj.topology.select(chain2selsubtraj)
    #parwise combinations of 2 chains
    chain1_chain2 = np.array(np.meshgrid(chain1_atoms, chain2_atoms)).T.reshape(-1, 2)
    # all_distances = md.compute_contacts(subtraj, contacts='all')
    all_distances = md.compute_distances(subtraj, atom_pairs=chain1_chain2)

    nframes = all_distances.shape[0]
    npair = all_distances.shape[1]

    # Create a dictionary to store the shortest distances between residues for each frame
    shortest_distances = defaultdict(lambda: defaultdict(lambda: float('inf')))

    for i in range(nframes):
        for j in range(npair):
            dist = all_distances[i, j]
            pair = chain1_chain2[j]
            res1 = get_residue_label(subtraj, pair[0], add=add1)
            res2 = get_residue_label(subtraj, pair[1], add=add2)
            
            # Check if the current distance is shorter than the stored shortest distance for the current frame
            if dist < shortest_distances[i][(res1, res2)]:
                shortest_distances[i][(res1, res2)] = dist


    # Convert the dictionary to a DataFrame with one column per frame
    df_shortest_distances = pd.DataFrame(shortest_distances)

    # set multilevel index names to "res1" and "res2"
    df_shortest_distances.index.names = ["res1", "res2"]

    return df_shortest_distances

def plot_contacts(df, cutoff, outputname="contacts.png", minimum_contacts=3, xaxis_label="selection 1", yaxis_label="selection 2", title="", vmax=None):

    ncol = len(df.columns)
    def count_values_below_threshold(row, threshold=0.4):
        return (row < threshold).sum()
    count_table  = df.apply(lambda x: (x < cutoff).sum(), axis=1).unstack(fill_value=0)

    # Trier les colonnes
    numeric_part_columns = count_table.columns.to_series().str.extract('(\d+)').astype(int)
    sorted_columns = numeric_part_columns[0].argsort()
    count_table = count_table.iloc[:, sorted_columns]

    # Trier l'index
    numeric_part_index = count_table.index.to_series().str.extract('(\d+)').astype(int)
    sorted_index = numeric_part_index[0].argsort()
    count_table = count_table.iloc[sorted_index]

    #keep only the values where the number of contacts is > 4
    count_table = count_table[count_table >= minimum_contacts].dropna(how='all', axis=0).dropna(how='all', axis=1)
    

    fig, ax = plt.subplots(figsize=(12,10))
    if vmax=="col":
        vmax = ncol
    else:
        vmax = count_table.max().max()
    g = sns.heatmap(count_table, cmap="Blues", annot=False,xticklabels=True, yticklabels=True, ax=ax, vmin=0, vmax=vmax)
    #Set xticklabels to display EVERY LABELS
    
    g.set_xticks(np.arange(0.5, len(count_table.columns), 1))
    g.set_xlabel(xaxis_label)
    g.set_ylabel(yaxis_label)
    #I need all xtick labels on the xaxis, turn by 45°
    g.set_xticklabels(g.get_xticklabels(), rotation=90, horizontalalignment='center', fontsize=8)
    #reduce de size of the label 
    
    g.set_title(title)
    plt.tight_layout()
    g.figure.savefig(outputname)


    return(g)


if __name__=="__main__":

    args = parse_args()
    folder = args.Folder
    os.chdir(folder)
    pdbs = glob(r"*.pdb")
    pdbs.remove("subtraj.pdb")
    cutoff = 0.5
    selection1 = args.selection1
    selection2 = args.selection2
    add1 = int(args.add1)
    add2 = int(args.add2)
    cutoff = float(args.cutoff)
    minimum_contacts = int(args.minimumContacts)
    df_shortest_distances = compute_contacts(pdbs, selection1, selection2, add1, add2, cutoff)
    

    outputname = args.outputname
    xaxis_label = args.xlabel
    yaxis_label = args.ylabel
    plot_contacts(df_shortest_distances, 
              cutoff,
              xaxis_label=xaxis_label,
              yaxis_label=yaxis_label, 
              minimum_contacts=minimum_contacts,
              title=f"Number of contacts among all models (minimum contact = {minimum_contacts})",
              outputname=outputname,)

