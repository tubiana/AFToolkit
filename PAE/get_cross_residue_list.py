import pandas as pd
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def parseArg():
    """
    
    """
    arguments = argparse.ArgumentParser(description="Tool to extract inter-complex distance information from AF2 models output jsons ")
    
    arguments.add_argument('-f', "--file", help="Json file", required=True, type=str)
    arguments.add_argument('-d', '--distance', help="distance threeshold", default=None, type=int)
    arguments.add_argument('-s', '--separation', help="Position of the NEXT sequence.", default=None, type=int)
    arguments.add_argument('-o', '--output', help="topfile", default=None, type=str)
  
    args = vars(arguments.parse_args())

    return args


def read_json(file):
    f = open(file)
    data = json.load(f)
    return data

def create_and_fill_matrix(data):
    ndata = len((data[0]["residue1"]))
    numres = int(np.sqrt(ndata))
    matrix = np.zeros((numres,numres), dtype=np.float32)

    for i in range(ndata):
        res1 = data[0]["residue1"][i] -1 # -1 to be 0 based
        res2 = data[0]["residue2"][i] -1 # -1 to be 0 based
        value = data[0]["distance"][i]
        matrix[res1,res2] = value
    return(matrix)

def show_matrix(matrix):
    sns.heatmap(matrix)
    plt.show()


def search_inter_complex(matrix, separation, distance):
    underDistance = np.argwhere(matrix <= distance)

    separation0 = separation-1 # (0based)
    keep = []
    for pair in underDistance:
        res1 = pair[0]
        res2=pair[1]
        if res1 <= separation0 and res2 > separation0:
            keep.append((res1+1,res2+1,matrix[res1,res2])) #+1 to be on 1 based index

    return(keep)


def save_data(data, matrix, output, separation):
    basename = Path(output).stem

    #Basic figure
    fig,ax = plt.subplots(figsize=(11,8))
    g = sns.heatmap(matrix, ax=ax)
    g.set(title=basename)
    plt.savefig(f"{basename}.png", dpi=300)

    #data in csv format.
    df = pd.DataFrame(data, columns=["Residue A","Residue B", "PAE"])
    #reset residue2 index
    df["Residue B (index corrected)"] = df["Residue B"] % separation
    df.to_csv(f"{basename}.csv", index=False, sep=";")



if __name__ == '__main__':
    args = parseArg()
    file = args["file"]
    distance = args["distance"]
    separation = args["separation"] 
    output = args["output"]

    data = read_json(file)
    matrix = create_and_fill_matrix(data)
    crossContact = search_inter_complex(matrix, separation, distance)
    save_data(crossContact,matrix, output, separation)




