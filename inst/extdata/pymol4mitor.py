from  pymol import cmd
import pandas as pd
import numpy as np
import sys

# get the protein name in the argument

df_path = sys.argv[2]
str_path =  sys.argv[3]
comp = sys.argv[4]

match comp:
    case "I":
        mt_chains={
            's': 'palegreen',
            'i': 'lightblue',
            'j': 'paleyellow',
            'r': 'lightpink',
            'k': 'palecyan',
            'l': 'lightorange',
            'm': 'bluewhite'
        }
        sel="chain 6+7"

    case "III":
        mt_chains={"J": "white"}
        sel="chain 3+4"
    case "IV":
        mt_chains={"A": "white",
           "B": "lightpink",
            "C": "lightorange"}
        sel="chain 1+2"
    case "V":
        mt_chains={"Q": "white", "N": "lightpink"}
        sel="chain p"

cmd.load(str_path)
cmd.spectrum()
cmd.select("nuc_sub", "not chain " + "+".join(mt_chains.keys()))
cmd.set("cartoon_transparency", 0.8, "nuc_sub")
# cmd.select("hetatm", "hetatm")
cmd.set("sphere_transparency", 0.8, "hetatm")
cmd.set("stick_transparency", 0.8, "hetatm")
cmd.select("membrane", sel)
cmd.hide("sphere", "membrane")
cmd.show("nb_spheres", "membrane")

for chain in mt_chains:
    cmd.color(mt_chains[chain], "chain " + chain)

df=pd.read_csv(df_path)

cmd.rotate("x", -90)
cmd.deselect()


for chain in np.unique(df["ChainID"]):
        ress="+".join(str(v) for v in df["Pos"][df["ChainID"]==chain].to_list())
        cmd.color("red", "chain " + chain + " and resi " + ress)




