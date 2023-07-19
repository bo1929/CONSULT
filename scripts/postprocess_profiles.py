import pandas as pd
import numpy as np
import sys

if __name__ == "__main__":
    for rank in ["species", "genus", "family", "class", "order", "phylum", "kingdom"]:
        df = pd.read_csv(sys.argv[1] + "-" + rank, sep="\t")
        df["FRACTION_TOTAL"] = df["FRACTION_TOTAL"] / df["FRACTION_TOTAL"].sum()
        df.to_csv(sys.argv[1] + "-" + rank + "-postp", sep="\t", index=False, float_format="%.15f")
