import pandas as pd 
import numpy as np

# Function that will perform DESeq2-like transformation on dataframe
# Features as rows and samples as columns 
def deseq2_norm(df):
    norm = (df.apply(np.log))                              # take a log of df
    norm['mean'] = norm.mean(numeric_only=True, axis=1)    # mean as pseudoreference
    norm = norm.loc[~norm['mean'].isin([np.inf, -np.inf])] # get rid of '-inf'
    norm = norm[norm.columns] - norm[['mean']].values      # substract mean
    norm = norm.drop('mean',axis=1)                        # drop mean
    medians = norm.median()                                # median values by samples
    factors = np.exp(medians)                              # convert to expotentials
    df = df/factors                                        # divide original counts by scaling factors
    return(df)
