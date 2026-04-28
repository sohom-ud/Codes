import numpy as np
import pandas as pd

def MVA(B, start_time, end_time):
    
    B_int = B.loc[start_time:end_time]
    
    M = np.cov(B_int.T, bias=True)
    
    eigvals, eigvecs = np.linalg.eig(M)
    
    eigvals = eigvals[::-1]
    eigvecs = eigvecs[:, ::-1]
    
    return eigvals, eigvecs

def rotate_timeseries(df, L, M, N):
    
    # Transformation matrix
    
    A = np.array([L, M, N])
    
    # Ainv = np.linalg.inv(A)
    
    df_rot = pd.DataFrame(df.values @ A.T, index=df.index, columns=['L', 'M', 'N'])
    
    return df_rot