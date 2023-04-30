import numpy as np
import pandas as pd

#t and tau are indices, could be necessary some conversion
def FC(t, tau, data):
  #Define the time window
  data_window = data[int((t-tau/2)):int((t+tau/2))]

  #.corr() method returns a Pandas DataFrame, .values makes it a Numpy Array (convenient for later consistency)
  #Lost information about names of DafaFrame columns, if any 
  FC_dataframe = data_window.corr(method="pearson")
  FC_array = FC_dataframe.values

  return FC_array


def dFC(FC_1, FC_2):
  #Get the upper triangular matrix of FC
  UpperTri_1 = FC_1[np.triu_indices(FC_1.shape[1], k=1)]
  UpperTri_2 = FC_2[np.triu_indices(FC_2.shape[1], k=1)]

  #Pearson correlation
  rho = np.corrcoef(UpperTri_1, UpperTri_2)[0, 1]
  return rho


#Note that it is being used the same tau both for FC window and speed
def v_dFC(t, tau, data):
  return  1 - dFC(FC(t, tau, data), FC(t+tau, tau, data))


def dFC_stream(data, tau):
  stream = np.zeros(shape=(1, data.shape[1], data.shape[1]))
  
  n = data.shape[0]
  #Iterate from 0 to n-tau (included) with step=tau
  for i in range(0, n-tau+1, tau):
    stream = np.append(stream, [FC(i+tau/2, tau, data)], axis=0)
  
  return stream[1:]


def v_stream(dFC_stream):
  v_array = np.array([])

  for i in range(dFC_stream.shape[0]-1):
    v = 1 - dFC(dFC_stream[i], dFC_stream[i+1])
    v_array = np.append(v_array, v)
  
  return v_array


def dFC_matrix(dFC_stream):
  #Number of FC in dFC_stream and number of elements in UpperTri matrices
  num_FC = dFC_stream.shape[0]
  n = dFC_stream.shape[1]
  num_UpperTri = int(n*(n-1)/2)
  
  #Save UpperTri vectors in a Pandas dataframe
  UpperTri_vectors = pd.DataFrame(np.zeros((num_UpperTri, num_FC)))

  for i in range(num_FC):
    FC_i = dFC_stream[i]
    UpperTri_vectors[i] = FC_i[np.triu_indices(n, k=1)]
  
  #Get correlation matrix
  dFC_dataframe = UpperTri_vectors.corr(method="pearson")
  matrix = dFC_dataframe.values

  return matrix


#Old version of dFC_matrix which can be useful for debugging
def dFC_matrix_2(dFC_stream):
  #Size of dFC_matrix: number of FC in dFC_stream
  num_FC = dFC_stream.shape[0]
  matrix = np.zeros(shape=(num_FC, num_FC))

  #Fill lower triangula matrix
  for i in range(num_FC):
    for j in range(i+1):
      matrix[i, j] = dFC(dFC_stream[i], dFC_stream[j])
  
  #Get the symmetric matrix
  matrix = matrix + matrix.T - np.diag(matrix)

  return matrix


def pooled_v_stream(data, tau_min, tau_max):
  v_array = np.array([])

  for tau in np.arange(start=tau_min, stop=tau_max+1):
    #Get dFC stream and v distribution at fixed tau
    stream_dFC = dFC_stream(data, tau)
    stream_v = v_stream(stream_dFC)

    v_array = np.append(v_array, stream_v)
  
  return v_array


def tSNE_evolution(dFC_stream, TSNE, n):
  #Create empty matrix of UpperTri: shape = (number of FC in dFC stream, len(UpperTri))
  dim_UpperTri = len(dFC_stream[0][np.triu_indices(dFC_stream[0].shape[1], k=1)])
  UpperTri_matrix = np.zeros(shape=(dFC_stream.shape[0], dim_UpperTri))

  #Fill UpperTri_matrix  
  for i in range(dFC_stream.shape[0]):
    #Get FC(t_i) and UpperTri(FC(t_i)) vector
    FC = dFC_stream[i]
    UpperTri = FC[np.triu_indices(FC.shape[1], k=1)]
    
    UpperTri_matrix[i] = UpperTri

  #Points in nD space
  embedded_points = TSNE(n_components=n, perplexity=30.0, early_exaggeration=4.0, method="exact").fit_transform(UpperTri_matrix)
  
  #Isolate feauture
  out = [embedded_points[:, i] for i in range(n)]
  
  return out