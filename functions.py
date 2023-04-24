import numpy as np
import pandas as pd

#t and tau are indices, could be necessary some conversion
def FC(t, tau, data):
  #Define the time window
  data_window = data[int((t-tau/2)):int((t+tau/2))]

  #.corr() method returns a Pandas DataFrame, .values makes it a Numpy Array (convenient for later consistency)
  #Lost information about names of DafaFrame columns 
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
  
  #range(100) needs to be updated
  for i in range(100):
    stream = np.append(stream, [FC(i+tau/2, tau, data)], axis=0)
  
  return stream[1:]


def v_stream(dFC_stream):
  v_array = np.array([])

  for i in range(dFC_stream.shape[0]-1):
    v = 1 - dFC(dFC_stream[i], dFC_stream[i+1])
    v_array = np.append(v_array, v)
  
  return v_array


def dFC_matrix(dFC_stream):
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

  for tau in range(tau_max - tau_min + 1):
    #Get dFC stream and v distribution at fixed tau
    stream_dFC = dFC_stream(data, tau)
    stream_v = v_stream(stream_dFC)

    v_array = np.append(v_array, stream_v)
  
  return v_array


def tSNE_evolution(dFC_stream, TSNE):
  #Create empty matrix of UpperTri: shape = (number of FC in dFC stream, len(UpperTri))
  dim_UpperTri = len(dFC_stream[0][np.triu_indices(dFC_stream[0].shape[1], k=1)])
  UpperTri_matrix = np.zeros(shape=(dFC_stream.shape[0], dim_UpperTri))

  #Fill UpperTri_matrix  
  for i in range(dFC_stream.shape[0]):
    #Get FC(t_i) and UpperTri(FC(t_i)) vector
    FC = dFC_stream[i]
    UpperTri = FC[np.triu_indices(FC.shape[1], k=1)]
    
    UpperTri_matrix[i] = UpperTri

  #Points in 2D space
  embedded_points = TSNE(n_components=2, perplexity=30.0, early_exaggeration=4.0, method="exact").fit_transform(UpperTri_matrix)
  
  #Isolate x, y
  feature1 = embedded_points[:, 0]
  feature2 = embedded_points[:, 1]
  
  return feature1, feature2