import numpy as np
import pandas as pd

#t and tau are indices (TRs), could be necessary some conversion
def FC(t, tau, data):
  #Define the time window
  data_window = data[int((t-tau/2)):int((t+tau/2))]

  #.corr() method returns a Pandas DataFrame, .values makes it a Numpy Array (convenient for later consistency)
  #Lost information about names of DafaFrame columns, if any 
  FC_dataframe = data_window.corr(method="pearson")
  FC_array = FC_dataframe.values

  return FC_array

#FC of the full scan session
full_FC = lambda df : df.corr("pearson").values

#Dynamic functional connectivity given two FC matrices
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

#Get the dFC stream
def dFC_stream(data, tau):
  stream = np.zeros(shape=(1, data.shape[1], data.shape[1]))
  
  n = data.shape[0]
  #Iterate from 0 to n-tau (included) with step=tau
  for i in range(0, n-tau+1, tau):
    stream = np.append(stream, [FC(i+tau/2, tau, data)], axis=0)
  
  return stream[1:]

#dFC speed from the dFC stream
def v_stream(dFC_stream):
  v_array = np.array([])

  for i in range(dFC_stream.shape[0]-1):
    v = 1 - dFC(dFC_stream[i], dFC_stream[i+1])
    v_array = np.append(v_array, v)
  
  return v_array

#dFC matrix from the dFC stream
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

#Old version of dFC_matrix which can be useful for debugging: equivalent but slower
def dFC_matrix_2(dFC_stream):
  #Size of dFC_matrix: number of FC in dFC_stream
  num_FC = dFC_stream.shape[0]
  matrix = np.zeros(shape=(num_FC, num_FC))

  #Fill lower triangula matrix
  for i in range(num_FC):
    for j in range(i+1):
      matrix[i, j] = dFC(dFC_stream[i], dFC_stream[j])
  
  #Get the symmetric matrix
  matrix = matrix + matrix.T - np.diag(np.diag(matrix))

  return matrix

#dFC speed from data by indices
#Can be used also using tau_min=tau_max to get dFC speed from data
def pooled_v_stream(data, tau_min, tau_max):
  v_array = np.array([])

  for tau in np.arange(start=tau_min, stop=tau_max+1):
    #Get dFC stream and v distribution at fixed tau
    stream_dFC = dFC_stream(data, tau)
    stream_v = v_stream(stream_dFC)

    v_array = np.append(v_array, stream_v)
  
  return v_array

#tSNE projection
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

#Use MINDy model and parameters to propagate signals
#hidden hyperparameters: TR and id_max
def propagation(W, alpha, D, sigma):
    id_max = 1200
    Xt = np.zeros((119, id_max))
    x1 = np.random.random((119))*1.5 - 1       #x1 random
    Xt[:, 0] = x1

    Psi = np.zeros((119, id_max))

    TR = 0.7
    b = 6.6667
    for i in range(id_max-1):

        t1 = (alpha**2 + (b*x1+0.5)**2)**0.5
        t2 = (alpha**2 + (b*x1-0.5)**2)**0.5
        psi = t1 - t2

        eps = np.random.normal(0, sigma, x1.size)
        x2 = x1 + (np.dot(W, psi) - D*x1)*TR + eps

        Psi[:, i+1] = psi

        Xt[:, i+1] = x2
        x1 = x2


    dfXt = pd.DataFrame(Xt)
    simul = dfXt.T
    return simul

#Get fluctuation strengths to perform DFA
def get_fluctuation(data, k):
    #Cumulative unbounded process
    integrated_data = np.cumsum(data) - np.mean(data)

    #Split into M segments of length k
    num_segments = len(integrated_data) // k
    segments = np.array_split(integrated_data[:num_segments * k], num_segments)

    #Local fit for each segment
    def fit_polynomial(segment, order=1):
        x = np.arange(len(segment))
        coeffs = np.polyfit(x, segment, order)
        fitted_values = np.polyval(coeffs, x)
        return fitted_values

    def fit_polynomials_to_segments(segments, order=1):
        fitted_segments = [fit_polynomial(segment, order) for segment in segments]
        return fitted_segments

    pol_order = 1   #Linear fit
    fitted_segments = fit_polynomials_to_segments(segments, pol_order)

    #Compute the fluctuation strength for each segment
    f_strengths = []
    for q in range(len(segments)):
        squared_fluctuations = (segments[q] - fitted_segments[q])**2
        f_strengths.append(np.mean(squared_fluctuations)**0.5)  
    
    #Return the fluctuation for the chosen k
    return np.mean(f_strengths)

#Actual Detrended fluctuation analysis (DFA)
def do_DFA(data):
    lengths = np.logspace(2, np.log(len(data)), base=np.e, num=20, dtype=int)   #spaced evenly on the log scale
    log_fluctuations = [np.log(get_fluctuation(data, k)) for k in lengths]
    log_lengths = np.log(lengths)

    #Lineared exponential fit
    A, C = np.polyfit(log_lengths, log_fluctuations, 1)
    return A

