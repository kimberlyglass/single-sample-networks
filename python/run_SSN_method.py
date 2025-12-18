import numpy as np


####################################
###########LIONESS##################
####################################
def LIONESS(LocExp):

    # determine size of input data
    NumNode, NumSamp = LocExp.shape

    # create filter to extract upper diagonal since these are all symmetric matrices / saves memory; 
    #NodeIdx are the indices of the elements returned
    uF = np.triu(np.ones((NumNode, NumNode), dtype=int), 1)
    i_idx, j_idx = np.where(uF)
    NodeIdx = np.vstack([i_idx, j_idx]).T
    uF_vec = uF.T.flatten()

    #pre-computation step
    AgNet_full = np.corrcoef(LocExp)
    AgNet = AgNet_full.flatten()[list(map(bool,uF_vec))]

    # allocate matrices
    #matrices to store the computed networks
    LION  = np.zeros((len(AgNet), NumSamp))

    
    for idx in range(NumSamp):

        #vector to extract all but sample idx
        idxvec = np.concatenate([np.arange(idx), np.arange(idx+1, NumSamp)])

        #compute Pearson correlation with all but sample idx;
        LocNet_full = np.corrcoef(LocExp[:, idxvec])
        LocNet = LocNet_full.flatten()[list(map(bool,uF_vec))]
#
        #LIONESS
        LION[:, idx] = NumSamp * (AgNet - LocNet) + LocNet

    return LION, AgNet, NodeIdx



####################################
###############SSN##################
####################################
def SSN(LocExp):

    # determine size of input data
    NumNode, NumSamp = LocExp.shape

    # create filter to extract upper diagonal since these are all symmetric matrices / saves memory; 
    #NodeIdx are the indices of the elements returned
    uF = np.triu(np.ones((NumNode, NumNode), dtype=int), 1)
    i_idx, j_idx = np.where(uF)
    NodeIdx = np.vstack([i_idx, j_idx]).T
    uF_vec = uF.T.flatten()

    #pre-computation step
    AgNet_full = np.corrcoef(LocExp)
    AgNet = AgNet_full.flatten()[list(map(bool,uF_vec))]


    # allocate matrices
    #matrices to store the computed networks
    SSN   = np.zeros((len(AgNet), NumSamp))

    
    for idx in range(NumSamp):

        #vector to extract all but sample idx
        idxvec = np.concatenate([np.arange(idx), np.arange(idx+1, NumSamp)])

        #compute Pearson correlation with all but sample idx;
        LocNet_full = np.corrcoef(LocExp[:, idxvec])
        LocNet = LocNet_full.flatten()[list(map(bool,uF_vec))]
#
        #SSN
        SSN[:, idx] = (NumSamp - 2) * (AgNet - LocNet) / (1 - LocNet * LocNet)

    return SSN, AgNet, NodeIdx


####################################
###############SWEET################
####################################

def SWEET(LocExp):

    # determine size of input data
    NumNode, NumSamp = LocExp.shape

    # create filter to extract upper diagonal since these are all symmetric matrices / saves memory; 
    #NodeIdx are the indices of the elements returned
    uF = np.triu(np.ones((NumNode, NumNode), dtype=int), 1)
    i_idx, j_idx = np.where(uF)
    NodeIdx = np.vstack([i_idx, j_idx]).T
    uF_vec = uF.T.flatten()

    #pre-computation step
    AgNet_full = np.corrcoef(LocExp)
    AgNet = AgNet_full.flatten()[list(map(bool,uF_vec))]

    #additional pre-computation step for SWEET
    PCCmat = np.corrcoef(LocExp, rowvar=False)  # MATLAB corr(LocExp) treats columns as variables
    muPCC = (np.sum(PCCmat, axis=0) - 1) / (NumSamp - 1)
    minPCC = np.min(muPCC)
    maxPCC = np.max(muPCC)

    # allocate matrices
    #matrices to store the computed networks
    SWEET = np.zeros((len(AgNet), NumSamp))

    
    for idx in range(NumSamp):

        #vector to extract all but sample idx
        idxvec = np.concatenate([np.arange(idx), np.arange(idx+1, NumSamp)])

        #compute Pearson correlation with all but sample idx;
        LocNet_full = np.corrcoef(LocExp[:, idxvec])
        LocNet = LocNet_full.flatten()[list(map(bool,uF_vec))]
#
        #SWEET
        Sp = (muPCC[idx] - minPCC + 0.01) / (maxPCC - minPCC + 0.01)
        SWEET[:, idx] = Sp * 0.1 * (NumSamp - 1) * (AgNet - LocNet) + LocNet

    return SWEET, AgNet, NodeIdx


####################################
###########BONOBO###################
####################################

def BONOBO(LocExp):

    # determine size of input data
    NumNode, NumSamp = LocExp.shape

    # create filter to extract upper diagonal since these are all symmetric matrices / saves memory; 
    #NodeIdx are the indices of the elements returned
    uF = np.triu(np.ones((NumNode, NumNode), dtype=int), 1)
    i_idx, j_idx = np.where(uF)
    NodeIdx = np.vstack([i_idx, j_idx]).T
    uF_vec = uF.T.flatten()

    #pre-computation step for LIONESS / SSN / SWEET
    AgNet_full = np.corrcoef(LocExp)
    AgNet = AgNet_full.flatten()[list(map(bool,uF_vec))]

    #pre-computation step for BONOBO (simplified for large N)
    AgCov = np.cov(LocExp)

    # allocate matrices
    #matrices to store the computed networks
    BONO  = np.zeros((len(AgNet), NumSamp))

    
    for idx in range(NumSamp):

        #vector to extract all but sample idx
        idxvec = np.concatenate([np.arange(idx), np.arange(idx+1, NumSamp)])
#
        #BONOBO
        LocCov = np.cov(LocExp[:, idxvec]) #BONOBO: covariance_matrix
        cov_diag = np.diag(LocCov)
        delta = 1 / (3 + 2 * np.mean(np.sqrt(cov_diag)) / np.var(cov_diag,ddof=1)) #BONOBO: delta

        sscov = delta * NumSamp * (AgCov - LocCov) + LocCov

        # sample-specific SD matrix
        sds = np.zeros((NumNode, NumNode))
        sds[np.diag_indices(NumNode)] = sscov.diagonal() ** (-0.5) #sample-specific standard deviations

        BON = sds @ sscov @ sds
        BONO[:, idx] = BON.flatten()[list(map(bool,uF_vec))]

    return BONO, AgNet, NodeIdx



####################################
##########ALL METHODS TOGETHER######
####################################

def GenerateLinearSSN(LocExp):

    # determine size of input data
    NumNode, NumSamp = LocExp.shape

    # create filter to extract upper diagonal since these are all symmetric matrices / saves memory; 
    #NodeIdx are the indices of the elements returned
    uF = np.triu(np.ones((NumNode, NumNode), dtype=int), 1)
    i_idx, j_idx = np.where(uF)
    NodeIdx = np.vstack([i_idx, j_idx]).T
    uF_vec = uF.T.flatten()

    #pre-computation step for LIONESS / SSN / SWEET
    AgNet_full = np.corrcoef(LocExp)
    AgNet = AgNet_full.flatten()[list(map(bool,uF_vec))]

    #pre-computation step for BONOBO (simplified for large N)
    AgCov = np.cov(LocExp)

    #additional pre-computation step for SWEET
    PCCmat = np.corrcoef(LocExp, rowvar=False)  # MATLAB corr(LocExp) treats columns as variables
    muPCC = (np.sum(PCCmat, axis=0) - 1) / (NumSamp - 1)
    minPCC = np.min(muPCC)
    maxPCC = np.max(muPCC)

    # allocate matrices
    #matrices to store the computed networks
    LION  = np.zeros((len(AgNet), NumSamp))
    SSN   = np.zeros((len(AgNet), NumSamp))
    SWEET = np.zeros((len(AgNet), NumSamp))
    BONO  = np.zeros((len(AgNet), NumSamp))

    
    for idx in range(NumSamp):

        #vector to extract all but sample idx
        idxvec = np.concatenate([np.arange(idx), np.arange(idx+1, NumSamp)])

        #compute Pearson correlation with all but sample idx; used for LIONESS / SSN / SWEET
        LocNet_full = np.corrcoef(LocExp[:, idxvec])
        LocNet = LocNet_full.flatten()[list(map(bool,uF_vec))]
#
        #LIONESS
        LION[:, idx] = NumSamp * (AgNet - LocNet) + LocNet
#
        #SSN
        SSN[:, idx] = (NumSamp - 2) * (AgNet - LocNet) / (1 - LocNet * LocNet)
#
        #SWEET
        Sp = (muPCC[idx] - minPCC + 0.01) / (maxPCC - minPCC + 0.01)
        SWEET[:, idx] = Sp * 0.1 * (NumSamp - 1) * (AgNet - LocNet) + LocNet
#
        #BONOBO
        LocCov = np.cov(LocExp[:, idxvec]) #BONOBO: covariance_matrix
        cov_diag = np.diag(LocCov)
        delta = 1 / (3 + 2 * np.mean(np.sqrt(cov_diag)) / np.var(cov_diag,ddof=1)) #BONOBO: delta

        sscov = delta * NumSamp * (AgCov - LocCov) + LocCov

        # sample-specific SD matrix
        sds = np.zeros((NumNode, NumNode))
        sds[np.diag_indices(NumNode)] = sscov.diagonal() ** (-0.5) #sample-specific standard deviations

        BON = sds @ sscov @ sds
        BONO[:, idx] = BON.flatten()[list(map(bool,uF_vec))]

    AllSSN = [LION, SSN, SWEET, BONO]
    return AllSSN, AgNet, NodeIdx




####################################
##########RUN METHOD WRAPPER########
####################################

#method_ID=1:LIONESS, #method_ID=2:SSN, #method_ID=3:SWEET, #method_ID=4:BONOBO, #method_ID=5:ALL,
def run_SSN_method(LocExp,method_ID):

	if method_ID==1:
		LIONESS_mat, AgNet, NodeIdx = LIONESS(LocExp)
		return LIONESS_mat, AgNet, NodeIdx

	elif method_ID==2:
		SSN_mat, AgNet, NodeIdx = SSN(LocExp)
		return SSN_mat, AgNet, NodeIdx

	elif method_ID==3:
		SWEET_mat, AgNet, NodeIdx = SWEET(LocExp)
		return SWEET_mat, AgNet, NodeIdx

	elif method_ID==4:
		BONO_mat, AgNet, NodeIdx = BONOBO(LocExp)
		return BONO_mat, AgNet, NodeIdx

	else:
		AllSSN, AgNet, NodeIdx = GenerateLinearSSN(LocExp)
		return AllSSN, AgNet, NodeIdx



if __name__ == "__main__":
	matr_0 = np.loadtxt('your_matrix.csv',delimiter=',')
	SSN_meth, AgNet, NodeIdx = run_SSN_method(matr_0,method_ID=1) #Change preferred method_ID



