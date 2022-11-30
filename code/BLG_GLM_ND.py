import numpy as np

# for testing (loads the appropriate data)
# data = io.loadmat("/Users/nadou/Documents/Documents - Nadège’s MacBook Air/GitHub/GitHub/bhack_td/results/bhack_task_03_output_distances.mat")
# Ds = data['Ds'].squeeze()
# X = Ds[0]
# Y = Ds[-1]
# Y = np.transpose(np.expand_dims(Y, axis=1), (0, 1, 3, 2))  # iEEG data, of size [ntimepairs 1 ntimechunks nsensors]


def BLG_GLM_ND(X,Y,dozscore=0,parloopflag=0):
    """
    X = predictors matrix of size[ndata npredictors dim3x, dim4x, ... dimnx]
    Y = dependent variable matrix of size[ndata npredictors dim3y, dim4y, ... dimny]
    dozscore = zscore
    parloopflag = compute inverse for different pages within parfor loop

    outputs:
    Beta = GLM coefficients
    Pred = prediction
    RSQ = R squared
    Xinv = inverse

    Handy features:
    - input X doesn't need to include an intercept term.
    - singleton dimensions>2 allowed and automatically expanded to match
      non-singleton dimensions in the other variable
     e.g.:
    X=rand([100,10,1,4]);
    Y=rand([100,1,10]);
    [Beta,Pred,RSQ,Xinv]=BLG_GLM_ND(X,Y)
    produces:
    Betas of size 11 1 10 4

    Bruno L. Giordano
    brungio@gmail.com
    November 2018
    June 2021: - implement automatic expansion of singleton dimensions>2 of X/Y
                 for corresponding non-singleton dimensions>2 in Y/X
               - input matrices can have more than 3 dimensions
    """

    sX = X.shape
    sY = Y.shape
    if sX[0] != sY[0]:
        raise ValueError("n datapoints for X and Y differs")
    if sY[1] != 1:
        raise ValueError('only one predicted Y variable allowed')

    nXdims = len(sX)
    nYdims = len(sY)
    if nXdims < nYdims:
        sX += tuple(np.ones([1, nYdims - nXdims], int)[0])
    elif nYdims < nXdims:
        sY += tuple(np.ones([1, nXdims - nYdims], int)[0])

    if len(sX) == 2:
        sX += (1,)
    if len(sY) == 2:
        sY += (1,)

    # if any((sX(3:end)-sY(3:end))~=0)
    #     warning('expanding variables along singleton pages')
    # end

    nXdims = len(sX)
    nYdims = len(sY)
    Xrep = np.ones(nXdims, int)
    Yrep = np.ones(nYdims, int)

    for i in range(2, nXdims):
        if sX[i] == 1 and sY[i] != 1:
            Xrep[i] = sY[i]
        elif sX[i] != 1 and sY[i] == 1:
            Yrep[i] = sX[i]
        elif sX[i] !=1 and sY[i] != 1 and sX[i] != sY[i]:
            raise ValueError('check consistency of X and Y pages')

    if any(Xrep != 1):
        redX = X

    if any(Xrep != 1) and len(X.shape) < len(Xrep):
        # first 'if' to avoid error in first line, second 'if' to avoid unnecessary computation
        missing_dim = int(np.where(Xrep != 1)[0])
        X = np.expand_dims(X, axis=missing_dim)
        X = np.repeat(X, Xrep[missing_dim], axis=missing_dim)

    if any(Yrep != 1) and len(Y.shape) < len(Yrep):
        missing_dim = int(np.where(Yrep != 1)[0])
        Y = np.expand_dims(Y, axis=missing_dim)
        Y = np.repeat(Y, Yrep[missing_dim], axis=missing_dim)

    sX = X.shape
    sY = Y.shape

    def fint(x):  # add intercept
        return np.concatenate((np.ones((x.shape[0], 1, x.shape[2], x.shape[3])), x), axis=1)

    X = fint(X)
    Beta = np.zeros((X.shape[1], *Y.shape[1:]))

    # THIS IS SUPER LONG (doing one loop instead of two nested loops didn't help)
    # look at sklearn ?
    for timechunk in range(X.shape[2]):
        for sensor in range(X.shape[3]):
            beta, _, _, _ = np.linalg.lstsq(X[:, :, timechunk, sensor], Y[:, :, timechunk, sensor])
            Beta[:, 0, timechunk, sensor] = beta.squeeze()


    # THE REST OF THE FUNCTION IS COMMENTED OUT FOR NOW, BECAUSE THE BETA COEFFICIENTS ARE COMPUTED ABOVE
    # AND ONLY THE Beta COEFFICIENTS (first original output) WERE USED IN THE BRAINHACK CODE
    #
    # X = np.reshape(X, [*sX[:2], np.prod(sX[2:])])
    # Y = np.reshape(Y, [*sY[:2], np.prod(sY[2:])])
    #
    # # sanity checks
    # # X = np.random.uniform(size=(100, 10))
    # # Y = np.random.uniform(size=(100,1))
    # # Beta, Pred, RSQ = BLG_GLM_ND(X, Y)
    # # Xint = np.concatenate((np.ones((X.shape[0], 1, X.shape[2])), fdemean(X)), axis=1)
    #
    # # [B, BINT, R, RINT, STATS] = regress(Y, Xint)
    # # corr(B, Beta)
    # # corr(Xint * B, Pred)
    # # disp(num2str([RSQ STATS(1)]))
    #
    # def fdemean(x):
    #     return x - x.mean(axis=0)
    #
    # def fzscore(x):
    #     return stats.zscore(x)
    #
    # # fpinv = @(x)  BLG_pinv3D(x, 1);
    #
    # def feye(x, npages):
    #     y = np.expand_dims(np.eye(x.shape[0]), axis=2)
    #     return np.repeat(y, npages, axis=2)
    #
    # # fint = @(s) ones(s); % add intercept
    #
    # def fresmatr(x, xpinv):
    #     #TODO: feye(x, size(x, 3)) - mtimesx(x, xpinv)  # residual forming matrix
    #     return None
    #
    # def fbeta(xinv, y):
    #     #TODO: mtimesx(xinv, 'n', y, 'n')  # betas
    #     return None
    #
    # def fregrdem(x):  # demean and add intercept
    #     return np.concatenate((np.ones((x.shape[0], 1, x.shape[2])), fdemean(x)), axis=1)
    #
    # def fregrzsc(x):  # zscore and add intercept
    #     return np.concatenate((np.ones((x.shape[0], 1, x.shape[2])), fzscore(x)), axis=1)
    #
    # def fpred(x, beta):
    #     #TODO: mtimesx(x, beta)  # predictions
    #     return None
    #
    # def fres(xresmatr, y):
    #     #TODO: mtimesx(xresmatr, y)  # residuals
    #     return None
    #
    # # def frsq(y, ypred):
    # #     # BLG_vectorformND(BLGmx_corr(cat(2, y, ypred))). ^ 2
    # #     return None
    #
    # def frsq(y, ypred):
    #     #TODO: BLGmx_corr2(y, ypred). ^ 2
    #     return None
    #
    # if dozscore:
    #     X = fregrzsc(X)
    #     Y = fzscore(Y)
    # else:  # else just demean predictors
    #     X = fregrdem(X)
    #
    # if any(Xrep != 1):
    #     if dozscore:
    #         redX = fregrzsc(redX)
    #     else:  # else just demean predictors
    #         redX = fregrdem(redX)
    #     sredx = redX.shape
    #     if len(sredx) == 2:
    #         sredx = sredx + tuple([1])
    #     redX = np.reshape(redX, [*sredx[:2], np.prod(sredx[2:])])
    #
    # #TODO:
    # #     # tmpXinv = zeros([sredx([2 1]) prod(sredx(3:end))])
    # #     tmpXinv = BLG_pinv3D(redX, parloopflag)
    # #     tmpXinv=reshape(tmpXinv,[sredx([2 1]) sredx(3:end)])
    # #     tmpXinv=repmat(tmpXinv,Xrep)
    # #     stmp=size(tmpXinv)
    # #     Xinv=reshape(tmpXinv,[stmp(1:2) prod(stmp(3:end))])
    # # else:
    # #     Xinv=BLG_pinv3D(X,parloopflag)
    #
    # Beta = fbeta(Xinv, Y)
    #
    # #TODO: if nargout > 1:
    #     Pred = fpred(X, Beta)
    # #TODO: if nargout > 2:
    #     RSQ = frsq(Y, Pred)
    #
    # if len(sX) > 2:
    #     Beta = np.reshape(Beta, [Beta.shape[0], Beta.shape[1], *sX[2:]])
    # #TODO:
    # #     if nargout > 1:
    # #         Pred = reshape(Pred, [size(Pred, 1) size(Pred, 2) sX(3:end)])
    # #     if nargout > 2:
    # #         RSQ = reshape(RSQ, [size(RSQ, 1) size(RSQ, 2) sX(3:end)])
    #
    # return Beta, Pred, RSQ, Xinv

    return Beta
