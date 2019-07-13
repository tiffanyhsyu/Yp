""" Module for PCA code"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import numpy as np
from matplotlib import pyplot as plt

import pdb as debugger

from pypeit.core import qa

def func_vander(x, func, deg, minv=None, maxv=None):
    if func == "polynomial":
        return np.polynomial.polynomial.polyvander(x, deg)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legvander(xv, deg)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebvander(xv, deg)
    else:
        print("Fitting function '{0:s}' is not implemented yet. Please choose from 'polynomial', 'legendre', 'chebyshev'.")

def robust_polyfit(xarray, yarray, order, weights=None, maxone=True, sigma=3.0,
                   function="polynomial", initialmask=None, forceimask=False,
                   minv=None, maxv=None, guesses=None, bspline_par=None, verbose=True):
    """
    A robust (equally weighted) polynomial fit is performed to the xarray, yarray pairs
    mask[i] = 1 are masked values

    :param xarray: independent variable values
    :param yarray: dependent variable values
    :param order: the order of the polynomial to be used in the fitting
    :param weights: weights to be used in the fitting (weights = 1/sigma)
    :param maxone: If True, only the most deviant point in a given iteration will be removed
    :param sigma: confidence interval for rejection
    :param function: which function should be used in the fitting (valid inputs: 'polynomial', 'legendre', 'chebyshev', 'bspline')
    :param initialmask: a mask can be supplied as input, these values will be masked for the first iteration. 1 = value masked
    :param forceimask: if True, the initialmask will be forced for all iterations
    :param minv: minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    :param maxv: maximum value in the array (or the right limit for a legendre/chebyshev polynomial)
    :return: mask, ct -- mask is an array of the masked values, ct is the coefficients of the robust polyfit.
    """

    # Setup the initial mask
    if initialmask is None:
        mask = np.zeros(xarray.size, dtype=np.int)
        if forceimask:
            print("Initial mask cannot be enforced -- no initital mask supplied")
            forceimask = False
    else:
        mask = initialmask.copy()
    mskcnt = np.sum(mask)
    # Iterate, and mask out new values on each iteration
    ct = guesses
    while True:
        w = np.where(mask == 0)
        xfit = xarray[w]
        yfit = yarray[w]
        if weights is not None:
            wfit = weights[w]
        else:
            wfit = None
        ct = func_fit(xfit, yfit, function, order, w=wfit,
                      guesses=ct, minv=minv, maxv=maxv, bspline_par=bspline_par)
        yrng = func_val(ct, xarray, function, minv=minv, maxv=maxv)
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        if xarray.size-np.sum(mask) <= order+2:
            if verbose:
                print("More parameters than data points - fit might be undesirable")
            break  # More data was masked than allowed by order
        if maxone:  # Only remove the most deviant point
            tst = np.abs(yarray[w]-yrng[w])
            m = np.argmax(tst)
            if tst[m] > sigma*sigmed:
                mask[w[0][m]] = 1
        else:
            if forceimask:
                w = np.where((np.abs(yarray-yrng) > sigma*sigmed) | (initialmask==1))
            else:
                w = np.where(np.abs(yarray-yrng) > sigma*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break  # No new values have been included in the mask
        mskcnt = np.sum(mask)
    # Final fit
    w = np.where(mask == 0)
    xfit = xarray[w]
    yfit = yarray[w]
    if weights is not None:
        wfit = weights[w]
    else:
        wfit = None
    ct = func_fit(xfit, yfit, function, order, w=wfit, minv=minv, maxv=maxv, bspline_par=bspline_par)
    return mask, ct

def func_val(coeff, x, func, minv=None, maxv=None):
    """ Generic routine to return an evaluated function
    Functional forms include:
      polynomial, legendre, chebyshev, bspline, gauss

    Parameters
    ----------
    coeff : ndarray, coefficients
    x : 
    func : desired functional form
    minv
    maxv

    Returns
    -------
    values : ndarray

    """
    if func == "polynomial":
        return np.polynomial.polynomial.polyval(x, coeff)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legval(xv, coeff)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebval(xv, coeff)
    elif func == "bspline":
        return interpolate.splev(x, coeff, ext=1)
    elif func == "gaussian":
        if len(coeff) == 2:
            return gauss_2deg(x, coeff[0], coeff[1])
        elif len(coeff) == 3:
            return gauss_3deg(x, coeff[0], coeff[1], coeff[2])
        else:
            print("Not ready for this type of gaussian")
    elif func == "moffat":
        if len(coeff) == 3:
            return moffat(x, coeff[0], coeff[1], coeff[2])
        else:
            print("Not ready for this type of Moffat")
    else:
        print("Fitting function '{0:s}' is not implemented yet. Please choose from 'polynomial', 'legendre', 'chebyshev', 'bspline'.")

def func_fit(x, y, func, deg, minv=None, maxv=None, w=None, guesses=None,
             bspline_par=None, return_errors=False):
    """ General routine to fit a function to a given set of x,y points

    Parameters
    ----------
    x : ndarray
    y : ndarray
    func : str
      polynomial, legendre, chebyshev, bspline, gaussian
    deg : int
      degree of the fit
    minv : float, optional
    maxv
    w
    guesses : tuple
    bspline_par : dict
      Passed to bspline_fit()

    Returns
    -------
    coeff : ndarray or tuple
      ndarray for standard function fits
      tuple for bspline

    """
    if func == "polynomial":
        return np.polynomial.polynomial.polyfit(x, y, deg, w=w)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legfit(xv, y, deg, w=w)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebfit(xv, y, deg, w=w)
    elif func == "bspline":
        if bspline_par is None:
            bspline_par = {}
        # TODO -- Deal with this kwargs-like kludge
        return bspline_fit(x, y, order=deg, w=w, **bspline_par)
    elif func in ["gaussian"]:
        # Guesses
        if guesses is None:
            ampl, cent, sigma = guess_gauss(x, y)
            # As first guess choose slope and intercept to be zero
            b = 0
            m = 0
        else:
            if deg == 2:
                ampl, sigma = guesses
            elif deg == 3:
                ampl, cent, sigma = guesses
            elif deg == 4:
                b, ampl, cent, sigma = guesses
            elif deg == 5:
                m, b, ampl, cent, sigma = guesses
        # Error
        if w is not None:
            sig_y = 1./w
        else:
            sig_y = None
        if deg == 2:  # 2 parameter fit
            popt, pcov = curve_fit(gauss_2deg, x, y, p0=[ampl, sigma], sigma=sig_y)
        elif deg == 3:  # Standard 3 parameters
            popt, pcov = curve_fit(gauss_3deg, x, y, p0=[ampl, cent, sigma],
                                   sigma=sig_y)
        elif deg == 4:  # 4 parameters
            popt, pcov = curve_fit(gauss_4deg, x, y, p0=[b, ampl, cent, sigma],sigma=sig_y)
        elif deg == 5:  # 5 parameters
            popt, pcov = curve_fit(gauss_5deg, x, y, p0=[m, b, ampl, cent, sigma],sigma=sig_y)
        else:
            print("Not prepared for deg={:d} for Gaussian fit".format(deg))
        # Return
        if return_errors:
            return popt, pcov
        else:
            return popt
    elif func in ["moffat"]:
        # Guesses
        if guesses is None:
            ampl, cent, sigma = guess_gauss(x, y)
            p0 = ampl
            p2 = 3. # Standard guess
            p1 = (2.355*sigma)/(2*np.sqrt(2**(1./p2)-1))
        else:
            p0,p1,p2 = guesses
        # Error
        if w is not None:
            sig_y = 1./w
        else:
            sig_y = None
        if deg == 3:  # Standard 3 parameters
            popt, pcov = curve_fit(moffat, x, y, p0=[p0,p1,p2], sigma=sig_y)
        else:
            print("Not prepared for deg={:d} for Moffat fit".format(deg))
        # Return
        return popt
    else:
        print("Fitting function '{0:s}' is not implemented yet. Please choose from 'polynomial', 'legendre', 'chebyshev','bspline'.")


def basis(xfit, yfit, coeff, npc, pnpc, weights=None, skipx0=True, x0in=None, mask=None,
          function='polynomial'):
    nrow = xfit.shape[0]
    ntrace = xfit.shape[1]
    if x0in is None:
        x0in = np.arange(float(ntrace))

    # Mask out some orders if they are bad
    if mask is None or mask.size == 0:
        usetrace = np.arange(ntrace)
        outmask = np.ones((nrow, ntrace))
    else:
        usetrace = np.where(np.in1d(np.arange(ntrace), mask) == False)[0]
        outmask = np.ones((nrow, ntrace))
        outmask[:,mask] = 0.0

    # Do the PCA analysis
    eigc, hidden = get_pc(coeff[1:npc+1, usetrace], npc)

    modl = func_vander(xfit[:,0], function, npc)
    eigv = np.dot(modl[:,1:], eigc)

    med_hidden = np.median(hidden, axis=1)
    med_highorder = med_hidden.copy()
    med_highorder[0] = 0

    high_order_matrix = med_highorder.T[np.newaxis,:].repeat(ntrace, axis=0)

    coeffstr = []
    for i in range(1, npc+1):
        if weights is not None:
            tmask, coeff0 = robust_polyfit(x0in[usetrace], hidden[i-1, :], pnpc[i],
                                                   weights=weights[usetrace], sigma=2.0, function=function,
                                                   minv=x0in[0], maxv=x0in[-1])
        else:
            tmask, coeff0 = robust_polyfit(x0in[usetrace], hidden[i-1, :], pnpc[i],
                                                   sigma=2.0, function=function,
                                                   minv=x0in[0], maxv=x0in[-1])
        coeffstr.append(coeff0)
        high_order_matrix[:, i-1] = func_val(coeff0, x0in, function, minv=x0in[0], maxv=x0in[-1])
    high_fit = high_order_matrix.copy()

    high_order_fit = np.dot(eigv, high_order_matrix.T)
    sub = (yfit - high_order_fit) * outmask

    numer = np.sum(sub, axis=0)
    denom = np.sum(outmask, axis=0)
    x0 = np.zeros(ntrace, dtype=np.float)
    fitmask = np.zeros(ntrace, dtype=np.float)
    x0fit = np.zeros(ntrace, dtype=np.float)
    chisqnu = 0.0
    chisqold = 0.0
    robust = True
    if not skipx0:
        fitmask = (np.abs(denom) > 5).astype(np.int)
        if robust:
            good = np.where(fitmask != 0)[0]
            bad = np.where(fitmask == 0)[0]
            x0[good] = numer[good]/denom[good]
            imask = np.zeros(ntrace, dtype=np.float)
            imask[bad] = 1.0
            ttmask, x0res = robust_polyfit(x0in, x0, pnpc[0], weights=weights, sigma=2.0,
                                                   function=function, minv=x0in[0], maxv=x0in[-1], initialmask=imask)
            x0fit = func_val(x0res, x0in, function, minv=x0in[0], maxv=x0in[-1])
            good = np.where(ttmask == 0)[0]
            xstd = 1.0  # This should represent the dispersion in the fit
            chisq = ((x0[good]-x0fit[good])/xstd)**2.0
            chisqnu = np.sum(chisq)/np.sum(1-ttmask)
            fitmask = 1.0-ttmask
            print("Reduced chi-squared = {0:E}".format(chisqnu))
        else:
            for i in range(1, 5):
                good = np.where(fitmask != 0)[0]
                x0[good] = numer[good]/denom[good]
                x0res = func_fit(x0in[good], x0[good], function, pnpc[0],
                                         weights=weights, minv=x0in[0], maxv=x0in[-1])
                x0fit = func_val(x0res, x0in, function, minv=x0in[0], maxv=x0in[-1])
                chisq = (x0[good]-x0fit[good])**2.0
                fitmask[good] *= (chisq < np.sum(chisq)/2.0).astype(np.int)
                chisqnu = np.sum(chisq)/np.sum(fitmask)
                print("  Reduced chi-squared = {0:E}".format(chisqnu))
                if chisqnu == chisqold:
                    break
                else:
                    chisqold = chisqnu
        if chisqnu > 2.0:
            print("PCA has very large residuals")
        elif chisqnu > 0.5:
            print("PCA has fairly large residuals")

    else:
        x0res = 0.0
    x3fit = np.dot(eigv,high_order_matrix.T) + np.outer(x0fit,np.ones(nrow)).T
    outpar = dict({'high_fit': high_fit, 'x0': x0, 'x0in': x0in, 'x0fit': x0fit, 'x0res': x0res, 'x0mask': fitmask,
                   'hidden': hidden, 'usetrc': usetrace, 'eigv': eigv, 'npc': npc, 'coeffstr': coeffstr})
    return x3fit, outpar


def do_pca(data, cov=False):
    tolerance = 1.0E-5
    Nobj, Mattr = data.shape

    if cov:
        colmean = (np.sum(data,0)/Nobj)[:,np.newaxis]
        temp = np.ones((Nobj,1))
        X = data - np.dot(temp,colmean.T)
    else:
        print("PCA without cov=True is not implemented. Unable to continue")
    A = np.dot(X.T,X)
    eigva, eigve = np.linalg.eig(A)
    indx = np.where(np.abs(eigva) <= tolerance*np.max(eigva))[0]
    if np.size(indx) != 0: eigva[indx] = 0.0
    indx = np.where(np.abs(eigve) <= tolerance*np.max(eigve))[0]
    if np.size(indx) != 0: eigve[indx] = 0.0

    # Sort by increasing eigenvalue
    indx = np.argsort(eigva)
    eigva = eigva[indx]
    eigve = eigve[:,indx]
    eigva = eigva[::-1]
    eigve = eigve.T[::-1,:]

    return eigva, eigve


def extrapolate(outpar, ords, function='polynomial'):
    nords = ords.size

    x0ex = func_val(outpar['x0res'], ords, function,
                            minv=outpar['x0in'][0], maxv=outpar['x0in'][-1])

    # Order centre
    high_matr = np.zeros((nords, outpar['npc']))
    for i in range(1, outpar['npc']+1):
        if outpar['coeffstr'][i-1][0] == -9.99E9:
            high_matr[:,i-1] = np.ones(nords)*outpar['high_fit'][0,i-1]
            continue
        high_matr[:,i-1] = func_val(outpar['coeffstr'][i-1], ords, function,
                                            minv=outpar['x0in'][0], maxv=outpar['x0in'][-1])
    extfit = np.dot(outpar['eigv'], high_matr.T) + np.outer(x0ex, np.ones(outpar['eigv'].shape[0])).T
    outpar['high_matr'] = high_matr
    return extfit, outpar


def refine_iter(outpar, orders, mask, irshft, relshift, fitord, function='polynomial'):
    fail = False
    x0ex = func_val(outpar['x0res'], orders, function,  minv=outpar['x0in'][0], maxv=outpar['x0in'][-1])
    # Make the refinement
    x0ex[irshft] += relshift
    # Refit the data to improve the refinement
    good = np.where(mask != 0.0)[0]
    null, x0res = robust_polyfit(orders[good], x0ex[good], fitord, sigma=2.0, function=function,
                                         minv=outpar['x0in'][0], maxv=outpar['x0in'][-1])
    x0fit = func_val(x0res, orders, function, minv=outpar['x0in'][0], maxv=outpar['x0in'][-1])
    chisq = (x0ex[good]-x0fit[good])**2.0
    chisqnu = np.sum(chisq)/np.sum(mask)
    print("  Reduced chi-squared = {0:E}".format(chisqnu))
    if chisqnu > 0.5: # The refinement went wrong, ignore the refinement
        fail = True
    outpar['x0res'] = x0res
    extfit = np.dot(outpar['eigv'], outpar['high_matr'].T) + np.outer(x0fit, np.ones(outpar['eigv'].shape[0])).T
    return extfit, outpar, fail


def get_pc(data, k, tol=0.0, maxiter=20, nofix=False, noortho=False):

    p = data.shape[0]
    if p == 0:
        print("You need to supply more components in the PCA")
    #n = np.size(data)/p
    if k > p:
        debugger.set_trace()
        print("The number of principal components must be less than or equal to the order of the fitting function")

    # Set the initial conditions
    eigv = np.zeros((p, k))
    eigv[:k, :k] = np.identity(k)

    niter = 0
    diff = tol*2.0 + 1.0

    while (niter < maxiter) and (diff > tol):
        hidden = np.dot(np.dot(np.linalg.inv(np.dot(eigv.T, eigv)), eigv.T), data)
        oldeigv = eigv.copy()
        eigv = np.dot(data, np.dot(hidden.T, np.linalg.inv(np.dot(hidden, hidden.T))))
        if tol > 0.0:
            diff = 0.0
            for i in range(k):
                diff += np.abs( 1.0 - np.sum(oldeigv[:,i]*eigv[:,i])/np.sqrt(np.sum(oldeigv[:,i]**2)*np.sum(eigv[:,i]**2)) )
        niter += 1

    # Orthonormalize?
    if not noortho:
        for b in range(k):
            # Orthogonalize
            for bp in range(b):
                dot = np.sum(eigv[:,b]*eigv[:,bp])
                eigv[:,b] -= dot*eigv[:,bp]
            # Normalize
            dot = np.sum(eigv[:,b]**2)
            dot = 1.0/np.sqrt(dot)
            eigv[:,b] *= dot
        # Project variables onto new coordinates?
        if not nofix:
            hidden = np.dot(eigv.T, data)
            eval_hidden, evec_hidden = do_pca(hidden.T, cov=True)
            eigv = np.dot(eigv, evec_hidden.T)
        hidden = np.dot(eigv.T, data)
    return eigv, hidden

################################################
# 2D Image PCA

def image_basis(img, numpc=0):
    """
    img is a masked array
    numpc is the number of principal components
    """
    # Compute the eigenvalues and eigenvectors of the covariance matrix
    # First, subtract the median (along the spatial direction)
    #imgmed = (img.data - img.mean(axis=1).data.reshape((img.data.shape[0],1))).T
    imgmed = (img.data - np.median(img.data,axis=1).reshape((img.data.shape[0],1))).T
    imgmed[np.where(img.mask.T)] = 0.0
    eigval, eigvec = np.linalg.eig(np.cov(imgmed))
    p = np.size(eigvec, axis=1)
    # Sort the eigenvalues in ascending order
    idx = np.argsort(eigval,kind='mergesort')
    idx = idx[::-1]
    # Sort eigenvectors according to the sorted eigenvalues
    eigvec = eigvec[:, idx]
    eigval = eigval[idx]
    # Select the first few principal components
    if (numpc < p) and (numpc >= 0):
        eigvec = eigvec[:,range(numpc)]
    # Project the data
    imgmed[np.where(img.mask.T)] = 0.0
    score = np.dot(eigvec.T, imgmed)
    return eigvec, eigval, score


def pca2d(img, numpc):
    # Compute eigenvalues and eigenvectors of covariance matrix
    imgsub = (img-np.mean(img.T, axis=1)).T  # subtract the mean (along a column)
    [latent, coeff] = np.linalg.eig(np.cov(imgsub))
    p = coeff.shape[1]
    # Sort the eigenvalues/vectors in ascending order
    idx = np.argsort(latent)
    idx = idx[::-1]
    coeff = coeff[:, idx]
    if numpc < p:
        coeff = coeff[:, range(numpc)]
    # projection of the data in the new space
    proj = np.dot(coeff.T, imgsub)
    # Reconstruct the image
    imgpca = np.dot(coeff, proj).T + np.mean(img, axis=0)
    return imgpca.astype(np.float)


def pca_plot(inpar, ofit, outroot, maxp=25, pcadesc="", addOne=True,
             show=False):
    """ Saves quality control plots for a PCA analysis
    Parameters
    ----------
    inpar
    ofit
    prefix : str
      prefix for the filenames
    maxp
    pcadesc
    addOne

    Returns
    -------

    """

    plt.rcdefaults()
    #plt.rcParams['font.family']= 'times new roman'

    npc = inpar['npc']+1
    pages, npp = qa.get_dimen(npc, maxp=maxp)
    #
    x0 = inpar['x0']
    ordernum = inpar['x0in']
    x0fit = inpar['x0fit']
    usetrc = inpar['usetrc']
    hidden = inpar['hidden']
    high_fit = inpar['high_fit']
    nc = np.max(ordernum[usetrc])
    # Loop through all pages and plot the results
    ndone = 0
    for i in range(len(pages)):
        plt.clf()
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        if i == 0:
            if pages[i][1] == 1: ind = (0,)
            elif pages[i][0] == 1: ind = (0,)
            else: ind = (0,0)
            axes[ind].plot(ordernum[usetrc], x0[usetrc], 'bx')
            axes[ind].plot(ordernum, x0fit, 'k-')
            amn, amx = np.min(x0fit), np.max(x0fit)
            diff = x0[usetrc]-x0fit[usetrc]
            tdiffv = np.median(diff)
            mdiffv = 1.4826*np.median(np.abs(tdiffv-diff))
            amn -= 2.0*mdiffv
            amx += 2.0*mdiffv
            mval = amn-0.15*(amx-amn)
            dmin, dmax = tdiffv-2.0*mdiffv, tdiffv+2.0*mdiffv
            diff = mval + diff*0.20*(amx-amn)/(dmax-dmin)
            wign = np.where(np.abs(diff-np.median(diff))<4.0*1.4826*np.median(np.abs(diff-np.median(diff))))[0]
            dmin, dmax = np.min(diff[wign]), np.max(diff[wign])
            axes[ind].plot(ordernum[usetrc], diff, 'rx')
            if addOne:
                axes[ind].plot([0, nc+1], [mval,mval], 'k-')
                axes[ind].axis([0, nc+1, dmin-0.5*(dmax-dmin), amx + 0.05*(amx-amn)])
            else:
                axes[ind].plot([0, nc], [mval, mval], 'k-')
                axes[ind].axis([0, nc, dmin-0.5*(dmax-dmin), amx + 0.05*(amx-amn)])
            axes[ind].set_title("Mean Value")
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
            npp[0] -= 1
        for j in range(npp[i]):
            if pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            axes[ind].plot(ordernum[usetrc], hidden[j+ndone,:], 'bx')
            axes[ind].plot(ordernum, high_fit[:,j+ndone], 'k-')
            vmin, vmax = np.min(hidden[j+ndone,:]), np.max(hidden[j+ndone,:])
            if ofit[1+j+ndone] != -1:
                cmn, cmx = np.min(high_fit[:,j+ndone]), np.max(high_fit[:,j+ndone])
                diff = hidden[j+ndone,:]-high_fit[:,j+ndone][usetrc]
                tdiffv = np.median(diff)
                mdiffv = 1.4826*np.median(np.abs(tdiffv-diff))
                cmn -= 2.0*mdiffv
                cmx += 2.0*mdiffv
                mval = cmn-0.15*(cmx-cmn)
                dmin, dmax = tdiffv-2.0*mdiffv, tdiffv+2.0*mdiffv
                #dmin, dmax = np.min(diff), np.max(diff)
                diff = mval + diff*0.20*(cmx-cmn)/(dmax-dmin)
                wign = np.where(np.abs(diff-np.median(diff))<4.0*1.4826*np.median(np.abs(diff-np.median(diff))))[0]
                dmin, dmax = np.min(diff[wign]), np.max(diff[wign])
                #vmin, vmax = np.min(hidden[j+ndone,:][wign]), np.max(hidden[j+ndone,:][wign])
                axes[ind].plot(ordernum[usetrc], diff, 'rx')
                axes[ind].plot([0, 1+nc], [mval, mval], 'k-')
#				ymin = np.min([(3.0*dmin-dmax)/2.0,vmin-0.1*(vmax-dmin),dmin-0.1*(vmax-dmin)])
#				ymax = np.max([np.max(high_fit[:,j+ndone]),vmax+0.1*(vmax-dmin),dmax+0.1*(vmax-dmin)])
                ymin = dmin-0.5*(dmax-dmin)
                ymax = cmx + 0.05*(cmx-cmn)
                if addOne: axes[ind].axis([0, nc+1, ymin, ymax])
                else: axes[ind].axis([0, nc, ymin, ymax])
            else:
                if addOne: axes[ind].axis([0, nc+1, vmin-0.1*(vmax-vmin), vmax+0.1*(vmax-vmin)])
                else: axes[ind].axis([0, nc, vmin-0.1*(vmax-vmin), vmax+0.1*(vmax-vmin)])
            axes[ind].set_title("PC {0:d}".format(j+ndone))
            axes[ind].tick_params(labelsize=8)
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        if i == 0: npp[0] = npp[0] + 1
        # Delete the unnecessary axes
        for j in range(npp[i], axes.size):
            if pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            f.delaxes(axes[ind])
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        ndone += npp[i]
        # Save the figure
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axes.size
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if pcadesc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(pcadesc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        if show:
            plt.show()
        else:
            outfile = outroot+'.pdf'.format(i)
            f.savefig(outfile, dpi=200)
        plt.close()
        f.clf()
    del f

    plt.rcdefaults()

    return