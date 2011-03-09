import numpy as np
from numpy.fft import rfft2, irfft2

def upsample(arr, factor):
    if arr.shape[0] != arr.shape[1]:
        raise ValueError("array argument must be square shape")
    arr_k = rfft2(arr)
    nkx, nky = arr_k.shape
    unkx = nkx*factor
    unky = (nky-1) * factor + 1
    upsample_fft = np.zeros((unkx, unky), dtype=np.complex128)
    upsample_fft[:nkx/2+1, :nky] = arr_k[:nkx/2+1, :]
    upsample_fft[-(nkx/2-1):, :nky] = arr_k[-(nkx/2-1):, :]
    upsampled = irfft2(upsample_fft)
    return upsampled

def _upsample(arr, factor):
    if arr.shape[0] != arr.shape[1]:
        raise ValueError("array argument must be square shape")
    arr_k = rfft2(arr)
    nkx, nky = arr_k.shape
    unkx = nkx*factor
    unky = (nky-1) * factor + 1
    upsampled = irfft2(arr_k, s=[arr.shape[0]*factor, arr.shape[1]*factor])
    return upsampled
