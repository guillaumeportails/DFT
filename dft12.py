import numpy as np

def fft_radix2(x):
    N = len(x)
    if N <= 1:
        return x
    even = fft_radix2(x[0::2])
    odd = fft_radix2(x[1::2])
    T = [np.exp(-2j * np.pi * k / N) * odd[k] for k in range(N // 2)]
    return [even[k] + T[k] for k in range(N // 2)] + \
           [even[k] - T[k] for k in range(N // 2)]

def fft_radix3(x):
    N = len(x)
    if N == 1:
        return x
    elif N % 3 != 0:
        raise ValueError("Size of x must be a power of 3")
    
    x0 = x[::3]
    x1 = x[1::3]
    x2 = x[2::3]
    
    X0 = fft_radix3(x0)
    X1 = fft_radix3(x1)
    X2 = fft_radix3(x2)
    
    W_N = np.exp(-2j * np.pi / N)
    W = np.array([W_N ** k for k in range(N)])
    X = np.zeros(N, dtype=complex)
    
    for k in range(N // 3):
        X[k] = X0[k] + W[k] * X1[k] + W[2*k] * X2[k]
        X[k + N // 3] = X0[k] + W[k + N // 3] * X1[k] + W[2 * (k + N // 3)] * X2[k]
        X[k + 2 * (N // 3)] = X0[k] + W[k + 2 * (N // 3)] * X1[k] + W[2 * (k + 2 * (N // 3))] * X2[k]
    
    return X

def mixed_radix_fft(x):
    N = len(x)
    
    # Assumes N is 12, which is 2^2 * 3
    if N != 12:
        raise ValueError("This implementation only supports N = 12")
    
    # Step 1: Split into 3 sequences of length 4 (Radix-3)
    x0 = x[0::3]       # start:end:step
    x1 = x[1::3]
    x2 = x[2::3]
    
    # Step 2: Apply Radix-2 FFT on each subsequence of length 4
    X0 = fft_radix2(x0)
    X1 = fft_radix2(x1)
    X2 = fft_radix2(x2)
    
    # Step 3: Combine results using twiddle factors
    W_N = np.exp(-2j * np.pi / N)
    X = np.zeros(N, dtype=complex)
    
    for k in range(4):  # 0..3
        X[k + 4*0] = X0[k] + W_N**(k + 4*0) * X1[k] + W_N**(2*(k + 4*0)) * X2[k]
        X[k + 4*1] = X0[k] + W_N**(k + 4*1) * X1[k] + W_N**(2*(k + 4*1)) * X2[k]
        X[k + 4*2] = X0[k] + W_N**(k + 4*2) * X1[k] + W_N**(2*(k + 4*2)) * X2[k]
    
    return X

# Example usage
X = np.random.random(12) + 1j * np.random.random(12)
Y = mixed_radix_fft(X)
print("X = [ ", *X, end="];\n")
print("Y = [ ", *Y, end="];\n")
print("e = max(abs(Y-fft(X)))")
