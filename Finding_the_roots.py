import numpy as np

def pagerank(M, num_iterations: int = 100, d: float = 0.85):
    N = M.shape[1]
    v = np.ones(N) / N
    M_hat = (d * M + (1 - d) / N)
    for i in range(num_iterations):
        v = v @ M_hat
    return v

if __name__ == '__main__':
    # M = np.array([[0, 0, 0, 0, 1],
    #           [0.5, 0, 0, 0, 0],
    #           [0.5, 0, 0, 0, 0],
    #           [0, 1, 0.5, 0, 0],
    #           [0, 0, 0.5, 1, 0]])
    M = np.array([[0, 1, 0, 1/3],
                  [1/2, 0, 1/2, 1/3],
                  [1/2, 0, 0, 1/3],
                  [0, 0, 1/3, 0]])
    v = pagerank(M, 100, 0.85)
    print(v)