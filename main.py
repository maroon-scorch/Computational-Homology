import numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

def normal_basis(n, k):
    """ Returns e_k of the rank n free Z-module """
    output = [0 for _ in range(n)]
    output[k] = 1
    return np.asarray(output)

def find_f_vector(complex):
    """ Given a Delta complex, computes its f-vector """
    vector = [len(elt) for elt in complex]
    return vector

def boundary_map(prev, next):
    """ Given C_n(X) and C_{n-1}(X) - represented as
    prev and next, finds the boundary map between the chain complexes. """
    # Chain Complex is 0
    if prev == []:
        # In this case the boundary map is the zero map, which we represent as 0
        return 0
    elif next == []:
        return np.zeros((len(prev), 1))
    else:
        next_rank = len(next)
        next_dict = [normal_basis(next_rank, k) for k in range(next_rank)]
        
        boundary_map = []
        for elt in prev:
            column_vector = np.zeros(next_rank)
            for list_idx, complex_idx in enumerate(elt):
                if list_idx % 2 == 0:
                    column_vector = column_vector + next_dict[complex_idx]
                else:
                    column_vector = column_vector - next_dict[complex_idx]
            boundary_map.append(column_vector)
        
        boundary_map = np.column_stack(boundary_map)
        print(boundary_map)
        print(np.shape(boundary_map))
        
        return boundary_map

def find_elementary_divisors(matrix, rank):
    """ Given a matrix and its rank, find its elementary divisors """
    mat = matrix.tolist()
    m = Matrix(mat)
    smf = smith_normal_form(m, domain=ZZ)
    smf = np.asarray(smf).tolist()
    
    elementary_divisors = []
    for i in range(rank):
        elementary_divisors.append(smf[i][i])
    
    return elementary_divisors

def compute_homology(complex):
    complex.append([])
    rank = find_f_vector(complex)
    
    homology = []
    for i in range(0, len(complex) - 1):
        if i == 0:
            # Computes the zero-th homology group
            A = boundary_map(complex[1], complex[0])
            B = boundary_map(complex[0], [])
            
            rank_a = np.linalg.matrix_rank(A)
            rank_b = np.linalg.matrix_rank(B)
            rank_h = rank[i] - rank_a - rank_b
            elementary_divisors = find_elementary_divisors(A, rank_a)
            homology.append((rank_h, elementary_divisors))
        else:
            # Computes the zero-th homology group
            A = boundary_map(complex[i+1], complex[i])
            B = boundary_map(complex[i], complex[i-1])
            
            if type(A) == type(0):
                print("Test")
                rank_a = 0
                rank_b = np.linalg.matrix_rank(B)
                rank_h = rank[i] - rank_a - rank_b
                elementary_divisors = []
            else:
                rank_a = np.linalg.matrix_rank(A)
                rank_b = np.linalg.matrix_rank(B)
                rank_h = rank[i] - rank_a - rank_b
                elementary_divisors = find_elementary_divisors(A, rank_a)
            homology.append((rank_h, elementary_divisors))
    print(homology)  
    return homology

def pretty_print(homology):
    print("Homology of X: -------------------------------------")
    for i in range(len(homology)):
        rank, divisors = homology[i]
        text = "H_" + str(i) + "(X) = Z^" + str(rank)
        for d in divisors:
            if d != 1 and d != -1:
                text = text + " X Z/" + str(d) + "Z"
        print(text)
    print("For all i greater than " + str(len(homology)) + ", H_i(X) is 0.")
    print("----------------------------------------------------")
        

if __name__ == "__main__":
    
    # Real Projective Plane
    input = [[(), ()],
             [(1,0), (1,0), (0,0)],
             [(1,0,2), (0, 1, 2)]]
    
    # Torus
    input = [[()],
             [(0, 0), (0,0), (0,0)],
             [(0,1,2), (2, 1, 0)]]
    
    dimension = len(input)
    f_vector = find_f_vector(input)
    
    print("Dimension of Complex: ", dimension)
    print("F-vector: ", f_vector)
    
    # boundary_map(input[2], input[1])
    # boundary_map(input[1], input[0])
    # boundary_map([], input[0])
    
    homology = compute_homology(input)
    pretty_print(homology)
    