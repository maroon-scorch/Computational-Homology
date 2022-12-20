import numpy as np, sys
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

def read_input(inputFile):
    """ Read and parse the input file """
    with open(inputFile, "r") as f:
        number_of_points = int(f.readline())
        points = [() for _ in range(number_of_points)]
        
        input = [points]
        for line in f.readlines():
            tokens = line.strip().split(";")
            dim_list = []
            for elt in tokens:
                t = [int(e) for e in elt.strip().split()]
                dim_list.append(t)
            input.append(dim_list) 
    return input

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
        # print(boundary_map)
        # print(np.shape(boundary_map))
        
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

def compute_simplicial_homology(complex):
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

# For the specific delta complex, please see Hatcher 2.1.8
def generate_lens_space(n):
    if n == 1:
        input = [[(), ()],
                  [(1, 0), (1, 1), (0, 0)],
                  [(0, 0, 1), (0, 2, 0)],
                  [(0, 1, 1, 0)]]
    elif n > 1:
        input = []
        # There are only 2 points in the Delta Complex
        # 0 - the point at the center
        # 1 - the point out at the polygon boundary
        input.append([(), ()])
        
        # Number of edges should be -
        # n from the boundary connecting to the top
        # 1 from bottom to top at the center
        # 1 on the boundary
        # We will make the list of edges in that order
        
        edge_list = []
        for i in range(n):
            # Appending n edges from boundary to center
            edge_list.append((1, 0))
        # the edge from bottom to top is a loop from 0 to 0
        # and the edge from boundary to itself
        edge_list.append((0,0))
        edge_list.append((1, 1))
        input.append(edge_list)
        # Number of Faces should be -
        # n on the top of the complex (ordered 0, 1, 1 vertex counter-clock wise)
        # n on the neighbors of the Tetrahedron (ordered 0, 0, 1 vertex counter-clock wise)
        # We will make the list of faces in that order
        
        face_list = []
        # The vertex of face will be ordered as such:
        for i in range(n):
            # Append each face on the top
            # start at face with left edge - 0, right edge - 1
            next = int((i + 1) % n)
            face_list.append((n+1, next, i))
        for i in range(n):
            # Append each face in neighbors
            next = int((i + 1) % n)
            face_list.append((next, i, n))
        input.append(face_list)
        
        # Number of Tetrahedrons should be n, vertex ordered 0, 0, 1, 1
        # We will make the list of Tetrahedrons in that order
        tetrahedron_list = []
        for i in range(n):
            next = int((i + 1) % n)
            tetrahedron_list.append((next, i, n + next, n + i))
        input.append(tetrahedron_list)   
    
    return input
    
if __name__ == "__main__":
    input_file = sys.argv[1]
    if input_file == "lens":
        num = int(sys.argv[2])
        input = generate_lens_space(num)
    else: 
        input = read_input(input_file)

    dimension = len(input)
    f_vector = find_f_vector(input)
    
    print("Dimension of Complex: ", dimension)
    print("F-vector: ", f_vector)
    
    homology = compute_simplicial_homology(input)
    pretty_print(homology)
    
    # generate_lens_space(100)
    