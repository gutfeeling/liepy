from sympy import zeros, sqrt


def list_product(list1,list2):
    '''
    Returns the dot product of two lists
    '''
    return sum([a*b for a,b in zip(list1,list2)])

def gram_schmidt_rotation(scalar_product_matrix):
    '''
    Given a matrix of scalar products of linearly independent
    vectors, returns a rotation matrix that can convert this
    set of vectors into an orthonormal set.
    '''
    dimension = scalar_product_matrix.shape[0]
    rotation_matrix = zeros(dimension)
    for i in range(dimension):
        rotation_matrix[i,i] = 1
        for j in range(i):
            scalar_product = sum([rotation_matrix[j,k]*scalar_product_matrix[i,k]
                for k in range(j+1)])
            for k in range(j+1):
                rotation_matrix[i,k] += -rotation_matrix[j,k]*scalar_product
        norm_factor = sum([rotation_matrix[i,j]*rotation_matrix[i,k]*
            scalar_product_matrix[j,k] for j in range(i+1)
            for k in range(i+1)])
        for j in range(i+1):
            rotation_matrix[i,j] = rotation_matrix[i,j]/sqrt(norm_factor)
    return rotation_matrix
