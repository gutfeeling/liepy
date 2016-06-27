from sympy import Rational, sqrt, Matrix, eye, zeros, together

from groups import LieGroup
from utils import list_product, gram_schmidt_rotation


class Representation(object):

    def __init__(self, lie_group, highest_weight):
        '''
        Initialize a Representation object as follows:
        r = Representation(lg, hw)
        Here lg is a LieGroup object and hw is a list denoting
        the highest weight.

        This function checks errors in the input parameters.
        '''

        # check for errors in lie_group

        if not isinstance(lie_group, LieGroup):
            raise TypeError("The first argument of the Representation"
                " class must be a LieGroup object, which can be imported"
                " from the file groups.py")

        self.lie_group = lie_group

        # check for errors in highest weight

        if not (isinstance(highest_weight, list)):
            raise TypeError("The second argument of the Representation"
            " class must be a list denoting the highest weight. Example:"
            " [1,1]")

        if len(highest_weight) != self.lie_group.dimension:
            raise ValueError("Highest weight must have dimension %s" %
                self.lie_group.dimension)

        if not all(isinstance(i, int) for i in highest_weight):
            raise TypeError("Entries in the highest weight must be integers")

        if not all(i>=0 for i in highest_weight):
            raise ValueError("Entries in the highest weight must be non negative")

        self.highest_weight = highest_weight

    def scalar_product(self, state1, state2, simple_root_list,
                       highest_weight_vector):
        '''
        Computes the scalar product between two states. The states are
        represented by lists containing intergers. The sequence of integers
        are basically the sequence of lowering operator that must act
        on the highest weight state to obtain this state.

        For example, the list [1,2,3] gives us a state which is obtained by
        lowering the highest weight state first by the simple root 3,
        then lowering the resulting state by simple root 2, and finally
        by simple root 1.
        '''

        ## this computation is very expensive.

        ## we need to compute many scalar products to resolve degenracies in
        ## the weight diagram.

        ## we store every computed scalar product to a dictionary in order to
        ## reuse results.

        # check if result is stored

        try:
            return self.scalar_product_dict[len(state1)][
                (tuple(state1),tuple(state2))]
        except KeyError:
            pass

        try:
            return self.scalar_product_dict[len(state1)][
                (tuple(state2),tuple(state1))]
        except KeyError:
            pass

        # if result is not stored, compute it

        if len(state1) ==0 and len(state2) == 0:
            return 1
        else:
            result = 0
            moving_operator = state1[0]
            for i in range(len(state2)):
                if moving_operator == state2[i]:
                    new_vector = highest_weight_vector
                    for j in range(i+1,len(state2)):
                        new_vector = [a-b for a,b in zip(
                            new_vector, simple_root_list[state2[j]])]
                    product = list_product(new_vector, simple_root_list[state1[0]])
                    new_state1 = state1[1:]
                    new_state2 = state2[:i] + state2[i+1:]
                    result+=product*self.scalar_product(new_state1, new_state2,
                        simple_root_list, highest_weight_vector)
            result = together(result)
            try:
                self.scalar_product_dict[len(state1)][
                    (tuple(state1),tuple(state2))] = result
            except KeyError:
                self.scalar_product_dict[len(state1)] = {
                    (tuple(state1),tuple(state2)) : result}
            return result

    def commutator_to_matrix(self, commutator, cartan_matrix, matrix_dict):
        '''
        Positive roots can be expressed as commutators of simple roots.
        For example, a positive root P may be equal to [1,[1,2]], which means
        commutator of the simple root 1 with the commutator of simple
        root 1 and simple root 2.

        After finding the representation matrices of the simple roots,
        we use this function to replace the integers (denoting simple roots)
        into their corresponding matrix expressions. We finally return a
        matrix which is equal to the commutator.
        '''

        root1 = self.lie_group.simple_root_pq(commutator[0], cartan_matrix)
        matrix1 = matrix_dict[tuple(root1)]
        if not isinstance(commutator[1],list):
            root2 = self.lie_group.simple_root_pq(commutator[1], cartan_matrix)
            matrix2 = matrix_dict[tuple(root2)]
        else:
            matrix2 = self.commutator_to_matrix(commutator[1],
                cartan_matrix, matrix_dict)
        result = matrix1*matrix2 - matrix2*matrix1

        # simplify converts the expressions into the simplest form.
        # have tried other functions like together, ratsimp etc. but
        # they dont work as well

        result.simplify()
        return result

    def __weights(self):
        '''
        A private function that returns the weights of the representation
        and the states belonging to it.

        The output is a dictionary, where the keys are the weights and the
        value is a list of two entries.

            - The first entry in the list contain the states in that weight space.
              These states are not orthogonal to each other.
            - The second entry provides the rotation matrix needed to turn the
              states into an orthogonal set.

        The first entry has the following form:
        [id, normalization factor, sequence of lowering operators, connection with
        parents]

            - 'id' is an unique id assigned to the state
            - 'normalization factor' is the factor that must be multiplied to the
              state obtained by the action of the lowering operators to get the
              normalized state
            - 'sequence of lowering operators' is the sequence of lowering
              operators that must act on the highest weight to form this state
            - 'connection with parents' is private data which are used by
              other methods of this class

        '''

        simple_root_list = self.lie_group.simple_roots()
        cartan_matrix = self.lie_group.cartan_matrix()
        fundamental_weight_list = self.lie_group.fundamental_weights()
        highest_weight_vector = [sum([fundamental_weight_list[i][j]*
            self.highest_weight[i] for i in range(self.lie_group.dimension)])
            for j in range(self.lie_group.dimension)]

        self.scalar_product_dict = {}

        weights_dict = {}
        weights_this_step = {}
        weights_last_step = {}

        state_index = 0

        # initialization

        p_value = [0 for i in range(self.lie_group.dimension)]
        q_value = [a+b for a,b in zip(self.highest_weight,p_value)]
        norm = 1
        lowering_chain = []
        parent = []
        rotation_to_ob = Matrix([[1,]])
        weights_last_step[tuple(self.highest_weight)] = [[[p_value,
            q_value, norm, lowering_chain, parent],], rotation_to_ob]

        # lower one step to generate new weights

        while len(weights_last_step) > 0:

            for (weights, value) in weights_last_step.items():
                start_index = state_index
                state_list = []
                for state_num in range(len(value[0])):
                    state_list.append([state_index,] + value[0][state_num][2:])
                    state_index+=1
                end_index = state_index -1

                weights_dict[weights] = [state_list,
                    [value[1],start_index, end_index]]

            print("Computed %s states so far..." % len(weights_dict))

            for weight in weights_last_step:
                for state_num in range(len(weights_last_step[weight][0])):
                    state = weights_last_step[weight][0][state_num]
                    q_value = state[1]
                    p_value = state[0]
                    norm = state[2]
                    for i in range(self.lie_group.dimension):
                        if q_value[i] > 0:
                            lowering_chain = list(state[3])
                            simple_root_i = [cartan_matrix[i,j]
                                        for j in range(self.lie_group.dimension)]
                            new_weight = [a-b for a,b in zip(weight,simple_root_i)]
                            j_value = (p_value[i] + q_value[i])/Rational(2,1)
                            m_value = -j_value + q_value[i]
                            new_norm = sqrt(2)/sqrt(list_product(simple_root_list[i],
                                simple_root_list[i])*(j_value - m_value +1)*
                                (j_value + m_value))
                            lowering_chain.insert(0,i)
                            new_p_value = [0 for j in range(self.lie_group.dimension)]
                            new_p_value[i]+=p_value[i]+1
                            new_q_value = None
                            new_state = [new_p_value, new_q_value, new_norm*norm,
                            lowering_chain,[[weight, state_num, i, 1/new_norm],]]
                            try:
                                weights_this_step[tuple(new_weight)][0].append(new_state)
                            except KeyError:
                                weights_this_step[tuple(new_weight)] = [
                                    [new_state,], None]

            # resolve degeneracy

            for weight in weights_this_step:
                all_states = weights_this_step[weight][0]
                degeneracy = len(all_states)
                scalar_product_matrix = Matrix(degeneracy, degeneracy,
                    lambda i,j : all_states[i][2]*all_states[j][2]*
                    self.scalar_product(all_states[i][3], all_states[j][3],
                    simple_root_list, highest_weight_vector)
                    if i > j else 0)
                scalar_product_matrix += scalar_product_matrix.T
                scalar_product_matrix += eye(degeneracy)
                rref = scalar_product_matrix.rref()
                dependent = [index for index in range(degeneracy)
                    if index not in rref[1]]
                for index in dependent:
                    direction = all_states[index][3][0]
                    parent = all_states[index][4][0][:2]
                    parent_p_value = weights_last_step[parent[0]][0][parent[1]][0]
                    for i in rref[1]:
                        if rref[0][i, index] != 0:
                            norm_independent = all_states[i][2]
                            matrix_element = sum([
                                (Rational(1,1)/norm_independent)*
                                rref[0][j,index]*scalar_product_matrix[i,j]
                                for j in rref[1] if rref[0][j,index] !=0])
                            all_states[i][0][direction] = parent_p_value[direction]+1
                            all_states[i][4].append([parent[0], parent[1],
                                direction, matrix_element])

                weights_this_step[weight][0] = [all_states[i] for i in rref[1]]
                for state in weights_this_step[weight][0]:
                    state[1] = [a+b for a,b in zip(state[0],weight)]

                # orthonormalize

                norm_matrix = scalar_product_matrix.extract(rref[1],rref[1])
                rotation_to_ob = gram_schmidt_rotation(norm_matrix)
                weights_this_step[weight][1] = rotation_to_ob

            weights_last_step = weights_this_step
            weights_this_step = {}

        self.scalar_product_dict = {}

        self.weights_dict = weights_dict
        return weights_dict

    def weights(self):
        '''
        This is the public method corresponding to the private method
        __weights(). The computation of weights in expensive. Therefore,
        we save the result after computing it once. Further calls to this
        method will only return the stored result.
        '''
        if hasattr(self, 'weights_dict'):
            return self.weights_dict
        else:
            return self.__weights()

    def dimension(self):
        '''
        Returns the dimension of the representation
        '''
        dimension = 0
        weights_dict = self.weights()
        for weight in weights_dict:
            degeneracy = len(weights_dict[weight][0])
            dimension += degeneracy
        return dimension


    def matrices(self):
        '''
        Returns the representation matrices in an orthonormal basis.

        The return value is a dictionary. The key represents the generators
        and the value is a matrix.

        The key may be an integer or a tuple. If it is an integer i, it
        denotes the ith Cartan generator. If it is a tuple, it denotes
        the root given by that tuple (in the q-p notation).
        '''

        cartan_matrix = self.lie_group.cartan_matrix()
        weights_dict = self.weights()
        representation_dimension = self.dimension()

        representation_matrix_dict = {}

        # calculate representation matrices for simple roots

        for i in range(self.lie_group.dimension):
            simple_root_i_matrix = zeros(representation_dimension)
            for weight in weights_dict:
                states = weights_dict[weight][0]
                for state in states:
                    state_index1 = state[0]
                    parents = state[3]
                    for parent in parents:
                        if parent[2] == i:
                            matrix_element = parent[3]
                            state_index2 = weights_dict[parent[0]][0][parent[1]][0]
                            simple_root_i_matrix[state_index2, state_index1] = parent[3]
            representation_matrix_dict[tuple(
                self.lie_group.simple_root_pq(i, cartan_matrix))] = simple_root_i_matrix

        # calculate representation matrices for cartan generators

        fundamental_weights = self.lie_group.fundamental_weights()

        cartan_matrices = [zeros(representation_dimension) for i in range(
            self.lie_group.dimension)]
        for weight in weights_dict:
            weight_vector = [sum([weight[i]*fundamental_weights[i][j]
                for i in range(self.lie_group.dimension)])
                for j in range(self.lie_group.dimension)]
            states = weights_dict[weight][0]
            for state in states:
                state_index = state[0]
                for i in range(self.lie_group.dimension):
                    cartan_matrices[i][state_index,state_index] = weight_vector[i]

        for i in range(self.lie_group.dimension):
            representation_matrix_dict[i] = cartan_matrices[i]

        # convert to orthonormal basis

        rotation_to_ob = zeros(representation_dimension)
        rotation_to_ob_inverse = zeros(representation_dimension)
        for weight in weights_dict:
            rotation_info = weights_dict[weight][1]
            rotation_matrix = rotation_info[0]
            rotation_matrix_inverse = rotation_matrix**-1
            start_index = rotation_info[1]
            end_index  = rotation_info[2]
            this_rotation = Matrix(representation_dimension,
                representation_dimension, lambda i,j :
                rotation_matrix[i-start_index,j-start_index]
                if i in range(start_index,end_index+1)
                and j in range(start_index,end_index+1) else 0)
            rotation_to_ob += this_rotation
            this_rotation_inverse = Matrix(representation_dimension,
                representation_dimension, lambda i,j :
                rotation_matrix_inverse[i-start_index,j-start_index]
                if i in range(start_index,end_index+1)
                and j in range(start_index,end_index+1) else 0)
            rotation_to_ob_inverse += this_rotation_inverse

        for key in representation_matrix_dict:
            rotated_matrix = rotation_to_ob* \
                representation_matrix_dict[key]*rotation_to_ob_inverse
            rotated_matrix.simplify()
            representation_matrix_dict[key] = rotated_matrix


        # calculate representation matrices for other positive roots

        positive_roots = self.lie_group.positive_roots()[0]


        for root in positive_roots:
            if isinstance(positive_roots[root][1], list):
                positive_root_matrix = self.commutator_to_matrix(
                    positive_roots[root][1], cartan_matrix,
                    representation_matrix_dict)
                factor = positive_roots[root][0]
                representation_matrix_dict[root] =  factor*positive_root_matrix

        # generate matrices for negative roots

        keys = list(representation_matrix_dict.keys())

        for key in keys:
            if isinstance(key,tuple):
                new_key = tuple([-i for i in key])
                representation_matrix_dict[new_key] = representation_matrix_dict[key].T

        return representation_matrix_dict
