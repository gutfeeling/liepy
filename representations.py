from sympy import Rational, sqrt, Matrix, eye, zeros, together, SparseMatrix

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

    def path_exists(self, lowering_chain, weights_dict, highest_weight):

        try:
            if lowering_chain in self.non_existent_lowering_chains[len(lowering_chain)]:
                return False
        except KeyError:
            pass

        current_weight = highest_weight

        for i in reversed(lowering_chain):
            children_list = [child for state in weights_dict[tuple(current_weight)][0]
                for child in state[3]]
            if all(child == [] for child in children_list):
                return True
            else:
                weight_list = [child[0] for child in children_list if child[2] == i]
                if weight_list == []:
                    try:
                        self.non_existent_lowering_chains[
                            len(lowering_chain)].append(lowering_chain)
                    except:
                        self.non_existent_lowering_chains[
                            len(lowering_chain)] = [lowering_chain]
                    return False
                else:
                    current_weight = weight_list[0]
        return True



    def scalar_product(self, state1, state2, simple_root_length_squared_list,
                       cartan_matrix, weights_dict, highest_weight):
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

        if len(state1) ==0 and len(state2) == 0:
            return 1

        if not self.path_exists(state1, weights_dict, highest_weight):
            return 0

        if not self.path_exists(state2, weights_dict, highest_weight):
            return 0


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

        result = 0
        moving_operator = state1[0]
        for i in range(len(state2)):
            if moving_operator == state2[i]:
                simple_root_length = simple_root_length_squared_list[moving_operator]
                product = self.highest_weight[moving_operator]
                for j in range(i+1,len(state2)):
                    product -= cartan_matrix[state2[j], moving_operator]
                product *= Rational(1,2)*simple_root_length
                new_state1 = state1[1:]
                new_state2 = state2[:i] + state2[i+1:]
                result+=product*self.scalar_product(new_state1, new_state2,
                    simple_root_length_squared_list, cartan_matrix, weights_dict,
                    highest_weight)
        result = together(result)
        try:
            self.scalar_product_dict[len(state1)][
                (tuple(state1),tuple(state2))] = result
        except KeyError:
            self.scalar_product_dict[len(state1)] = {
                (tuple(state1),tuple(state2)) : result}
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

        cartan_matrix = self.lie_group.cartan_matrix()
        simple_root_length_squared_list = self.lie_group.simple_root_length_squared_list()

        self.scalar_product_dict = {}
        self.non_existent_lowering_chains = {}

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
                    parents = value[0][state_num][4]
                    for parent in parents:
                        parent_weight = parent[0]
                        parent_state_num = parent[1]
                        weights_dict[parent_weight][0][parent_state_num][3].append(
                            [weights, state_num, parent[2], parent[3]]
                        )
                for state_num in range(len(value[0])):
                    state_list.append([state_index,] + value[0][state_num][2:4]
                        + [[],])
                    state_index+=1
                end_index = state_index -1

                weights_dict[weights] = [state_list,
                    [value[1],start_index, end_index]]

            all_states = [state for weight in weights_dict for state in
                weights_dict[weight][0]]

            print("Computed {0} weights and {1} states so far...".format(
                len(weights_dict), len(all_states)))

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
                            lowering_chain.insert(0,i)
                            new_norm = 1/sqrt(self.scalar_product(lowering_chain,
                                lowering_chain, simple_root_length_squared_list,
                                cartan_matrix, weights_dict, self.highest_weight))
                            new_p_value = [0 for j in range(self.lie_group.dimension)]
                            new_p_value[i]+=p_value[i]+1
                            new_q_value = None
                            matrix_element = Rational(1,1)*norm/new_norm
                            new_state = [new_p_value, new_q_value, new_norm,
                            lowering_chain,[[weight, state_num, i, matrix_element],]]
                            try:
                                weights_this_step[tuple(new_weight)][0].append(new_state)
                            except KeyError:
                                weights_this_step[tuple(new_weight)] = [
                                    [new_state,], None]

            # resolve degeneracy

            for weight in weights_this_step:
                all_states = weights_this_step[weight][0]
                degeneracy = len(all_states)
                if degeneracy > 1:
                    scalar_product_matrix = Matrix(degeneracy, degeneracy,
                        lambda i,j : all_states[i][2]*all_states[j][2]*
                        self.scalar_product(all_states[i][3], all_states[j][3],
                        simple_root_length_squared_list, cartan_matrix,
                        weights_dict, self.highest_weight)
                        if i > j else 0)
                    scalar_product_matrix += scalar_product_matrix.T
                    scalar_product_matrix += eye(degeneracy)
                    rref = scalar_product_matrix.rref()
                    dependent = [index for index in range(degeneracy)
                        if index not in rref[1]]

                    for i in rref[1]:
                        for j in range(degeneracy):
                            if not i==j:
                                direction = all_states[j][4][0][2]
                                parent = all_states[j][4][0][:2]
                                parent_norm = weights_last_step[
                                    parent[0]][0][parent[1]][2]
                                this_state_norm = all_states[j][2]
                                matrix_element = (Rational(1,1)*scalar_product_matrix[i,j]*
                                    parent_norm/this_state_norm)
                                all_states[i][4].append([parent[0],parent[1],direction,
                                    matrix_element])
                                if not matrix_element == 0:
                                    parent_p_value = weights_last_step[
                                        parent[0]][0][parent[1]][0]
                                    all_states[i][0][direction] = parent_p_value[direction]+1

                    weights_this_step[weight][0] = [all_states[i] for i in rref[1]]

                    # orthonormalize

                    norm_matrix = scalar_product_matrix.extract(rref[1],rref[1])
                    rotation_to_ob = gram_schmidt_rotation(norm_matrix)
                    weights_this_step[weight][1] = rotation_to_ob

                else:
                    weights_this_step[weight][1] = Matrix([[1]])

                for state in weights_this_step[weight][0]:
                    state[1] = [a+b for a,b in zip(state[0],weight)]

            weights_last_step = weights_this_step
            weights_this_step = {}

        self.scalar_product_dict = {}
        self.non_existent_lowering_chains = {}

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
            simple_root_i_nonzero = {}
            for weight in weights_dict:
                states = weights_dict[weight][0]
                for state in states:
                    state_index1 = state[0]
                    children = state[3]
                    for child in children:
                        if child[2] == i:
                            state_index2 = weights_dict[child[0]][0][child[1]][0]
                            simple_root_i_nonzero[(state_index1, state_index2)] = child[3]
            simple_root_i_matrix = SparseMatrix(representation_dimension,
                                                representation_dimension,
                                                simple_root_i_nonzero)
            representation_matrix_dict[tuple(
                self.lie_group.simple_root_pq(i, cartan_matrix))] = simple_root_i_matrix

        # convert to orthonormal basis

        rotation_to_ob_nonzero = {}
        for weight in weights_dict:
            rotation_info = weights_dict[weight][1]
            rotation_matrix = rotation_info[0]
            start_index = rotation_info[1]
            end_index = rotation_info[2]
            rotation_matrix_nonzero = dict(
                [((i,j),rotation_matrix[i-start_index,j-start_index]) for
                i in range(start_index, end_index+1) for
                j in range(start_index, end_index+1)]
            )
            rotation_to_ob_nonzero.update(rotation_matrix_nonzero)

        rotation_to_ob = SparseMatrix(representation_dimension,
                                      representation_dimension,
                                      rotation_to_ob_nonzero)

        for key in representation_matrix_dict:
            rotated_matrix = rotation_to_ob.multiply(
                representation_matrix_dict[key].multiply(rotation_to_ob.T))
            rotated_matrix.simplify()
            representation_matrix_dict[key] = rotated_matrix

        # calculate representation matrices for cartan generators

        fundamental_weights = self.lie_group.fundamental_weights()

        cartan_matrices_nonzero = [{} for i in range(
            self.lie_group.dimension)]
        for weight in weights_dict:
            weight_vector = [sum([weight[i]*fundamental_weights[i][j]
                for i in range(self.lie_group.dimension)])
                for j in range(self.lie_group.dimension)]
            states = weights_dict[weight][0]
            for state in states:
                state_index = state[0]
                for i in range(self.lie_group.dimension):
                    cartan_matrices_nonzero[i][(state_index,state_index)] = weight_vector[i]

        for i in range(self.lie_group.dimension):
            representation_matrix_dict[i] = SparseMatrix(representation_dimension,
                                                         representation_dimension,
                                                         cartan_matrices_nonzero[i])

        # calculate representation matrices for other positive roots

        positive_roots = self.lie_group.positive_roots()[0]
        positive_roots_list = [(key,) +  positive_roots[key] for key in positive_roots
                               if isinstance(positive_roots[key][1],list)]
        positive_roots_list = sorted(positive_roots_list,
                                     key = lambda item : item[3])

        for i in range(len(positive_roots_list)):
                matrix1 = representation_matrix_dict[positive_roots_list[i][2][0]]
                matrix2 = representation_matrix_dict[positive_roots_list[i][2][1]]
                positive_root_matrix = matrix1.multiply(matrix2).add(
                    - matrix2.multiply(matrix1))
                factor = positive_roots_list[i][1]
                root = positive_roots_list[i][0]

                # sympy cannot multiply a Pow object and a SparseMatrix object
                # to produce a SparseMatrix object. It produces a Mul object
                # instead. That's why factor*positive_root_matrix doesn't work.

                representation_matrix_dict[root] =  positive_root_matrix.applyfunc(
                    lambda i : factor*i)

        # generate matrices for negative roots

        keys = list(representation_matrix_dict.keys())

        for key in keys:
            if isinstance(key,tuple):
                new_key = tuple([-i for i in key])
                representation_matrix_dict[new_key] = representation_matrix_dict[key].T

        return representation_matrix_dict
