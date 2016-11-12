import copy

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
            raise ValueError(
                "Entries in the highest weight must be non negative")

        self.highest_weight = highest_weight

    def path_exists(self, lowering_chain, weights):
        '''
        This function checks if a certain lowering chain e.g. [1,4,2]
        is actually present in the weight diagram. Returns True if present,
        returns False otherwise.
        '''

        if tuple(lowering_chain) in self.non_existent_paths:
            return False

        current_weight = self.highest_weight
        current_level = 0
        for root in reversed(lowering_chain):
            try:
                current_weight = weights[current_level][tuple(current_weight)][
                    "connections"][root]
                current_level += 1
            except KeyError:
                self.non_existent_paths.add(tuple(lowering_chain))
                return False
        return True



    def scalar_product(self, state1, state2, simple_root_length_squared_list,
                       cartan_matrix, weights):
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

        if len(state1) ==0 and len(state2) == 0:
            return 1

        if not self.path_exists(state1, weights):
            return 0

        if not self.path_exists(state2, weights):
            return 0

        # check if result is stored

        try:
            return self.scalar_products[(tuple(state1),tuple(state2))]
        except KeyError:
            pass

        try:
            return self.scalar_products[(tuple(state2),tuple(state1))]
        except KeyError:
            pass

        # if result is not stored, compute it

        result = 0
        moving_operator = state1[0]
        for i in range(len(state2)):
            if moving_operator == state2[i]:
                simple_root_length = simple_root_length_squared_list[
                    moving_operator]
                product = self.highest_weight[moving_operator]
                for j in range(i+1,len(state2)):
                    product -= cartan_matrix[state2[j], moving_operator]
                product *= Rational(1,2)*simple_root_length
                new_state1 = state1[1:]
                new_state2 = state2[:i] + state2[i+1:]
                result+=product*self.scalar_product(new_state1, new_state2,
                    simple_root_length_squared_list, cartan_matrix, weights)
        result = together(result)
        self.scalar_products[(tuple(state1),tuple(state2))] = result
        return result

    def weights(self):
        '''
        A function that returns the weights of the representation
        and the connections between them.

        The output is a dictionary that looks as follows:

        { 0 : {                     # level in the weight diagram
            (1, 1) : {              # weights in this level
                "connections" : {
                    0 : (2, -1),    # (1,1) is connected to (2, 1) via root 0
                    1 : (-1, 2)     # (1,1) is connected to (-1, 2) via root 1
                    }
                }
            },
          1 : {...},                # second level and so on...
        }
        '''

        cartan_matrix = self.lie_group.cartan_matrix()

        # initialization
        p_values = [0 for i in range(self.lie_group.dimension)]
        completely_lowered = [False for i in
                              range(self.lie_group.dimension)]
        weights = {0 : {
            tuple(self.highest_weight) : {
                "p_values" : p_values,
                "q_values" : None,
                "completely_lowered" : completely_lowered,
                "connections" : {},
                }
            }
        }
        level = 0

        while True:
            # condition for termination
            try:
                weights_this_level = weights[level]
            except KeyError:
                break

            # lower weights completely for this level in any direction
            # in which it has not already been lowered

            for weight in weights_this_level:
                p_values = weights_this_level[weight]["p_values"]
                completely_lowered = weights_this_level[weight][
                    "completely_lowered"]
                q_values = [a+b for a,b in zip(weight, p_values)]
                for root_num in range(self.lie_group.dimension):
                    if (q_values[root_num] > 0 and
                    not completely_lowered[root_num]):
                        # lower completely
                        q_value = q_values[root_num]
                        current_level = level
                        current_weight = weight
                        current_p_value = p_values[root_num]
                        for num in range(q_value):
                            new_level = current_level + 1
                            new_p_value  = current_p_value + 1
                            this_simple_root = [cartan_matrix[root_num,i]
                                for i in range(self.lie_group.dimension)]
                            new_weight = [a-b for a,b in
                                zip(current_weight, this_simple_root)]
                            first_time = {
                                "p_values" : [0 for i in
                                    range(self.lie_group.dimension)],
                                "q_values" : None,
                                "completely_lowered" : [0 for i in
                                    range(self.lie_group.dimension)],
                                "connections" : {},
                                }
                            try:
                                raise_keyerror_if_level_absent = (
                                    weights[new_level])
                                try:
                                    raise_keyerror_if_weight_absent = (
                                        weights[new_level][
                                            tuple(new_weight)])
                                except KeyError:
                                    weights[new_level][
                                        tuple(new_weight)] = first_time
                            except KeyError:
                                weights[new_level] = {
                                    tuple(new_weight) : first_time
                                    }
                            weights[new_level][tuple(new_weight)][
                                "p_values"][root_num] = current_p_value
                            weights[new_level][tuple(new_weight)][
                                "completely_lowered"][root_num] = True

                            weights[current_level][tuple(current_weight)][
                                    "connections"][root_num] = new_weight

                            current_level = new_level
                            current_weight = new_weight
                            current_p_value = new_p_value
                        weights_this_level[weight]["completely_lowered"][
                            root_num] = True

            level += 1

        for level in weights:
            for weight in weights[level]:
                weights[level][weight] = {
                    "connections" : weights[level][weight]["connections"]}
        return weights

    def __states(self):
        '''
        A private function that returns the states in the representation.

        The return value is a dictionary that looks as follows:

        { 2 : {    # level in the weight diagram
            (0, 0) : {    # weights in this level
                "states" : [{   # independent states in this weight space
                    "norm" : sqrt(2),    # norm of the state
                    "lowering_chain" : [0,1],    # lowering chain for the state
                    "matrix_element_information" : [...]    # private
                    },{...}],               # other states in this weight
                "rotation_to_ob" : Matrix(...)    # rotation matrix to
                }                                 # go to orthonormal basis
            }
         0 : {...},
        }
        '''
        cartan_matrix = self.lie_group.cartan_matrix()
        simple_root_length_squared_list = (
            self.lie_group.simple_root_length_squared_list())
        weights = self.weights()

        self.scalar_products = {}
        self.non_existent_paths = set()

        # initialization
        states = {0 : {
            tuple(self.highest_weight) : {
                "states" : [{"lowering_chain" : [],
                             "norm" : 1,
                             "matrix_element_information" : [],
                             }],
                "rotation_to_ob" : Matrix([[1]]),
                }
            }
        }
        deepest_level = max([level for level in weights])

        # start generating states
        for level in range(deepest_level):
            dimension = 0
            for level in states:
                for weight in states[level]:
                    degeneracy = len(states[level][weight]["states"])
                    dimension += degeneracy
            print("generated {0} states so far".format(dimension))
            for weight in weights[level]:
                connections_for_this_weight = weights[level][weight][
                    "connections"]
                for root_num in connections_for_this_weight:
                    new_weight = connections_for_this_weight[root_num]
                    # lower each state
                    for state_num in range(len(
                    states[level][weight]["states"])):
                        state = states[level][weight]["states"][state_num]
                        norm = state["norm"]
                        lowering_chain = state["lowering_chain"]
                        new_lowering_chain = [root_num] + lowering_chain
                        new_scalar_product = sqrt(self.scalar_product(
                            new_lowering_chain, new_lowering_chain,
                            simple_root_length_squared_list, cartan_matrix,
                            weights))
                        if new_scalar_product == 0:
                            # this is not a valid state
                            continue
                        new_norm = 1/new_scalar_product
                        new_level = level + 1
                        new_state = {"lowering_chain" : new_lowering_chain,
                                     "norm" : new_norm,
                                     "matrix_element_information" :
                                        [{"level" : level,
                                          "weight" : weight,
                                          "state_num" : state_num,
                                          "direction" : root_num,
                                          "matrix_element" :
                                            Rational(1,1)*norm/new_norm
                                        }],
                                    }
                        try:
                            raise_keyerror_if_level_absent = (
                                states[new_level])
                            try:
                                raise_keyerror_if_weight_absent = (
                                    states[new_level][tuple(new_weight)])
                                states[new_level][tuple(new_weight)][
                                    "states"].append(new_state)
                            except KeyError:
                                states[new_level][tuple(new_weight)] = {
                                    "states" : [new_state],
                                    "rotation_to_ob" : Matrix([[1]]),
                                }
                        except KeyError:
                            states[new_level] = {
                                tuple(new_weight) : {
                                    "states" : [new_state],
                                    "rotation_to_ob" : Matrix([[1]]),
                                }
                            }
            #resolve degeneracy
            new_level = level + 1
            for weight in weights[new_level]:
                states_for_this_weight = states[new_level][weight]["states"]
                degeneracy = len(states_for_this_weight)
                if degeneracy == 1:
                    continue
                scalar_product_matrix = Matrix(degeneracy, degeneracy,
                    lambda i,j : states_for_this_weight[i]["norm"]*
                    states_for_this_weight[j]["norm"]*self.scalar_product(
                        states_for_this_weight[i]["lowering_chain"],
                        states_for_this_weight[j]["lowering_chain"],
                        simple_root_length_squared_list, cartan_matrix,
                        weights)
                    if i > j else 0)
                scalar_product_matrix += scalar_product_matrix.T
                scalar_product_matrix += eye(degeneracy)
                rref = scalar_product_matrix.rref()
                dependents = [index for index in range(degeneracy)
                              if index not in rref[1]]
                # calculate additional matrix elements (if any)
                for independent in rref[1]:
                    for state_num in range(degeneracy):
                        if state_num != independent:
                            # state_num's parent might be linked to independent
                            state_num_norm = states_for_this_weight[state_num][
                                "norm"]
                            state_num_matrix_element_information = (
                                states_for_this_weight[state_num][
                                    "matrix_element_information"])
                            parent_level = (
                                state_num_matrix_element_information[0][
                                    "level"])
                            parent_weight = (
                                state_num_matrix_element_information[0][
                                    "weight"])
                            parent_state_num = (
                                state_num_matrix_element_information[0][
                                    "state_num"])
                            direction = (
                                state_num_matrix_element_information[0][
                                    "direction"])
                            parent_norm = states[parent_level][parent_weight][
                                "states"][parent_state_num]["norm"]
                            matrix_element = (
                                Rational(1,1)*scalar_product_matrix[
                                    state_num, independent]*parent_norm/
                                state_num_norm)
                            if matrix_element == 0:
                                continue
                            states[new_level][weight]["states"][independent][
                                "matrix_element_information"].append(
                                    {"level" : parent_level,
                                     "weight" : parent_weight,
                                     "state_num" : parent_state_num,
                                     "direction" : direction,
                                     "matrix_element" : matrix_element
                                     }
                                )
                # orthonormalize
                if degeneracy == 1:
                    states[new_level][weight]["rotation_to_ob"] = Matrix([[1]])
                else:
                    norm_matrix = scalar_product_matrix.extract(
                        rref[1], rref[1])
                    rotation_to_ob = gram_schmidt_rotation(norm_matrix)
                    states[new_level][weight]["rotation_to_ob"] = rotation_to_ob

                # keep only independent states
                states[new_level][weight]["states"] = [
                    states_for_this_weight[i] for i in rref[1]]

        # add an unique index to all states
        state_index = 0
        for level in states:
            for weight in states[level]:
                states[level][weight]["start_index"] = state_index
                for state_num in range(len(states[level][weight]["states"])):
                    states[level][weight]["states"][state_num]["index"] = (
                        state_index)
                    state_index += 1
                states[level][weight]["end_index"] = state_index - 1
        print("total number of states is {0}".format(state_index))
        del self.scalar_products
        del self.non_existent_paths
        return states

    def states(self):
        '''
        Public function corresponding to the private function __states()
        Implements caching of the output of the private function.
        '''
        if not hasattr(self, "stored_states"):
            self.stored_states = self.__states()
        return self.stored_states

    def dimension(self):
        '''
        Returns the dimension of the representation
        '''
        dimension = 0
        states = self.states()
        for level in states:
            for weight in states[level]:
                degeneracy = len(states[level][weight]["states"])
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
        states = self.states()
        representation_dimension = self.dimension()

        representation_matrices = {}

        # calculate representation matrices for simple roots

        for i in range(self.lie_group.dimension):
            simple_root_i_nonzero = {}
            for level in states:
                for weight in states[level]:
                    states_for_this_weight = states[level][weight]["states"]
                    for state in states_for_this_weight:
                        state_index1 = state["index"]
                        matrix_element_information = state[
                            "matrix_element_information"]
                        for entry in matrix_element_information:
                            if entry["direction"] == i:
                                parent_state = states[entry["level"]][
                                    entry["weight"]]["states"][
                                        entry["state_num"]]
                                state_index2 = parent_state["index"]
                                simple_root_i_nonzero[(
                                    state_index2, state_index1)] = (
                                        entry["matrix_element"])
            simple_root_i_matrix = SparseMatrix(representation_dimension,
                                                representation_dimension,
                                                simple_root_i_nonzero)
            representation_matrices[tuple(
                self.lie_group.simple_root_pq(i, cartan_matrix))] = (
                    simple_root_i_matrix)

        # convert to orthonormal basis

        rotation_to_ob_nonzero = {}
        for level in states:
            for weight in states[level]:
                rotation_matrix = states[level][weight]["rotation_to_ob"]
                start_index = states[level][weight]["start_index"]
                end_index = states[level][weight]["end_index"]
                rotation_matrix_nonzero = dict(
                    [((i,j),rotation_matrix[i-start_index, j-start_index]) for
                    i in range(start_index, end_index+1) for
                    j in range(start_index, end_index+1)]
                    )
                rotation_to_ob_nonzero.update(rotation_matrix_nonzero)

        rotation_to_ob = SparseMatrix(representation_dimension,
                                      representation_dimension,
                                      rotation_to_ob_nonzero)

        for key in representation_matrices:
            rotated_matrix = rotation_to_ob.multiply(
                representation_matrices[key].multiply(rotation_to_ob.T))
            rotated_matrix.simplify()
            representation_matrices[key] = rotated_matrix

        # calculate representation matrices for cartan generators

        fundamental_weights = self.lie_group.fundamental_weights()

        cartan_matrices_nonzero = [{} for i in range(
            self.lie_group.dimension)]
        for level in states:
            for weight in states[level]:
                weight_vector = [sum([weight[i]*fundamental_weights[i][j]
                    for i in range(self.lie_group.dimension)])
                    for j in range(self.lie_group.dimension)]
                states_for_this_weight = states[level][weight]["states"]
                for state in states_for_this_weight:
                    state_index = state["index"]
                    for i in range(self.lie_group.dimension):
                        cartan_matrices_nonzero[i][(state_index,
                            state_index)] = weight_vector[i]

        for i in range(self.lie_group.dimension):
            representation_matrices[i] = SparseMatrix(representation_dimension,
                representation_dimension, cartan_matrices_nonzero[i])

        # calculate representation matrices for other positive roots

        positive_roots = self.lie_group.positive_roots()[0]
        positive_roots_list = [(key,) +  positive_roots[key]
                               for key in positive_roots
                               if isinstance(positive_roots[key][1],list)]
        positive_roots_list = sorted(positive_roots_list,
                                     key = lambda item : item[3])

        for i in range(len(positive_roots_list)):
                matrix1 = representation_matrices[positive_roots_list[i][2][0]]
                matrix2 = representation_matrices[positive_roots_list[i][2][1]]
                positive_root_matrix = matrix1.multiply(matrix2).add(
                    - matrix2.multiply(matrix1))
                factor = positive_roots_list[i][1]
                root = positive_roots_list[i][0]

                # sympy cannot multiply a Pow object and a SparseMatrix object
                # to produce a SparseMatrix object. It produces a Mul object
                # instead. That's why factor*positive_root_matrix doesn't work.

                representation_matrices[root] =  positive_root_matrix.applyfunc(
                    lambda i : factor*i)

        # generate matrices for negative roots

        keys = list(representation_matrices.keys())

        for key in keys:
            if isinstance(key,tuple):
                new_key = tuple([-i for i in key])
                representation_matrices[new_key] = (
                    representation_matrices[key].T)

        return representation_matrices
