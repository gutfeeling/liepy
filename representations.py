from sympy import Matrix, Rational, sqrt, together

from groups import LieGroup
from utils import list_product


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

    def states(self):

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

        # initialization

        p_value = [0 for i in range(self.lie_group.dimension)]
        q_value = [a+b for a,b in zip(self.highest_weight,p_value)]
        factor = 1
        lowering_chain = []
        parent = None
        weights_last_step[tuple(self.highest_weight)] = [[p_value,
            q_value, factor, lowering_chain, parent],]

        # lower one step to generate new weights

        while len(weights_last_step) > 0:

            for (weights, value) in weights_last_step.items():
                weights_dict[weights] = [l[2:] for l in value]

            print("Computed %s states so far..." % len(weights_dict))

            for weight in weights_last_step:
                for state_num in range(len(weights_last_step[weight])):
                    state = weights_last_step[weight][state_num]
                    q_value = state[1]
                    p_value = state[0]
                    # ensures that we don't accidentally modify the original list
                    for i in range(self.lie_group.dimension):
                        if q_value[i] > 0:
                            lowering_chain = list(state[3])
                            simple_root_i = [cartan_matrix[i,j]
                                        for j in range(self.lie_group.dimension)]
                            new_weight = [a-b for a,b in zip(weight,simple_root_i)]
                            j_value = (p_value[i] + q_value[i])/Rational(2,1)
                            m_value = -j_value + q_value[i]
                            factor = sqrt(2)/sqrt(list_product(simple_root_list[i],
                                simple_root_list[i])*(j_value - m_value +1)*
                                (j_value + m_value))
                            lowering_chain.insert(0,i)
                            new_p_value = [0 for j in range(self.lie_group.dimension)]
                            new_p_value[i]+=p_value[i]+1
                            new_q_value = None
                            new_state = [new_p_value, new_q_value, factor,
                            lowering_chain,[[weight, state_num, i, 1/factor],]]
                            try:
                                weights_this_step[tuple(new_weight)].append(new_state)
                            except KeyError:
                                weights_this_step[tuple(new_weight)] = [new_state,]

            # resolve degeneracy

            for weight in weights_this_step:
                all_states = weights_this_step[weight]
                degeneracy = len(all_states)
                scalar_product_matrix = Matrix(degeneracy, degeneracy,
                    lambda i,j : self.scalar_product(all_states[i][3],
                    all_states[j][3], simple_root_list, highest_weight_vector))
                rref = scalar_product_matrix.rref()
                dependent = [index for index in range(degeneracy)
                    if index not in rref[1]]
                for index in dependent:
                    direction = all_states[index][3][0]
                    parent = all_states[index][4][0][:2]
                    parent_p_value = weights_last_step[parent[0]][parent[1]][0]
                    for i in rref[1]:
                        if rref[0][i, index] != 0:
                            factor_independent = all_states[i][2]
                            new_factor = sum([
                                (Rational(1,1)/factor_independent)*
                                rref[0][j,index]*scalar_product_matrix[i,j]
                                for j in rref[1] if rref[0][j,index] !=0])
                            all_states[i][0][direction] = parent_p_value[direction]+1
                            all_states[i][4].append([parent[0], parent[1],
                                direction, new_factor])
                weights_this_step[weight] = [all_states[i] for i in rref[1]]
                for state in weights_this_step[weight]:
                    state[1] = [a+b for a,b in zip(state[0],weight)]


            weights_last_step = weights_this_step
            weights_this_step = {}

        self.scalar_product_dict = {}
        return weights_dict
