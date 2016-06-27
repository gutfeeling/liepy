from sympy import Rational, sqrt, Matrix

from utils import list_product


class LieGroup(object):

    def __init__(self,family,dimension):
        '''
        This function performs two tasks :

        1. Checks for possible errors in the family and dimension
        information provided by the user.

        2. Converts the family and dimension information
        to the standard form. The standard form uses alphabetic
        family names 'a', 'b', 'c' 'd', 'e', 'f' and 'g'.

        Example : su(4) -> a(3), so(10) -> d(5) etc.
        '''

        # check errors in user input

        if family not in ['a','b','c','d','e','f','g','su','so','sp']:
            raise ValueError("Sorry, I can't recognize the Lie group family."
                " Allowed options are 'a','b','c','d','e','f','g','su', 'so'"
                " and 'sp'.")

        if not isinstance(dimension,int):
            raise TypeError("The dimension of the Lie group must be an integer")

        # dimension checks

        if family == 'su':
            if dimension < 2:
                raise ValueError("Dimension of the su family cannot be less than 2")
        elif family == 'so':
            if dimension % 2 == 0 and dimension <8:
                raise ValueError("Dimension of even dimensional so groups cannot"
                " be less than 8")
            elif dimension % 2 != 0 and dimension <5:
                raise ValueError("Dimension of odd dimensional so groups cannot"
                " be less than 5")
        elif family == 'sp':
            if dimension % 2 != 0:
                raise ValueError("The family sp can only have even dimensions")
            elif dimension < 6:
                raise ValueError("Dimension of the sp family cannot be less than 6")
        elif family == 'a':
            if dimension < 1:
                raise ValueError("Dimension of the a family cannot be less than 1")
        elif family == 'b':
            if dimension < 2:
                raise ValueError("Dimension of the b family cannot be less than 2")
        elif family == 'c':
            if dimension < 3:
                raise ValueError("Dimension of the c family cannot be less than 3")
        elif family == 'd':
            if dimension < 4:
                raise ValueError("Dimension of the d family cannot be less than 4")
        elif family == 'e' and dimension not in [6,7,8]:
            raise ValueError("The family e can only have dimensions 6,7 and 8")
        elif family == 'f' and dimension != 4:
            raise ValueError("The family f can only have dimension 4")
        elif family == 'g' and dimension != 2:
            raise ValueError("The family g can only have dimension 2")


        # convert to standard form

        if family == 'su':
            self.family, self.dimension = 'a', dimension-1
        elif family == 'so':
            if dimension % 2 == 0:
                self.family, self.dimension = 'd', dimension/2
            else:
                self.family, self.dimension = 'b', (dimension -1)/2
        elif family == 'sp':
            self.family, self.dimension = 'c', dimension/2
        else:
            self.family, self.dimension = family, dimension





    def simple_roots(self):
        '''
        Returns a list of simple roots
        '''

        simple_root_list = []

        if self.family == 'a':
            for i in range(self.dimension):
                vector1=[0 for j in range(i-1)]   \
                        +[-sqrt(i)/sqrt(2*(i+1)) for j in range(int(i>0))]  \
                        +[1/sqrt(2*j*(j+1)) for j in range(i+1,self.dimension+1)]
                vector2=[0 for j in range(i)]  \
                        +[-sqrt(i+1)/sqrt(2*(i+2))]  \
                        +[1/sqrt(2*j*(j+1)) for j in range(i+2,self.dimension+1)]
                simple_root_list.append([a -b for a,b in zip(vector1,vector2)])

        elif self.family=='b':
            for i in range(self.dimension-1):
                simple_root_list.append([0 for j in range(i)]+[1,-1]
                                        +[0 for j in range(self.dimension-i-2)])
            simple_root_list.append([0 for j in range(self.dimension-1)]+[1,])

        elif self.family=='c':
            for i in range(self.dimension-1):
                vector1=[0 for j in range(i-1)]  \
                        +[-sqrt(i)/sqrt(2*(i+1)) for j in range(int(i>0))]  \
                        +[1/sqrt(2*j*(j+1)) for j in range(i+1,self.dimension)]  \
                        +[0,]
                vector2=[0 for j in range(i)]  \
                        +[-sqrt(i+1)/sqrt(2*(i+2))]  \
                        +[1/sqrt(2*j*(j+1)) for j in range(i+2,self.dimension)]  \
                        +[0,]
                simple_root_list.append([a-b for a,b in zip(vector1,vector2)])

            vector3=[0 for j in range(self.dimension-2)]  \
                    +[-sqrt(2*(self.dimension-1))/sqrt(self.dimension),0]
            vector4=[0 for j in range(self.dimension-1)]  \
                    +[-sqrt(2)/sqrt(self.dimension)]

            simple_root_list.append([a-b for a,b in zip(vector3,vector4)])

        elif self.family=='d':
            for i in range(self.dimension-1):
                simple_root_list.append([0 for j in range(i)]
                                       +[1,-1]
                                       +[0 for j in range(self.dimension-i-2)])
            simple_root_list.append([0 for j in range(self.dimension-2)]
                                    +[1,1])

        elif self.family=='e':

            if self.dimension==6:
                simple_root_list=[[1,-1,0,0,0,0],
                                  [0,1,-1,0,0,0],
                                  [0,0,1,-1,0,0],
                                  [0,0,0,1,-1,0],
                                  [-Rational(1,2),-Rational(1,2),-Rational(1,2),
                                  -Rational(1,2),Rational(1,2),sqrt(3)/2],
                                  [0,0,0,1,1,0]]
            elif self.dimension==7:
                simple_root_list=[[0, 0, 0, 0, 0, -1, 1],
                                  [0, 0, 0, 0, -1, 1, 0],
                                  [0, 0, 0, -1, 1, 0, 0],
                                  [0, 0, -1, 1, 0, 0, 0],
                                  [0, -1, 1, 0, 0, 0, 0],
                                  [0, 1, 1, 0, 0, 0, 0],
                                  [1/sqrt(2),Rational(1,2),-Rational(1,2),
                                  -Rational(1,2),-Rational(1,2),-Rational(1,2),
                                 -Rational(1,2)]]
            else:
                simple_root_list=[[0, 0, 0, 0, 0, -1, 1, 0],
                                  [0, 0, 0, 0, -1, 1, 0, 0],
                                  [0, 0, 0, -1, 1, 0, 0, 0],
                                  [0, 0, -1, 1, 0, 0, 0, 0],
                                  [0, -1, 1, 0, 0, 0, 0, 0],
                                  [-1, 1, 0, 0, 0, 0, 0, 0],
                                  [Rational(1,2),-Rational(1,2),-Rational(1,2),
                                  -Rational(1,2),-Rational(1,2),-Rational(1,2),
                                  -Rational(1,2),Rational(1,2)],
                                  [1,1,0,0,0,0,0,0]]
        elif self.family=='f':
            simple_root_list=[[Rational(1,2),-Rational(1,2),
                              -Rational(1,2),-Rational(1,2)],
                              [0,0,0,1],
                              [0,0,1,-1],
                              [0,1,-1,0]]

        elif self.family=='g':
            simple_root_list=[[0,1], [sqrt(3)/2,-Rational(3,2)]]

        return simple_root_list

    def cartan_matrix(self):
        '''
        Computes the Cartan matrix from the list of simple roots
        '''
        simple_root_list = self.simple_roots()
        dimension = self.dimension
        # Multiplying by Rational(1,2) ensures that the resulting expression
        # is a SymPy expression.
        return Matrix(dimension, dimension, lambda i,j :
            Rational(2,1)*list_product(simple_root_list[i],simple_root_list[j])/
            list_product(simple_root_list[j],simple_root_list[j]))

    @staticmethod
    def simple_root_pq(i, cartan_matrix):
        '''
        Returns the ith simple root in the q-p notation by returning
        the relevant row of the Cartan matrix
        '''
        return [cartan_matrix[i,k] for k in range(cartan_matrix.shape[0])]


    def fundamental_weights(self):
        '''
        Returns the list of fundamental weights
        '''
        simple_root_list = self.simple_roots()
        fundamental_weight_list = []
        A = Matrix(self.dimension, self.dimension,
            lambda i,j : Rational(2,1)*simple_root_list[i][j]/
            list_product(simple_root_list[i],simple_root_list[i]))
        for num in range(self.dimension):
            Y = Matrix(self.dimension, 1,
                 lambda i,j : 1 if i == num else 0)
            X = A.LUsolve(Y)
            fundamental_weight_list.append(
                 [X[i,0] for i in range(self.dimension)])
        return fundamental_weight_list

    def positive_roots(self):
        '''
        Computes all positive roots of the lie algebra
        Expresses the positive roots as commutators of simple roots
        The output is a dictionary where the key is the positive root (expressed
        in the q-p notation)and the value is a list of two items. The first item is a
        number and the second is a commutator of simple roots. This root can be
        expressed as the product of the number and the commutator.
        '''
        simple_root_list = self.simple_roots()
        cartan_matrix = self.cartan_matrix()

        roots_dict = {}
        roots_last_step = {}
        roots_this_step = {}

        for i in range(self.dimension):
            simple_root_i = [cartan_matrix[i,j] for j in range(self.dimension)]
            q_value = [0 for j in range(i)] + [2,] +  \
                       [0 for j in range(self.dimension -i -1)]
            p_value = [a-b for a,b in zip(q_value, simple_root_i)]
            roots_last_step[tuple(simple_root_i)] = [p_value, q_value, 1 , i]

        while True:

            for (roots, value) in roots_last_step.items():
                roots_dict[roots] = (value[2],value[3])

            for roots in roots_last_step:
                p_value = roots_last_step[roots][0]
                q_value = roots_last_step[roots][1]
                for i in range(self.dimension):
                    if p_value[i] > 0:
                        simple_root_i = [cartan_matrix[i,j]
                                        for j in range(self.dimension)]
                        new_root = [a+b for a,b in zip(roots,simple_root_i)]
                        try:
                            roots_this_step[tuple(new_root)][1][i]=q_value[i]+1
                        except KeyError:
                            j_value = (p_value[i] + q_value[i])/Rational(2,1)
                            m_value = -j_value + q_value[i]
                            factor = sqrt(2)/sqrt(list_product(simple_root_list[i],
                                simple_root_list[i])*(j_value + m_value +1)*
                                (j_value - m_value))
                            last_factor = roots_last_step[roots][2]
                            new_factor = factor*last_factor
                            last_commutator = roots_last_step[roots][3]
                            new_commutator = [i, last_commutator]
                            new_q_value = [0 for j in range(self.dimension)]
                            new_q_value[i]+=q_value[i]+1
                            roots_this_step[tuple(new_root)] = [None,
                                new_q_value, new_factor, new_commutator]

            for roots in roots_this_step:
                roots_this_step[roots][0] = [
                    a-b for a,b in zip(roots_this_step[roots][1], roots)]

            if len(roots_this_step) == 0:
                adjoint_rep = list([key for key in roots_last_step][0])
                return (roots_dict, adjoint_rep)
            else:
                roots_last_step = roots_this_step
                roots_this_step = {}
