from sympy import Rational, sqrt

def standard_name(family,dimension):
    '''
    This function converts the family and dimension information
    to the standard form. The standard form uses alphabetic
    family names 'a', 'b', 'c' 'd', 'e', 'f' and 'g'.

    Example : su(4) -> a(3), so(10) -> d(5) etc.
    '''
    if family == 'su':
        return ('a', dimension-1)
    elif family == 'so':
        if dimension % 2 == 0:
            return ('d', dimension/2)
        else:
            return ('b', (dimension -1)/2)
    elif family == 'sp':
        if dimension % 2 == 0:
            return ('c', dimension/2)
        else:
            raise ValueError('The family sp can only have even dimensions')
    else:
        return (family, dimension)

def simple_roots(family,dimension):
    '''
    Returns a list of vectors that have the same lengths and
    scalar products as the simple roots of the given lie group.

    These vectors might not be the simple roots themselves. See
    e.g e(7), for which we use 8 dimensional vectors.
    '''

    simple_root_list = []

    if family == 'a':
        for i in range(dimension):
			vector1=[0 for j in range(i-1)]
                    +[-sqrt(i)/sqrt(2*(i+1)) for j in range(int(i>0))]
                    +[1/sqrt(2*j*(j+1)) for j in range(i+1,dimension+1)]
			vector2=[0 for j in range(i)]
                    +[-sqrt(i+1)/sqrt(2*(i+2))]
                    +[1/sqrt(2*j*(j+1)) for j in range(i+2,dimension+1)]
			simple_root_list.append([a -b for a,b in zip(vector1,vector2)])

    elif family=='b':
		for i in range(dimension-1):
			simple_root_list.append([0 for j in range(i)]
                                    +[1,-1]
                                    +[0 for j in range(dimension-i-2)])
		simple_root_list.append([0 for j in range(dimension-1)]+[1,])

    elif family=='c':
		for i in range(dimension-1):
			vector1=[0 for j in range(i-1)]
                   +[-sqrt(i)/sqrt(2*(i+1)) for j in range(int(i>0))]
                   +[1/sqrt(2*j*(j+1)) for j in range(i+1,dimension)]
                   +[0,]
			vector2=[0 for j in range(i)]
                   +[-sqrt(i+1)/sqrt(2*(i+2))]
                   +[1/sqrt(2*j*(j+1)) for j in range(i+2,dimension)]
                   +[0,]
			simple_root_list.append([a-b for a,b in zip(vector1,vector2)])

		vector3=[0 for j in range(dimension-2)]
                +[-sqrt(2*(dimension-1))/sqrt(dimension),0]
		vector4=[0 for j in range(dimension-1)]
                +[-sqrt(2)/sqrt(dimension)]

		simple_root_list.append([a-b for a,b in zip(vector3,vector4)])

    elif family=='d':
		for i in range(dimension-1):
			simple_root_list.append([0 for j in range(i)]
                                    +[1,-1]
                                    +[0 for j in range(dimension-i-2)])
		simple_root_list.append([0 for j in range(dimension-2)]
                                +[1,1])

    elif family=='e':

		if dimension==6:
			simple_root_list=[[1,-1,0,0,0,0],
                              [0,1,-1,0,0,0],
                              [0,0,1,-1,0,0],
                              [0,0,0,1,-1,0],
                              [-Rational(1,2),-Rational(1,2),-Rational(1,2),
                              -Rational(1,2),Rational(1,2),sqrt(3)/2],
                              [0,0,0,1,1,0]]
		elif dimension==7:
			simple_root_list=[[0, 0, 0, 0, 0, 0, -1, 1],
                              [0, 0, 0, 0, 0, -1, 1, 0],
                              [0, 0, 0, 0, -1, 1, 0, 0],
                              [0, 0, 0, -1, 1, 0, 0, 0],
                              [0, 0, -1, 1, 0, 0, 0, 0],
                              [0, -1, 1, 0, 0, 0, 0, 0],
                              [Rational(1,2),Rational(1,2),Rational(1,2),
                              Rational(1,2),-Rational(1,2),-Rational(1,2),
                              -Rational(1,2),-Rational(1,2)]]
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
    elif family=='f':
        simple_root_list=[[Rational(1,2),-Rational(1,2),
                         -Rational(1,2),-Rational(1,2)],
                         [0,0,0,1],
                         [0,0,1,-1],
                         [0,1,-1,0]]

    elif family=='g':
		simple_root_list=[[0,1], [sqrt(3)/2,-Rational(3,2)]]

	return simple_root_list
