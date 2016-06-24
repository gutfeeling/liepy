class Representation(object):

    def __init__(self, lie_group, highest_weight):
        '''
        Checks for errors in highest weight
        '''
        self.lie_group = lie_group

        # check for errors in highest weight

        if len(highest_weight) != self.lie_group.dimension:
            raise ValueError("Highest weight must have dimension %s" %
                self.lie_group.dimension)

        if not all(isinstance(i, int) for i in highest_weight):
            raise TypeError("Entries in the highest weight must be integers")

        if not all(i>=0 for i in highest_weight):
            raise ValueError("Entries in the highest weight must be non negative")

        self.highest_weight = highest_weight
