import math
from math import sqrt
import numbers

def zeroes(height, width):
        """
        Creates a matrix of zeroes.
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
 
    def determinant(self):
        """
        Calculates the determinant of a 1x1 or 2x2 matrix.
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate determinant of non-square matrix.")
        if self.h > 2:
            raise(NotImplementedError, "Calculating determinant not implemented for matrices largerer than 2x2.")
        
        """
        Returns the determinant of a 1x1 matrix
        
        Example:

        > my_matrix = Matrix([[1]])
        > my_matrix.determinant()
          1
        """
        if self.h == 1:
            return self.g[0][0]
        
        """
        Returns the determinant of a 2x2 matrix
        Example:

        > my_matrix = Matrix([[1, 2], [3, 4])
        > my_matrix.determinant()
          -2
        """
        a, b, c, d = self.g[0][0], self.g[0][1], self.g[1][0], self.g[1][1]
        res = (a * d) - (c * b)
            
        return res
            
        

    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries)
  
        Returns the determinant of a 2x2 matrix
        Example:

        > my_matrix = Matrix([[1, 2], [3, 4])
        > my_matrix.trace()
          5                      
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate the trace of a non-square matrix.")
        
        res = 0
        
        for i in range(self.h):
            res += self.g[i][i]
            
        return res
        

    def inverse(self):
        """
        Calculates the inverse of a 1x1 or 2x2 Matrix.
        """
        if not self.is_square():
            raise(ValueError, "Non-square Matrix does not have an inverse.")
        if self.h > 2:
            raise(NotImplementedError, "inversion not implemented for matrices larger than 2x2.")
        
        res = zeroes(self.h, self.w)
        
        if self.h == 1:
            res[0][0] = 1 / self.g[0][0]
            return res
               
        normalizer = 1 / self.determinant()
        a, b, c, d = self.g[0][0], self.g[0][1], self.g[1][0], self.g[1][1]
        
        
        res[0][0] = d * normalizer
        res[0][1] = -1 * b * normalizer
        res[1][0] = -1 * c * normalizer
        res[1][1] = a * normalizer
        
        return res
        

    def T(self):
        """
        Returns a transposed copy of this Matrix.
        """
        res = zeroes(self.w, self.h)
        
        for c in range(self.w):
            for r in range(self.h):
                res.g[c][r] = self.g[r][c]
        return res

    def is_square(self):
        return self.h == self.w

    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be added if the dimensions are the same") 
        
        res = zeroes(self.h, self.w)
        
        for r in range(self.h):
            for c in range(self.w):
                res.g[r][c] = self.g[r][c] + other.g[r][c]
            
        
        return res

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """   
        res = zeroes(self.h, self.w)
        
        for r in range(self.h):
            for c in range(self.w):
                res.g[r][c] = self.g[r][c] * -1
                      
        return res

    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be subtracted if the dimensions are the same") 
        
        res = zeroes(self.h, self.w)
        
        for r in range(self.h):
            for c in range(self.w):
                res.g[r][c] = self.g[r][c] - other.g[r][c]
                
        return res

    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        def dot_product(vector1, vector2):
            """
            gets the dot product of two vectors
            """
            res = 0         
            for v1, v2 in zip(vector1, vector2):
                res += v1 * v2              
            return res
            
        transposed_other = other.T()
        
        res = zeroes(self.h, other.w)
        for r, row_A in enumerate(self.g):
            for c, row_B in enumerate(transposed_other):
                res.g[r][c] = dot_product(row_A, row_B)
        
        return res
        

    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.

        Example:

        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        if isinstance(other, numbers.Number):
            res = zeroes(self.h, self.w)
        
            for r in range(self.h):
                for c in range(self.w):
                    res.g[r][c] = self.g[r][c] * other
            return res
            