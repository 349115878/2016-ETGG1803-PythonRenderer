from vector import *


class MatrixN:
    """The matrixN class will be used to manipulate and create matrices."""
    sStrPrecision = None
    def __init__(self, rows, cols, *args):
        """This method will create a new matrix instance, will error if it does not get the right amount of arguments"""
        tempvect = []
        self.mData = []
        self.mRows = rows
        self.mCols = cols
        total = rows * cols
        if len(args) != 0 and len(args[0]) == total:
            if len(args[0]) == total:
                for i in range(total):
                    tempvect.append(args[0][i])
                    if len(tempvect) >= self.mCols:
                        a = VectorN(*tempvect)
                        self.mData.append(a)
                        tempvect = []

        elif len(args) == 0:
            zerolist = []
            for c in range(self.mCols):
                zerolist.append(0)
            for r in range(self.mRows):
                b = VectorN(*zerolist)
                self.mData.append(b)
        else:
            raise ValueError("You must pass exactly " + str(total) + " values in the data array to populate this " +
                             str(self.mRows) + " X " + str(self.mCols) + " MatrixN")

    def __mul__(self, other):
        """Mul can be used to multiply a matrix by a vector, matrix, or scalar"""
        if isinstance(other, MatrixN):
            if other.mRows != self.mCols:
                raise ValueError("These matrices cannot be multipled, The Rows and Cols must be equal.")
            else:
                r = MatrixN(self.mRows, other.mCols)
                for i in range(self.mRows):
                    for c in range(other.mCols):
                        z = dot(self.getRow(i), other.getColumn(c))
                        r[i, c] = z
                return r
        elif isinstance(other, VectorN):
            if other.mDim != self.mCols:
                raise ValueError("This vector has incorrect amount of rows to be multiplied with this matrix")
            else:
                r = MatrixN(self.mRows, 1)
                for i in range(self.mRows):
                    z = dot(self.getRow(i), other)
                    r[i, 0] = z
                return r.getColumn(0)
        elif isinstance(other, int) or isinstance(other, float):
            r = MatrixN(self.mRows, self.mCols)
            for i in range(self.mRows):
                for c in range(self.mCols):
                    z = float(self.mData[i][c] * other)
                    r[i, c] = z
            return r
        else:
            raise ValueError("You must pass a MatrixN or a VectorN, you passed a" + str(other))

    def __rmul__(self, other):
        """Rmul will be used for right handed systems"""
        if isinstance(other, VectorN):
            if other.mDim != self.mRows:
                raise ValueError("This vector has incorrect amount of rows to be multiplied with this matrix")
            else:
                r = MatrixN(1, self.mCols)
                for i in range(self.mCols):
                    z = dot(self.getColumn(i), other)
                    r[0, i] = z
                return r.getRow(0)
        elif isinstance(other, int) or isinstance(other, float):
            r = MatrixN(self.mRows, self.mCols)
            for i in range(self.mRows):
                for c in range(self.mCols):
                    z = float(self.mData[i][c] * other)
                    r[i, c] = z
            return r
        else:
            raise ValueError("You must pass a VectorN or MatrixN, you passed a" + str(other))

    def __getitem__(self, index):
        """This methods will get the value of an item of a particular index"""
        r = index[0]
        c = index[1]
        if r > self.mRows or c > self.mCols:
            raise IndexError("List index out of range.")
        return self.mData[r][c]

    def __setitem__(self, index, value):
        """This method will set the item of a particular index"""
        r = index[0]
        c = index[1]
        self.mData[r][c] = float(value)


    def __str__(self):
        """This method will return a string form of the matrix"""
        total = []
        temp_a = ""
        for i in range(len(self.mData[0])):
            if i == 0:
                if MatrixN.sStrPrecision is not None:
                    temp_a += str(round(self.mData[0][i], MatrixN.sStrPrecision))
                else:
                    temp_a += str(self.mData[0][i])
            else:
                if MatrixN.sStrPrecision is not None:
                    temp_a += " " + str(round(self.mData[0][i], MatrixN.sStrPrecision))
                else:
                    temp_a += "  " + str(self.mData[0][i])
        top = "/" + str(temp_a) + "\\" + "\n"
        total.append(top)
        temp_a = ""
        if len(self.mData) > 2:
            for i in range(self.mRows - 1):
                if i != 0 and i != self.mRows:
                    for y in range(len(self.mData[0])):
                        if y == 0:
                            if MatrixN.sStrPrecision is not None:
                                temp_a += str(round(self.mData[i][y], MatrixN.sStrPrecision))
                            else:
                                temp_a += str(self.mData[i][y])
                        else:
                            if MatrixN.sStrPrecision is not None:
                                temp_a += "  " + str(round(self.mData[i][y], MatrixN.sStrPrecision))
                            else:
                                temp_a += "  " + str(self.mData[i][y])

                    mid = "|" + str(temp_a) + "|" + "\n"
                    total.append(mid)
                    temp_a = ""

        for i in range(len(self.mData[0])):
            if i == 0:
                if MatrixN.sStrPrecision is not None:
                    temp_a += str(round(self.mData[-1][i], MatrixN.sStrPrecision))
                else:
                    temp_a += str(self.mData[-1][i])
            else:
                if MatrixN.sStrPrecision is not None:
                    temp_a += "  " + str(round(self.mData[-1][i], MatrixN.sStrPrecision))
                else:
                    temp_a += "  " + str(self.mData[-1][i])
        bottom = "\\" + str(temp_a) + "/"
        total.append(bottom)
        string = ""
        for t in range(len(total)):
            string += total[t]
        return string

    def transpose(self):
        """This will take a matrix and flip its dimensions, 1x3 becomes 3x1, etc."""
        row = self.mCols
        cols = self.mRows
        temp = MatrixN(row, cols)
        for r in range(row):
            for c in range(cols):
                temp[r, c] = self.mData[c][r]
        return temp

    def inverse(self):
        """Gets the inverse of the matrix and returns it"""
        inv = identity(self.mRows)
        m = self.copy()
        if self.mRows == self.mCols:
            for i in range(self.mCols):
                a = m.findPivot(i)
                if a is None:
                    return None
                else:
                    m = m.RowSwap(a[0], i)
                    inv = inv.RowSwap(a[0], i)
                    m = m.ScaleRow(i, 1 / a[1])
                    inv = inv.ScaleRow(i, 1 / a[1])
                    for r in range(self.mRows):
                        if r != i:
                            z = -m.mData[r][i]
                            m = m.MultRowAdd(i, r, z)
                            inv = inv.MultRowAdd(i, r, z)
            return inv
        else:
            return None


    def RowSwap(self, i, j):
        """Swaps two rows, i and j."""
        m = self.copy()
        m.setRow(i, self.mData[j])
        m.setRow(j, self.mData[i])
        return m

    def ScaleRow(self, i, s):
        """Once the rows have been swapped, the pivot will now need to equal 1, do this by doing 1 / pivotval"""
        m = self.copy()
        m.mData[i] *= s
        return m

    def MultRowAdd(self, i, s, ni):
        """This method will then sweep all the columns above and below the pivot point"""
        m = self.copy()
        m.setRow(s, m.getRow(s) + m.getRow(i) * ni)
        return m

    def findPivot(self, c):
        """This method will find the largest value in a particular column"""
        maxi = 0
        for r in range(c, self.mRows):
            curval = self.mData[r][c]
            if curval > maxi or maxi == 0:
                maxi = curval
                row = r
        if max == 0:
            return None
        return[row, maxi]

    def copy(self):
        """This method makes a deep copy of the matrix"""
        rows = self.mRows
        cols = self.mCols
        z = MatrixN(rows, cols)
        for r in range(rows):
            for c in range(cols):
                z[r, c] = self.mData[r][c]
        return z

    def getRow(self, row):
        """This method will acquire the value of a single row"""
        m = self.copy()
        return m.mData[row]

    def getColumn(self, col):
        """This method will acquire the value of a single column"""
        tempdata = []
        for r in range(self.mRows):
            tempdata.append(self.mData[r][col])
        return VectorN(*tempdata)

    def setRow(self, row, vector):
        """This method will set the value of a single row"""
        if vector.mDim == self.mCols:
            self.mData[row] = vector
        else:
            raise ValueError("Vector must be size " + str(self.mCols) +
                             " you passed a vector with length" + str(vector.mDim))

    def setColumn(self, col, vector):
        """This method will set the value of a single column"""
        if vector.mDim == self.mRows:
            for r in range(self.mRows):
                self.mData[r][col] = vector[r]
        else:
            raise ValueError("Vector must be size " + str(self.mRows) +
                             " you passed a vector with length" + str(vector.mDim))


def identity(num):
        """This function makes an identity matrix that is num x num dimensions."""
        row = num
        col = num
        temp = MatrixN(row, col)
        for r in range(row):
            for c in range(col):
                if c == r:
                    temp[r, c] = 1
                else:
                    temp[r, c] = 0
        return temp

def Scale(x, y, z, righthand):
    """This function will allow scaling of objs"""
    m = identity(4)
    m[0, 0] = x
    m[1, 1] = y
    m[2, 2] = z
    if righthand:
        m.transpose()
    return m

def Translate(x, y, z, righthand):
    """This function will translate a 3d obj with a 4d matrix"""
    m = identity(4)
    v = VectorN(x, y, z, 1)
    m.setRow(3, v)
    if righthand:
        m.transpose()
    return m

def RotateX(rad, righthand):
    """Rotatation matrix around the X axis"""
    m = identity(4)
    x = VectorN(1, 0, 0, 0)
    m.setRow(0, x)
    y = VectorN(0, math.cos(rad), math.sin(rad), 0)
    m.setRow(1, y)
    z = VectorN(0, -math.sin(rad), math.cos(rad), 0)
    m.setRow(2, z)
    if righthand:
        m.transpose()
    return m

def RotateY(rad, righthand):
    """Rotation matrix around the Y axis"""
    m = identity(4)
    x = VectorN(math.cos(rad), 0, -math.sin(rad), 0)
    m.setRow(0, x)
    y = VectorN(0, 1, 0, 0)
    m.setRow(1, y)
    z = VectorN(math.sin(rad), 0,  math.cos(rad), 0)
    m.setRow(2, z)
    if righthand:
        m.transpose()
    return m

def RotateZ(rad, righthand):
    """Rotation matrix around the Z axis"""
    m = identity(4)
    x = VectorN(math.cos(rad), math.sin(rad), 0, 0)
    m.setRow(0, x)
    y = VectorN(-math.sin(rad), math.cos(rad), 0, 0)
    m.setRow(1, y)
    z = VectorN(0, 0, 1, 0)
    m.setRow(2, z)
    if righthand:
        m.transpose()
    return m
