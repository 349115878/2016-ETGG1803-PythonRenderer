import math


class VectorN(object):
    """This is the base class that will work for N dimensional vectors."""
    def __init__(self, *args):
        """This is the constructor, obviously. i set it to crunch through all of the numbers
        it gets passed and makes both its dimension and its list size
        dependent on how many things are passed."""
        self.mData = []
        for i in range(len(args)):
            z = float(args[i])
            self.mData.append(z)
        self.mDim = len(self.mData)

        if self.mDim == 2:
            self.__class__ = Vector2
        if self.mDim == 3:
            self.__class__ = Vector3

    def pairwise_mult(self, other):
        """This method will take two vectors and make a pairwise vector out of them."""
        if isinstance(other, VectorN):
            if len(other) != self.mDim:
                raise ValueError("You can only pairwise another " + str(self.__class__.__name__) +
                                 " to this " + str(self.__class__.__name__) +
                                 " (you passed " + "'" + str(self) + "').")
        if not isinstance(other, VectorN):
            raise ValueError("You can only pairwise another " + str(self.__class__.__name__) + " to this " + str(
                self.__class__.__name__) +
                             " (you passed " + "'" + str(self) + "').")
        i = self.copy()
        for z in range(self.mDim):
            i[z] *= other[z]
        return i

    def clamp(self):
        """Goes through the vector and clamps any part that is above 1.0, or below 0.0."""
        for i in range(self.mDim):
            if self.mData[i] > 1:
                self.mData[i] = float(1)
            if self.mData[i] < 0:
                self.mData[i] = float(0)
        return Vector3(self.mData[0], self.mData[1], self.mData[2])

    def __str__(self):
        """This method will return one main string from taking every one number out and adding it to the string."""
        s = "<Vector" + str(self.mDim) + ": "
        for i in range(len(self.mData)):
            s += str(self.mData[i])
            if i < len(self.mData) - 1:
                s += ", "
        s += ">"
        return s

    def __len__(self):
        """this method will give back the total dimensions (length) of the vector"""
        return self.mDim

    def __getitem__(self, index):
        """This method gets items from the vector and returns them to be printed."""
        return self.mData[index]

    def __setitem__(self, index, value):
        """This method will allow you to change a specific point on the vector with a new attribute,
        provided you give it the right index and don't put 20 as the index in a 2d vector or something."""
        x = float(value)
        self.mData[index] = x

    def __eq__(self, value):
        """Will check two vectors, if they are the same, returns true, if not returns false."""
        if value == self.mData:
            return True
        else:
            return False

    def copy(self):
        """Makes a deep copy of VectorN as opposed to a shallow copy."""
        x = VectorN(*self.mData)
        x.__class__ = self.__class__
        return x

    def __neg__(self):
        """Will be called by a vector and return the negative value of the vector."""
        r = self.copy()
        for i in range(self.mDim):
            r[i] = -r[i]
        return r

    def __add__(self, vector):
        """Takes two vectors of the same dimensions and adds them together."""
        if isinstance(vector, VectorN):
            if len(vector) != self.mDim:
                raise ValueError("You can only add another " + str(self.__class__.__name__) +
                                 " to this " + str(self.__class__.__name__) +
                                 " (you passed " + "'" + str(vector) + "').")
        if not isinstance(vector, VectorN):
            raise ValueError("You can only add another " + str(self.__class__.__name__) + " to this " + str(
                self.__class__.__name__) +
                             " (you passed " + "'" + str(vector) + "').")
        i = self.copy()
        for z in range(self.mDim):
            i[z] += vector[z]
        return i

    def __sub__(self, vector):
        """Takes two vectors of the same dimensions and subtracts them."""
        if isinstance(vector, VectorN):
            if len(vector) != self.mDim:
                raise ValueError("You can only subtract another " + str(self.__class__.__name__) +
                                 " from this " + str(self.__class__.__name__) +
                                 " (you passed " + "'" + str(vector) + "').")
        if not isinstance(vector, VectorN):
            raise ValueError("You can only subtract another " + str(self.__class__.__name__) + " from this "
                             + str(self.__class__.__name__) + " (you passed " + "'" + str(vector) + "').")
        i = self.copy()
        for z in range(self.mDim):
            i[z] -= vector[z]
        return i

    def __mul__(self, scalar):
        """multiplies a vector with a scalar and then returns the result."""
        if not isinstance(scalar, int) and not isinstance(scalar, float):
                return NotImplemented

        i = self.copy()
        for z in range(self.mDim):
            i[z] *= float(scalar)
        return i

    def __rmul__(self, scalar):
        """Multiples a vector with a scalar when the order is reversed."""
        if isinstance(scalar, VectorN) or isinstance(scalar, str) or not isinstance(scalar, float):
                raise ValueError("You can only multiply this " + str(self.__class__.__name__) + " and a scalar " +
                                 " You attempted to multiply by " + "'" + str(scalar) + "'.")
        i = self.copy()
        for z in range(self.mDim):
            i[z] *= float(scalar)
        return i

    def __truediv__(self, scalar):
        """Takes a vector that you pass to it and divides it by a scalar"""
        if isinstance(scalar, VectorN) or isinstance(scalar, str):
                raise ValueError("You can only divide this " + str(self.__class__.__name__) + " by a scalar " +
                                 " you attempted to divide by " + " '" + str(scalar) + "'.")
        i = self.copy()
        for z in range(self.mDim):
            i[z] /= float(scalar)
        return i

    def magnitude(self):
        """Will determine the magnitude of a vector with scalars. Will use the pythagorean theorem to determine
         the hypotenuse."""
        i = self.copy()
        total = 0
        for m in range(self.mDim):
            i[m] **= 2
            total += i[m]
        return total ** (1/2)

    def magnitudeSquared(self):
        """Takes the magnitude function and squares it to get a new power. useful for something i imagine."""
        return self.magnitude() ** 2

    def normalized(self):
        """Will normalize a vector so that its magnitude will always equal 1.0 ."""
        normal = self.copy()
        if normal.isZero():
            return normal
        else:
            for i in range(self.mDim):
                normal[i] /= self.magnitude()
            return normal

    def isZero(self):
        """This method will check if an instance is equal to zero. If it is, returns True, otherwise returns False."""
        totalzero = 0
        for i in range(self.mDim):
            if self[i] == float(0):
                totalzero += 1
        if totalzero == self.mDim:
            return True
        else:
            return False

    @property
    def x(self):
        """Getter method for the first data point in a vector list."""
        return self.mData[0]

    @property
    def i(self):
        """will convert all the floats in the vector to integers."""
        tempdata = []
        for i in range(len(self.mData)):
            x = int(self.mData[i])
            tempdata.append(x)
        return tuple(tempdata)

    @x.setter
    def x(self, new):
        """Setter method for the first data point in a vector list, will change the value of it."""
        self.mData[0] = float(new)

    @property
    def y(self):
        """Getter method for second data point in a vector list"""
        return self.mData[1]

    @y.setter
    def y(self, new):
        """Setter method for the second data point in a vector list, will change the value of it."""
        self.mData[1] = float(new)

    @property
    def z(self):
        """Getter method for z, only useful for 3d vectors, so its not going to be in 2d vectors."""
        return self.mData[2]

    @z.setter
    def z(self, new):
        """Setter method for z, will take whatever value is passed and change the existing value for z."""
        self.mData[2] = new


class Vector2(VectorN):
    """this class will be for 2d vector"""

    @property
    def y(self):
        """Getter method for second data point in a vector list"""
        return self.mData[1]

    @y.setter
    def y(self, new):
        """Setter method for the second data point in a vector list, will change the value of it."""
        self.mData[1] = float(new)

    @property
    def radians(self):
        """Getter method for radians, it will find how many radians there are with the first and second data points."""
        x = math.atan2(self.mData[1], self.mData[0])
        return x
    @property
    def radians_inv(self):
        """Getter method for radians, it will find how many radians there are with the first and second data points."""
        x = math.atan2(-self.mData[1], self.mData[0])
        return x


    @property
    def degrees(self):
        """Getter method for degrees, same as radians, but converts it to degrees."""
        x = math.atan2(self.mData[1], self.mData[0])
        x = x * 180/3.141592653589793
        return x

    @property
    def perpendicular(self):
        """Gets the perpendicular for a 2D vector by doing stuff"""
        Vperp = Vector2(-self.mData[1], self.mData[0])
        return Vperp


class Vector3(VectorN):
    """This class is for 3d vectors."""

    @property
    def y(self):
        """Getter method for second data point in a vector list"""
        return self.mData[1]

    @y.setter
    def y(self, new):
        """Setter method for the second data point in a vector list, will change the value of it."""
        self.mData[1] = float(new)

    @property
    def z(self):
        """Getter method for z, only useful for 3d vectors, so its not going to be in 2d vectors."""
        return self.mData[2]

    @z.setter
    def z(self, new):
        """Setter method for z, will take whatever value is passed and change the existing value for z."""
        self.mData[2] = new

    def rotateX(self, angle):
        y = self.y * math.cos(angle) - self.z * math.sin(angle)
        z = self.y * math.sin(angle) + self.z * math.cos(angle)
        return Vector3(self.x, y, z).normalized()

    def rotateY(self, angle):
        x = self.x * math.cos(angle) - self.z * math.sin(angle)
        z = self.x * math.sin(angle) + self.z * math.cos(angle)
        return Vector3(x, self.y, z).normalized()

    def rotateZ(self, angle):
        x = self.x * math.cos(angle) - self.y * math.sin(angle)
        y = self.x * math.sin(angle) + self.y * math.cos(angle)
        return Vector3(x, y, self.z).normalized()

def polar_to_vector2(radian, hypo, invertY):
    """Takes a radian and sends it through sin and cosine
     and times it by the hypotenuse to convert it to vector coordinates."""
    if invertY:
        o = -hypo * math.sin(radian)
        a = hypo * math.cos(radian)
    else:
        o = hypo * math.sin(radian)
        a = hypo * math.cos(radian)
    return Vector2(a, o)

def dot(vector, other):
    """Takes two vectors and gives back the dot product, errors if the other is not a vector"""
    if isinstance(other, VectorN):
        output = 0
        if len(vector) == len(other):
            for i in range(len(vector)):
                raw = vector[i] * other[i]
                output += raw
            return output
    else:
        raise ValueError("You can only take the dot product of two Vectors, you input: " + str(other))

def cross(vector, other):
    """Takes two Vector3's and crosses them, will error if you input any other kind of vector."""
    if isinstance(other, Vector3):
        if len(vector) == len(other) and len(vector) == 3:
            output_x = vector[1] * other[2] - vector[2] * other[1]
            output_y = vector[2] * other[0] - vector[0] * other[2]
            output_z = vector[0] * other[1] - vector[1] * other[0]
            cross_product = Vector3(output_x, output_y, output_z)
            return cross_product
    else:
        raise ValueError("You can only take the cross product of two Vector 3's, you input: " + str(other))

