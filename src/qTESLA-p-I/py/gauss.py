class GaussianInteger:
    re = 0
    imag = 0
    p = 343576577

    def __init__(self, re, imag = 0, p = None):
        assert type(re) == int
        assert type(imag) == int

        self.re = re
        self.imag = imag

        if type(p) == int:
            self.p = p
        elif p is not None:
            raise Exception(
                "Bad value of p received." + str(p)
                )

    def conjugate(self):
        return GaussianInteger(self.re, -self.imag)

    def __valid_operand(self, b):
        if type(b) == int:
            b = GaussianInteger(b)
        elif type(b) == float:
            raise Exception(
                "Error! We shouldn't be dealing with floats. Found " + str(b)
                )
        return b

    def __repr__(self):
        return "%d + i%d" % (self.re, self.imag)

    def __eq__(self, b):
        b = self.__valid_operand(b)

        return self.re == b.re and self.imag == b.imag

    # Galois add
    def __add__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re + b.re) % self.p,
            (self.imag + b.imag) % self.p
            )

    def __radd__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re + b.re) % self.p,
            (self.imag + b.imag) % self.p
            )

    # Galois sub
    def __sub__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re - b.re) % self.p,
            (self.imag - b.imag) % self.p
            )

    def __rsub__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re - b.re) % self.p,
            (self.imag - b.imag) % self.p
            )
    # Galois mul
    def __mul__(self, b):
        b = self.__valid_operand(b)
        # https://stackoverflow.com/questions/19621686/complex-numbers-product-using-only-three-multiplications
        # 
        # S1=ac,S2=bd, and S3=(a+b)(c+d). Now you can compute the results as 
        # A=S1−S2 and B=S3−S1−S2.
        # 
        s1 = self.re * b.re
        s2 = self.imag * b.imag
        s3 = (self.re + self.imag) * (b.re + b.imag) 

        return GaussianInteger(
            (s1 - s2) % self.p,
            (s3 - s1 - s2) % self.p 
            )

    def __rmul__(self, b):
        b = self.__valid_operand(b)

        s1 = self.re * b.re
        s2 = self.imag * b.imag
        s3 = (self.re + self.imag) * (b.re + b.imag) 

        return GaussianInteger(
            (s1 - s2) % self.p,
            (s3 - s1 - s2) % self.p
            )

    def __pow__(self, b):
        # Square and multiply
        if(b == 0):
            return GaussianInteger(1)

        exp = bin(b)[3:]
        value = self

        for i in range(len(exp)):
            value = value * value
            if(exp[i:i+1] == '1'):
                value = value * self
        return value

    # About divisions and remainders:
    # https://math.stackexchange.com/questions/889809/calculating-the-reminder-when-dividing-complex-numbers
    def __floordiv__(self, b):
        if b.imag == 0:
            b = b.re
        if type(b) == int:
            assert b != 0

            return GaussianInteger(
                self.re // b,
                self.imag // b
                )
        else:
            assert isinstance(b, GaussianInteger)

            # We don't want to reduce before the rounding
            s1 = self.re * b.conjugate().re
            s2 = self.imag * b.conjugate().imag
            s3 = (self.re + self.imag) * (b.conjugate().re + b.conjugate().imag) 

            num = GaussianInteger(
                (s1 - s2),
                (s3 - s1 - s2)
                )

            s1 = b.re * b.conjugate().re
            s2 = b.imag * b.conjugate().imag

            den = (s1 - s2)

            return GaussianInteger(
                int(round(num.re / den)),
                int(round(num.imag / den))
                )

    def __mod__(self, b):
        if type(b) == int:
            assert b != 0

            return GaussianInteger(
                int(self.re % b),
                int(self.imag % b)
                )
        else:
            assert isinstance(b, GaussianInteger)

            return self - (self // b) * b