class beyond():
    def __init__(self, num, exp):
        self.num = num ## number for example 6.434E#, the 6.434 is the number, and should be converted to a interger
        self.exp = exp ## After the decimal. which is how percise I want the number ot be

    def __str__(self):
        return f"{self.num}e{self.exp}"
    def add(self, other):
        if self.exp == other.exp: ## Makes sure the exponent is the same
            num1 = int(self.num * 10**self.exp)
            num2 = int(other.num * 10**other.exp)
            result = (num1 + num2) / (10**self.exp)
            return beyond(result, self.exp)
        else:
            a = [self.exp, other.exp]
            a = max(a)
            num1 = int(self.num * 10**self.exp)
            num2 = int(other.num * 10**other.exp)
            result = (num1 + num2) / (10 **a)
            return beyond(result, a)
        return
    ## defining how multipliation works
    def times(self, other):
        a = self.exp + other.exp ## The resulting exp
        num1 = int(self.num * 10 **self.exp)
        num2 = int(other.num * 10 ** other.exp)
        result = (num1 * num2) / (10 ** a)
        return beyond(result, a)

    def div(self, other):
        a = self.exp - other.exp ##
        num1 = str(self.num)
        num2 = str(other.num)
        res = str(self.num / other.num)
        tmp = [len(num1), len(num2)]
        tmp = min(tmp)
        res = float(res[:tmp]) ## Truncating the float
        return beyond(res,a)


if __name__ == '__main__':
    x1 = beyond(2.29,4)
    x2 = beyond(1.59,10)
    x = x1.div(x2)
    print(f"Divide: {x}")
    x1 = beyond(2.29,4)
    x2 = beyond(-1.59,10)
    x = x1.add(x2)
    print(f"sub: {x}")
    x1 = beyond(2.29,4)
    x2 = beyond(1.59,10)
    x = x1.add(x2)
    print(f"add: {x}")
    x1 = beyond(2.29,4)
    x2 = beyond(1.59,10)
    x = x1.times(x2)
    print(f"times: {x}")

    ## Testing the accuracy
    num = beyond(4.32, 10)
    num1 = beyond(1,-4)
    num2 = beyond(1,2)
    a = num1.times(num2)
    final_res = num.add(a)
    print(f"Result before small increment: {final_res}")
    N = 100
    for i in range(N):
        num1 = num1.add(num2)
    print(f"Result after small increments: {num1}")