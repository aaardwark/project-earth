class Variable:
    def __init__(self, identifier):
        # alphanumeric case-sensitive, leading character from a-zA-Z
        self.name = str(identifier)
    def __repr__(self):
        return f'Variable({self.name})'
    def __str__(self):
        return self.name
    def __lt__(self, other):
        return self.name < other.name
    
class Product:
    def __init__(self, coefficient, *factor_list):
        # sorted list of variables and a real coefficient
        # repeated variable is equivalent to exponentiation
        self.factors = sorted(list(factor_list))
        self.coeff = float(coefficient)
    def __repr__(self):
        return f'SimpleProduct({self.coeff},*{self.factors})'
    def __str__(self):
        output = f'{self.coeff}*'
        for f in self.factors:
            output += str(f) + '*'
        return output[:-1]
    
class Sum:
    def __init__(self, *term_list):
        # list of simpleproducts
        # single variables should be SimpleProduct(var)
        # repeated terms should be collected into a simpleproduct with integer coefficient
        self.terms = list(term_list)
    def __repr__(self):
        return f'Sum(*{self.terms})'
    def __str__(self):
        output = ''
        for f in self.terms:
            output += str(f) + ' + '
        return output[:-3]


def pvar(var):
    return Product(1,var)
def svar(var):
    return Sum(Product(1,var))

def scale_prod(product, real):
    return Product(float(real)*product.coeff, *product.factors)
def scale_sum(sum, real):
    return Sum(*[scale_prod(term, real) for term in sum.terms])
def mult_sum_prod(sum, product):
    return Sum(*[mult_2prods(term, product) for term in sum.terms])


def mult_2prods(product1, product2):
    comb_factors = sorted(product1.factors + product2.factors)
    comb_coeff = product1.coeff * product2.coeff
    return Product(comb_coeff, *comb_factors)
def mult_prods(*products):
    if len(products) == 0:
        raise ValueError('Cannot multiply zero Products')
    final_product = products[0]
    for p in products[1:]:
        final_product = mult_2prods(final_product,p)
    return final_product



def add_2sums(sum1, sum2):
    terms1_factors = [term.factors for term in sum1.terms]
    terms2_factors = [term.factors for term in sum2.terms]
    combined_terms = []

    for index1, factor_list in enumerate(terms1_factors):
        try:
            index2 = terms2_factors.index(factor_list)
            terms2_factors[index2] = None

            collected_product = Product(sum1.terms[index1].coeff + sum2.terms[index2].coeff, *factor_list)
            combined_terms.append(collected_product)
        except ValueError:
            combined_terms.append(sum1.terms[index1])
    
    for index2, factor_list in enumerate(terms2_factors):
        if factor_list:
            combined_terms.append(sum2.terms[index2])

    return Sum(*combined_terms)
def add_sums(*sums):
    if len(sums) == 0:
        raise ValueError('Cannot add zero Sums')
    final_sum = sums[0]
    for p in sums[1:]:
        final_sum = add_2sums(final_sum,p)
    return final_sum


def mult_2sums(sum1, sum2):
    sum_list = []
    for term_from1 in sum1.terms:
        multiplied = [mult_2prods(term_from1,term_from2) for term_from2 in sum2.terms]
        sum_list.append(Sum(*multiplied))

    final_sum = sum_list[0]
    for n in range(1, len(sum_list)):
        final_sum = add_2sums(final_sum, sum_list[n])
    return final_sum
def mult_sums(*sums):
    if len(sums) == 0:
        raise ValueError('Cannot multiply zero Sums')
    final_sum = sums[0]
    for p in sums[1:]:
        final_sum = mult_2sums(final_sum,p)
    return final_sum

####
def save(expression, filename):
    with open(f'/Users/arundhati/Downloads/HW_Parijatha/ProjectEarth/{filename}.txt', 'w') as writeto:
        import re
        def strip_dot0s(match):
            return match[1] + '*'
        
        expr_str = str(expression)
        expr_str = re.sub('\+ \-','- ',expr_str)
        expr_str = re.sub('1\.0\*', '', expr_str)
        expr_str = re.sub('(\d+)\.0\*', strip_dot0s, expr_str)

        writeto.write(expr_str)


##########
x,y,z = Variable('x'), Variable('y'), Variable('z')
a1,a2 = Variable('a1'), Variable('a2')
b1,b2 = Variable('b1'), Variable('b2')
c1,c2 = Variable('c1'), Variable('c2')
K1,K2 = Variable('K1'), Variable('K2')
s1,s2 = Variable('s1'), Variable('s2')

p1 = Sum(Product(1,a1,x), Product(1,b1,y), Product(1,c1,z))
p2 = Sum(Product(1,a2,x), Product(1,b2,y), Product(1,c2,z))
# a1x + b1y + c1z | sum

q12 = Sum(Product(1,a1,a2), Product(1,b1,b2), Product(1,c1,c2))
# a1a2 + b1b2 + c1c2 | sum

print(mult_2sums(q12,q12))