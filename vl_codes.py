from math import log2, ceil
from random import random
from trees import xtree2newick

def shannon_fano(p):
    """
    Shannon-Fano algorithm for constructing a prefix-free code from a probability distribution

    Prameters
    ---------
    p : dict [label : str, probability : float]
        Probability distribution

    Returns
    -------
    code : dict [label : str, codeword : list [0,1]]
        Code table

    """

    p = dict(sorted([(a,p[a]) for a in p if p[a]>0.0], key = lambda el: el[1], reverse = True))
    first_key = list(p.keys())[0]
    if p[first_key] == 0:
        p.pop(first_key)

    del first_key

    temp = 0
    f = {}
    for key, value in p.items():
        f[key] = temp
        temp += value

    code = {} 
    for symbol, prob in p.items(): # for each probability
        length  = ceil(-log2(prob))

        codeword = [] # initialise current codeword
        myf = f[symbol]
        for _ in range(length):
            myf *= 2
            if myf >= 1 :
                codeword.append(1)
                myf -=1
            else:
                codeword.append(0)

        code[symbol] = codeword # assign the codeword
        
    return code # return the code table


def huffman(p):
    """
    Huffman's algorithm for constructing a prefix-free code from a probability distribution
    
    Prameters
    ---------
    p : dict {[label : str, probability : float]}
        Probability distribution
    
    Returns
    -------
    code : list [parent(index): int , [ childeren(index): int ] , label: str]
        coding tree in xtree format
    """

    xt = [[-1,[], a] for a in p]
    p = [(k,p[a]) for k,a in zip(range(len(p)),p)]

    nodelabel = len(p)

    while len(p) > 1:
        p = sorted(p, key = lambda el:el[1])

        xt.append([-1,[], str(nodelabel)])

        # we incrase the variable nodelabel by 1 for its next use
        nodelabel += 1

        xt[p[0][0]][0] = len(xt)-1
        xt[p[1][0]][0] = len(xt)-1
        
        xt[-1][1] = [p[0][0], p[1][0]]

        p.append((len(xt)-1, p[0][1]+p[1][1]))

        p.pop(0)
        p.pop(0)
        
    return(xt)



def bits2bytes(x):
    n = len(x)+3
    r = (8 - n % 8) % 8
    prefix = format(r, '03b')
    x = ''.join(str(a) for a in x)
    suffix = '0'*r
    x = prefix + x + suffix
    x = [x[k:k+8] for k in range(0,len(x),8)]
    y = []
    for a in x:
        y.append(int(a,2))

    return y

def bytes2bits(y):
    x = [format(a, '08b') for a in y]
    r = int(x[0][0:3],2)
    x = ''.join(x)
    x = [int(a) for a in x]
    for k in range(3):
        x.pop(0)
    for k in range(r):
        x.pop()
    return x


def vl_encode(x, c):
    y = []
    for a in x:
        y.extend(c[a])
    return y


def vl_decode(y, xt):
    x = []
    root = [k for k in range(len(xt)) if xt[k][0]==-1]
    if len(root) != 1:
        raise NameError('Tree with no or multiple roots!')
    root = root[0]
    leaves = [k for k in range(len(xt)) if len(xt[k][1]) == 0]

    n = root
    for k in y:
        if len(xt[n][1]) < k:
            raise NameError('Symbol exceeds alphabet size in tree node')
        if xt[n][1][k] == -1:
            raise NameError('Symbol not assigned in tree node')
        n = xt[n][1][k]
        if len(xt[n][1]) == 0: # it's a leaf!
            x.append(xt[n][2])
            n = root
    return x


if __name__ == "__main__":
    p = [random() for k in range(16)]
    p = dict([(chr(k+ord('a')),p[k]/sum(p)) for k in range(len(p))])
    print(f'Probability distribution: {p}\n')
    order_p = sorted(p, key = lambda x: p[x], reverse = True)
    print(order_p)
    # code_f = shannon_fano(p)
    code_h = huffman(p)
    print(xtree2newick(code_h))
