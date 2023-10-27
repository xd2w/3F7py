from math import log2, ceil
from random import random

def shannon_fano(p):

    # Begin by sorting the probabilities in decreasing order, as required
    # in Shannon's paper.
    p = dict(sorted([(a,p[a]) for a in p if p[a]>0.0], key = lambda el: el[1], reverse = True))
    first_key = list(p.keys())[0]
    if p[first_key] == 0:
        p.pop(first_key)

    del first_key

    # Compute the cumulative probability distribution
    temp = 0
    f = {}
    for key, value in p.items():
        f[key] = temp
        temp += value

    # assign the codewords
    code = {} # initialise as an empty dictionary
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
    # create an xtree with all the source symbols (to be the leaves) initially orphaned
    xt = [[-1,[], a] for a in p]
    # label the probabilities with a "pointer" to their corresponding nodes in the tree
    # in the process, we convert the probability vector from a Python dictionary to a
    # list of tuples (so we can modify it)
    p = [(k,p[a]) for k,a in zip(range(len(p)),p)]

    # the leaves are labeled according to the symbols in the probability vector. We will
    # label the remaining tree nodes (mainly for visualisation purposes) with numbers
    # starting at len(p)
    nodelabel = len(p)

    # this loop will gradually increase the tree and reduce the probability list.
    # It will run until there is only one probability left in the list (which probability
    # will by then be 1.0)
    while len(p) > 1:
        # sort probabilities in ascending order (so [0] and [1] are smallest)
        # using the command "sorted" and, as its second element, the expression
        #  "key = lambda el:el[1]" (which is an inline function to retrieve the
        # second entry from a tuple.) Note that the natural order of "sorted" is
        # increasing, which is as you want it in this case.
        # The output of sorted() can be written back to p (i.e. p = sorted(p, ....))
        #...

        # Now append a new node to the tree xt with no parent (parent = -1),
        # no children (children = []) and label str(nodelabel)
        #...

        # we incrase the variable nodelabel by 1 for its next use
        nodelabel += 1

        # assign parent of the nodes pointed to by the smallest probabilities 
        # Note that the smallest probabilities are now p[0] and p[1], so their
        # pointers to nodes in xt are p[0][0] and p[1][0]. The corresponding
        # xt nodes should be assigned the new node you created as a parent,
        # whose index is len(xt)-1 since it has been appended at the end of xt
        #...
        #...
        
        # assign the children of new node to be the nodes pointed to by
        # p[0] and p[1]. Note that the new node can be addressed as xt[-1]
        #...

        # create a new entry pointing to the new node in the list of probabilities
        # This new entry should be a tuple with len(xt)-1 as its first element,
        # and the sum of the probabilities in p[0] and p[1] as its second element
        #...

        # remove the two nodes with the smallest probability
        p.pop(0)
        p.pop(0)
        # (using pop(0) twice removes the first element of the list twice
        # and hence removes the first 2 elements.)
        
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
    code = shannon_fano(p)
    print(code)
