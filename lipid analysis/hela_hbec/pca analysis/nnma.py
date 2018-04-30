import numpy as NP
def cf(q, v):
    """ the cost function """
    qv = (q - v)**2
    return NP.sum(NP.sum(qv, axis=0))


def nnma(d, max_iter=1000):
    x, y = d.shape
    z = y
    w = NP.random.rand(x, y)
    h = NP.random.rand(y, z)
    for i in range(max_iter):
        wh = NP.dot(w, h)
        cost = cf(d, wh)
        if cost == 0: 
            break
        hn = NP.dot(w.T, d)
        hd = NP.dot(NP.dot(w.T, w), h)
        h *= hn/hd
        wn = NP.dot(d, h.T)
        wd = NP.dot(NP.dot(w, h), h.T)
        w *= wn/wd
    return NP.dot(w, h)

def main():
    file = open('species_mt.txt', 'rb')
    table = [row.strip().split('\t') for row in file]
    table = [[float(y) for y in x] for x in table]
    
    d = NP.array(table)
    d1 = nnma(d)
    #nmf_mdl = pymf.NMF(d, num_bases=2, niter=1000)
    #nmf_mdl.initialization()
    #nmf_mdl.factorize()
    #d1 = NP.dot(nmf_mdl.W, nmf_mdl.H)
    NP.savetxt('mt_result.txt', d1)
    #print(d1)

main()
