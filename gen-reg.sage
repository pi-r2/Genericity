def complete_basis(B, E) :
    """
    This function takes as input a list of linearily independant vectors B in E and 
    completes them into a basis naively.
    If the list is empty, a random basis of E is returned
    """
    res = list(B) #Avoid modifying the input.
    if len(res) == 0 :
        b = E.random_element()
        while b.is_zero() :
            b = E.random_element()
        res.append(b)

    while len(res) != E.dimension() :
        b = E.random_element()
        if not(b in span(res)) :
            res.append(b)
    return res

    


def KeyGen(q,k,v, m=None, verbose = False, triangular = False):
    """
    This function generates a key pair for UOV parameters k,v where k is the dimension of the oil subspace and v the dimension of the vinegar subspace.
    To do this, we sample random matrices with a block of zero of size to produce the private key, and impose that they must be invertible.
    We sample a random invertible change of variables A, compute G = F\circ A, and return the pair (A,F), G.
    (A,F) is the UOV private key, G is the UOV public key. The corresponding systems are obtained by evaluating the quadratic forms over FF[x1,...,xn].
    Important notice: This code is a demonstration tool and should not used for applications where security matters. 
    """
    if m is None :
        m = k
    n = k+v 
    FF = GF(q)
    RR = None 
    if verbose :
        RR = PolynomialRing(FF, 'x', n)
    #Define A

    A = matrix.identity(FF,n)
    for i in range(n-k) :
        for j in range(k) :

            A[n-1-i,j] = FF.random_element()
    #print(A)
    #Define F, G
    F = []
    G = []
    for e in range(m):
        Fe = matrix(FF, n, n)
        while Fe.determinant() == 0 :
            B0 = matrix(FF,k,k)
            B1 = random_matrix(FF,k,v)
            B2 = random_matrix(FF,v,k)
            B3 = random_matrix(FF,v,v)
            Fe = block_matrix([[B0,B1],[B2,B3]])
        
        if FF.characteristic() != 2 : 
            Fe = FF(2)**(-1)*(Fe+Fe.transpose())
        Ge = A.transpose()*Fe*A
        F.append(Fe)
        G.append(Ge)
        
    if verbose :
        print("The public key system is:")
        x = vector(RR.gens())
        for Ge in G :
            print(x*Ge*x)
    return ((A,F), G) 

    
def ToMSolve(F, finput="/tmp/in.ms"): #From msolve library interfaces
    """Convert a system of sage polynomials into a msolve input file.

    Inputs :
    F (list of polynomials): system of polynomial to solve
    finput (string): name of the msolve input file.

    """
    A = F[0].parent()
    assert all(A1 == A for A1 in map(parent,F)),\
            "The polynomials in the system must belong to the same polynomial ring."
    variables, char = A.variable_names(), A.characteristic()
    s = (", ".join(variables) + " \n"
            + str(char) + "\n")

    B = A.change_ring(order = 'degrevlex') 
    F2 = [ str(B(f)).replace(" ", "") for f in F ]
    if "0" in F2:
        F2.remove("0")
    s += ",\n".join(F2) + "\n"

    fd = open(finput, 'w')
    fd.write(s)
    fd.close()

def FromMsolve(output,R):
    with open(output, 'r') as o :
        sols = []
        for l2 in o.readlines()[2:] :
            l = ''
            for c in l2 :
                if c not in ['[', ']', '\n', ',' , ':'] :
                    l+=c    
            if l[0] == '#' :
                continue
            sols.append(R(l))
    return sols


### Define a collection of primes to test.

primes =  [2] + [31, 251, 1021, 4093]

### Define families of UOV parameters with different expected dimension of the singular locus, and of various rank defects.

parameters = [
    #n, m, o. In all cases, n > 2o.
    [10, 4, 4], #Regular UOV               
    [10, 5, 4], #Regular UOV               
    [10, 5, 5], #O has the same dim as V                  
    [12, 5, 5], #Regular UOV                
    [12, 6, 5], #Regular UOV                
    [12, 7, 5], #O has the same dim as V    
    [12, 6, 6], #O has the same dim as V    
    [14, 6, 6], #Regular UOV                
    [14, 7, 6], #Regular UOV                
    [14, 7, 7], #O has the same dim as V    
    [16, 7, 7], #Regular UOV               
    [16, 8, 7], #Regular UOV
]


count=1000
with open("reg.log",'w') as fout:
    for n,m,o in parameters:
        print("Parameter set:", n,m,o, file=fout)
        print("Parameter set:", n,m,o)        
        for p in primes : 
            print("Prime:", p, file=fout)
            print("Prime:", p)
            success = 0
            for _ in range(count):
                _, G = KeyGen(p,o,n-o, m)
                
                R = PolynomialRing(GF(p), 'X', n)
                X = vector(R.gens())
                eqs = [X*g*X for g in G]
                file = "/tmp/reg"+str(m)+"_"+str(r)+"_"+str(p)
                ToMSolve(eqs, file+".ms")
                os.system("./msolve -t16 -g1 -f "+file+".ms -o "+file)#+" > "+file+".log")
                gb= FromMsolve(file, R)
                I = Ideal([i.lt() for i in gb])
                hs_e = I.hilbert_series() 
                z = hs_e.parent().gens()[0]
                hs_t = (1-z^2)^m / (1-z)^n 
                if hs_e==hs_t:
                    success+=1
            print("Successes:", success, file=fout)
            print("Ratio:", round((success/count)*100, 1), file=fout)
        print( )
        
