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

    

K.<x> = GF(16)
field_dict = {};
field_dict[0] = K(0)
field_dict[1] = K(1)
field_dict[2] = K(x)
field_dict[3] = K(x+1)
field_dict[4] = K(x*x)
field_dict[5] = K(x*x+1)
field_dict[6] = K(x*x+x)
field_dict[7] = K(x*x+x+1)
field_dict[8] = K(x**3)
field_dict[9] = K(x**3 + 1)
field_dict[10] = K(x**3 + x)
field_dict[11] = K(x**3 + x+1)
field_dict[12] = K(x**3 + x*x)
field_dict[13] = K(x**3 + x*x+1)
field_dict[14] = K(x**3 + x*x+x)
field_dict[15] = K(x**3 + x*x+x+1)
inverse_field_dict = { field_dict[x]:x for x in field_dict }


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
    if q == 16 :
        FF = K
    else:
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

    
def ToMSolve(F, finput="/tmp/in.ms"): #From msolve library interfaces https://github.com/algebraic-solving/msolve
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



def ToM4GB(F, finput="/tmp/in.ms"): 
    """
    Interface found in https://github.com/WardBeullens/BreakingSNOVA/blob/master/Attack_demo.sage
    
    Convert a system of sage polynomials into a m4gb input file.

    Inputs :
    F (list of polynomials): system of polynomial to solve
    finput (string): name of the input file.

    """
    A = F[0].parent()
    q = A.base_ring().cardinality()
    vars = vector(A.gens())
    n = len(A.gens())
    # write file for m4gb
    with open(finput, 'w') as f:
        f.write(f'$fieldsize {q}\n')
        f.write('$vars '+" ".join(A.variable_names())+'\n')    
        for equation in F:
            line = ""
            for c,M in equation:
                if q ==16:
                    line += str(inverse_field_dict[c])+"*"+str(M)+" + "
                else:
                    line += str(c)+"*"+str(M)+" + "
            line = line[:-2]
            f.write(line+"\n")


def FromM4GB_F16(output, R, LT = False):
    """
    Interface adapted from https://github.com/WardBeullens/BreakingSNOVA/blob/master/Attack_demo.sage
    """
    with open(output, 'r') as file:
        sols = []
        
        q = R.base_ring().cardinality()
        for line in file:
            monomials = line.split("+")
            poly = 0
            for m in monomials :
                m = m.split("*")
                try:
                    poly +=field_dict[int(m[0])]*product(map(R,m[1:]))
                except Exception as e:
                    poly += product(map(R, m))
            if not LT:  
                sols.append(poly) 
            else :
                sols.append(poly.lt()) 
        return sols

 


### Define a collection of fields to test.

fields =  [16] + [31, 251, 1021, 4093]

### Define families of UOV parameters with different expected dimension of the singular locus, and of various rank defects.

parameters = [
    #n, m, o. In all cases, n > 2o.
    [10, 4, 4], #Regular UOV               ; r=1 dim 1
    [10, 5, 5], #Regular UOV               ; r=1 dim 1
    [12, 5, 5], #O has the same dim as V   ; r=1 dim 3 
    [12, 7, 5], #O has the same dim as V   ; r=1 dim 3 
    [12, 6, 6], #O has the same dim as V   ; r=1 dim 3 
    [14, 6, 6], #O has the same dim as V   ; r=1 dim 3 
    [14, 8, 6], #O has the same dim as V   ; r=1 dim 3 
    [14, 7, 7], #O has the same dim as V   ; r=1 dim 3 
    [16, 8, 8],  
]

count=10
with open("testgen-det-f16.log",'w') as fout:
    for n,m,o in parameters:
        print("Parameter set:", n,m,o, file=fout)
        print("Parameter set:", n,m,o)        
        r = min(1, n-o)
        delta = r*(m - r) + o - r*(n - o)
        while delta > 0 and r <= min(n-o, m):
            print("rank defect:",r, "dimension", delta, file=fout )
            print("rank defect:",r, "dimension", delta)
            
            for p in fields : 
                print("Field size:", p, file=fout)
                print("Field size:", p)
                success = 0
                for _ in range(count):
                    _, G = KeyGen(p,o,n-o, m)
                    
                    vars_X = ['X'+str(i) for i in range(n)]
                    vars_Y = ['Y'+str(i) for i in range((m-r)*r)]
                    variables = vars_X + vars_Y
                    if p == 16:
                        R = PolynomialRing(K, variables, n + (m-r)*r)
                    else:
                        R = PolynomialRing(GF(p), variables, n + (m-r)*r)
                    X = vector(R.gens()[:n])
                    Y = block_matrix(R, 1, 2, [
                        matrix.identity(r), matrix(R, r, m-r, R.gens()[n:])
                    ])
                    eqs = [X*g*X for g in G]
                    J = matrix([(g+g.T)*X for g in G])
                    if p == 16:
                        H = [K.random_element() + sum([K.random_element()*xi for xi in X]) for _ in range(delta)]
                    else:
                        H = [GF(p).random_element() + sum([GF(p).random_element()*xi for xi in X]) for _ in range(delta)]
                    
                    sys = eqs + (Y*J).coefficients() + H
                    file = "/tmp/sing"+str(m)+"_"+str(r)+"_"+str(p)
                    if p.is_prime():
                        ToMSolve(sys, file+".ms")
                        os.system("msolve -t4 -g1 -f "+file+".ms -o "+file)#+" > "+file+".log")
                        gb= FromMsolve(file, R)
                        I = Ideal(gb)
                    else:                    
                        ToM4GB(sys,file+".in")
                        os.system(f"cd ~/m4gb-master && ./solver.sh -v 0 -t 8 -n {n + (m-r)*r} {file}.in {file}.ans")
                        gbM4 = FromM4GB_F16(file+".ans",R, LT=True)
                        I = Ideal(gbM4)
                    if I.dimension() == 0 :
                        success+=1
                    else:
                        print(I.dimension())    
                print("Successes:", success, file=fout)
                print("Ratio:", round((success/count)*100, 1), file=fout)
            r+=1
            delta = r*(m - r) + o - r*(n - o)
        print( )
        
