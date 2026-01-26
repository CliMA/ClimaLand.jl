from sympy import *

n = 6
r = 0

def w(r, x):
    return (x**r) * log(1/x)

half = Rational(1,2)
a = [0, half]
b = [0, 1]
for l in range(2, 2*n+1):
    a.append(half)
    b.append( (l-1) / (2*sqrt((2*l-1)*(2*l-3))) ) 

x, y = symbols('x y')

phat = [1, (x - a[1])/b[2]]
for l in range(2,2*n):
    phat.append( ((x-a[l])*phat[l-1] - b[l]*phat[l-2])/b[l+1] )

p = [1, x - a[1]]
for l in range(2,2*n):
    p.append((x-a[l])*p[l-1] - (b[l]**2)*p[l-2])

q = [1]
nrmq = [ sqrt(integrate(q[0]**2 * w(r,x), (x,0,1))) ]
qhat = [ 1/nrmq[0] ]
alpha = [ 0, integrate(x *qhat[0]**2 * w(r,x), (x,0,1)) ]
beta = [ 0, sqrt(integrate(w(r,x),(x,0,1))) ]
q.append( x-alpha[1] )
nrmq.append( sqrt(integrate(q[1]**2 * w(r,x), (x,0,1))) )
qhat.append( q[1]/nrmq[1] )
alpha.append( integrate(x * qhat[1]**2 * w(r,x), (x,0,1)) )
beta.append( nrmq[1]/nrmq[0] )
for l in range(2,n):
    q.append( ((x-alpha[l])*q[l-1] - beta[l]**2 * q[l-2]) )
    nrmq.append( sqrt( integrate(q[l]**2 * w(r,x), (x,0,1)) ) )
    qhat.append( q[l] / nrmq[l] )
    alpha.append( integrate(x * qhat[l]**2 * w(r,x), (x,0,1)) )
    beta.append( nrmq[l] / nrmq[l-1] )

def chk_orthp():
    for k in range(0,2*n):
        s = integrate(phat[k]*phat[k], (x, 0, 1))
        assert(s==1)
        s = integrate(p[k]*p[k], (x, 0, 1))
        assert(s*(2*k+1)*binomial(2*k,k)**2==1)
        for l in range(k+1,2*n):
            s = integrate(phat[k]*phat[l], (x, 0, 1))
            assert(s==0)
            s = integrate(p[k]*p[l], (x, 0, 1))
            assert(s==0)

def chk_orthq(m):
    for k in range(0,m):
        s = integrate(qhat[k]*qhat[k]*w(r,x), (x, 0, 1))
        assert(s==1)
        for l in range(k+1,m):
            s = integrate(q[k]*q[l]*w(r,x), (x, 0, 1))
            assert(s==0)
            s = integrate(qhat[k]*qhat[l]*w(r,x), (x, 0, 1))
            assert(s==0)
    
         

nu = [0]
for l in range(1,2*n):
    nu.append( integrate(phat[l-1]*w(r,x), (x, 0, 1)) )

def showsigma(k):
    for l in range(k, 2*n-k):
        sigma = integrate(phat[l-1]*qhat[k-1]*w(r,x), (x,0,1))
        print("sigma_{0},{1} = {2:10.6f}".format(l, k, float(sigma)))

def show_coefs():
    print("{0:2s}  {1:10s}  {2:10s}  {3:10s}  {4:10s}".format(
          "k", "alpha_k", "beta_k", "a_k", "b_k"))
    for k in range(1, n):
        print("{0:2d}  {1:10.6f}  {2:10.6f}  {3:10.6f}  {4:10.6f}".format(
              k, float(alpha[k]), float(beta[k]), float(a[k]), float(b[k])))
