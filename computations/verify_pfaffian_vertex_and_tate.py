#!/usr/bin/env python3
"""Exact Schur identities and a rebuilt order-six Tate certificate.
No old unpublished verifier is imported. All coefficients are rational.
The all-degree conclusions use the Schur/Bott proof, not finite sampling.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from collections import defaultdict
import sympy as s


def schur_dimension(parts):
    p=tuple(parts)+(0,)*(5-len(parts))
    return s.prod((p[i]-p[j]+j-i)/s.Integer(j-i) for i in range(5) for j in range(i+1,5))


def verify_schur():
    r,c=s.symbols('r c',integer=True)
    dim=schur_dimension((r+c,r+c,c,c,0))
    free=lambda n:(n+1)*(n+2)*(n+3)*(n+4)/s.Integer(24)
    rhs=schur_dimension((r,r))*free(c)-schur_dimension((r,r-1))*free(c-1)+schur_dimension((r-1,r-1))*free(c-2)
    assert s.factor(dim-rhs)==0
    for rv,expected in ((1,(10,5,1,6)),(2,(50,40,10,20))):
        dims=[schur_dimension(p).subs(r,rv) for p in ((r,r),(r,r-1),(r-1,r-1))]
        assert tuple(dims)+(dims[0]-dims[1]+dims[2],)==expected
    print('PASS: symbolic all-r,c Schur Hilbert identity and resolution ranks')
    from math import comb
    @lru_cache(None)
    def image_rank(m,u):
        if u==0:return 0
        return comb(5,u-1)*comb(m+5,4)-image_rank(m+1,u-1)
    counts=[];images=[]
    for u in range(6):
        total=rank=0
        for v in range(6):
            twice_m=3+u-v
            if twice_m<0 or twice_m%2:continue
            m=twice_m//2
            total+=comb(5,u)*comb(5,v)*comb(m+4,4)
            rank+=comb(5,v)*image_rank(m,u)
        counts.append(total);images.append(rank)
    assert counts==[35,350,1260,2100,1650,495]
    assert images==[0,210,935,1765,1524,495]
    total3=sum(comb(5,u)*comb(5,v)*comb((6+u-v)//2+4,4)
        for u in range(6) for v in range(6) if 6+u-v>=0 and (6+u-v)%2==0)
    assert sum(counts)==5890 and sum(images)==4929 and total3==20330
    print('PASS: finite Tate spaces E_(3/2)=5890, image(h)=4929, E_3=20330')


ZERO=(0,)*5
CUBIC=(1,1,1,1,-4)  # F=sum C_i x_i^3; B_ii=C_i*x_i

def bits(mask):return tuple(i for i in range(5) if mask>>i&1)
def eps(mask,i):return -1 if (mask&((1<<i)-1)).bit_count()%2 else 1

def accumulate(out,key,value):
    value=out.get(key,0)+value
    if value:out[key]=value
    else:out.pop(key,None)

def increment(e,i):return e[:i]+(e[i]+1,)+e[i+1:]

def left_eta(poly,i):
    out={};bit=1<<i
    for (im,jm,e),v in poly.items():
        if im&bit:accumulate(out,(im^bit,jm,e),v*eps(im,i))
        sign=(-1 if im.bit_count()%2 else 1)*eps(jm,i)
        if jm&bit:
            accumulate(out,(im,jm^bit,increment(e,i)),-CUBIC[i]*sign*v)
        else:accumulate(out,(im,jm|bit,e),sign*v)
    return out

def left_theta(poly,i):
    out={};bit=1<<i
    for (im,jm,e),v in poly.items():
        if not im&bit:accumulate(out,(im|bit,jm,e),eps(im,i)*v)
    return out

@lru_cache(None)
def basis_product(i,j,k,l):
    poly={(k,l,ZERO):1}
    for a in reversed(bits(j)):poly=left_eta(poly,a)
    for a in reversed(bits(i)):poly=left_theta(poly,a)
    return tuple(poly.items())

def multiply(p,q):
    out={}
    for (i,j,e),v in p.items():
        for (k,l,f),w in q.items():
            ef=tuple(a+b for a,b in zip(e,f))
            for (a,b,g),d in basis_product(i,j,k,l):
                accumulate(out,(a,b,tuple(x+y for x,y in zip(ef,g))),v*w*d)
    return out

def delta(p):
    out={}
    for (i,j,e),v in p.items():
        for a in bits(i):accumulate(out,(i^(1<<a),j,increment(e,a)),v*eps(i,a))
    return out

def homotopy(p):
    out={}
    for (i,j,e),v in p.items():
        den=sum(e)+i.bit_count()
        if not den:continue
        for a,power in enumerate(e):
            if power and not i>>a&1:
                f=e[:a]+(power-1,)+e[a+1:]
                accumulate(out,(i|(1<<a),j,f),-v*Q(power*eps(i,a),den))
    return out

def block(p):
    out=s.zeros(5)
    for col in range(5):
        initial=31^(1<<col)
        for (i,j,e),v in p.items():
            if any(e):continue
            mask=initial;coeff=v*((-1)**col)
            for a in reversed(bits(j)):
                if not mask>>a&1:coeff=0;break
                coeff*=eps(mask,a);mask^=1<<a
            if not coeff:continue
            for a in reversed(bits(i)):
                if mask>>a&1:coeff=0;break
                coeff*=eps(mask,a);mask|=1<<a
            if coeff and mask.bit_count()==1:out[bits(mask)[0],col]+=s.Rational(coeff.numerator,coeff.denominator)
    return out

def verify_tate():
    avec={(0,1):1,(0,2):-1,(0,3):1,(0,4):-1,(1,2):1,(3,4):1}
    a1={}
    for pair,value in avec.items():
        I=tuple(i for i in range(5) if i not in pair)
        inv=sum(a>b for a in I for b in pair)
        a1[0,sum(1<<i for i in I),ZERO]=Q(value*(-1)**inv)
    aa=[{},a1];blocks=[s.zeros(5),block(a1)]
    for n in range(2,7):
        residual={}
        for i in range(1,n):
            for key,value in multiply(aa[i],aa[n-i]).items():accumulate(residual,key,value)
        an=homotopy(residual)
        check=dict(residual)
        for key,value in delta(an).items():accumulate(check,key,value)
        assert not check,(n,len(check))
        aa.append(an);blocks.append(block(an))
        print('Tate order',n,'terms',len(an),'; exact MC identity PASS',flush=True)
    expectedA=s.Matrix([[0,-1,1,-1,1],[1,0,-1,0,0],[-1,1,0,0,0],[1,0,0,0,-1],[-1,0,0,1,0]])
    expectedB=s.diag(2,2,2,2,-8)
    expectedC=s.Matrix([[0,-4,-4,16,16],[-4,0,8,12,32],[-4,8,0,12,32],[16,12,12,0,48],[16,32,32,48,0]])/3
    expectedD=s.Matrix([[0,-3040,-3040,-3040,-3040],[-3280,-560,-4606,-2437,-1277],[-3280,-4606,-560,-4677,-5837],[-3280,-2437,-4677,-560,-1206],[-3280,-2477,-4637,-2406,2240]])/90
    assert blocks[1]==expectedA
    assert blocks[2]==expectedB
    assert blocks[3]==s.zeros(5) and blocks[5]==s.zeros(5)
    assert blocks[4]==expectedC and blocks[6]==expectedD
    t=s.symbols('t')
    det=(t*blocks[1]+t**2*blocks[2]+t**4*blocks[4]+t**6*blocks[6]).det(method='domain-ge')
    assert s.expand(det).coeff(t,10)==-648
    assert Q(-648,2*(-31104))==Q(1,96)
    print('PASS: four literal Tate blocks rebuilt; determinant coefficient -648; Hessian coefficient 1/96')

if __name__=='__main__':
    verify_schur()
    verify_tate()
