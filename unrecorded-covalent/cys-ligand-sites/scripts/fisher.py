from math import lgamma, exp, log
def _lp(a,b,c,d):
    n=a+b+c+d
    return (lgamma(a+b+1)+lgamma(c+d+1)+lgamma(a+c+1)+lgamma(b+d+1)-lgamma(n+1)-lgamma(a+1)-lgamma(b+1)-lgamma(c+1)-lgamma(d+1))
def fisher(a,b,c,d):
    """two-sided Fisher exact p for [[a,b],[c,d]]; returns (odds_ratio, p)"""
    r1,c1,n=a+b,a+c,a+b+c+d
    p0=_lp(a,b,c,d); p=0.0
    for x in range(max(0,r1+c1-n),min(r1,c1)+1):
        lp=_lp(x,r1-x,c1-x,n-r1-c1+x)
        if lp<=p0+1e-9: p+=exp(lp)
    orr=(a*d)/(b*c) if b*c else float('inf')
    return orr,min(p,1.0)
