from sympy import *
from sympy.utilities.codegen import codegen
import re
t  = Symbol("t")
pre,rho,pot = Symbol("pre"),Symbol("rho"),Symbol("pot")
n1, n2, n3, n4  = symbols("n1 n2 n3 n4")
#ILU = Function('ILU')(n1, n2, n3, n4)
x = Function('a1')(t)
y = Function('a2')(t)
z = Function('a3')(t)
th = Function('th')(t)
ph = Function('ph')(t)
r = Function('r')(t)
vr, vth, vph = symbols("vr vth vph")


x = r*sin(th)*cos(ph)
y = r*sin(th)*sin(ph)
z = r*cos(th)
print(latex(x,mode='plain'))
print(latex(y,mode='plain'))
print(latex(z,mode='plain'))
xpot = sin(th)*cos(ph)*pot
ypot = sin(th)*sin(ph)*pot
zpot = cos(th)*pot


#vx = diff(x,t)
#vy = diff(y,t)
#vz = diff(z,t)
vx = sin(th)*cos(ph)*vr+cos(th)*cos(ph)*vth - sin(ph)*vph
vy = sin(th)*sin(ph)*vr+cos(th)*sin(ph)*vth + cos(ph)*vph
vz = cos(th)*vr-sin(th)*vth
print(vx)
print(vy)
print(vz)

Q11 = rho*((vx*vx-x*xpot))
Q22 = rho*((vy*vy-y*ypot))
Q33 = rho*((vz*vz-z*zpot))
Q12 = rho*(2.0*vx*vy-x*ypot-y*xpot)
Q13 = rho*(2.0*vx*vz-x*zpot-z*xpot)
Q23 = rho*(2.0*vy*vz-y*zpot-z*ypot)


#print(latex(Q11,mode='plain'))
#print(latex(Q22,mode='plain'))
#print(latex(Q33,mode='plain'))
#print(latex(Q12,mode='plain'))
#print(latex(Q23,mode='plain'))
#print(latex(Q23,mode='plain'))


Qs = (Q11,Q22,Q33,Q12,Q13,Q23)
print(Q23)
f = open('out.tst','w')
for x in Qs:
   x = expand(x)
   qfs = 0
   for i in range(3):
       for j in range(3):
           for k in range(3):
               for q in range(3):
                   if((k != 1)):
                       Qtt = x.coeff(sin(th),i).coeff(cos(th),j).coeff(sin(ph),k).coeff(cos(ph),q)
                       qfs =qfs + Qtt*sympify("ILU(l1,m1,l2,m2,l3,m3,"+str(j)+","+str(i)+","+str(k)+","+str(q)+")")
               
   rg = re.sub(r'(\(t\))','',str(expand(qfs)))
   rg = re.sub(r'Derivative\(([a-z0-9]*)[^,]*, t\)',r'v\1',rg)
   rg = re.sub(r'([\+\-])',r'\1&\n',rg)
   rg = re.sub(r'(0\.666666666666667)',r'trd',rg)
   rg = re.sub(r'(1\.33333333333333)',r'frd',rg)
   rg = re.sub(r'(2\.66666666666667)',r'erd',rg)
   
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 2, 0, 0, 2\))',r'tmp1',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 0, 2, 0\))',r'tmp2',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 1, 0, 2\))',r'tmp3',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 2, 0, 2\))',r'tmp4',rg)

   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 2, 0, 2, 0\))',r'tmp5',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 0, 0, 2\))',r'tmp6',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 1, 2, 0\))',r'tmp7',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 2, 2, 0\))',r'tmp8',rg)
   
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 2, 0, 0\))',r'tmp9',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 1, 0, 0\))',r'tmp10',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 2, 0, 0, 0\))',r'tmp11',rg)

   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 0, 0, 2\))',r'tmp12',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 0, 2, 0\))',r'tmp13',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 1, 0, 2\))',r'tmp14',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 1, 2, 0\))',r'tmp15',rg)

   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 1, 0, 1\))',r'tmp16',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 2, 0, 1\))',r'tmp17',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 2, 0, 0, 1\))',r'tmp18',rg)

   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 0, 1, 0, 1\))',r'tmp19',rg)
   rg = re.sub(r'(ILU\(l1, m1, l2, m2, l3, m3, 1, 0, 0, 1\))',r'tmp20',rg)







   
   print(rg,"\n")
   print(rg,file=f)
   print("\n",file=f)
  
print((Q23))

