from math import pi, cos, cosh, sin, sinh

def f(y,z,b,c):
    sum = 0.0
    pb2b = pi/(2*b)
    for n in range(0,5):
        tnp1 = 2.0*n+1
        term = (cosh(pb2b*tnp1*z)*cos(tnp1*pb2b*y))/cosh(pb2b*tnp1*c)
        term = (((-1.0)**n)/(tnp1**3))*term
        sum = sum + term
    
    return 0.5*b**2-0.5*y**2-16.0*b*b*sum/(pi**3)

def dfdy(y,z,b,c):
    sum = 0
    pb2b = pi/(2*b)
    for n in range(0,5):
        tnp1 = 2*n+1
        term = (cosh(pb2b*tnp1*z)*sin(tnp1*pb2b*y))/cosh(pb2b*tnp1*c)
        term = (((-1.0)**n)/(tnp1**2))*term
        sum = sum + term
    
    return -y + 8.0*b*sum/(pi**2)

def dfdz(y,z,b,c):
    sum = 0
    pb2b = pi/(2*b)
    for n in range(0,5):
        tnp1 = 2*n+1
        term = (sinh(pb2b*tnp1*z)*cos(tnp1*pb2b*y))/cosh(pb2b*tnp1*c)
        term = (((-1.0)**n)/(tnp1**2))*term
        sum = sum + term
    
    return -8.0*b*sum/(pi**2)


n=40
x1=[]


for i in range(0,n+1):
    x1.append(cos(i*pi/n))

scale = f(x1[n/2],x1[n/2],1.0,1.0);
print scale
x2 = x1
#print x1
thread = ''
w = open("w.txt",'w')
for i in x1:
    for j in x2:
        if abs(i)==1 or abs(j)==1:
            w.write('0.0'+'\n')
#            thread = thread +'\t'+ "0.0"
        else:
            w.write(str(1.0*f(i,j,1.0,1.0)/scale)+'\n')
#            thread = thread +'\t'+ str(f(i,j,1.0,1.0))
#    w.write(thread)
#    w.write('\n')
#    thread = ''
w.close()

#thread = ''
w = open("dwdx.txt",'w')
for i in x1:
    for j in x2:
        w.write(str(1.0*dfdy(i,j,1.0,1.0)/scale)+'\n')
#        thread = thread +'\t'+ str(dfdy(i,j,1.0,1.0))
#    w.write(thread)
#    w.write('\n')
#    thread = ''
w.close()

thread = ''
w = open("dwdy.txt",'w')
for i in x1:
    for j in x2:
        w.write(str(1.0*dfdz(i,j,1.0,1.0)/scale)+'\n')
#        thread = thread +'\t'+ str(dfdz(i,j,1.0,1.0))
#    w.write(thread)
#    w.write('\n')
#    thread = ''
w.close()
