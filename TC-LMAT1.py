# -*- coding: utf-8 -*- 


 
from __future__ import division
from __future__ import print_function
from LabMat1_Lib import *
import pylab as pl

clc()

#%% expressões gerais
a,c,E,R,C,dv,dt,V,tau=sym('a,c,varepsilon,R,C,dV,dt,V,tau')
b=symbols('b', positive=True)
dvt=dv/dt                                         
display(sym('Express\~{o}es\;Gerais\:'))
f1=inline(dvt*R*C+V,t)                            
v1=inline(a*e**(-t/b)+c,t)                                        
display(sym('fun{c}\c\~ao\;(1)\:'))
display(Eq(f1(t),0))                              
display(sym('fun{c}\c\~ao\;(2)\:'))
display(Eq(sym('V(t)'),v1(t)))                    

#%% alínea 1.)
display(sym('Subdivis\~ao\;1.)'))
v0t=Eq(E-v1(0),0)                       
display(Eq(sym('V(0)'),E))                                    
display(Eq(sym('V(0)'),v1(0)))                    
display(Eq(v1(0),E))                               
display(v0t)

#%% alínea 2.)
display(sym('Subdivis\~ao\;2.)'))
edvltinf=Limit(v1(t),t,inf)                       
fvltinf=limit(v1(t),t,inf)                        
display(Eq(edvltinf,0))                           
fvltinfe=Eq(fvltinf,0)
display(fvltinfe)                                 
                                     
#%% alínea 3.)
display(sym('Subdivis\~ao\;3.)'))
display(sym('Derivada\;de\;V(t)\;em\;fun{c}\c\~ao\;a\;t\:'))
fderiva=diff(v1(t),t)                             
display(Eq(dv/dt,fderiva))                        
display(sym('Sistema\;de\;equa{c}\c\~oes\:'))
fure=R*C*fderiva+v1(t)                              
display(Eq(fure,0))                               
display(fvltinfe)                                 
display(v0t)                                      
S=solve([v0t,fvltinf,fure],[a,b,c])              
display(sym('Solu{c}\c\~ao\;do\;Sistema\:'))
display([a,b,c],S)                                

#%% alínea 4.)
display(sym('Subdivis\~ao\;4.)'))
v2=inline(subs(v1(t),a,S[0][0]),t)                
v3=inline(subs(v2(t),b,S[0][1]),t)               
v4=inline(subs(v3(t),c,S[0][2]),t)               
display(Eq(sym('V_D'),v4(t)))

#%% alínea 5.)
display(sym('Subdivis\~ao\;5.)'))

#%% a.)
display(sym('Declive\;da\;tangente\;\`a\;curva\;de\;carga\;na\;origem\:'))
dervd=diff(v4(t),t)           
declv=subs(dervd,t,0)                          
display(declv)
ftan=inline(declv*t+E,t)                      
display(sym('Equa{c}\c\~ao\;da\;recta\:'))
vtau=Eq(ftan(t),0)
display(vtau)
Stau=solve(vtau,t)
display(sym('Valor\;de\;\mathcal{T}\:'),Stau)

#%% b.) regra 3 simples nao esta a dar o valor em e
descrg=inline(subs(v4(t),t,Stau[0]),t)
display(Eq(sym('V_(RC)'),descrg(t)))
dsct=(E-descrg(t))
value=dsct/E
prc=value*100
display(prc)

#%% alínea 6.)
display(sym('Subdivis\~ao\;6.)'))

#%% a.)
FEM=50; RL=400; CL=0.010
vdesc=inline(FEM*exp(-t/(RL*CL)),t) 
g6a=ezplot(vdesc(t),[0,16])
axis([0,16,-2,55])
grid(1)
xlabel('tempo (s)')
ylabel(u'Tensão (V)')
setp(g6a,c='b',lw=3,ls='-')
title(u'Gráfico')

#%% b.)
dvdesc=diff(vdesc(t),t)     
m=subs(dvdesc,t,0)        
display(Eq(sym('Declive\;da\;recta\;tangente'),m))
y=inline(m*t+50,t)
g6b=ezplot(y(t),[0,16])
axis([0,16,-2,55])
grid(1)                                   
setp(g6b,c='c',lw=2,ls='--')            

#%% c.)
ah=limit(vdesc(t),t,inf)
g6c=ezplot(ah,[0,16])
display(Eq(Limit(vdesc(t),t,inf),ah))
setp(g6c,c='m',lw=2,ls='--')

#%% d.) assintota → y=0, para isso temos equação da tangente = 0
eqpx=Eq(m*t+50,0)
display(eqpx)
xint=solve(eqpx,t)
display(sym('Interse{c}\c\~ao\;da\;tangente\;e\;assintota\:'))
display([xint,0])

#%% e.) 
eq2=FEM*exp(-xint[0]/(RL*CL))          
display(sym('Interse{c}\c\~ao\;de\;A\;e\;B\:'))     
display([xint[0],eq2])                
g6d=ezplot(eq2,[0,4])                                  
setp(g6d,c='r',lw=2,ls='-')
g6e=ezplot(4,t,[0,18.4]) 
setp(g6e,c='r',lw=2,ls='-')
legend('ABCD')                   
show()                           
