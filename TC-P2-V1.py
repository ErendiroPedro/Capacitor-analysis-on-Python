# -*- coding: utf-8 -*-
"""
Laboratórios de Matemática 2
Parte 2 - Versão 1.2.teste

Data: 2017/04/08

Autores
1. <1160555> <Erendiro Sangueve Njunjuvili Pedro>
2. <1160775> <João Pedro Freitas Neves>
3. <1161053> <Tiago Miguel Pereira Ribeiro>
4. <1161312> <João Pedro Gomes da Costa>
"""
from __future__ import division
from __future__ import print_function
from LabMat1_Lib import *
import pylab as pl


clc()

#%% Bloco de identificação do trabalho
print('\n\nTrabalho curricular\nPARTE 1\nVERSÃO 1\n\n')
print('AUTORES:')
print('\t<1160555>\t<Erendiro Pedro>')
print('\t<1161053>\t<Tiago Ribeiro>')
print('\t<1160775>\t<João Neves>')
print('\t<1161312>\t<João Costa>')

#%% Leitura do ficheiro de dados
# Estrutura de dados: dados=array(Nlinhas,3colunas)
#  t/s |  Vin/V  | Vc/V
#
Ficheiro="Dados_TC_P1_V1.txt"
print('\n\nLeitura do ficheiro de dados:',Ficheiro,'...\n\n')
dados=loadtxt(Ficheiro)
print('      t/s      Vin/V     Vc/V ')
print(35*'-')
for i in range(0,5):
    print((3*"%10.3f") % tuple(dados[i,:]))
print('     .....     .....     .....')
for i in range(-5,0):
    print((3*"%10.3f") % tuple(dados[i,:]))
print(35*'-','\n\n')
    
#%% Resolução do trabalho

# Dados

 #Resistencia
R=4700#'-+ 5% = T_R
erro_R=4700*0.05

C=330*10**-6 # 330 microfardas

# Em t=0 (s), K está na posição 1
# Ao longo do tempo foram registados os valores de V_in e V_c  (volt)

 #Erros de leitura
T_r=0.001 #(s)
V_r=0.005 #(v)


#%% Pergunta 1

print("Pergunta 1\n\n")
print("a)\n")
    
tt=dados[:,0]
Vin=dados[:,1]
Vc=dados[:,2]


td=10     #inicio da primeira descarga
tc=25  #inicio da primeira carga

for i in range(0,len(dados)):
    if (tt[i]==td): break; #valor obtido atravez da primeira parte do trabalho
for j in range(0,len(dados)):
    if (tt[j]==tc): break; #valor obtido atravez da primeira parte do trabalho

Vd=Vc[i]    #tencao no momento da descarga
Vcarg=Vc[j]

print("Valores extremos da tensão:")
print("Tensão no momento de descarga: ", Vd)
print("Tensão no momento de carga: ", Vcarg)

print("\nE como tal, os valores para VM e V0, serão respectivamente: ")
VM=round(Vd); V0=round(Vcarg)
print(VM,",",V0," (V)")

t,T,tau,n=sym('t T tau n')
n,N=symbols('n,N',positive=True,integer=True)
m=inline(t-T*floor((t-tau)/T),(t,tau,T))  # definir função modulus

tau1=round(R*C,1)  #(tau=RC=~1.551=~1.5) aproximacao de tau para nao dar erro nos harmonicos
T=25    #periodo da primeira descarga completa
tau=10       #inteiros mais proximos dos valores resgistrados

for w in range(0,len(tt)): #Procurar posição do em que o tempo assum valor de 10
    if(tt[w]==td):
        break

LVc=VM*(1-e**(-(t-tc)/tau1)); #Defenicao de Vc no momento de carga
LVd=Vc[w]*e**(-(t-td)/tau1);     #Defenicao de Vc no momento de descarga

VC=inline(Piecewise((LVc,t>25),(LVd,t<=25)),t) # definir os 2 ramos da funcao Vc        
Vcp=inline(VC(m(t,tau,T)),t)    #defenicao da funcao periodicamente

          
Vi=inline(Piecewise((VM,t>25),(V0,t<=25)),t) # definir os 2 ramos da Funcao Vin
Vip=inline(Vi(m(t,tau,T)),t)    #defenicao da funcao periodicamente

print("\n\nRepresentação gráfica dos sinais Vin e Vc:")
          
h1=ezplot(Vcp(t),[0,100])
h2=ezplot(Vip(t),[0,100])
title("Vin e Vc");xlabel("t (s)");ylabel("Vin / Vc (v)")
setp(h2,color='b')
setp(h1, color='r')
axis([0,100,0,8])
grid(1)
show()

print("\t\t--------//--------\n\n\n")
#%% b)
print("b)\n")
#Calculo da potencia media
print("Aplicando a formula que nos é apresentada no formulário, segue-se:\n ")
PVin=intc(Vi(t)**2,t,tau,(tau+T))/T
print('PVin=', FIX(PVin,2),'(W)\n')

PVc=intc(VC(t)**2,t,tau,(tau+T))/T
print('PVc=', FIX(PVc,2),'(W)')
print("\n\n\t\t--------//--------\n\n\n")
#%% c) 
print("c);d)\n")
#series na forma tringonometricas dos sinais

W0=2*pi/T #frequencia angular do sinal 

#para vin

a0=2*intc(Vi(t),t,tau,(tau+T))/T
an=2*intc(Vi(t)*cos(n*W0*t),t,tau,(tau+T))/T
bn=2*intc(Vi(t)*sin(n*W0*t),t,tau,(tau+T))/T

#Harmonicos de ordem n
Hvin=inline(an*cos(n*W0*t)+bn*sin(n*W0*t),(n,t))

#para vin

A0=2*intc(VC(t),t,tau,(tau+T))/T
An=2*intc(VC(t)*cos(n*W0*t),t ,tau,(tau+T))/T
Bn=2*intc(VC(t)*sin(n*W0*t),t,tau,(tau+T))/T

#Harmonicos de ordem n
Hvc=inline(An*cos(n*W0*t)+Bn*sin(n*W0*t), (n,t))

#%%d)

#Grafico das aproximacoes

SFVin=a0/2+symsum(Hvin(n,t),n,1,10)
h1=ezplot(SFVin,[0,100])
h2=ezplot(Vip(t),[0,100])
setp(h2,color='b')
setp(h1, color='r')
xlabel('t (s)'); ylabel('Vin (v)')
title(u'Vin + Aproximação pela Serie de fourier')
grid('on')
show()

SFVc=a0/2+symsum(Hvc(n,t),n,1,10)
h1=ezplot(SFVc,[0,100])
h2=ezplot(Vcp(t),[0,100])
setp(h2,color='b')
setp(h1, color='r')
axis([0,100,0,8])
xlabel('t (s)'); ylabel('Vc (v)')
title(u'Vc + Aproximação pela Serie de fourier')
grid('on')
show()

print("\n\n\t\t--------//--------\n\n\n")
#%% e)
#Espectros de Vin
print("e)\n")
Ani=(an**2+bn**2)
AVi=sqrt(Ani) # calculo das amplitudes de Vi
phin=atan2(-bn,an)  #calculo das fases de Vi

Pn=(a0/2)**2+(1/2)*symsum(Ani,n,1,26) #Valor de potencia media Pn
RMSE=100*sqrt((PVin-Pn)/PVin)       #Calculo da percentagem do erro
print('Valor de Pn', FIX(Pn,2),'Watts')
print('Valor de RMSE:', FIX(RMSE,2),'%')



W=[i*W0 for i in range(27)] # calculo da frequencia angular
HVi1=[abs(a0/2)]+[subs(AVi,n,i) for i in range(1,27)]
PHI=[subs(phin,n,i) for i in range(1,27)]
subplot(1,2,1);stem(W,HVi1);axis([0,9,0,3.5]);grid(1)
title(u'Amplitude de Vi')
subplot(1,2,2);stem(W[1:],PHI);axis([0,9,-4,4]);grid(1)
title(u'Fases de Vi')
show()



#%% Espectros de Vc

Anc=An**2+Bn**2
AVc=sqrt(Anc) # calculo das amplitudes de Vc
phin1=atan2(-Bn,An) #calculo das fases de Vc

Pnvc=(A0/2)**2+(1/2)*symsum(Anc,n,1,6) #Valor de potencia media Pn
RMSE1=100*sqrt((PVc-Pnvc)/PVc)         #Calculo da percentgem do erro
print('Valor de Pn', FIX(Pnvc,2),'Watts') 
print('Valor de RMSE:', FIX(RMSE1,2),'%')


W=[i*W0 for i in range(6)] # calculo da frequencia angular
HVc1=[abs(A0/2)]+[subs(AVc,n,i) for i in range(1,6)]
PHI1=[subs(phin1,n,i) for i in range(1,6)]
subplot(1,2,1);stem(W,HVc1);axis([0,2,0,3.5]);grid(1)
title(u'Amplitude de Vc')
subplot(1,2,2);stem(W[1:],PHI1);axis([0,2,-4,4]);grid(1)
title(u'Fases de Vc')
show()


print("\n\n\t\t------------------//------------------\n\n\n")

#%% Resolucao numerica da Serie de fourier para Vc

#Descobrir os indices correrpondentes ao tempo a analizar
taui=10
tf=35

for i in range(0,len(dados)):
    if (tt[i]==taui): break; #Valor de indice do tempo de inicio
for j in range(0,len(dados)):
    if (tt[j]==tf): break; #valor do tempo final
    
    
dt=round(tt[len(dados)-1]/len(dados),4) #calculo do intrevalo de tempo dt

Vcdt=[Vc[k]*dt for k in range(i,j)] #calculo numerico do integral de Vc para o calculo de a0
Vcdta=[Vc[k]*cos(n*W0*tt[k])*dt for k in range(i,j)] # integral para an
Vcdtb=[Vc[k]*sin(n*W0*tt[k])*dt for k in range(i,j)] #integral para bn

Na0=(2/T)*sum(Vcdt) #defenicao de a0

Nan=(2/T)*sum(Vcdta) #defenicao de an
    
Nbn=(2/T)*sum(Vcdtb) #defenicao de bn
    
NHN=inline(Nan*cos(n*W0*t)+Nbn*sin(n*W0*t),(n,t)) # defenicao da serie de harmonicos

#%% Grafico de arpoximacao numeirca da serie de fourier

SFNVc=Na0/2+symsum(NHN(n,t),n,1,10) #serie de fourier
h1=ezplot(SFNVc,[0,100])
h2=ezplot(Vcp(t),[0,100])
setp(h2,color='b')
setp(h1, color='r')

xlabel('t (s)'); ylabel('Vc (v)')
title(u'Vc + Aproximação pela Serie de fourier')
grid('on')
show()

#%% 

NAnc=(Nan**2+Nbn**2)
NAVc=sqrt(NAnc)
Nphin1=atan2(-Nbn,Nan)
    
    
NPnvc=(Na0/2)**2+1/2*symsum(NAnc,n,1,6)
RMSE2=100*sqrt((PVc-NPnvc)/PVc)         #Calculo da percentgem do erro
print('Valor de Pn', FIX(Pnvc,2),'Watts') 
print('Valor de RMSE:', FIX(RMSE1,2),'%')     
    
W=[i*W0 for i in range(6)] # calculo da frequencia angular
NHVc1=[abs(Na0/2)]+[subs(NAVc,n,i) for i in range(1,6)]
NPHI1=[subs(Nphin1,n,i) for i in range(1,6)]
subplot(1,2,1);stem(W,NHVc1);axis([0,2,0,3.5]);grid(1)
title(u'Amplitude de Vc')
subplot(1,2,2);stem(W[1:],NPHI1);axis([0,2,-4,4]);grid(1)
title(u'Fases de Vc')
show()
    
    
     
print("\n\n\t\t------------------//------------------\n\n\n")    
    
#%% 3.
print("/n  3.")  
    
# estudo de Vin

# THD = sqrt(h2**2+h3**2+...+h10**2)/h1
# SFVin=a0/2+symsum(Hvin(n,t),n,1,10)
# SFVc=a0/2+symsum(Hvc(n,t),n,1,10)
x_3_vin=0
for i in range(1,10):
    x_3_vin=x_3_vin+(a0/2+symsum(Hvin(n,t),n,i,i))**2
     
THD_vin=sqrt(x_3_vin)/(a0/2+symsum(Hvin(n,t),n,i,i))
print(THD_vin)
print("/noi/n")
#disp(simplify(THD_vin))
   
    
    
    
    