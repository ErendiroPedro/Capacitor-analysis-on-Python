# -*- coding: utf-8 -*-
"""
Laboratórios de Matemática 2
Parte 1 - Versão 1

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

 #Condensador
# C=??? (um objectivo 'conhecer seu valor')


# Em t=0 (s), K está na posição 1
# Ao longo do tempo foram registados os valores de V_in e V_c  (volt)

 #Erros de leitura
T_r=0.001 #(s)
V_r=0.005 #(v)

#%% Pergunta 1

print('Pergunta 1\n\n')

# nº de amostras
ndados=len(dados)
print(u'- Foram feitas ', ndados,'amostras')


# intervalo de tempo entre amostras
T=100/ndados
print (u'- Intervalo de tempo entre amotras = ',format(T,'.3f'),' +- ',T_r,' (s) ')


# frequncia de amostragens
f=1/T  #valor da frequência

t_=sym('t')
f_=1/t
f_r=diff(f_,t_)*T_r  #erro associado ao calculo da frequência
f_r=abs(subs(f_r,t_,T))

print (u'- Frequência de amostragem = ', format(f,'.7f'),' +- ', format(f_r,'.7f'), '(Hz) ')


# tempo durante ocorreu aquisição de  V(in) e V(c)
tempaq=dados[:,0][-1]-dados[:,0][0]
print (u'- Tempo de aquisição = ', tempaq,' +- ', T_r,' (s)')

print('\n\t -----//-----')



#%% Pergunta 2

print('\n\nPergunta 2\n\n')

t=dados[:,0]
Vin=dados[:,1]
Vc=dados[:,2]

plot(t,Vin,t,Vc)
xlabel('t (s)'); ylabel('V (V)')
title(u'Deteção de ocorrências')  
show()

print('\n\t -----//-----')

#%% Pergunta 3

print('\n\nPergunta 3\n\n')

l=[]#Lista de derivadas
l1=[]#lista de indices onde a derivada é negativa
l2=[]
tempdes=[]
#definicao de uma lista de derivadas
for i in range(0,len(dados)):
    if(i<len(dados)-1):
        l.append((Vin[i+1]-Vin[i])/T)
    else:
        l.append((Vin[i]-Vin[i-1])/T)

#Pesquisar na lista de derivadas, as derivadas negativas       
#Lista de indices e instantes onde a derivade é negativa



for i in range(0,len(l)): 
    if (l[i]<-1):  #neste caso temos <-1 para desprezar as derivadas de muito pequeno valor que representam erros de leitura.
        l1.append(i)
        tempdes=tempdes+[t[i]]
    if(l[i]>1):
        l2.append(i)
        tempdes=tempdes+[t[i]]
        
tempdes=tempdes+[t[-1]]  



print(u'As descargas do condesador ocorrem nos seguintes intervalos: \n')
for i in range(0,len(tempdes)-1):
    if ((i%2)==0):
        print('[',tempdes[i],',', tempdes[i+1],']')
print(u'\nTodas com um erro = ', '+- ', T_r, ' (s)')


print('\n\t -----//-----')


#%% Pergunta 4

print('\n\nPergunta 4\n\n')

#%% A)
Vdc=Vc[l1[0]] #Vc no momento da descarga
Vtau=0.369*Vdc #Vc para calcular o tau

for i in range(l1[0],l2[0]):    #ciclo para descobrir o indice para calcular tau
    if(Vc[i]<Vtau):
        itau=i
        break

tau=t[i]-t[l1[0]]   #calculo do tau

print(u'O Valor de tau e: \n')
print(tau, '(s) +-', T_r)

print('\n\t -----//-----')

#%% B)

#Equacao para calcular C
C1=tau/R
C=FIX(C1,5)
#derivada em ordem a R
dCr=-(tau/(R*R))
#derivada em ordem a tau
dCtau=1/R

merro=abs(dCtau*T_r)+abs(dCr*erro_R)
merro=2*10**(-5)

print(u'O Valor de C e: \n')
print(C, 'f  +-', merro)

#%% C)

#Definicao dos arrays
dQt=[]
Qt=[]

#definir de Q(t)
for i in range(0,len(dados)):
    Qt.append(C1*Vc[i])

#Definir um array com as derivadas de Vc(t)
dQt.append((Qt[0+1]-Qt[0])/T)#Primeira derivada

for i in range(1,len(dados)-1): #defenicao das derivadas, ou TVMC
    if(i>0):
        if(i<len(dados)):
            dQt.append((Qt[i+1]-Qt[i-1])/2*T)
            
dQt.append((Qt[i]-Qt[i-1])/T)

plot(t,Qt)
axis([10,25,0,0.0020])
xlabel('t (s)'); ylabel('Q (t)')
title(u'Grafico de Q(t)')  
show()

plot(t,dQt)
axis([10,25,-0.00005,0.00005])
xlabel('t (s)'); ylabel('I (t)')
title(u'Grafico de i(t)')  
grid(0.5)
show()

#%% D)

# Valor de inicio da descarga.
z1=sym('z1');z1=0;
z2=sym('z2');z2=0;

for i in range(0,len(dados)):
    if(dados[i,0]==10.000):
        z1=i
for i in range(0,len(dados)):
    if(dados[i,0]==25.000):
        z2=i
# Como a função está em ordem ao tempo, a largura dos retangulos será o intervalo de tempo entre amostras.
# E as alturas dos retangulos serão os valores de i**2*R, logo a area sera (i**2*R)*T 
#Defínir as variaveis   
z=[]
zM=sym('zM');zM=0;
zm=sym('zm');zm=0;
#obter os valores da função i**2*R
for i in range(z1,z2):
    z.append((dQt[i]**2)*R)
# Calculo da soma de Riemann a direita (Majorante)
tz=len(z)
for i in range(1,tz):
    zM=(z[i]*T)+zM

print("O valor do Majorante é: ",zM)
# Calculo da soma de Riemann a esquerda (Minorante)
for i in range(0,tz-1):
    zm=(z[i]*T)+zm
print("O valor do minorante é: ",zm)
# Calculo do Erro
erroz=abs(zM-zm)/2
print("O valor do erro é: ",erroz)
# Calculo do Valor Medio
areaz=(zM+zm)/2
print("O valor médio é: ",areaz)