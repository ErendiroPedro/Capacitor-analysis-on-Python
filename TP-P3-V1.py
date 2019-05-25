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
from LabMat2_Lib import *
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
data_p3=loadtxt('data_TC_P3_versao_1.txt')
tt=data_p3[:,0] #discretização do tempo no periodo de definição T.
SFVct=data_p3[:,1] #Aproximação de Fourier de SFVc(t) no periodo de definição T.
print(35*'-','\n\n')
    
#%% Resolução do trabalho

#Constantes:
R= 1.2 #Ohm
L=2 #H
C_=1 #F


print("Trabalho curricular-Lmat2-Parte 3\n\n")

#%%
print('1)\n')
print(u"*Sem output necessário. Verificar código*")

#%%
print('\n\n2)\n')

ti=tt[0]; tf=tt[-1]

for i in range(0,len(data_p3)):
    if(SFVct[i]==min(SFVct)) :
        break;
icum=i     #mais tarde será utilizado
tcum=tt[i]

print('O instante em que o comutador muda de posição é: ', tcum,' (s)' )
print('O inicio do periodo de definição = ', ti,' (s)' )
print('O final do periodo de definição = ', tf,' (s)' )

#%%
print('\n\n3)\n')
    
        
t=sym('t')				# definir t como v. ind.
i=Function('i')  			# definir i como v. dep
C=symbols('C{}'.format(1))  	# definir a constante C1
EDO=Eq(L*Diff(i(t),t),-R*i(t))	# definir a EDO
    
#Antes da comutação
#Antes da comutação é como se o condensador não estivesse no circuito, e assim sendo:
#SFVc(t)=Vl+Vr      
print("-Antes da comutação o circuito é modelado pela seguinte equação diferencial:" )
display(Eq(sym('SFVct'),L*Diff(i(t),t)+R*i(t)))


#Após a comutação
#Após a comutação o condensador já tem efeito no circuito e portanto:
#SFVc(t)=Vl+Vr+Vc  
print("\n\n-Depois da comutação o circuito é modelado pela seguinte equação diferencial:" )
display(Eq(sym('SFVct'),L*Diff(i(t),t)+R*i(t)+(1/C_)*Intc(i(t),t,tcum,tf)))

#%%
print('\n\n4)\n')

print("Condições iniciais para a equação diferencial antes da comutação: " )
print("i(ti)=0")

print("\nCondições iniciais para a equação diferencial depois da comutação: " )
print("i(tcum)=in e q(tcum)=0")
print("\nLegenda:")

print("ti: instante de inicio do periodo de definição ")
print("tf: intante final do periodo de definição")
print("in: último valor da corrente antes do interruptor abrir ")

#%%
print('\n\n5)\n')

    #Antes do instante de comutação
ii=[]
ii.append(0) #Condição inicial
h=tt[1]-tt[0] #Diferença entre duas amostras consecutivas

for k in range(1,icum):
    ii.append(ii[k-1]+h*(1/L)*(SFVct[k-1]-R*ii[k-1])) #método de Euler

#gráfico
plot(tt[0:len(ii)],ii)
title(u'i(t) antes do instante de comutação')
grid('on');show()

    #Depois do instante de comutação

#De modo a facilitar o calculo pretendido e evitar o calculo de um integral numericamente 
#fizemos a seguinte substituição: (dq/dt)=i

qq=[]
qq.append(0)    #Condição inicial
iid=[]
iid.append(ii[-1]) #Condição inicial
ttd=[]
ttd.append(tcum)

j=0;

i1,q=sym('i1 q')

#resolução
for z in range(icum+1,len(tt)):
    
    j+=1
    #método de Euler  
    f=inline(i1,i1)
    g=inline((1/L)*(SFVct[z]-R*i1-(1/C_)*q),(i1,q))    
    
    iid.append(iid[j-1]+h*g(iid[j-1],qq[j-1])) #1 equação do sistema de EDOs
    qq.append(qq[j-1]+h*f(iid[j-1]))    #2 equação do sistema de EDOs

#Deste sistema obtemos duas listas, uma com os valores da carga ao longo do tempo e outra com os valores da corrente

plot(tt[icum:len(tt)],iid)
title(u'i(t) depois do instante de comutação')
grid('on');show()

plot(tt[icum:len(tt)],qq)
title(u'q(t) depois do instante de comutação')
grid('on');show()



#%%
print('\n\n6)\n')
    

I=ii+iid #Concatenação dos vetores de corrente 

Vr=[]   #Vetor com as tensões na resistência
Vl=[]   #Vetor com as tensões na bobina
Vc=[]   #Vetor com as tensões no condensador
di_dt=[]    #Vetor a derivada da corrente

#Tensão na resistência
k=0
while (k<len(tt)):
    Vr.append(R*(I[k])) 
    k+=1      

#Tensão no condensador
k=0
for k in range(0,len(ii)):
    Vc.append(0)    #Inicialmente a tensão no condensador é nula

p=0
for k in range(icum,len(tt)):
#Segundo leis fisicas sabemos que a tensão no condesador é dada por: Vc=q*C logo
    Vc.append(qq[p]*C_) #Tensão no condensador após a comutação
    p+=1
    
#Tensão na bobina
k=0

#Aplicando a lei de Kirchoff para tensões em malhas fechadas, sabemos que a soma das 
#quedas de tensão nos componentes do circuito tem que ser igual a tensão debitada pela fonte, 
#Assim sendo, após a comutação teremos que Vc=SFVct-Vr-Vl

for k in range(0,len(tt)):
    Vl.append(SFVct[k]-Vr[k]-Vc[k])
    
print(u"*Sem output necessário. Verificar código*")

#%%
print('\n\n7)\n')


#Representação gráfica da corrente durante o periodo de definição
plot(tt[0:len(tt)],I)
title(u'I(t) durante o periodo de definição do sinal')
grid('on');show()

#Tensão na resistência duarnte o periodo de definição
plot(tt[0:len(tt)],Vr)
title(u'Vr(t) durante o periodo de definição do sinal')
grid('on');show()

#Tensão na bobina duarnte o periodo de definição
plot(tt[0:len(tt)],Vl)
title(u'Vl(t) durante o periodo de definição do sinal')
grid('on');show()

#Tensão no condensador duarnte o periodo de definição
plot(tt[0:len(tt)],Vc)
title(u'Vc(t) durante o periodo de definição do sinal')
grid('on');show()

#%%
print('\n\n\n\t\t',32*'-')
print("\t\t\tFim do trabalho!")