import numpy as np
import matplotlib.pyplot as mpl
from sympy import *
from matplotlib import interactive
import matplotlib.patheffects as pe

global P
P = 106*0.593  # Persistence length of DNA in nm*kcal/(mol*rad^2)
global C
C = 1.8*1.4393 # Persistence length of melted strands in kcal/(mol*rad^2)

global twists

twists = [0.0,-20.0,-11.2,19.0,-2.0,-11.6,-15.0,0.0]

sl=['c','g','c','a','t','g','a','a','c','t','g','c','a','g','t','t','a','t','a','t','g','g','a','c','c','t','c','g','a','t','g','c','g','g','c','g','t','a','c','a','g','t','a','c','g','c']
sl =np.asarray(sl)
global Fs

x=Symbol('x')
y=Symbol('y')
r=Symbol('r')
beta=Symbol('beta')

ks=C*np.pi**2/180**2
k_con=0.0244
Fs = 3.76

def string_split(seq):
    return [char.lower() for char in seq]




def make_energy(seq):
    e, e_bp = [], []
    sl = string_split(seq)
    ll=0
    while ll<len(sl):
          if sl[ll]=='c' or sl[ll]=='g':
             e_bp.append(0.12)
          else:
             e_bp.append(0.64)
          if ll<len(sl)-1:
             if (sl[ll]=='a' and sl[ll+1]=='g') or (sl[ll]=='c' and sl[ll+1]=='t'):
                e.append(-1.44)
             elif (sl[ll]=='a' and sl[ll+1]=='c') or (sl[ll]=='g' and sl[ll+1]=='t'):
                e.append(-2.19)
             elif (sl[ll]=='a' and sl[ll+1]=='t'):
                e.append(-1.72)
             elif (sl[ll]=='a' and sl[ll+1]=='a') or (sl[ll]=='t' and sl[ll+1]=='t'):
                e.append(-1.49)
             elif sl[ll]=='t' and sl[ll+1]=='a':
                e.append(-0.57)
             elif (sl[ll]=='t' and sl[ll+1]=='c') or (sl[ll]=='g' and sl[ll+1]=='a'):
                e.append(-1.81)
             elif (sl[ll]=='t' and sl[ll+1]=='g') or (sl[ll]=='c' and sl[ll+1]=='a'):
                e.append(-0.93)
             elif (sl[ll]=='g' and sl[ll+1]=='g') or (sl[ll]=='c' and sl[ll+1]=='c'):
                e.append(-1.82)
             elif sl[ll]=='g' and sl[ll+1]=='c':
                e.append(-2.55)
             elif sl[ll]=='c' and sl[ll+1]=='g':
               e.append(-1.29)
          ll+=1
    return e, e_bp

def duplex_energy(nbp, tw_real):
    P_opt=100.0
    f_dna = 0.593*P_opt/(0.334*nbp)*np.pi**2/180.0**2
    E_dupl= 0.5*f_dna*tw_real**2
    return E_dupl

global rg_factor

rg_factor=4.1

global f_bubble

f_bubble=C*np.pi**2/180**2    


global lengx

lengx=[60,10,10,10,10,10,10,60]


class dna:
      def __init__(self,seq):
          self.e, self.e_bp = make_energy(seq)
          
      def E_bubble(self,n,pos,tw): #Calculate energy of bubble, n=bubblesize,pos=position index of bubble-initiation (base-pair)
          l=1
          E_stacking = 0
          E_st       = 0
          while l<=n:
                E_stacking = E_stacking - self.e[pos+l-1] - self.e_bp[pos+l-1]
                l+=1
          E_st       = rg_factor + E_stacking - self.e[pos-1]  # Stacking energy depends on stacking bef&aft flipped base + base pairing + ring contribution
          E_dangl    = 0.5*f_bubble*(tw-34.5*(n+1))**2/n             #Dangling/Twisting of bubble
          E_bub = E_st + E_dangl
          return E_bub




      def loop_exclude(self,idx, seq_idx):
          n=1
          pos=0
          Z_bub = []
          Z_dead= []
          while pos<len(self.e_bp):
                n=1
                while n<=15 and (pos+n-1)<len(self.e_bp)-1:
                      E_bub = self.E_bubble(n,pos, twists[seq_idx])
                      Z_bub.append(np.exp(-E_bub/0.593))
                      if pos<=idx and idx <= (pos+n):
                         Z_dead.append(np.exp(-E_bub/0.593))
                      n+=1
                pos+=1
          E_duplex=duplex_energy(lengx[seq_idx],twists[seq_idx])
          Z_reg = np.exp(-E_duplex/0.593)
          p = np.sum(Z_dead)/(Z_reg+np.sum(Z_bub))
          return p





def melting_prob(seq):
    interactive(True)
    outer_ii=1
    indices=[]
    t=[]
    probs=[]
    ii=1
    while ii<=len(seq)-1:
              if ii<=60:
                 lp=0
              elif ii>60 and ii<=70:
                 lp=1
              elif ii>70 and ii<=80:
                 lp=2
              elif ii>80 and ii<=90:
                 lp=3
              elif ii>90 and ii<=100:
                 lp=4
              elif ii>100 and ii<=110:
                 lp=5
              elif ii>110 and ii<=120:
                 lp=6
              else:
                 lp=7
              x=dna(seq)
              indices.append(outer_ii)
              p = x.loop_exclude(ii,lp)
              probs.append(p)
              ii+=1
              outer_ii = outer_ii+1
    mpl.style.use('classic')
    mpl.rc('text', usetex=True)
    mpl.rc('font', family='Arial')
    mpl.rcParams['xtick.major.size']=10
    mpl.rcParams['ytick.major.size']=10
    mpl.rc('xtick.major', pad=10)
    mpl.rc('ytick.major', pad=10)
    mpl.xticks(fontsize=20)
    mpl.yticks(fontsize=20)
    xt = np.linspace(0,180,181)
    xtick_labels=seq
    options = {'A':'red','T':'red','G':'blue','C':'blue'}
    ax=mpl.gca()
    ax.set_xticks(xt)
    ax.set_xticklabels(xtick_labels,fontsize=7)
    colo= [options[i] for i in xtick_labels]
    for xtick,colors in zip(ax.get_xticklabels(),colo):
            xtick.set_color(colors)
    #options={ "gc":"cyan","at":"lightgreen"}
    mpl.plot(indices,probs,'-',color='black',lw=2,)
    mpl.xlabel(r'Base-pair position ', fontsize=40)
    mpl.ylabel("Melting Probability", fontsize=40)
    mpl.legend(prop={'size':35},loc=1)
    return probs


