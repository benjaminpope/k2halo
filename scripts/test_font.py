import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['figure.figsize']=(8.0,6.0)    #(6.0,4.0)
mpl.rcParams['font.size']=18               #10 
mpl.rcParams['savefig.dpi']= 200             #72 
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams["font.family"] = "Times New Roman"

x = np.linspace(0,1,100)

y = x**2

plt.plot(x,y,'.')

plt.xlabel('x label')
plt.ylabel('y label')
plt.title('font test')

plt.savefig('font_test.png')