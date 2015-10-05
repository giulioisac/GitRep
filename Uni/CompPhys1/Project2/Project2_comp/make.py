import commands,os
import matplotlib.pyplot as plt
import numpy as np

def getarray(a):
	"""

	function analogous to the one in c++
	a=name of the program
	
	"""
	
	os.system("".join(["c++ -Wall -c ",a,".cpp"]))
	os.system("".join(["c++ -o", " " , a , " " , a , ".o -larmadillo -llapack -lblas"]))
	out = commands.getoutput("".join(["./",a]))
	
	b = out.split()
	i=0
	array= dict()
	j=0
	while (b[i]!="break"):
		array[j] = []
		while (b[i]!='change'):
			array[j].append(float(b[i]))
			i+=1
		i+=1
		j+=1
	for x in range(j):
		yield array[x]

ro, psi1, psi2, psi3,ro2, psi11, psi22, psi33 = getarray("main")
# a=0
# b=0
# c=0
# for i in range(0,len(ro)):
# 	a+=psi1[i]*psi1[i]
# 	b+=psi2[i]*psi2[i]
# 	c+=psi3[i]*psi3[i]
# print a*4*np.pi
# print b*4*np.pi
# print c*4*np.pi
dr=2.5/200
psi1=psi1/sum(np.absolute(psi1)*dr)
psi2=psi2/sum(np.absolute(psi2)*dr)
psi3=psi3/sum(np.absolute(psi3)*dr)
psi11=psi11/sum(np.absolute(psi11)*dr)
psi22=psi22/sum(np.absolute(psi22)*dr)
psi33=psi33/sum(np.absolute(psi33)*dr)

plt.plot(ro, psi1, '-b',
		ro, psi2,'-r',
		ro, psi3,'-g',
		ro, psi11, '--b',
		ro, psi22,'--r',
		ro, psi33,'--g')

plt.legend(['E1 INT', 'E2 INT',"E3 INT",'E1 NO-INT', 'E2 NO-INT',"E3 NO-INT"], loc='upper right')
plt.grid()
plt.xlabel("ro")
plt.ylabel("psi^2")
plt.title("Probability")
plt.savefig('OMEGA_001_ARMA.png')
plt.show()
