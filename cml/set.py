import matplotlib.pyplot as plt
import CML
import maps
import sys
import numpy as np
import math
from scipy.sparse.linalg import eigs, ArpackError

def plot(mat):
	plt.clf()
	plt.imshow(mat,cmap=plt.get_cmap('gray'))
	plt.colorbar()


def setMap(mapping):
	if(mapping=='logistic'):
		cmlMap = maps.logisticMap
		par = [4.0]
	elif(mapping=='doubling'):
		cmlMap = maps.doublingMap
		par = []
	elif(mapping=='SOM'):
		cmlMap = maps.somMap
		par = [0.6,0.2]
	else:
		raise Exception('Unsupported map',mapping)
	return cmlMap, par

if __name__ == "__main__":
	if('-h' in sys.argv) or ('--help' in sys.argv):
		print('================================')
		print('Syntax:')
		print('python main.py commands')
		print('================================')
		print('Commands:')
		print('================================')
		print('-d: display the last iteration of CML')
		print('-mat <type> <matlen>: initial condition type(random, gaussian or bessel), matlen(>0)')
		print('-alpha <min> <max> <nsteps>: logistic alpha')
		print('-nit <nit>: number of iterations')
		print('-c <min> <max> <nsteps>: coupling factor(float)')
		print('-csv: if set, saves each iteration csv files')
		print('-h or --help: display help')
		print('================================')
		exit()
	
	initialMat,matLen,nit,mapping,neigh,output,nlyap = 'gaussian',128,50,'logistic',[(1,0),(0,1),(-1,0),(0,-1)],False,100
	amin,amax,adelta = 3.0, 4.01, 0.025
	cmin,cmax,cdelta = 0.05, 1.01, 0.025
	#amin,amax,adelta = 3.0, 4.01, 0.5
	#cmin,cmax,cdelta = 0.01, 1.01, 0.5
	for it in range(len(sys.argv)):
		if(sys.argv[it] == '-mat'):
			if(it+2 < len(sys.argv)):
				initialMat = sys.argv[it+1]
				matLen = int(sys.argv[it+2])
			else:
				raise Exception('Wrong Syntax','-mat')
		elif(sys.argv[it] == '-alpha'):
			if(it+3 < len(sys.argv)):
				amin = sys.argv[it+1]
				amax = sys.argv[it+2]
				adelta = sys.argv[it+3]
			else:
				raise Exception('Wrong Syntax','-map')
		elif(sys.argv[it] == '-nit'):
			if(it+1 < len(sys.argv)):
				nit=int(sys.argv[it+1])
			else:
				raise Exception('Wrong Syntax','-nit')
		elif(sys.argv[it] == '-c'):
			if(it+3 < len(sys.argv)):
				cmin=float(sys.argv[it+1])
				cmax=float(sys.argv[it+2])
				cdelta=float(sys.argv[it+3])
			else:
				raise Exception('Wrong Syntax','-c')
		elif(sys.argv[it] == '-lyap'):
			if(it+1 < len(sys.argv)):
				nlyap=float(sys.argv[it+1])
			else:
				raise Exception('Wrong Syntax','-lyap')
		#elif(sys.argv[it] != '-d') and (sys.argv[it] != '-o'):
			#print("Warning! Non used argument: "+(sys.argv[it]))

	mapCML,par=setMap(mapping)
	grad = []
	if(initialMat=='gaussian') or (initialMat=='bessel') or (initialMat=='random'):
		c = CML.CML(initialMat,matLen)
	else:
		raise Exception('Unknown matrix type')
	mat = []
	matMax = []
	xticks = []
	yticks = []
	x,y = 0,0
	
	for coupling in np.arange(cmin,cmax,cdelta):
		mat.append([])
		matMax.append([])
		coupling = round(coupling,ndigits=4)
		for alpha in np.arange(amin,amax,adelta):
			alpha = round(alpha,ndigits=4)
			if(alpha>4.0):
				break
			lyap = []
			lyapMax = []
			print(coupling,alpha)
			c = CML.CML(initialMat,matLen)
			for i in range(nit):
				par = [alpha]
				c.getCML(neigh,mapCML,coupling,par)
				if i<nit-nlyap:
					continue
				if not('-lyap' in sys.argv):
					continue
				try:
					jac = c.getJacobian(maps.dlogisticMap,coupling,par)
					vals, vecs = eigs(jac,k=1,which = 'LM')
					maxEigen = np.fabs(np.real(vals))[0]
					if (maxEigen > 1e-6):
						lyap.append(np.log(maxEigen))
				except ArpackError:
					print("Not convergent for ",i," (LM)")
					
				try:
					vals,vecs = eigs(jac,k=1,which='LR')
					maxEigen = np.fabs(np.real(vals))[0]
					if (maxEigen > 1e-6) :
						lyapMax.append(np.log(maxEigen))
				except ArpackError:
					print("Not convergent for ",i," (LR)")
				#print("Lyapunov: ",np.mean(lyap),"(LM)",np.mean(lyapMax),"(LR)")
			if '-d' in sys.argv:
				np.savetxt("output/"+str(alpha)+"_"+str(coupling)+".txt",c.mat)
				plt.clf()
				plot(c.mat)
				plt.savefig("output/"+str(alpha)+"_"+str(coupling)+".png")
			print("Lyapunov: ",np.mean(lyap),"(LM)",np.mean(lyapMax),"(LR)")
			mat[x].append(np.mean(lyap))
			matMax[x].append(np.mean(lyapMax))
			y=y+1
		x = x+1 
	plt.clf()
	plt.imshow(np.transpose(mat),extent=[cmin,cmax,amin,amax],origin="lower",interpolation="bilinear",cmap=plt.get_cmap('seismic'))
	plt.colorbar()
	plt.savefig("Lyap.png")
	np.savetxt("lyapMat.txt",mat)
	np.savetxt("lyapMaxEigen.txt",matMax)
