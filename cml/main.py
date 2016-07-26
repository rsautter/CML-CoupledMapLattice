import matplotlib.pyplot as plt
import CML
import maps
import sys

def plot(mat):
    plt.clf()
    plt.imshow(mat,cmap=plt.get_cmap('gray'))
    plt.colorbar()


def setMap(mapping):
    if(mapping=='logistic'):
        cmlMap = maps.logisticMap
    elif(mapping=='doubling'):
        cmlMap = maps.doublingMap
    elif(mapping=='SOM'):
        cmlMap = maps.somMap
    else:
        raise Exception('Unsupported map',mapping)
    return cmlMap

if __name__ == "__main__":
    if('-h' in sys.argv) or ('--help' in sys.argv):
        print('================================')
        print('Syntax:')
        print('python main.py commands')
        print('================================')
        print('Commands:')
        print('================================')
        print('-d: display the last iteration of CML')
        print('-mat type matlen: initial condition type(gaussian or bessel), matlen(>0)')
        print('-map map: type of map(logistic, som, doubling)')
        print('-nit nit: number of iterations')
        print('-c coupling: coupling factor(float)')
        print('-o: if set, saves each iteration png files')
        print('-h or --help: display help')
        print('================================')
        exit()
    
    initialMat,matLen,nit,mapping,coupling,neigh,output = 'gaussian',128,50,'logistic',0.5,[(1,0),(0,1),(-1,0),(0,-1)],False
    for it in range(len(sys.argv)):
        if(sys.argv[it] == '-mat'):
            if(it+2 < len(sys.argv)):
                initialMat = sys.argv[it+1]
                matLen = int(sys.argv[it+2])
            else:
                raise Exception('Wrong Syntax','-mat')
        elif(sys.argv[it] == '-map'):
            if(it+1 < len(sys.argv)):
                mapping=sys.argv[it+1]
            else:
                raise Exception('Wrong Syntax','-map')
        elif(sys.argv[it] == '-nit'):
            if(it+1 < len(sys.argv)):
                nit=int(sys.argv[it+1])
            else:
                raise Exception('Wrong Syntax','-nit')
        elif(sys.argv[it] == '-c'):
            if(it+1 < len(sys.argv)):
                coupling=float(sys.argv[it+1])
            else:
                raise Exception('Wrong Syntax','-c')
        #elif(sys.argv[it] != '-d') and (sys.argv[it] != '-o'):
            #print("Warning! Non used argument: "+(sys.argv[it]))

    mapCML=setMap(mapping)
    if(initialMat=='gaussian') or (initialMat=='bessel') :
        c = CML.CML(initialMat,matLen)
    else:
        raise Exception('Unknown matrix type')
    for i in range(nit):
        if ('-o' in sys.argv):
            plot(c.mat)
            plt.savefig('output/it'+str(i)+'.png')
        c.getCML(neigh,mapCML,coupling)
    if ('-d' in sys.argv):
        plot(c.mat)
        plt.show()
