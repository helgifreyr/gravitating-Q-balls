from os import listdir,chdir
import scipy

def get_data(dir):
    chdir(dir)
    datas = scipy.zeros((len(listdir('.')),15))
    i=0
    for folder in sorted(listdir('.')):
        if 'w' in folder:
            data = ''.join([line.replace('\n','').replace(r'{','').replace(r'}','') for line in open(folder+'/tmp.txt')]).split(',')
            a = data[0]
            w = data[1]
            c1 = data[2]
            c2 = data[3]
            c3 = data[4]
            massINF = data[5]
            jINF = data[6]
            Mint = data[7]
            Jint = data[8]
            minf0 = data[9]
            f0H = data[10]
            f1H = data[11]
            maxW = data[12]
            Zm = data[13]
            rm = data[14]
            datas[i] = [a,w,c1,c2,c3,massINF,jINF,Mint,Jint,minf0,f0H,f1H,maxW,Zm,rm]
            i+=1
    chdir('../../')
    return datas
