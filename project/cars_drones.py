# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\潘sir\.spyder2\.temp.py
"""
#[565.0,575.0],[25.0,185.0],[345.0,750.0],[945.0,685.0],[845.0,655.0]
import numpy as np
np.set_printoptions(threshold='nan')
import math
import matplotlib.pyplot as plt
import dongtai

import datetime

numant = 10 #蚂蚁个数

w_weight = 0#无人车初始重量
alpha = 1   #信息素重要程度因子
#beta = 5    #启发函数重要程度因子
beta = 2
rho = 0.1   #信息素的挥发速度
Q = 1
speed_car=1 #无人车速度
speed_drone=2#无人机速度
itermax = 500

lengthaver = np.zeros(itermax) #各代路径的平均长度
lengthbest = np.zeros(itermax) #各代及其之前遇到的最佳路径长度


def getdistmat(coordinates):
    num = coordinates.shape[0]
    distmat = np.zeros((num,num))
    for i in range(num):
        for j in range(i,num):
            distmat[i][j] = distmat[j][i]=np.linalg.norm(coordinates[i]-coordinates[j])#norm 两点距离函数
    return distmat






def Ant_algorithm(distmat,weight,numcity,maxcarweight):##蚁群算法进行初始无人车路径规划

    pheromonetable  = np.ones((numcity,numcity)) # 信息素矩阵  ones返回初始值为1的矩阵
    pathtable = np.zeros((numant,numcity*2-1)).astype(int) #路径记录表    
    pathbest = np.zeros((itermax,numcity*2-1)) # 各代及其之前遇到的最佳路径长度
    etatable = 1.0/(distmat+np.diag([1e10]*numcity)) #启发函数矩阵，表示蚂蚁从城市i转移到矩阵j的期望程度  diag 形成对角矩阵  
    iter = 0
    
    for i in range(numant):
        for j in range(numcity*2-1):
            pathtable[i][j] = -1    
    
    while iter < itermax:
    
        '''
        # 随机产生各个蚂蚁的起点城市
        if numant <= numcity:#城市数比蚂蚁数多
            pathtable[:,0] = np.random.permutation(range(0,numcity))[:numant]#产生numant个数，放在第0列上
        else: #蚂蚁数比城市数多，需要补足
            pathtable[:numcity,0] = np.random.permutation(range(0,numcity))[:]
            pathtable[numcity:,0] = np.random.permutation(range(0,numcity))[:numant-numcity]
        '''

        pathtable[:,0] = 0
        length = np.zeros(numant) #计算各个蚂蚁的路径距离

        for i in range(numant):
    
            w_weight = 0
            visiting = pathtable[i,0] # 当前所在的城市
    
            #visited = set() #已访问过的城市，防止重复
            #visited.add(visiting) #增加元素
            unvisited = set(range(numcity))#未访问的城市
            unvisited.remove(visiting) #删除元素
            j=1
            while(len(unvisited) != 0):

                listunvisited = list(unvisited)
                flag=0
                probtrans = np.zeros(len(listunvisited))
                can_unvisited=[0]*len(listunvisited)
                
                for l in range (len(listunvisited)):
                    if(w_weight + weight[listunvisited[l]] <= maxcarweight):#判断是否超出无人车最大载重，若没有，则加入无人车
                        can_unvisited[l] = listunvisited[l]

                for k in range(len(can_unvisited)):
                    if(can_unvisited[k]>0):#计算信息素
    
                        probtrans[k] = np.power(pheromonetable[visiting][listunvisited[k]],alpha)\
                                *np.power(etatable[visiting][listunvisited[k]],beta)
                        flag = 1
                if(flag == 1):
                    cumsumprobtrans = (probtrans/sum(probtrans)).cumsum()#使最大概率变为1，轮盘必会找到城市
                    
                    cumsumprobtrans -= np.random.rand()

                    for cc in range(len(cumsumprobtrans)):
                        if(cumsumprobtrans[cc]>0):
                            kw=cc
                            break
                    k = listunvisited[[kw][0]]
                    pathtable[i,j] = k

                    w_weight = w_weight + weight[k]
                    unvisited.remove(k)

                    length[i] += distmat[visiting][k]
        
                    visiting = k
                    j += 1

                else:
                    length[i] += distmat[visiting][0]
                    w_weight = 0
                    visiting = 0
                    pathtable[i,j] = 0
                    

                    j += 1

            length[i] += distmat[visiting][0]
            pathtable[i,j] = 0      #添加结束点0
    
    

        # 包含所有蚂蚁的一个迭代结束后，统计本次迭代的若干统计参数
        lengthaver[iter] = length.mean()
    
        if iter == 0:#寻找总路程最短
            lengthbest[iter] = length.min()
            pathbest[iter] = pathtable[length.argmin()].copy()#argmin 最小值下标      
        else:
            if length.min() > lengthbest[iter-1]:
                lengthbest[iter] = lengthbest[iter-1]
                pathbest[iter] = pathbest[iter-1].copy()
    
            else:
                lengthbest[iter] = length.min()
                pathbest[iter] = pathtable[length.argmin()].copy()    
    
    
        # 更新信息素
        changepheromonetable = np.zeros((numcity,numcity))
        for i in range(numant):
            #nnnum=0
            for j in range(numcity*2-2):
                if((pathtable[i,j] != -1) and (pathtable[i,j] != -1) ):
                    
                    #if(distmat[pathtable[i,j]][pathtable[i,j+1]] == 0):
                        #continue
                    
                    changepheromonetable[pathtable[i,j]][pathtable[i,j+1]] += Q/distmat[pathtable[i,j]][pathtable[i,j+1]]

    
            #changepheromonetable[pathtable[i,nnnum]][pathtable[i,0]] += Q/distmat[pathtable[i,nnnum]][pathtable[i,0]]
    
        pheromonetable = (1-rho)*pheromonetable + changepheromonetable
    
    
        iter += 1 #迭代次数指示器+1

    #print "outyiqun"
    return pathbest


def print_path(pathbest,numcity):#建立多维数组，每行代表一条路径
    bestpath = pathbest[-1]

    num = 0
    m=1
    x = 0
    y = 1
    flag = 1
    for j in range(len(bestpath)):
        if(bestpath[j] == 0):
            num += 1
    finalpathtable = np.zeros((num-1,numcity+1)).astype(int)
    for i in range(num-1):
        for j in range(numcity+1):
            finalpathtable[i][j] = -1 
    finalpathtable[:,0] = 0
    while(flag == 1):
        if(bestpath[m] == 0 and bestpath[m+1] == -1):
            break
        if(bestpath[m] == 0):
            finalpathtable[x][y] = 0
            x += 1
            y = 1

        else:
            finalpathtable[x][y] = bestpath[m]
            y += 1
        m += 1

    finalpathtable[x][y] = 0
    return finalpathtable        


def outqingjian(path1,weight,numcity,coordinates,weightflag,drone_weight):

    jilubiao = [[[0 for a in range(2)] for b in range(len(path1[0]))] for c in range(len(path1))]

    m=[]
    zonglu=[]
    zhongwe=[]
    for i in range (len(weight)):
        if(weight[i] < drone_weight and weight[i] > 0):
             weightflag[i] = 0                   
        else:
             weightflag[i] = 1
    weightflag[0] = 1
    for i in range (len(path1)):      
        for k in path1[i]:#将所有点加入有个分组，和所有重件点加入一个分组
            if(k!=-1):
                zonglu.append(coordinates[k])
            if(weightflag[k] == 1 and k !=-1):
                zhongwe.append(coordinates[k])
    for i in range (len(path1)):
        count=0
        l=[]       
        for k in path1[i]:
            if(weightflag[k] == 1 and k !=-1):
                l.append(coordinates[k])
                jilubiao[i][count][0]=k
                jilubiao[i][count][1]=1
            else:
                jilubiao[i][count][0]=k
                jilubiao[i][count][1]=0
            count+=1

        m.append(l)
    #print "outqingjian"
    return m,weightflag,jilubiao,zonglu,zhongwe
def outqingjian_all(path1,weight,numcity,coordinates,weightflag):

    jilubiao = [[[0 for a in range(2)] for b in range(len(path1[0]))] for c in range(len(path1))]

    m=[]
    zonglu=[]
    zhongwe=[]
    '''
    for i in range (len(weight)):
        if(weight[i] < drone_weight and weight[i] > 0):
             weightflag[i] = 0                   
        else:
             weightflag[i] = 1
    '''
    #for i in range(len(weightflag))
    for i in range (len(path1)):
        path_len=path1[i].tolist().index(-1)-1
        if(path_len >= 20):
            weightflag[path1[i][6]] = 1
            weightflag[path1[i][12]] = 1
            weightflag[path1[i][18]] = 1
        elif(path_len >= 15):
            weightflag[path1[i][6]] = 1
            weightflag[path1[i][12]] = 1
        elif (path_len >= 10):
            weightflag[path1[i][3]] = 1
            weightflag[path1[i][8]] = 1
        elif (path_len >= 5):
            weightflag[path1[i][3]] = 1

    weightflag[0] = 1
    for i in range (len(path1)):
        for k in path1[i]:#将所有点加入有个分组，和所有重件点加入一个分组
            if(k!=-1):
                zonglu.append(coordinates[k])
            if(weightflag[k] == 1 and k !=-1):
                zhongwe.append(coordinates[k])
    for i in range (len(path1)):
        count=0
        l=[]
        for k in path1[i]:
            if(weightflag[k] == 1 and k !=-1):
                l.append(coordinates[k])
                jilubiao[i][count][0]=k
                jilubiao[i][count][1]=1
            else:
                jilubiao[i][count][0]=k
                jilubiao[i][count][1]=0
            count+=1

        m.append(l)
    #print "outqingjian"
    return m,weightflag,jilubiao,zonglu,zhongwe

def wurenji(jilu,coordinates,weight,num_drone,drone_weight):#无人机路径

     if(num_drone==1):     
         drflag=[]#将每个闭合路径的轻件点放在一起
         dr_connect=[]
         jiaohui=[]
         first_drone=[]
         for i in range (len(jilu)):
             drone_road=[]
             maxweight=0#各点重量加起来是否超过无人机载重
             count=0
             while (count <= len(jilu[i]) ):

                 if(count > 0 and count == len(jilu[i])-1):
                     #print "tiaochu"
                     break
                 elif(jilu[i][count][1]==1 and jilu[i][count+1][1]==0 and jilu[i][count+1][0]==-1):

                     break
                 elif(jilu[i][count][1]==1 and jilu[i][count+1][1]==0 and jilu[i][count+1][0]!=-1):#遇到第一个轻件点起飞
                     #print jilu[i][count]," up"
                     start=jilu[i][count][0]
                     dr_connect.append(coordinates[jilu[i][count][0]][0])
                     drone_road.append(coordinates[jilu[i][count][0]])
                 elif(jilu[i][count][1]==0):
                     maxweight += weight[jilu[i][count][0]]
                     if(maxweight>drone_weight):
                         maxweight = 0
                         end=count

                         while jilu[i][end][1]!=1:
                             end+=1

                         end=jilu[i][end][0]
                         #print start,end,jilu[i][count-1][0]
                         temp=[]
                         temp=jiangluo(start,end,jilu[i][count-1][0],coordinates)
                         drone_road.append(temp)#尽管下一个点为轻件点，但超过无人机重量，降落
                         #dr_connect.append(1)
                         #print temp," down**"
                         #print temp," up**"
                         jiaohui.append(temp[0])
                         '''
                         temp=[]
                         temp=jiangluo(start,end,jilu[i][count-1][0],coordinates)
                         if(temp[0]<=coordinates[start][0]):
                             drone_road.append(coordinates[start])
                         elif(temp[0]>=coordinates[end][0]):
                             drone_road.append(coordinates[end])
                         else:
                             drone_road.append(temp)#尽管下一个点为轻件点，但超过无人机重量，降落
                             dr_connect.append(1)                     
                             print temp," down**"
                             print temp," up**"
                             jiaohui.append(jiangluo(start,end,jilu[i][count-1][0],coordinates)[0])
                         '''
                         count=count-1
                     else:
                         drone_road.append(coordinates[jilu[i][count][0]])
                         #dr_connect.append(0)
                         if(jilu[i][count+1][1]==1):
                             #print "count*",count
                             #print jilu[i][count+1]," down"
                             maxweight=0
                             dr_connect.append(coordinates[jilu[i][count+1][0]][0])
                             drone_road.append(coordinates[jilu[i][count+1][0]])#遇到下一个点为重件点是降落

                 count = count + 1
             drflag.append(drone_road)#求每个环路的轻件点飞行路径
             first_drone.append(drone_road)
     else:
         drflag=[]#将每个闭合路径的轻件点放在一起
         dr_connect=[]
         jiaohui=[]
         first_drone=[] #通过比较，第一辆无人机飞行时间最长，用于和无人机比较
         for i in range (len(jilu)):
             temp=[]
             for pp in (jilu[i]):
                 temp.append(pp)
             #print "temp=",temp
             drone_road=[]
                 
             for j in range (num_drone):
                 maxweight=0#各点重量加起来是否超过无人机载重
                 count=0
                 

                 while (count <= len(temp) ):

                     if(temp[count][0]==-1):
                         break
                     elif(temp[count][1]==1 and temp[count+1][1]==0 and temp[count+1][0]==-1):

                         break

                     elif(temp[count][1]==1 and temp[count+1][1]==0 and temp[count+1][0]!=-1):#遇到第一个轻件点起飞
                         #print temp[count]," up"
                         start=temp[count][0]
                         dr_connect.append(coordinates[temp[count][0]][0])
                         drone_road.append(coordinates[temp[count][0]])
                     elif(temp[count][1]==0 and temp[count][0]!=-1):
                         maxweight += weight[temp[count][0]]

                         if(maxweight>drone_weight):
                             maxweight = 0
                             end=count

                             while temp[end][1]!=1:
                                 end+=1

                             end=temp[end][0]

                             tempp=[]
                             tempp=jiangluo(start,end,temp[count-1][0],coordinates)
                             drone_road.append(tempp)#尽管下一个点为轻件点，但超过无人机重量，降落
                             #dr_connect.append(1)
                             #print tempp," down**"
                             #print tempp," up**"
                             jiaohui.append(tempp[0])
                             count=count-1
                         else:

                             drone_road.append(coordinates[temp[count][0]])
                             #dr_connect.append(0)
                             del temp[count]

                             flag_heavy=0
                             for k in range(num_drone-j):
                                 if(temp[count+k][1]==1):
                                     flag_heavy=1

                                     #print temp[count+k]," down"
                                     maxweight=0
                                     dr_connect.append(coordinates[temp[count+k][0]][0])
                                     drone_road.append(coordinates[temp[count+k][0]])#遇到下一个点为重件点是降落
                                     count=count+k-1
                                     break
                                     drflag.append(coordinates[jilu[i][k+1][0]][0])

                             if(flag_heavy==0):
                                 count=count+num_drone-2-j

                     count = count + 1
                 if(j==0):
                     tttemp=[]

                     for aa in drone_road:
                         tttemp.append(aa)
                    
                     first_drone.append(tttemp)

             drflag.append(drone_road)#求每个环路的轻件点飞行路径

     return drflag,dr_connect,jiaohui,first_drone


def jiangluo(start,end,jiang,coordinates):#无人机降落坐标
    temp=[0]*2
    kk=(coordinates[end][1]-coordinates[start][1])/(coordinates[end][0]-coordinates[start][0])
    if(kk==0):
        kk=0.000001
    bb=coordinates[end][1]-kk*coordinates[end][0]
    p=-1.0/kk #垂直线斜率
    q=coordinates[jiang][1]-p*coordinates[jiang][0] #垂直线方程中的y=kx+b
    u=(bb-q)/(p-kk) #交点横坐标
    v=kk*u+bb #交点纵坐标
    temp[0]=u
    temp[1]=v     
    return temp

def picture(zhongjian,drflag,zonglu,dr_connect,jiaohui):
    connect=np.array(dr_connect)
    fflag=0
    cc = np.array(zhongjian)
    
    zz=np.array(zonglu)
    #print "cc"
    #print cc
    
    numcity = len(cc)
    
    numcity2=len(zz)
    bestpath=[]
    bestpath1=[]
    bestpath2=[]



        
    for k in range(len(drflag)) : #对轻件点画图
        aa= np.array(drflag[k])
        #print "aa"
        #print aa
        numcity1 = len(aa)
        if(numcity1==0):
            break
        '''
        for kk in aa:#打印无人车路途中回车装货
            for ll in jiaohui:
                if(kk[0]==ll):
                    plt.plot(kk[0],kk[1],'bs')
                    #print "jiaohui***",kk[0]
                    break
        '''

        for j in range (numcity1):
            bestpath1.append(j)
        plt.plot(aa[:,0],aa[:,1],'r.')
        for j in range(numcity1-1):
        
            m,n = bestpath1[j],bestpath1[j+1]
      #  for k in range(len(drflag)):
    #        if(aa[m][0]==)

            for xx in range(len(connect)-1):
                if(aa[m][0]==connect[xx] and aa[n][0]==connect[xx+1] and aa[n][0] != 0 and aa[m][0] != 0):
                    #print "connect=",connect[xx]
                    fflag=1
                    break
                
            if(fflag == 1):
                fflag=0
                continue

            plt.plot([aa[m][0],aa[n][0]],[aa[m][1],aa[n][1]],'k--',linewidth=1.0)
    
    for i in range (numcity):
        bestpath.append(i)

    plt.plot(cc[:, 0], cc[:, 1], 'bs')

    plt.xlim([-150,150])
    plt.ylim([-150,150])
    for i in range(numcity - 1):
        m, n = bestpath[i], bestpath[i + 1]
        plt.plot([cc[m][0], cc[n][0]], [cc[m][1], cc[n][1]], 'k-',linewidth=1.0)

    plt.plot([cc[bestpath[0]][0], cc[n][0]], [cc[bestpath[0]][1], cc[n][1]], 'k-',linewidth=1.0)
    
    
     
    ax=plt.gca()
    ax.set_title("vrp3")
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y_axis')
    
    plt.savefig('vrp3.png',dpi=500,bbox_inches='tight')
    plt.close() 
    for k in range (numcity2):
        bestpath2.append(k) 
    plt.plot(zz[:,0],zz[:,1],'bs')
    plt.xlim([-150,150])
    plt.ylim([-150,150])
    for k in range(numcity2-1):
        
        m,n = bestpath2[k],bestpath2[k+1]
        plt.plot([zz[m][0],zz[n][0]],[zz[m][1],zz[n][1]],'k',linewidth=1.0)
    
    plt.plot([zz[bestpath2[0]][0],zz[n][0]],[zz[bestpath2[0]][1],zz[n][1]],'k',linewidth=1.0)
    ax=plt.gca()
    ax.set_title("vrp2")
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y_axis')
    
    plt.savefig('vrp2.png',dpi=500,bbox_inches='tight')
    plt.close() 
    
    plt.plot(zz[:,0],zz[:,1],'bs')
    
    plt.xlim([-150,150])
    plt.ylim([-150,150])
    
    ax=plt.gca()
    ax.set_title("vrp1")
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y_axis')
    
    plt.savefig('vrp1.png',dpi=500,bbox_inches='tight')
    plt.close() 



def main():

    maxcarweight = 400  # 无人车最大载重量

    drone_weight = 20# 无人机载重量
    #coordinates = np.array([[1.0,0.0],[2.0,5.0],[4.0,2.0],[5.0,4.0],[10.0,3.0],
                        #[20.0,96.0],[86.0,45.0],[-3.0,6.8],[-6.8,6.3],[9.8,-6.6]])
    #weight=[0,3,5,6,2,4,3,2,1,10]
    '''
    print "input the number of city:"
    num_city = int(raw_input())#输入快件点数量，包含轻件点和重件点
    #num_drone =1
    coordinates = np.random.uniform(-100, 100, size=(num_city,2))#在-100--100的坐标系内生成坐标
    '''
    print "input the number of drone:"



    num_drone = int(raw_input())#输入无人机数量

    starttime = datetime.datetime.now()

    coordinates =[[	35	,	35	]	,
[	41	,	49	]	,
[	35	,	17	]	,
[	55	,	45	]	,
[	55	,	20	]	,
[	15	,	30	]	,
[	25	,	30	]	,
[	20	,	50	]	,
[	10	,	43	]	,
[	55	,	60	]	,
[	30	,	60	]	,
[	20	,	65	]	,
[	50	,	35	]	,
[	30	,	25	]	,
[	15	,	10	]	,
[	30	,	5	]	,
[	10	,	20	]	,
[	5	,	30	]	,
[	20	,	40	]	,
[	15	,	60	]	,
[	45	,	65	]	,
[	45	,	20	]	,
[	45	,	10	]	,
[	55	,	5	]	,
[	65	,	35	]	,
[	65	,	20	]	,
[	45	,	30	]	,
[	35	,	40	]	,
[	41	,	37	]	,
[	64	,	42	]	,
[	40	,	60	]	,
[	31	,	52	]	,
[	35	,	69	]	,
[	53	,	52	]	,
[	65	,	55	]	,
[	63	,	65	]	,
[	2	,	60	]	,
[	20	,	20	]	,
[	5	,	5	]	,
[	60	,	12	]	,
[	40	,	25	]	,
[	42	,	7	]	,
[	24	,	12	]	,
[	23	,	3	]	,
[	11	,	14	]	,
[	6	,	38	]	,
[	2	,	48	]	,
[	8	,	56	]	,
[	13	,	52	]	,
[	6	,	68	]	,
[	47	,	47	]	,
[	49	,	58	]	,
[	27	,	43	]	,
[	37	,	31	]	,
[	57	,	29	]	,
[	63	,	23	]	,
[	53	,	12	]	,
[	32	,	12	]	,
[	36	,	26	]	,
[	21	,	24	]	,
[	17	,	34	]	,
[	12	,	24	]	,
[	24	,	58	]	,
[	27	,	69	]	,
[	15	,	77	]	,
[	62	,	77	]	,
[	49	,	73	]	,
[	67	,	5	]	,
[	56	,	39	]	,
[	37	,	47	]	,
[	37	,	56	]	,
[	57	,	68	]	,
[	47	,	16	]	,
[	44	,	17	]	,
[	46	,	13	]	,
[	49	,	11	]	,
[	49	,	42	]	,
[	53	,	43	]	,
[	61	,	52	]	,
[	57	,	48	]	,
[	56	,	37	]	,
[	55	,	54	]	,
[	15	,	47	]	,
[	14	,	37	]	,
[	11	,	31	]	,
[	16	,	22	]	,
[	4	,	18	]	,
[	28	,	18	]	,
[	26	,	52	]	,
[	26	,	35	]	,
[	31	,	67	]	,
[	15	,	19	]	,
[	22	,	22	]	,
[	18	,	24	]	,
[	26	,	27	]	,
[	25	,	24	]	,
[	22	,	27	]	,
[	25	,	21	]	,
[	19	,	21	]	,
[	20	,	26	]	,
[	18	,	18	]	]








    '''
    total_x=0
    total_y=0
    for i in range (len(coordinates)):
        total_x = coordinates[i][0] + total_x
        total_y = coordinates[i][1] + total_y
    total_x=total_x/len(coordinates)
    total_y=total_y/len(coordinates)
    start_pon=[]
    start_pon.append(total_x)
    start_pon.append(total_y)
    coordinates.insert(0,start_pon)
    

    
    #print start_pon
    coordinates[0][0] = 0.0
    coordinates[0][1] = 0.0
    '''
    coordinates = np.array(coordinates)
    num_city=len(coordinates)
    #weight = np.random.uniform(0, 5, size=(num_city))#随机生成快件点的重量

    weight =[	0	,
	10	,
	7	,
	13	,
	19	,
	26	,
	3	,
	5	,
	9	,
	16	,
	16	,
	12	,
	19	,
	23	,
	20	,
	8	,
	19	,
	2	,
	12	,
	17	,
	9	,
	11	,
	18	,
	29	,
	3	,
	6	,
	17	,
	16	,
	16	,
	9	,
	21	,
	27	,
	23	,
	11	,
	14	,
	8	,
	5	,
	8	,
	16	,
	31	,
	9	,
	5	,
	5	,
	7	,
	18	,
	16	,
	1	,
	27	,
	36	,
	30	,
	13	,
	10	,
	9	,
	14	,
	18	,
	2	,
	6	,
	7	,
	18	,
	28	,
	3	,
	13	,
	19	,
	10	,
	9	,
	20	,
	25	,
	25	,
	36	,
	6	,
	5	,
	15	,
	25	,
	9	,
	8	,
	18	,
	13	,
	14	,
	3	,
	23	,
	6	,
	26	,
	16	,
	11	,
	7	,
	41	,
	35	,
	26	,
	9	,
	15	,
	3	,
	1	,
	2	,
	22	,
	27	,
	20	,
	11	,
	12	,
	10	,
	9	,
	17	]











    weight[0] = 0.0 
    weight=np.array(weight)
    #print coordinates
    #print weight
    weightflag = [0 for x in range(0, num_city)] #重件点与轻件点标志
    numcity = coordinates.shape[0] #城市个数
    distmat = getdistmat(coordinates) #城市的距离矩阵
    all_qingjian=1
    for i in range(len(weight)):
        if(weight[i] <= drone_weight):
            continue
        else:
            all_qingjian=0

    if(all_qingjian==0):
        pathbest = Ant_algorithm(distmat,weight,numcity,maxcarweight)

        path1 = print_path(pathbest,numcity)
        if(num_drone == 0):
            #print path1
            total_dis=0
            zhenghe_path=[]
            for i in range (len(path1)):
                for j in range (len(path1[i])):
                    if(path1[i][j]==-1):
                        break
                    else:
                        zhenghe_path.append(path1[i][j])
            for i in range(len(zhenghe_path)-1):
                total_dis = total_dis + np.linalg.norm(coordinates[zhenghe_path[i+1]]-coordinates[zhenghe_path[i]])
            total_time=total_dis/speed_car
            endtime = datetime.datetime.now()

            print (endtime - starttime).seconds
            print "total_dis=",total_dis
            print "total_time=",total_time
        else:
            #print path1
            jilu = [[[0 for a in range(2)] for b in range(len(path1[0]))] for c in range(len(path1))]
            mostwei_arr,weightflag,jilu,zonglu,zhongwe=outqingjian(path1,weight,numcity,coordinates,weightflag,drone_weight)#mostwei_arr为分割好的每段重件点，
            #if (num_drone != 0):
            drflag,dr_connect,jiaohui,first_drone=wurenji(jilu,coordinates,weight,num_drone,drone_weight)

            #print "jilu=++++++++++"
            #print jilu
            #print dr_connect
            #print zonglu
            #print mostwei_arr
            ######算距离
            mostwei_arr1=[]#进行数组转化，判断是否两个挨着的为重件
            drflag1=[]
            dis_flag=[0 for a in range(len(drflag1))]
            for i in range (len(mostwei_arr)):
                for j in range (len(mostwei_arr[i])):
                    mostwei_arr1.append(mostwei_arr[i][j])
            for i in range (len(drflag)):
                for j in range (len(drflag[i])):
                    drflag1.append(drflag[i][j])
            mostwei_arr1=np.array(mostwei_arr1)
            drflag1=np.array(drflag1)
            dis_flag=[0 for a in range(len(drflag1))]
            #print "mostwei_arr1=",mostwei_arr1
            #print "dis_flag1=",drflag1
            #print "first_drone=",first_drone
            for i in range (len(drflag1)):
                #print drflag1[i] in mostwei_arr1
                if(drflag1[i] in mostwei_arr1):
                    #print (mostwei_arr1-drflag1[i]).any()
                    dis_flag[i]=1#记录无人机路径中的轻件和重件

            #print "dis_flag=",dis_flag
            drone_dis=0
            car_dis=0
            for i in range (len(drflag1)-1):
                if(dis_flag[i+1]==1 and dis_flag[i]==1):
                    continue
                else:
                    drone_dis = drone_dis + np.linalg.norm(drflag1[i+1]-drflag1[i])
            #print mostwei_arr1
            for i in range (len(mostwei_arr1)-1):
                car_dis= car_dis + np.linalg.norm(mostwei_arr1[i+1]-mostwei_arr1[i])
            total_dis = car_dis+drone_dis

            #算时间
            total_time=0
            #print mostwei_arr
            #print first_drone

            for i in range (len(jilu)):
                for j in range (len(jilu[i])):
                    if (jilu[i][j][1] == 1 and jilu[i][j+1][0] == -1):
                        break
                    elif (jilu[i][j][1] == 1 and jilu[i][j+1][1] == 1):
                        total_time = total_time + np.linalg.norm(coordinates[jilu[i][j+1][0]]-coordinates[jilu[i][j][0]])/speed_car
                    elif (jilu[i][j][1] == 1 and jilu[i][j+1][1] == 0):
                        time_tmpcar=0
                        time_tmpdrone=0
                        start_point=coordinates[jilu[i][j][0]]
                        k=j+1
                        while(jilu[i][k][1] != 1):
                            #time_tmpcar = time_tmpcar + np.linalg.norm(coordinates[jilu[i][k][0]]-coordinates[jilu[i][k-1][0]])/speed_car
                            k=k+1
                            #print "k***=",k
                        time_tmpcar = time_tmpcar + np.linalg.norm(
                            coordinates[jilu[i][k][0]] - coordinates[jilu[i][j][0]]) / speed_car
                        end_point=coordinates[jilu[i][k][0]]
                        first_drone1=np.array(first_drone[i])
                        index_s=first_drone1.tolist().index(start_point.tolist())
                        index_e=first_drone1.tolist().index(end_point.tolist())
                        #print index_s,index_e
                        for kk in range(index_s,index_e):
                            time_tmpdrone = time_tmpdrone + np.linalg.norm(first_drone1[kk+1]-first_drone1[kk])/speed_drone
                        if(time_tmpcar>time_tmpdrone):
                            total_time=total_time+time_tmpcar
                        else:
                            total_time=total_time+time_tmpdrone
                    elif (jilu[i][j][1] == 0):
                        continue
            endtime = datetime.datetime.now()

            print (endtime - starttime).seconds
            print "car_dis,drone_dis,total_dis=", car_dis, drone_dis, total_dis
            print "total_time=",total_time
            '''
            picture(zhongwe,drflag,zonglu,dr_connect,jiaohui)

            print "point or not:"  #动态揽收
            #ppp = int(raw_input())
            ppp=1
            if (ppp == 1):
                dongtai.dongtai(mostwei_arr,jilu,coordinates,speed_car,weight,maxcarweight,num_drone,drone_weight)
        #print drone_road
            '''

        #picture(zhongjian,drone_road,drflag,zonglu)
            '''
            print "jilu=：",jilu
            print "dr_connect=：",dr_connect
            print "zonglu=：",zonglu
            print "mostwei_arr=：",mostwei_arr
            print "drflag:=",drflag
            picture(zhongwe,drflag,zonglu,dr_connect,jiaohui)
            '''

        ##wurenji(path1,numcity)
        #print mostwei_arr,weightflag
        #zuixiaojuli(mostwei_arr,path1,coordinates,weightflag)
    else:##全是轻件场景
        #maxcarweight1=1.5*drone_weight
        pathbest = Ant_algorithm(distmat, weight, numcity, maxcarweight)

        path1 = print_path(pathbest, numcity)
        if (num_drone == 0):
            # print path1
            total_dis = 0
            zhenghe_path = []
            for i in range(len(path1)):
                for j in range(len(path1[i])):
                    if (path1[i][j] == -1):
                        break
                    else:
                        zhenghe_path.append(path1[i][j])
            for i in range(len(zhenghe_path) - 1):
                total_dis = total_dis + np.linalg.norm(coordinates[zhenghe_path[i + 1]] - coordinates[zhenghe_path[i]])
            total_time = total_dis / speed_car
            endtime = datetime.datetime.now()

            print (endtime - starttime).seconds
            print "total_dis=", total_dis
            print "total_time=", total_time
        else:
            drone_weight1=drone_weight*0.8
            mostwei_arr, weightflag, jilu, zonglu, zhongwe = outqingjian_all(path1,weight,numcity,coordinates,weightflag)
            #mostwei_arr, weightflag, jilu, zonglu, zhongwe = outqingjian(path1, weight, numcity, coordinates,weightflag, drone_weight1)# mostwei_arr为分割好的每段重件点，
            # if (num_drone != 0):
            drflag, dr_connect, jiaohui, first_drone = wurenji(jilu, coordinates, weight, num_drone, drone_weight)
            '''
            # print "jilu=++++++++++"
            print "jilu=",jilu
            # print dr_connect
            print "zonglu=",zonglu
            print "mostwei_arr=",mostwei_arr
            '''
            ######算距离
            mostwei_arr1 = []  # 进行数组转化，判断是否两个挨着的为重件
            drflag1 = []
            dis_flag = [0 for a in range(len(drflag1))]
            for i in range(len(mostwei_arr)):
                for j in range(len(mostwei_arr[i])):
                    mostwei_arr1.append(mostwei_arr[i][j])
            for i in range(len(drflag)):
                for j in range(len(drflag[i])):
                    drflag1.append(drflag[i][j])
            mostwei_arr1 = np.array(mostwei_arr1)
            #print "mostwei_arr1:  ",mostwei_arr1
            drflag1 = np.array(drflag1)
            dis_flag = [0 for a in range(len(drflag1))]
            # print "mostwei_arr1=",mostwei_arr1
            # print "dis_flag1=",drflag1
            # print "first_drone=",first_drone
            for i in range(len(drflag1)):
                # print drflag1[i] in mostwei_arr1
                if (drflag1[i] in mostwei_arr1):
                    # print (mostwei_arr1-drflag1[i]).any()
                    dis_flag[i] = 1  # 记录无人机路径中的轻件和重件

            # print "dis_flag=",dis_flag
            drone_dis = 0
            car_dis = 0
            for i in range(len(drflag1) - 1):
                if (dis_flag[i + 1] == 1 and dis_flag[i] == 1):
                    continue
                else:
                    drone_dis = drone_dis + np.linalg.norm(drflag1[i + 1] - drflag1[i])
            for i in range(len(mostwei_arr1) - 1):
                car_dis = car_dis + np.linalg.norm(mostwei_arr1[i + 1] - mostwei_arr1[i])

            total_dis = car_dis + drone_dis

            # 算时间
            total_time = 0
            # print mostwei_arr
            # print first_drone

            for i in range(len(jilu)):
                for j in range(len(jilu[i])):
                    if (jilu[i][j][1] == 1 and jilu[i][j + 1][0] == -1):
                        break
                    elif (jilu[i][j][1] == 1 and jilu[i][j + 1][1] == 1):
                        total_time = total_time + np.linalg.norm(
                            coordinates[jilu[i][j + 1][0]] - coordinates[jilu[i][j][0]]) / speed_car
                    elif (jilu[i][j][1] == 1 and jilu[i][j + 1][1] == 0):
                        time_tmpcar = 0
                        time_tmpdrone = 0
                        start_point = coordinates[jilu[i][j][0]]
                        k = j + 1
                        while (jilu[i][k][1] != 1):
                            #time_tmpcar = time_tmpcar + np.linalg.norm(
                                #coordinates[jilu[i][k][0]] - coordinates[jilu[i][k - 1][0]]) / speed_car
                            k = k + 1
                            # print "k***=",k
                        time_tmpcar = time_tmpcar + np.linalg.norm(
                                coordinates[jilu[i][k][0]] - coordinates[jilu[i][j][0]]) / speed_car
                        end_point = coordinates[jilu[i][k][0]]
                        first_drone1 = np.array(first_drone[i])
                        index_s = first_drone1.tolist().index(start_point.tolist())
                        index_e = first_drone1.tolist().index(end_point.tolist())
                        # print index_s,index_e
                        for kk in range(index_s, index_e):
                            time_tmpdrone = time_tmpdrone + np.linalg.norm(
                                first_drone1[kk + 1] - first_drone1[kk]) / speed_drone
                        if (time_tmpcar > time_tmpdrone):
                            total_time = total_time + time_tmpcar
                        else:
                            total_time = total_time + time_tmpdrone
                    elif (jilu[i][j][1] == 0):
                        continue
            endtime = datetime.datetime.now()

            print (endtime - starttime).seconds
            print "car_dis,drone_dis,total_dis=", car_dis, drone_dis, total_dis
            print "total_time=", total_time

    
    
if __name__ == "__main__":
    main()
    
    
    
