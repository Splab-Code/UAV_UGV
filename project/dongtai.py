# -*- coding: utf-8 -*-
"""
Created on Mon Apr 09 10:33:30 2018

@author: 潘sir
"""

import numpy as np
np.set_printoptions(threshold='nan')
import cars_drones
import matplotlib.pyplot as plt
numant = 10 #蚂蚁个数

w_weight = 0#无人车初始重量
#drone_weight = 1.8 #无人机载重量
alpha = 1   #信息素重要程度因子
#beta = 5    #启发函数重要程度因子
beta = 2
rho = 0.1   #信息素的挥发速度
Q = 1
speed_car=1 #无人车速度
speed_drone=2#无人机速度
itermax = 300
lengthaver = np.zeros(itermax) #各代路径的平均长度
lengthbest = np.zeros(itermax) #各代及其之前遇到的最佳路径长度

def Ant_algorithm(distmat, weight, numcity,maxcarweight):  ##蚁群算法进行初始无人车路径规划

    pheromonetable = np.ones((numcity, numcity))  # 信息素矩阵  ones返回初始值为1的矩阵
    pathtable = np.zeros((numant, numcity * 2 - 1)).astype(int)  # 路径记录表
    pathbest = np.zeros((itermax, numcity * 2 - 1))  # 各代及其之前遇到的最佳路径长度
    etatable = 1.0 / (distmat + np.diag([1e10] * numcity))  # 启发函数矩阵，表示蚂蚁从城市i转移到矩阵j的期望程度  diag 形成对角矩阵
    iter = 0

    for i in range(numant):
        for j in range(numcity * 2 - 1):
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
        # print iter
        pathtable[:, 0] = 0
        length = np.zeros(numant)  # 计算各个蚂蚁的路径距离
        # pathtable
        for i in range(numant):

            w_weight = 0
            visiting = pathtable[i, 0]  # 当前所在的城市

            # visited = set() #已访问过的城市，防止重复
            # visited.add(visiting) #增加元素
            unvisited = set(range(numcity))  # 未访问的城市
            unvisited.remove(visiting)  # 删除元素
            j = 1
            while (len(unvisited) != 0):
                # print "j=",j
                # print "visisting=",visiting
                listunvisited = list(unvisited)
                flag = 0
                probtrans = np.zeros(len(listunvisited))
                can_unvisited = [0] * len(listunvisited)

                for l in range(len(listunvisited)):
                    #if (w_weight + weight[listunvisited[l]] <= maxcarweight):  # 判断是否超出无人车最大载重，若没有，则加入无人车
                    can_unvisited[l] = listunvisited[l]
                # print "can_unvisited=",can_unvisited
                # print "listunvisited=",listunvisited
                for k in range(len(can_unvisited)):
                    if (can_unvisited[k] > 0):  # 计算信息素

                        probtrans[k] = np.power(pheromonetable[visiting][listunvisited[k]], alpha) \
                                       * np.power(etatable[visiting][listunvisited[k]], beta)
                        flag = 1
                if (flag == 1):
                    cumsumprobtrans = (probtrans / sum(probtrans)).cumsum()  # 使最大概率变为1，轮盘必会找到城市

                    cumsumprobtrans -= np.random.rand()
                    # print "listunvisited=",listunvisited
                    # print "cumsumprobtrans=",cumsumprobtrans
                    # k = listunvisited[find(cumsumprobtrans>0)[0]]
                    for cc in range(len(cumsumprobtrans)):
                        if (cumsumprobtrans[cc] > 0):
                            kw = cc
                            break
                    k = listunvisited[[kw][0]]
                    pathtable[i, j] = k
                    # print "k=",k
                    # print "len(listunvisited)=",len(listunvisited)
                    # print "len(weight)=",len(weight)

                    w_weight = w_weight + weight[k]
                    unvisited.remove(k)
                    # visited.add(k)

                    length[i] += distmat[visiting][k]

                    visiting = k
                    j += 1
                    # print "w_weight=",w_weight
                    # print "============="
                else:
                    length[i] += distmat[visiting][0]
                    w_weight = 0
                    visiting = 0
                    pathtable[i, j] = 0

                    # print "**********"
                    j += 1
                # print "length[",i,"]=",length[i]
                # print "pathtable[",i,",",j-1,"]=",pathtable[i,j-1]
            length[i] += distmat[visiting][0]
            pathtable[i, j] = 0  # 添加结束点0

        # print length
        # 包含所有蚂蚁的一个迭代结束后，统计本次迭代的若干统计参数
        lengthaver[iter] = length.mean()

        if iter == 0:  # 寻找总路程最短
            lengthbest[iter] = length.min()
            pathbest[iter] = pathtable[length.argmin()].copy()  # argmin 最小值下标
        else:
            if length.min() > lengthbest[iter - 1]:
                lengthbest[iter] = lengthbest[iter - 1]
                pathbest[iter] = pathbest[iter - 1].copy()

            else:
                lengthbest[iter] = length.min()
                pathbest[iter] = pathtable[length.argmin()].copy()

                # 更新信息素
        changepheromonetable = np.zeros((numcity, numcity))
        for i in range(numant):
            # nnnum=0
            for j in range(numcity * 2 - 2):
                if ((pathtable[i, j] != -1) and (pathtable[i, j] != -1)):
                    # if(distmat[pathtable[i,j]][pathtable[i,j+1]] == 0):
                    # continue

                    changepheromonetable[pathtable[i, j]][pathtable[i, j + 1]] += Q / distmat[pathtable[i, j]][
                        pathtable[i, j + 1]]
                    # nnnum += 1

            # changepheromonetable[pathtable[i,nnnum]][pathtable[i,0]] += Q/distmat[pathtable[i,nnnum]][pathtable[i,0]]

        pheromonetable = (1 - rho) * pheromonetable + changepheromonetable

        iter += 1  # 迭代次数指示器+1
    # print pathbest
    print "outyiqun"
    return pathbest


def picture_d(zhongjian, drflag, zonglu, dr_connect, jiaohui,point):
    connect = np.array(dr_connect)
    fflag = 0
    cc = np.array(zhongjian)

    zz = np.array(zonglu)
    # print "cc"
    # print cc

    numcity = len(cc)

    numcity2 = len(zz)
    bestpath = []
    bestpath1 = []
    bestpath2 = []

    for k in range(len(drflag)):  # 对轻件点画图
        aa = np.array(drflag[k])
        # print "aa"
        # print aa
        numcity1 = len(aa)
        if (numcity1 == 0):
            break
        '''
        for kk in aa:  # 打印无人车路途中回车装货
            for ll in jiaohui:
                if (kk[0] == ll):
                    plt.plot(kk[0], kk[1], 'bs')
                    # print "jiaohui***",kk[0]
                    break
        '''
        for j in range(numcity1):
            bestpath1.append(j)
        plt.plot(aa[:, 0], aa[:, 1], 'r.')

        for j in range(numcity1 - 1):

            m, n = bestpath1[j], bestpath1[j + 1]
            #  for k in range(len(drflag)):
            #        if(aa[m][0]==)

            for xx in range(len(connect) - 1):
                if (aa[m][0] == connect[xx] and aa[n][0] == connect[xx + 1] and aa[n][0] != 0 and aa[m][0] != 0):
                    # print "connect=",connect[xx]
                    fflag = 1
                    break

            if (fflag == 1):
                fflag = 0
                continue

            plt.plot([aa[m][0], aa[n][0]], [aa[m][1], aa[n][1]], 'k--',linewidth=1.0)

    for i in range(numcity):
        bestpath.append(i)

    plt.plot(cc[:, 0], cc[:, 1], 'bs')
    plt.plot(point[0], point[1], 'ws')
    plt.plot(point[0],point[1],'g*')
    plt.xlim([-150, 150])
    plt.ylim([-150, 150])
    for i in range(numcity - 1):
        m, n = bestpath[i], bestpath[i + 1]
        plt.plot([cc[m][0], cc[n][0]], [cc[m][1], cc[n][1]], 'k-',linewidth=1.0)

    #plt.plot([cc[bestpath[0]][0], cc[n][0]], [cc[bestpath[0]][1], cc[n][1]], 'k-',linewidth=1.0)

    ax = plt.gca()
    ax.set_title("vrp4")
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y_axis')

    plt.savefig('vrp4.png', dpi=500, bbox_inches='tight')
    plt.close()

def dongtai(mostwei_arr,jilu,coordinates,speed_car,weight,maxcarweight,num_drone,drone_weight):
    print "**********start dongtai***********"    
    time=60
    hight_weight=np.array(mostwei_arr)
    point=np.random.uniform(-100, 100, size=2)
    point_weight = np.random.uniform(0, 4)
    print "point_weight",point_weight
    jilu2=[]
    new_weight=[]
    road=[]
    
    flag=0#判断是否属于该段路的标志
    l_min = 100000
    len_hight=hight_weight.shape[0]
    for i in range(len_hight):
        total2_weight=0
        total_time=0 
        flag=0
        for j in range (len(mostwei_arr[i])-2):
            length=np.linalg.norm(mostwei_arr[i][j+1]-mostwei_arr[i][j])#计算路程中每点距离，算时间，是否达到动态快递出现时间
            t=length/speed_car
            total_time=total_time + t
            if(total_time > time):
                #for k in range (???? , len(jilu)): #目前无法确定是否加入动态点后是否超重
                    
                for k in range (j+1 , len(mostwei_arr[i])-1):
                    llength=np.linalg.norm(mostwei_arr[i][k]-point)#找最近的重件点进行归属
                    if(l_min > llength):
                        start_point=mostwei_arr[i][j+1]
                        l_min=llength
                        road_num=i
                        flag=1
            if(flag == 1):
                break
    coordinates1=coordinates.tolist()
    print coordinates1
    print start_point.tolist()
    indx1=coordinates1.index(start_point.tolist())#开始判断有动态点生成时初始点在coordinates中下标
    for i in range (len(jilu[road_num])):#求出jilu中的下标，并将剩下点加入新路径
        if(indx1==jilu[road_num][i][0]):
            j=i
            while(jilu[road_num][j][0] != -1):
                jilu2.append(jilu[road_num][j][0])
                road.append(coordinates[jilu[road_num][j][0]])
                new_weight.append(weight[jilu[road_num][j][0]])
                j=j+1
            road.insert(len(road)-2,point)#新生成的点加入末尾
            new_weight.insert(len(new_weight)-2,point_weight)
            break
    road1=road[:len(road)-1] #将终点切去，带入蚁群
    new_weight1=new_weight[:len(new_weight)-1]
    road1=np.array(road1)
    road=np.array(road)
    print "road1=",road1
    print "new_weight1=",new_weight1
    #road=np.array(road)
    #jilu2=np.array(jilu2)
    #new_weight=np.array(new_weight)
    #print jilu2
    #print road
    #print new_weight
    print "road=",road
    print "new_weight=",new_weight
    numcity=len(road1)
    new_weightflag = [-1 for x in range(0, numcity+1)] #重件点与轻件点标志
    distmat = cars_drones.getdistmat(road1) #城市的距离矩阵
    pathbest = Ant_algorithm(distmat,new_weight1,numcity,maxcarweight)
    path1 = cars_drones.print_path(pathbest,numcity)
    print "path1**=",path1    
    path1=path1.tolist()
    print "path1++=",path1 
    #path1.insert(len(path1)-2,len(road1))
    #path1=path1[:len(path1)-1]
    path1[0][len(path1[0])-1]=len(road1)#将新路径的回归初始点替换为最初仓库
    path1=np.array(path1)
    print "path1=",path1
    numcity = numcity + 1
    mostwei_arr1,new_weightflag,jilu3,zonglu2,zhongwe=cars_drones.outqingjian(path1,new_weight,numcity,road,new_weightflag,drone_weight)#mostwei_arr为分割好的每段重件点，
        #if (num_drone != 0):
    print "new_weightflag=",new_weightflag   
    print "jilu3=",jilu3
    drflag2,dr_connect2,jiaohui2,first_drone=cars_drones.wurenji(jilu3,road,new_weight,num_drone,drone_weight)
    picture_d(zhongwe, drflag2, zonglu2, dr_connect2, jiaohui2,point)
