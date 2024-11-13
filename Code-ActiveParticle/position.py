from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import os
import drawwithplt as Dplt


##第一步、加载文件，获取文件路径以及标签
train_path="C:/Users/jing/Desktop/position"
allpath=[]
lllables=[]
def get_lableandwav(path,dir):
    dirs = os.listdir(path)
    for a in dirs:
        print(a)
        print(os.path.isfile(path+"/"+a))
        if os.path.isfile(path+"/"+a):
            allpath.append(dirs)
            if dir!="":
                lllables.append(dir)
        else:
            get_lableandwav(str(path)+"/"+str(a),a)
         ##循环遍历这个文件夹
 
    return allpath,lllables
##第一步、加载文件，获取文件路径以及标签
[allpath,lllables]=get_lableandwav(train_path,"")
print(allpath)
print("----------")
print(lllables)
for path in allpath[0]:
    print(path)
    f=open('C:/Users/jing/Desktop/position/%s'%path)
    sentimentlist = []
    for line in f:
        s = line.strip().split()
        s = np.array(s).astype(float)
        sentimentlist.append(s)
    f.close()
    sentimentlist = np.array(sentimentlist)
    print(sentimentlist.shape)
    # df_train=pd.DataFrame(sentimentlist,columns=['s_no','deal_code','text']
    fig = plt.figure()
    ax = fig.subplots()
    ax.set_xlabel('X') 
    ax.set_xlim(0, 40) 
    ax.set_ylabel('Y') 
    ax.set_ylim(0, 120) 
    ax.scatter(sentimentlist[ :, 0], sentimentlist[ :, 1])

    Dplt.SaveFig(1,'%s.png'%path[2:-4])

# 当elevation=0时，视角为沿x1负方向看，当elevation=90时，视角沿x3负方向看。
# 当azimuth=0时，视角为沿x1负方向看，当azimuth=90时，视角沿x2负方向看。
# 随着azimuth的增加，从x3负方向看，x1x2平面是顺时针旋转的。
# 逆时针旋转，能把x1,x2的大小顺序调整为常规平面坐标系。