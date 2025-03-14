import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# 1. 从文件中读入实验数据，存放在data:{thread_num: [speedup1, speedup2...]}
data = {}
threads_num = range(2, 17, 2)
for i in threads_num:
    filename = "Result_T" + str(i) + ".log"
    with open(filename, "r") as fp:
        lines = fp.readlines()

        # 以空行隔开各个数据集
        speedups = []
        for line in lines:
            if len(line) > 7 and line[:7] == "Speedup":
                speedups.append(eval(line[9:]))
        
        data[i] = speedups[::]


# 2. 绘制一幅折线图
ax = plt.subplot()

# 处理刻度线
ax.minorticks_on()
ax.tick_params(axis='y', which='major', width = 1, length=5)
ax.tick_params(axis='y', which='minor', width = 1, length=3)
ax.tick_params(axis='x', which = 'minor', length=0)
ax.tick_params(axis='x', which = 'major', direction='in')

# 标题
plt.title("The Speedup rate on different Batch Size and Threads Number")
plt.xlabel("The Number of Threads")
plt.ylabel("The Speedup Rate")
batch_size = ['1K', '5K', '10K', '100K']

# 绘图
for i in range(len(batch_size)):
    speedup_data = [speedup[i] for speedup in data.values()]
    ax.plot(threads_num, speedup_data, label = batch_size[i])
    ax.scatter(threads_num, speedup_data)
    ax.legend()
plt.savefig("img/LineGraph.png")
plt.close("all")

# 3. 画出一幅三维曲面图
fig = plt.figure()
ax = fig.add_axes(Axes3D(fig))

# Y为数据量对应的索引值
X = threads_num
Y = range(0, 4, 1)
# 计算3维曲面分格线坐标
X,Y = np.meshgrid(X,Y)
# Z值为speedup_data
def f(X,Y):
    Z = []
    for i in range(len(X)):
        row = []
        for j in range(len(X[i])):
            row.append(data[X[i][j]][Y[i][j]])
        Z.append(np.array(row))
    return np.array(Z)



ax.set_xlabel("The Number of Threads")
ax.set_ylabel("The Batch Size")
ax.set_zlabel("The Speedup Rate")
ax.set_yticks(range(0, 4, 1), batch_size)

ax.plot_surface(X, Y, f(X,Y),rstride=1,cstride=1,cmap=plt.cm.rainbow, alpha=0.9)
plt.savefig('img/Surface.png')