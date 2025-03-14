# 多线程并行快速排序算法

## 项目简介
本项目使用OpenMP编程方法实现多线程并行快速排序。测试数据集为随机生成的int型数组，数据量包括1K、5K、10K、100K。  

## 使用方法  
提供两个shell脚本文件 `run.sh` 与 `draw.sh`，分别用于计算各数据集上的算法性能以及将结果可视化，二者都需要在 `QuickSort_OMP` 工作目录下执行。 
1. `run.sh`  
读取 `~/dataset/*.txt` 中的各数据集，文件名的后缀表示 **数据量大小**。分别运行并行快排和串行快排程序，测算时间并计算加速比，结果写入 `~/Results/*.log` 中各对应文件下，文件名的后缀表示 **并行执行的线程数量**。  

2. `draw.sh`  
从上一步得到的各log文件中读取测试结果并绘图，图像存放在 `~/Results/img` 下，包括二维的折线图与三维的曲面图。  