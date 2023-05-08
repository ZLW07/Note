## PQP代码分析

两种模型介绍：

<font color=blue size=4>**RSS模型：**</font>

![](.\image\RSS模型.png)



<font color=blue size=4>**OBB模型：**</font>

![](.\image\OBB模型.bmp)

------

$\color{blue}{Tri}$类分析：该类是存储单个三角面片信息，报告三角面片编号和三角形的顶点信息

其成员变量为：

1. 三角形的三个顶点（都是double类型）： P1[3]  P2[3]  P3[3]
2. 当前三角形的编号（为int类型）： id

$\color{blue}{BV}$类分析：

其成员变量有：

1. RSS模型或者是OBB模型的方向向量： $R[3][3]$

   对于RSS模型而言：

2. 模型的中心点位置：$Tr[3]$
3. 模型的边长：$l[2]$
4. 模型半径：$r$

​    对于构建OBB模型：

5. OBB模型的中心点：$To[3]$
6. OBB模型边的长度（一半的长度）： $d[3]$

$\color{blue}{PQP{\_Model}}$类分析：

其成员变量有：

有关于Tris类的：

1. Tri类型的指针： tris
2. 记录有多少个tris： num_tris
3. 记录给tris分配了个内存：num_tris_alloced 

有关于BV类的：

   

