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

------

$\color{blue}{BV}$类分析：

其成员变量有：

1. RSS模型或者是OBB模型的方向向量： $R[3][3]$

1. 第一个子编号： first_child（其中正值表示BV树的第一个）


对于RSS模型而言：

1. 模型的中心点位置：$Tr[3]$
2. 模型的边长：$l[2]$
3. 模型半径：$r$

对于构建OBB模型：

4. OBB模型的中心点：$To[3]$

5. OBB模型边的长度（一半的长度）： $d[3]$

针对RSS模型和OBB模型共有的变量：

​	

其成员函数有：

```C++
//判断是否还有子项  0:没有  1: 有
int Leaf()
{
    return first_child < 0;
}
```

```C++
// 获取模型的大小
inline PQP_REAL BV::GetSize()
{
#if PQP_BV_TYPE & RSS_TYPE // 对于RSS模型而言是获取对角的长度
    return (sqrt(l[0] * l[0] + l[1] * l[1]) + 2 * r);
#else //对于OBB是获取对角顶点一半的平方长度
    return (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
#endif
}
```

```C++
//判断是否两个模型是否相互之间的位置状态
int BV_Overlap(PQP_REAL R[3][3], PQP_REAL T[3], BV *b1, BV *b2)
{
#if PQP_BV_TYPE & OBB_TYPE
    return (obb_disjoint(R, T, b1->d, b2->d) == 0);//计算两个OBB模型的间距；返回值为0表示未发生碰撞
#else //计算两个RSS模型的间距；返回值为0表示未发生碰撞
    PQP_REAL dist = RectDist(R, T, b1->l, b2->l);
    if (dist <= (b1->r + b2->r))
        return 1;
    return 0;
#endif
}
```

```C++
#if PQP_BV_TYPE & RSS_TYPE
//获取两个RSS模型的实际距离，当二者相交时，距离为0
PQP_REAL BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], BV *b1, BV *b2)
{
    PQP_REAL dist = RectDist(R, T, b1->l, b2->l);
    dist -= (b1->r + b2->r);
    return (dist < (PQP_REAL)0.0) ? (PQP_REAL)0.0 : dist;
}
#endif
```



------

$\color{blue}{PQP{\_Model}}$类分析：

其成员变量有：

有关于Tris类的：

1. Tri类型的指针： tris
2. 记录有多少个tris： num_tris
3. 记录给tris分配了个内存：num_tris_alloced 
4. 记录最后一个Tri数据的指针：last_tri

有关于BV类的：

 	1. BV类型指针：b
 	2. 记录有多少个BV模型：num_bvs
 	3. 记录给BV分配了多少个类型的BV： num_bvs_alloced

其他成员变量：

 1. 模型建立的状态标志： build_state

    ```C++
    enum BUILD_STATE
    {
        PQP_BUILD_STATE_EMPTY,    // empty state, immediately after constructor
        PQP_BUILD_STATE_BEGUN,    // after BeginModel(), state for adding triangles
        PQP_BUILD_STATE_PROCESSED // after tree has been built, ready to use
    };
    ```

    

其成员函数：

1.  获取指定的BV模型的地址

```C++
// 获取第n个BV类的地址
BV *child(int n)
{
    return &b[n];
}
```

2. 开始构建模型

```C++
int PQP_Model::BeginModel(int n) //默认参数为n=8
{
    // reset to initial state if necessary
    if (build_state != PQP_BUILD_STATE_EMPTY)
    {
        delete[] b;
        delete[] tris;
        num_tris = num_bvs = num_tris_alloced = num_bvs_alloced = 0;
    }

    // prepare model for addition of triangles
    if (n <= 0)
        n = 8;
    num_tris_alloced = n;
    tris = new Tri[n];
    if (!tris)
    {
        ZLOG_ERR << "PQP Error!  Out of memory for tri array on "
                    "BeginModel() call!";
        return PQP_ERR_MODEL_OUT_OF_MEMORY;
    }

    // give a warning if called out of sequence

    if (build_state != PQP_BUILD_STATE_EMPTY)
    {
        ZLOG_ERR << "PQP Warning! Called BeginModel() on a PQP_Model that was not empty. This model was cleared and "
                    "previous triangle additions were lost.";
        build_state = PQP_BUILD_STATE_BEGUN;
        return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
    }

    build_state = PQP_BUILD_STATE_BEGUN;
    return PQP_OK;
}
```

3. 添加三角面片

```C++
int PQP_Model::AddTri(const PQP_REAL *p1, const PQP_REAL *p2, const PQP_REAL *p3, int id)
{
    if (build_state == PQP_BUILD_STATE_EMPTY)
    {
        BeginModel();
    }
    else if (build_state == PQP_BUILD_STATE_PROCESSED)
    {
        ZLOG_ERR << "PQP Warning! Called AddTri() on PQP_Model object that was already ended. AddTri function was "
                    "ignored.  Must "
                    "do a BeginModel() to clear. The model for addition of new triangles";
        return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
    }

    // allocate for new triangles
    if (num_tris >= num_tris_alloced)
    {
        Tri *temp = new Tri[num_tris_alloced * 2];
        memcpy(temp, tris, sizeof(Tri) * num_tris);
        delete[] tris;
        tris = temp;
        num_tris_alloced = num_tris_alloced * 2;
    }

    // initialize the new triangle

    tris[num_tris].p1[0] = p1[0];
    tris[num_tris].p1[1] = p1[1];
    tris[num_tris].p1[2] = p1[2];

    tris[num_tris].p2[0] = p2[0];
    tris[num_tris].p2[1] = p2[1];
    tris[num_tris].p2[2] = p2[2];

    tris[num_tris].p3[0] = p3[0];
    tris[num_tris].p3[1] = p3[1];
    tris[num_tris].p3[2] = p3[2];

    tris[num_tris].id = id;

    num_tris += 1;

    return PQP_OK;
}
```

4. 模型构建,包括将多分配的存储三角形数据的内存去除掉。然后根据存储的三角形数据来构建BV模型。

   构建BV模型的函数为：build_model

```C++
int PQP_Model::EndModel()
{
    if (build_state == PQP_BUILD_STATE_PROCESSED)
    {
        ZLOG_ERR << "PQP Warning! Called EndModel() on PQP_Model object that was already ended. EndModel() was "
                    "ignored.  Must do a BeginModel() to clear the model for addition of new triangles";
        return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
    }

    // report error is no tris
    if (num_tris == 0)
    {
        ZLOG_ERR << "PQP Error! EndModel() called on model with"
                    " no triangles";
        return PQP_ERR_BUILD_EMPTY_MODEL;
    }

    // shrink fit tris array
    if (num_tris_alloced > num_tris)
    {
        Tri *new_tris = new Tri[num_tris];
        memcpy(new_tris, tris, sizeof(Tri) * num_tris);
        delete[] tris;
        tris = new_tris;
        num_tris_alloced = num_tris;
    }

    // create an array of BVs for the model
    b = new BV[2 * num_tris - 1];
    if (!b)
    {
        ZLOG_ERR << "PQP Error! out of memory for BV array "
                    "in EndModel()";
        return PQP_ERR_MODEL_OUT_OF_MEMORY;
    }
    num_bvs_alloced = 2 * num_tris - 1;
    num_bvs = 0;

    // we should build the model now.
    build_model(this);
    build_state = PQP_BUILD_STATE_PROCESSED;
    last_tri = tris;
    return PQP_OK;
}
```

构建BV树模型的函数：

1. 求取协方差函数：

函数的依据的数学公式：
$$
平均值：{\rm{\bar X  =  }}\frac{{\sum\limits_{i = 1}^n {{X_i}} }}{n}\\
标准差：S = \sqrt {\frac{{\sum\limits_{i = 1}^n {{{({X_i} - {\rm{\bar X}})}^2}} }}{{n - 1}}}  \\
方差：{S^2} = \frac{{\sum\limits_{i = 1}^n {{{({X_i} - {\rm{\bar X}})}^2}} }}{{n - 1}} \\
协方差:\begin{array}{l}
{\mathop{\rm var}} (X) = \frac{{\sum\limits_{i = 1}^n {({X_i} - {\rm{\bar X}})({X_i} - {\rm{\bar X}})} }}{{n - 1}}\\
{\mathop{\rm cov}} (X,Y) = \frac{{\sum\limits_{i = 1}^n {({X_i} - {\rm{\bar X}})({Y_i} - {\rm{\bar Y}})} }}{{n - 1}}
\end{array} \\
\begin{array}{l} \\
且有:{\mathop{\rm cov}} (X,X) = {\mathop{\rm var}} (X) \\
{~~~~~~~~~~}{\mathop{\rm cov}} (X,Y) = {\mathop{\rm cov}} (Y,X)
\end{array}
$$

协方差矩阵：

数据为：${x_i} = \left[ {\begin{array}{*{20}{c}}
{{a_{i1}}}&{{a_{i2}}}& \cdots &{{a_{in}}}
\end{array}} \right]$
$$
{X_{m \times n}} = \left[ {\begin{array}{*{20}{c}}
{{a_{11}}}&{{a_{12}}}& \cdots &{{a_{1n}}}\\
{{a_{21}}}&{{a_{22}}}& \cdots &{{a_{2n}}}\\
 \vdots & \vdots & \vdots & \vdots \\
{{a_{m1}}}&{{a_{m2}}}&{}&{{a_{mn}}}
\end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
{{X_1}}&{{X_2}}& \cdots &{{X_n}}
\end{array}} \right]\\
{\mathop{\rm cov}} {[X_{m \times n}]} = \left[ {\begin{array}{*{20}{c}}
{{\mathop{\rm cov}} ({X_1},{X_1})}& \cdots &{{\mathop{\rm cov}} ({X_1},{X_n})}\\
 \vdots & \ddots & \vdots \\
{{\mathop{\rm cov}} ({X_n},{X_1})}& \cdots &{{\mathop{\rm cov}} ({X_n},{X_n})}
\end{array}} \right]
$$



```C++
void get_covariance_triverts(PQP_REAL M[3][3], Tri *tris, int num_tris)
{
    int i;
    PQP_REAL S1[3];
    PQP_REAL S2[3][3];

    S1[0] = S1[1] = S1[2] = 0.0;
    S2[0][0] = S2[1][0] = S2[2][0] = 0.0;
    S2[0][1] = S2[1][1] = S2[2][1] = 0.0;
    S2[0][2] = S2[1][2] = S2[2][2] = 0.0;

    // get center of mass
    for (i = 0; i < num_tris; i++)
    {
        PQP_REAL *p1 = tris[i].p1;
        PQP_REAL *p2 = tris[i].p2;
        PQP_REAL *p3 = tris[i].p3;

        S1[0] += p1[0] + p2[0] + p3[0];
        S1[1] += p1[1] + p2[1] + p3[1];
        S1[2] += p1[2] + p2[2] + p3[2];

        S2[0][0] += (p1[0] * p1[0] + p2[0] * p2[0] + p3[0] * p3[0]);
        S2[1][1] += (p1[1] * p1[1] + p2[1] * p2[1] + p3[1] * p3[1]);
        S2[2][2] += (p1[2] * p1[2] + p2[2] * p2[2] + p3[2] * p3[2]);
        S2[0][1] += (p1[0] * p1[1] + p2[0] * p2[1] + p3[0] * p3[1]);
        S2[0][2] += (p1[0] * p1[2] + p2[0] * p2[2] + p3[0] * p3[2]);
        S2[1][2] += (p1[1] * p1[2] + p2[1] * p2[2] + p3[1] * p3[2]);
    }

    PQP_REAL n = (PQP_REAL)(3 * num_tris);

    // now get covariances

    M[0][0] = (S2[0][0] - S1[0] * S1[0] / n) / (n - 1);
    M[1][1] = (S2[1][1] - S1[1] * S1[1] / n) / (n - 1);
    M[2][2] = (S2[2][2] - S1[2] * S1[2] / n) / (n - 1);
    M[0][1] = (S2[0][1] - S1[0] * S1[1] / n) / (n - 1);
    M[1][2] = (S2[1][2] - S1[1] * S1[2] / n) / (n - 1);
    M[0][2] = (S2[0][2] - S1[0] * S1[2] / n) / (n - 1);
    M[1][0] = M[0][1];
    M[2][0] = M[0][2];
    M[2][1] = M[1][2];
}
```

2. 求取特征值$dout[3]$和特征向量$vout[3][3]$

```C++
void inline Meigen(PQP_REAL vout[3][3], PQP_REAL dout[3], PQP_REAL a[3][3])
{
    int n = 3;
    int j, iq, ip, i;
    PQP_REAL tresh, theta, tau, t, sm, s, h, g, c;
    int nrot;
    PQP_REAL b[3];
    PQP_REAL z[3];
    PQP_REAL v[3][3];
    PQP_REAL d[3];

    Midentity(v);
    for (ip = 0; ip < n; ip++)
    {
        b[ip] = a[ip][ip];
        d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }

    nrot = 0;

    for (i = 0; i < 50; i++)
    {
       
        sm = 0.0;
        for (ip = 0; ip < n; ip++)
        {
            for (iq = ip + 1; iq < n; iq++)
            {
                //判断上三角矩阵对应值的和是否为0,如下
        		//x  0  0
                //x  x  0
                //x  x  x
                sm += fabs(a[ip][iq]);
            }
        }

        if (sm == 0.0)
        {
            McM(vout, v);//将v赋给vout
            VcV(dout, d);//将d赋给dout
            return;
        }

        if (i < 3)
            tresh = (PQP_REAL)0.2 * sm / (n * n);
        else
            tresh = 0.0;

        for (ip = 0; ip < n; ip++)
            for (iq = ip + 1; iq < n; iq++)
            {
                g = (PQP_REAL)100.0 * fabs(a[ip][iq]);
                if (i > 3 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(d[iq]) + g == fabs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh)
                {
                    h = d[iq] - d[ip];
                    if (fabs(h) + g == fabs(h))
                        t = (a[ip][iq]) / h;
                    else
                    {
                        theta = (PQP_REAL)0.5 * h / (a[ip][iq]);
                        t = (PQP_REAL)(1.0 / (fabs(theta) + sqrt(1.0 + theta * theta)));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = (PQP_REAL)1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / ((PQP_REAL)1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j < ip; j++)
                    {
                        ROTATE(a, j, ip, j, iq);
                    }
                    for (j = ip + 1; j < iq; j++)
                    {
                        ROTATE(a, ip, j, j, iq);
                    }
                    for (j = iq + 1; j < n; j++)
                    {
                        ROTATE(a, ip, j, iq, j);
                    }
                    for (j = 0; j < n; j++)
                    {
                        ROTATE(v, j, ip, j, iq);
                    }
                    nrot++;
                }
            }
        for (ip = 0; ip < n; ip++)
        {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    fprintf(stderr, "eigen: too many iterations in Jacobi transform.\n");

    return;
}
```



```C++
int build_recurse(PQP_Model *m, int bn, int first_tri, int num_tris)
{
    BV *b = m->child(bn); // （开始时，取PQP_Model的BV树中的第一个BV， bn = first_tri = 0）

    // compute a rotation matrix

    PQP_REAL C[3][3], E[3][3], R[3][3], s[3], axis[3], mean[3], coord;

#if RAPID2_FIT
    moment *tri_moment = new moment[num_tris];
    compute_moments(tri_moment, &(m->tris[first_tri]), num_tris);
    accum acc;
    clear_accum(acc);
    for (int i = 0; i < num_tris; i++) accum_moment(acc, tri_moment[i]);
    delete[] tri_moment;
    covariance_from_accum(C, acc);
#else
    get_covariance_triverts(C, &m->tris[first_tri], num_tris);
#endif

    Meigen(E, s, C);

    // place axes of E in order of increasing s

    int min, mid, max;
    if (s[0] > s[1])
    {
        max = 0;
        min = 1;
    }
    else
    {
        min = 0;
        max = 1;
    }
    if (s[2] < s[min])
    {
        mid = min;
        min = 2;
    }
    else if (s[2] > s[max])
    {
        mid = max;
        max = 2;
    }
    else
    {
        mid = 2;
    }
    McolcMcol(R, 0, E, max);
    McolcMcol(R, 1, E, mid);
    R[0][2] = E[1][max] * E[2][mid] - E[1][mid] * E[2][max];
    R[1][2] = E[0][mid] * E[2][max] - E[0][max] * E[2][mid];
    R[2][2] = E[0][max] * E[1][mid] - E[0][mid] * E[1][max];

    // fit the BV

    b->FitToTris(R, &m->tris[first_tri], num_tris);

    if (num_tris == 1)
    {
        // BV is a leaf BV - first_child will index a triangle

        b->first_child = -(first_tri + 1);
    }
    else if (num_tris > 1)
    {
        // BV not a leaf - first_child will index a BV

        b->first_child = m->num_bvs;
        m->num_bvs += 2;

        // choose splitting axis and splitting coord

        McolcV(axis, R, 0);

#if RAPID2_FIT
        mean_from_accum(mean, acc);
#else
        get_centroid_triverts(mean, &m->tris[first_tri], num_tris);
#endif
        coord = VdotV(axis, mean);

        // now split

        int num_first_half = split_tris(&m->tris[first_tri], num_tris, axis, coord);

        // recursively build the children

        build_recurse(m, m->child(bn)->first_child, first_tri, num_first_half);
        build_recurse(m, m->child(bn)->first_child + 1, first_tri + num_first_half, num_tris - num_first_half);
    }
    return PQP_OK;
}
```
