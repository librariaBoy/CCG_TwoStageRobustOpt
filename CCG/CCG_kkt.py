'''两阶段鲁棒优化CC&G算法'''
import numpy as np
import gurobipy as gp
import matplotlib.pyplot as plt

# 下界和上界集合
ub_set = []
lb_set = []

# i 处建造仓库的费用
f = np.array([400,414,326]) 

# i 处仓库的单位容量产生的费用
a = np.array([18,25,20])

# i->j 运输产生的费用
c = np.array([[22,33,24],
              [33,23,30],
              [20,25,27]])

'''CC&G算法模型'''

## 主问题
MP = gp.Model("CC&G主问题")
# i 处仓库是否建造
y = MP.addVars(range(3), vtype=gp.GRB.BINARY, name="y")

# i 处仓库的容量
z = MP.addVars(range(3), lb=0, vtype=gp.GRB.CONTINUOUS, name="z")

# i->j 的运输
x_SP = MP.addVars([0,1,2],[0,1,2], lb=0, vtype=gp.GRB.CONTINUOUS, name="x_SP")

# 第一阶段问题
yita = MP.addVar(vtype=gp.GRB.CONTINUOUS, name="yita")

# 目标函数
MP_Obj = gp.quicksum(f[i] * y[i] + a[i] * z[i] for i in range(3))
MP_Obj = MP_Obj + yita

MP.setObjective(MP_Obj, gp.GRB.MINIMIZE)

# MP 初始约束条件
MP.addConstrs((z[i] <= 800*y[i] for i in range(3)), "z(i) <= 800*y(i)")
MP.addConstr(sum(z[i] for i in range(3)) >= 772, "sum_z >= 772")
init_d = [206,274,220]
# (2)
for i in range(3):
    row_sum_x = sum(x_SP[i,j] for j in range(3))
    MP.addConstr(row_sum_x <= z[i], name=f"sum_x_{i} <= z{i}")
# (3)
for j in range(3):
    col_sum_x = sum(x_SP[i,j] for i in range(3))
    MP.addConstr(col_sum_x >= init_d[j], name=f"sum_x_{j} >= opt_d{j}")


## 子问题
SP = gp.Model("CC&G子问题")

# i->j 的运输
x = SP.addVars([0,1,2], [0,1,2], lb=0, ub=gp.GRB.INFINITY, vtype=gp.GRB.CONTINUOUS, name="x")

# 不确定集
g = SP.addVars([0,1,2], lb=0, ub=1, vtype=gp.GRB.CONTINUOUS, name="g" )
#g = SP.addVars([0,1,2], vtype=gp.GRB.BINARY, name="g" )
d = SP.addVars([0,1,2], lb=0, ub=gp.GRB.INFINITY, vtype=gp.GRB.CONTINUOUS, name="d" )
SP.addConstr(d[0] == 206 + 40 * g[0], name="d0 = 206 + 40*g0")
SP.addConstr(d[1] == 274 + 40 * g[1], name="d1 = 274 + 40*g1")
SP.addConstr(d[2] == 220 + 40 * g[2], name="d2 = 220 + 40*g2")
SP.addConstr(sum(g[i] for i in range(3)) <= 3)
SP.addConstr(g[0] + g[1] + g[2] <= 1.8, name="g0 + g1 + g2 <= 1.8")
SP.addConstr(g[0] + g[1] <= 1.2, name="g0 + g1 <= 1.2")
d_base = np.array([206,274,220])

# 上下界
UB = float("inf")
LB = float("-inf")

# 迭代次数
k = 0;

# 主循环
while (UB - LB > 1e-5):

    print("第",k,"次迭代")
    print("上界：",UB)
    print("下界：",LB)
    print("---------------------")
    print("       MP k          ")
    print("---------------------")

    # 主问题求解
    MP.optimize()

    # 获取 z 的最优值
    z_opt = np.array([z[i].X for i in range(3)])
    print(z_opt)
    # 更新下界
    LB = max(LB,MP.objVal)
    lb_set.append(LB)

    # 子问题根据主问题的解加约束
    if k > 0:
        # 删除之前的子问题约束
        SP.remove(sp1)
        SP.remove(sp2)
        SP.remove(sp3)
        SP.remove(sp4)
        SP.remove(sp5)
        SP.remove(sp6)
        SP.remove(sp7)
        SP.remove(sp8)
        SP.remove(sp9)
        SP.remove(v)
        SP.remove(w)
        SP.remove(pi)
        SP.remove(h)
        SP.update()
   
   # KKT 处理后的约束
    M = 1e5 # 大 M
    pi = SP.addVars([0,1,2], lb=-gp.GRB.INFINITY, ub=0, vtype=gp.GRB.CONTINUOUS, name="pi")
    theta = SP.addVars([0,1,2], lb=-gp.GRB.INFINITY, ub=0, vtype=gp.GRB.CONTINUOUS, name="theta")
    v = SP.addVars([0,1,2], vtype=gp.GRB.BINARY, name="v")
    w = SP.addVars([0,1,2], vtype=gp.GRB.BINARY, name="w")
    h = SP.addVars(range(3), range(3), vtype=gp.GRB.BINARY, name="h")


    sp1 = SP.addConstrs(sum(x[i,j] for j in range(3)) <= z_opt[i] for i in range(3))
    sp2 = SP.addConstrs(sum(x[i,j] for i in range(3)) >= d[j] for j in range(3))
    sp3 = SP.addConstrs(pi[i] - theta[j] <= c[i][j] for i in range(3) for j in range(3))
    sp4 = SP.addConstrs(-pi[i] <= M*v[i] for i in range(3))
    sp5 = SP.addConstrs(z_opt[i] - sum(x[i,j] for j in range(3)) <= M*(1-v[i]) for i in range(3) for j in range(3))
    sp6 = SP.addConstrs(-theta[j] <= M*w[j] for j in range(3))
    sp7 = SP.addConstrs(sum(x[i,j] for i in range(3)) - d[j] <= M*(1-w[j]) for j in range(3))
    sp8 = SP.addConstrs(x[i,j] <= M*h[i,j] for i in range(3) for j in range(3))
    sp9 = SP.addConstrs(c[i][j] - pi[i] + theta[j] <= M*(1-h[i,j]) for i in range(3) for j in range(3))

    SP_Obj = sum(c[i][j]*x[i,j] for i in range(3) for j in range(3))
    
    ## 对偶尝试
    #for i in range(3):
    #    for j in range(3):
    #        SP.addConstr((pi[i] - theta[j] <= c[i][j]), name=f"pi_{i}-theta_{j} <= c_{i}_{j}")

    #SP.addConstr(sum(g[i] for i in range(3)) <= 1)        
    #SP_Obj = sum(z_opt[i]*pi[i] for i in range(3)) - sum(d[j]*theta[j]  for j in range(3))
    
    # 子问题目标函数
    # SP_Obj = sum(c[i,j]*x[i,j] for i in range(3) for j in range(3))


    SP.setObjective(SP_Obj, gp.GRB.MAXIMIZE)

    print("---------------------")
    print("       SP k          ")
    print("---------------------")
    # 子问题求解
    SP.optimize()


    # 迭代次数累加
    k = k + 1
    
    if SP.Status != gp.GRB.INFEASIBLE:
        # 获取 d 的最优值
        # opt_d = np.array([d_base[i] + g_c[i].X/x_dual[i+3].X for i in range(3)])
        opt_d = np.array([d[i].X for i in range(3)])
        y_stage1 = np.array([y[i].X for i in range(3)])
        z_stage1 = np.array([z[i].X for i in range(3)])


    # 增加最优割和可行割
    if SP.status == gp.GRB.OPTIMAL:
        # 更新上界
        UB = min(UB,SP.ObjVal + sum(f[i]*y_stage1[i] + a[i]*z_stage1[i] for i in range(3)))
        ub_set.append(UB)

        # 根据论文公式 (12)、(13) 添加约束
        var_name = f'x{k}'
        new_x = MP.addVars([0,1,2],[0,1,2], lb=0, vtype=gp.GRB.CONTINUOUS, name=var_name)

        # (1)
        sum_cx = 0
        for i in range(3):
            for j in range(3):
                sum_cx += c[i,j]*new_x[i,j]
        MP.addConstr(yita >= sum_cx, name="yita >= sum_cx")

        # (2)
        for i in range(3):
            row_sum_x = sum(new_x[i,j] for j in range(3))
            MP.addConstr(row_sum_x <= z[i], name=f"sum_x_{i} <= z{i}")

        # (3)
        for j in range(3):
            col_sum_x = sum(new_x[i,j] for i in range(3))
            MP.addConstr(col_sum_x >= opt_d[j], name=f"sum_x_{j} >= opt_d{j}")
        
        MP.update()

    elif SP.Status == gp.GRB.UNBOUNDED:
        # 更新上界
        UB = min(UB,SP.ObjVal)
        ub_set.append(UB)

        # 根据论文公式 (14) 添加约束
        var_name = f'x{k}'
        MP.addVars([0,1,2],[0,1,2], lb=0, vtype=gp.GRB.CONTINUOUS, name=var_name)

        # (2)
        for i in range(3):
            row_sum_x = sum(new_x[i,j] for j in range(3))
            MP.addConstr(row_sum_x <= z[i], name=f"sum_x_{i} <= z{i}")

        # (3)
        for j in range(3):
            col_sum_x = sum(new_x[i,j] for i in range(3))
            MP.addConstr(col_sum_x >= opt_d[j], name=f"sum_x_{j} >= opt_d{j}")

    else:
        print("子问题求解失败,原问题无解")
        break
    
    if k > 10:
        break

print("")
print("-------------------------")
print("-------------------------")
print("求解结束")
print("上界：",UB)
print("下界：",LB)
print("迭代次数：",k)


plt.figure()
plt.plot(range(k), lb_set, label="lower bound")
plt.plot(range(k), ub_set, label="upper bound")
plt.xlabel("iteration")
plt.ylabel("value")
plt.legend()
plt.show()