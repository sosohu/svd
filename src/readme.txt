各版本名称以及改变
1 gesvd: 最初的串行单边jacobi算法
2 newgesvd: 使用了贪心的最大值匹配算法的串行单边jacobi算法
3 odd_even_gesvd:使用奇偶序列方式把版本1并行化
4 odd_even_gesvd_optim:对版本3每次迭代增加一个阀值以减少运算量
5 odd_even_gesvd_opt_a:对版本4每次取指作出改变，使用128位取指指令  另外改变数组方式，使得按列数据相邻
