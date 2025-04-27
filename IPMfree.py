import numpy as np

def gf2_k_nultiply(a,b,k,mod_poly):
    result = 0
    for i in range(k):
        if (b >> i) & 1:
            result ^= a << i
    for i in range(2*k-1,k-1,-1):
        if (result >> i) & 1:
            result ^= mod_poly << (i-k)
    return result & (2**k-1)

def Tolist(x):
    return [(x >> i) & 1 for i in range(3)][::-1]

def ToInt(bit_list):
    return sum(bit << i for i, bit in enumerate(bit_list[::-1]))


def compute_C0(x, y, k):  
    C0 = np.zeros((k, k), dtype=int)
    for i in range(k):
        for j in range(k - i):
            C0[i, i + j] = (x >> i)  * (y >> j) & 1
        C0[i] = C0[i][::-1]
    return C0

# def construct_Vp(q, k):
#     V = {}
#     for p in range(1, k):
#         # 矩阵 V^p 的行数为 k - p，列数为 k
#         Vp = np.zeros((k - p, k), dtype=int)
#         for i in range(k - p):
#             # 计算 x^{p-1} * q(x) 的系数
#             for j in range(k):
#                 Vp[i][k-j-1] += (q << p-1) >> j & 1
#         V[p] = Vp
    
#     return V

def compute_x_power_mod(m, q_int, k):
    result = 1  # 初始化为 x^0 = 1
    for _ in range(m):
        result <<= 1        # 乘以 x
        if result & (1 << k):  # 若次数 ≥k，则异或不可约多项式 q(x)
            result ^= q_int
    # 提取系数：索引 0 为 x^0 项，索引 k-1 为 x^{k-1} 项
    return [(result >> i) & 1 for i in range(k)]

def construct_Vp(q_int, k):
    V = {}
    for p in range(1, k):
        # 计算 x^{k+p-1} mod q(x) 的系数
        coeffs_low_to_high = compute_x_power_mod(k + p - 1, q_int, k)
        # 转换为高次到低次排列以匹配矩阵列顺序
        coeffs_high_to_low = coeffs_low_to_high[::-1]
        # 构造 V_p 矩阵 (行向量)
        V[p] = np.array([coeffs_high_to_low], dtype=int)
    return V

def compute_Cp(x, y, k, V):
    CP = {}
    for p in range(1, k):
        cp = np.zeros((0, k), dtype=int)
        for i in range(p, k):
            j = (k - 1 - i) + p
            # 计算值并填充第 i 行
            value = (x >> i & 1) * (y >> j & 1)
            m = np.full((1, k), value, dtype=int)
            cp = np.vstack((cp, m))
        cp = cp & V[p]
        CP[p] = cp
    return CP

def ipm_free(n_shares, k_bits, L , x, qx):
    # 生成n个随机共享
    shares = [np.random.randint(0, 2**k_bits) for _ in range(n_shares)]
    shares[0] = x
    sum_x = 0
    for i in range(1, n_shares):
        sum_x += gf2_k_nultiply(shares[i] ,L[i],k_bits,qx)
    shares[0] ^= sum_x
    return shares

def ipm_t(t, n_shares, k_bits, L, x, qx):
    overall_xor = [0] * n_shares  
    # 生成 t 个 ipm_free 异或
    for _ in range(t):
        shares = ipm_free(n_shares, k_bits, L, x, qx)
        overall_xor = [a ^ b for a, b in zip(overall_xor, shares)]
    return overall_xor

def gadget1(k,n,qx,L,t):
    # 生成Vp矩阵
    V = construct_Vp(qx, k)
    # print(f'V: {V}')
    s_init = np.array([ipm_t(t, n, k, L, 0, qx) for _ in range(n*k*(k+1)//2)])
    # print(f's_init: {s_init}')
    r_init = np.array([ipm_t(t, n, k, L, 0, qx) for _ in range(k*(k+1)//2)])
    # print(f'r_init: {r_init}')
    S = np.zeros((n, n, k*(k+1)//2), dtype=int)
    for i in range(n):
        for j in range(n):
            for l in range(k*(k+1)//2):
                S[i][j][l] = s_init[i*k*(k+1)//2  + l][j]
    # print(f'S: {S}')

    R = np.zeros((n, k*(k+1)//2), dtype=int)
    for i in range(n):
        for j in range(k*(k+1)//2):
            R[i][j] = r_init[j][i]
    print(f'R: {R}')
    for i in range(n):
        S[i][0] ^= R[i]
    # print(f"S_update:{S}")
    return S, V


def gadget2(x_i, y_j, S_ij, V, k):
    C0 = compute_C0(x_i, y_j, k)
    CP = compute_Cp(x_i, y_j, k, V)
    C = np.vstack((C0, *CP.values()))
    # print('generate C --------------------------')
    # print(f'C:\n{C}')
    S_ij_matrix = np.array([[elem >> i & 1 for i in range(k-1,-1,-1)]for elem in S_ij])
    # print(f'S_ij_matrix:\n{S_ij_matrix}')
    C_refresh = C ^ S_ij_matrix
    # print(f'C_refresh:\n{C_refresh}')
    z = np.zeros(k, dtype=int)
    for i in range(k*(k+1)//2):
        z ^= C_refresh[i]
    print(f'z_ij: {z}')
    return z

def gadget3(x_musk,y_musk, S, V, k, n, L , qx):
    z = np.zeros((n, n), dtype=int)
    Z = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(n):
            z[i][j] = ToInt(gadget2(x_musk[i], y_musk[j], S[i][j], V, k))
    for i in range(n):
        Z[i] = 0
        for j in range(n):
            Z[i] ^= gf2_k_nultiply(z[i][j],L[j],k,qx)
    return Z

def verify(z_musk, k, L , x,y ,n,qx):
    ans = sum([gf2_k_nultiply(z_musk[i],L[i],k,qx) for i in range(1, n)])
    ans ^= z_musk[0]
    print(f"z是：{ans}")
    return ans == gf2_k_nultiply(x,y,k,qx)

if __name__ == '__main__':
    x = 0b11011
    y = 0b10001
    k = 5 #k为多项式的次数
    n = 2
    t = 3
    L = [1, 6]
    qx = 0b101110 

    x_musk = ipm_free(n, k, L, x, qx)
    print(f"x是：{x},x的掩码是：{x_musk}")  # [7, 3]
    y_musk = ipm_free(n, k, L, y , qx)
    print(f"y是：{y},y的掩码是：{y_musk}")  # [7, 3]

    S, V = gadget1(k,n,qx,L,t)
    print(f'S: {S}')
    print(f'V: {V}')
    z_musk = gadget3(x_musk, y_musk, S, V, k, n, L, qx)
    print(f'z_musk: {z_musk}')
    print(f'verify: {verify(z_musk, k, L , x,y ,n,qx)}')

