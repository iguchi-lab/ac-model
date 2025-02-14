import ac_model.main as ac
from math import floor

C_af_H = 0.8
C_df_H = lambda Theta, h: 0.77 if (Theta < 5.0) and (h >= 80.0) else 1.0
C_af_C = 0.85
C_hm_C = 1.15
SHF_L_min_c = 0.4

table_3 = [
    (-0.00236,  0.01324,  0.08418, -0.47143, -1.16944,   6.54886),
    ( 0.00427, -0.02392, -0.19226,  0.94213,  2.58632, -12.85618),
    (-0.00275,  0.01542,  0.14947, -0.68303, -2.03594,  10.60561),
    ( 0.00063, -0.00351, -0.02865,  0.10522,  0.37336,  -1.09499),
    (-0.00005,  0.00028,  0.00184, -0.01090, -0.09609,   0.59229)
]
table_4_A = [
    (-0.000056,  0.000786,  0.071625),
    (-0.000145,  0.003337, -0.143643),
    (-0.000240, -0.029471,  1.954343),
    (-0.000035, -0.050909,  1.389751),
    ( 0.0,       0.0,       0.076800)
]
table_4_B = [
    ( 0.000108, -0.035658,  3.063873),
    (-0.000017,  0.062546, -5.471556),
    (-0.000245, -0.025126,  4.057590),
    ( 0.000323, -0.021166,  0.575459),
    ( 0.0,       0.000330,  0.047500)
]
table_4_C = [
    (-0.001465, -0.030500,  1.920431),
    ( 0.002824,  0.041081, -1.835302),
    (-0.001929, -0.009738,  1.582898),
    ( 0.000616, -0.014239,  0.546204),
    ( 0.0,      -0.000110,  0.023100)
]

table_5 = [
    (0.00000, 0.00000,  0.00000,  0.00000,  0.00000,  0.00000),
    (0.00000, 0.00000, -0.00036,  0.05080, -0.20346,  0.47765),
    (0.00000, 0.00000,  0.00227, -0.03952,  0.04115,  0.23099),
    (0.00000, 0.00000, -0.00911,  0.07102,  0.14950, -1.07335),
    (0.00000, 0.00000,  0.00044, -0.00214, -0.06250,  0.35150)
]
table_6_A = [
    (-0.0004078,  0.01035, -0.03248),
    ( 0.0,        0.04099, -0.818809),
    ( 0.0,       -0.04615,  2.10666),
    ( 0.0013382, -0.01179, -0.41778),
    ( 0.0000000, -0.00102,  0.09270)
]
table_6_B = [
    (-0.000056, -0.003539, -0.430566),
    ( 0.0,       0.015237,  1.188850),
    ( 0.0,       0.000527, -0.304645),
    (-0.000179,  0.020543,  0.130373),
    ( 0.0,       0.000240,  0.013500)
]
table_6_C = [
    (-0.0001598,  0.004848,  0.047097),
    ( 0.0,        0.016675,  0.362141),
    ( 0.0,       -0.008134, -0.023535),
    (-0.0000772,  0.012558,  0.056185),
    ( 0.0,       -0.000110,  0.010300)
]

def calc_Q_r_max_H(q_rtd_C, q_r_max_H, Theta_ex_d_t):
    a2, a1, a0 = calc_a_eq3(q_r_max_H, q_rtd_C)
    return a2 * (Theta_ex_d_t - 7) ** 2 + a1 * (Theta_ex_d_t - 7) + a0

def calc_a_eq3(q_r_max_H, q_rtd_C):
    b2 =  0.000181 * q_rtd_C * 10 ** (-3) - 0.000184
    b1 =  0.002322 * q_rtd_C * 10 ** (-3) + 0.013904
    b0 =  0.003556 * q_rtd_C * 10 ** (-3) + 0.993431
    c2 = -0.000173 * q_rtd_C * 10 ** (-3) + 0.000367
    c1 = -0.003980 * q_rtd_C * 10 ** (-3) + 0.003983
    c0 = -0.002870 * q_rtd_C * 10 ** (-3) + 0.006376
    a2 = b2 * q_r_max_H + c2
    a1 = b1 * q_r_max_H + c1
    a0 = b0 * q_r_max_H + c0
    return a2, a1, a0

def calc_f_H_Theta(x, q_rtd_C, Theta_ex, dualcompressor=False):
    a0, a1, a2, a3, a4 = calc_a_eq7(q_rtd_C, dualcompressor, Theta_ex)
    return a4 * x ** 4 + a3 * x ** 3 + a2 * x ** 2 + a1 * x + a0

def calc_a_eq7(q_rtd_C, dualcompressor, Theta_ex):
    if dualcompressor == False:     calc_p_i = calc_p_i_eq8     # dualcompressorでない
    else:                           calc_p_i = calc_p_i_eq9     # dualcompressor

    p_42 = calc_p_i(42, q_rtd_C)
    p_41 = calc_p_i(41, q_rtd_C)
    p_40 = calc_p_i(40, q_rtd_C)
    p_32 = calc_p_i(32, q_rtd_C)
    p_31 = calc_p_i(31, q_rtd_C)
    p_30 = calc_p_i(30, q_rtd_C)
    p_22 = calc_p_i(22, q_rtd_C)
    p_21 = calc_p_i(21, q_rtd_C)
    p_20 = calc_p_i(20, q_rtd_C)
    p_12 = calc_p_i(12, q_rtd_C)
    p_11 = calc_p_i(11, q_rtd_C)
    p_10 = calc_p_i(10, q_rtd_C)
    p_02 = calc_p_i(2,  q_rtd_C)
    p_01 = calc_p_i(1,  q_rtd_C)
    p_00 = calc_p_i(0,  q_rtd_C)

    a4 = p_42 * Theta_ex ** 2 + p_41 * Theta_ex + p_40 * 1
    a3 = p_32 * Theta_ex ** 2 + p_31 * Theta_ex + p_30 * 1
    a2 = p_22 * Theta_ex ** 2 + p_21 * Theta_ex + p_20 * 1
    a1 = p_12 * Theta_ex ** 2 + p_11 * Theta_ex + p_10 * 1
    a0 = p_02 * Theta_ex ** 2 + p_01 * Theta_ex + p_00 * 1

    return a0, a1, a2, a3, a4

# (容量可変型コンプレッサー搭載ルームエアコンディショナーでないルームエアコンディショナー)
def calc_p_i_eq8(i, q_rtd_C):
    s_i = table_3[4 - floor(i / 10)][(2 - (i % 10)) * 2]
    t_i = table_3[4 - floor(i / 10)][(2 - (i % 10)) * 2 + 1]
    return s_i * q_rtd_C * 10 ** (-3) + t_i

# (容量可変型コンプレッサー搭載ルームエアコンディショナー)
def calc_p_i_eq9(i, q_rtd_C):
    p_i_A = table_4_A[4 - floor(i / 10)][(2 - (i % 10))]
    p_i_B = table_4_B[4 - floor(i / 10)][(2 - (i % 10))]
    p_i_C = table_4_C[4 - floor(i / 10)][(2 - (i % 10))]

    if   q_rtd_C <= 2200:                       return p_i_A
    elif 2200 < q_rtd_C and q_rtd_C <= 4000:    return p_i_A * ((4000 - q_rtd_C) / (4000 - 2200)) +\
                                                       p_i_B * ((q_rtd_C - 2200) / (4000 - 2200))
    elif 4000 < q_rtd_C and q_rtd_C < 7100:     return p_i_B * ((7100 - q_rtd_C) / (7100 - 4000)) +\
                                                       p_i_C * ((q_rtd_C - 4000) / (7100 - 4000))
    elif 7100 <= q_rtd_C:                       return p_i_C

def calc_Q_r_max_C(q_r_max_C, q_rtd_C, Theta_ex_d_t):
    a2, a1, a0 = calc_a_eq13(q_r_max_C, q_rtd_C)
    return a2 * (Theta_ex_d_t - 35) ** 2 + a1 * (Theta_ex_d_t - 35) + a0

def calc_a_eq13(q_r_max_C, q_rtd_C):
    b2 =  0.000812 * q_rtd_C * 10 ** (-3) - 0.001480
    b1 =  0.003527 * q_rtd_C * 10 ** (-3) - 0.023000
    b0 = -0.011490 * q_rtd_C * 10 ** (-3) + 1.024328
    c2 = -0.000350 * q_rtd_C * 10 ** (-3) + 0.000800
    c1 = -0.001280 * q_rtd_C * 10 ** (-3) + 0.003621
    c0 =  0.004772 * q_rtd_C * 10 ** (-3) - 0.011170
    a2 = b2 * q_r_max_C + c2
    a1 = b1 * q_r_max_C + c1
    a0 = b0 * q_r_max_C + c0
    return a2, a1, a0

def calc_f_C_Theta(x, Theta_ex, q_rtd_C, dualcompressor=False):
    a0, a1, a2, a3, a4 = calc_a_eq22(Theta_ex, q_rtd_C, dualcompressor=dualcompressor)
    return a4 * x ** 4 + a3 * x ** 3 + a2 * x ** 2 + a1 * x + a0

def calc_a_eq22(Theta_ex, q_rtd_C, dualcompressor=False):
    if dualcompressor == False:     calc_p_i = calc_p_i_eq23        # dualcompressorでない
    else:                           calc_p_i = calc_p_i_eq24        # dualcompressor

    # 係数p_i
    p42 = calc_p_i(42, q_rtd_C)
    p41 = calc_p_i(41, q_rtd_C)
    p40 = calc_p_i(40, q_rtd_C)
    p32 = calc_p_i(32, q_rtd_C)
    p31 = calc_p_i(31, q_rtd_C)
    p30 = calc_p_i(30, q_rtd_C)
    p22 = calc_p_i(22, q_rtd_C)
    p21 = calc_p_i(21, q_rtd_C)
    p20 = calc_p_i(20, q_rtd_C)
    p12 = calc_p_i(12, q_rtd_C)
    p11 = calc_p_i(11, q_rtd_C)
    p10 = calc_p_i(10, q_rtd_C)
    p02 = calc_p_i(2, q_rtd_C)
    p01 = calc_p_i(1, q_rtd_C)
    p00 = calc_p_i(0, q_rtd_C)

    a4 = p42 * Theta_ex ** 2 + p41 * Theta_ex + p40 * 1
    a3 = p32 * Theta_ex ** 2 + p31 * Theta_ex + p30 * 1
    a2 = p22 * Theta_ex ** 2 + p21 * Theta_ex + p20 * 1
    a1 = p12 * Theta_ex ** 2 + p11 * Theta_ex + p10 * 1
    a0 = p02 * Theta_ex ** 2 + p01 * Theta_ex + p00 * 1

    return a0, a1, a2, a3, a4

def calc_p_i_eq23(i, q_rtd_C):
    q_rtd_C = min(5600, q_rtd_C)
    s_i = table_5[4 - floor(i / 10)][(2 - (i % 10)) * 2]
    t_i = table_5[4 - floor(i / 10)][(2 - (i % 10)) * 2 + 1]
    return s_i * q_rtd_C * 10 ** (-3) + t_i

def calc_p_i_eq24(i, q_rtd_C):
    p_i_A = table_6_A[4 - floor(i / 10)][(2 - (i % 10))]
    p_i_B = table_6_B[4 - floor(i / 10)][(2 - (i % 10))]
    p_i_C = table_6_C[4 - floor(i / 10)][(2 - (i % 10))]

    if q_rtd_C <= 2200:                         return p_i_A
    elif 2200 < q_rtd_C and q_rtd_C <= 4000:    return p_i_A * ((4000 - q_rtd_C) / (4000 - 2200)) +\
                                                       p_i_B * ((q_rtd_C - 2200) / (4000 - 2200))
    elif 4000 <= q_rtd_C and q_rtd_C < 7100:    return p_i_B * ((7100 - q_rtd_C) / (7100 - 4000)) +\
                                                       p_i_C * ((q_rtd_C - 4000) / (7100 - 4000))
    elif 7100 <= q_rtd_C:                       return p_i_C

class AirconModel_RAC(AirconSpec):
    """
    RACモデルを適用するエアコンの派生クラス。
    """
    def __init__(self, spec):
        super().__init__(spec)

    def estimate_COP(self, mode, inputdata):
        q_rtd_C        = self.spec['Q']['cooling']['rtd'] * 1000        #kW -> W
        q_max_C        = self.spec['Q']['cooling']['max'] * 1000        #kW -> W
        p_rtd_C        = self.spec['P']['cooling']['rtd'] * 1000        #kW -> W
        q_rtd_H        = self.spec['Q']['heating']['rtd'] * 1000        #kW -> W
        q_max_H        = self.spec['Q']['heating']['max'] * 1000        #kW -> W
        p_rtd_H        = self.spec['P']['heating']['rtd'] * 1000        #kW -> W
        Theta          = inputdata['T_ex']
        h              = inputdata['X_ex'] / ae.absolute_humidity_from_e(ae.saturation_vapor_pressure(Theta)) * 100

        if(inputdata['Q_S'] < 0.0):
            L_CS           = inputdata['Q_S'] * 3.6
            L_CL           = inputdata['Q_L'] * 3.6
            return self._estimate_cooling_COP(q_rtd_C, q_max_C, p_rtd_C, L_CS, L_CL, Theta)
        elif(inputdata['Q_S'] > 0.0):
            L_H            = inputdata['Q_S'] * 3.6
            return self._estimate_heating_COP(q_rtd_C, q_rtd_H, q_max_H, p_rtd_H, L_H, Theta, h)
        else:
            return 0.0

    def _estimate_cooling_COP(self, q_rtd_C, q_max_C, p_rtd_C, L_CS, L_CL, Theta):
        ratio      = C_af_C * C_hm_C

        q_r_max_C  = q_max_C / q_rtd_C                                              # 最大冷房能力比
        Q_r_max_C  = calc_Q_r_max_C(q_r_max_C, q_rtd_C, Theta)                      # 最大冷房出力比
        Q_max_C    = Q_r_max_C * q_rtd_C * ratio * 3600 * 10 ** (-6)                # 最大冷房出力
        L_max_CL   = L_CS * ((1.0 - SHF_L_min_c) / SHF_L_min_c) # 最大冷房潜熱負荷
        L_dash_CL  = np.clip(L_CL, 0, L_max_CL)                                     # 補正冷房潜熱負荷
        L_dash_C   = L_CS + L_dash_CL
        SHF_dash   = L_CS / L_dash_C if L_dash_C > 0 else 0                         # 冷房負荷補正顕熱比
        Q_max_CS   = Q_max_C * SHF_dash                                             # 最大冷房顕熱出力, 最大冷房潜熱出力
        Q_max_CL   = np.clip(Q_max_C * (1.0 - SHF_dash), 0, L_dash_CL)
        Q_T_CS     = min(Q_max_CS, L_CS)                                            # 処理冷房負荷
        Q_T_CL     = min(Q_max_CL, L_CL)
        Q_dash_T_C = (Q_T_CS + Q_T_CL) * 1.0 / C_af_C * C_hm_C  # 補正処理冷房負荷

        x1 = calc_f_C_Theta(Q_dash_T_C / (q_max_C * 3600 * 10 ** (-6)), Theta, q_rtd_C)
        x2 = calc_f_C_Theta(1.0 / q_r_max_C,                             35.0, q_rtd_C)
        E_E_C = x1 / x2 * (p_rtd_C) * 10 ** (-3)
        return (L_CS + L_CL) / 3.6 / E_E_C, E_E_C

    def _estimate_heating_COP(self, q_rtd_C, q_rtd_H, q_max_H, p_rtd_H, L_H, Theta, h):
        ratio = C_af_H * C_df_H(Theta, h)

        q_r_max_H  = q_max_H / q_rtd_H                                              # 最大暖房能力比
        Q_r_max_H  = calc_Q_r_max_H(q_rtd_C, q_r_max_H, Theta)                      # 最大暖房出力比
        Q_max_H    = Q_r_max_H * q_rtd_H * ratio * 3600 * 10 ** (-6)                # 最大暖房出力
        Q_T_H      = min(Q_max_H, L_H)                                              # 処理暖房負荷
        Q_dash_T_H = Q_T_H * 1.0 / ratio                                            # 補正処理暖房負荷

        x1 = calc_f_H_Theta(Q_dash_T_H / (q_max_H * 3600 * 10 ** (-6)), q_rtd_C, Theta)
        x2 = calc_f_H_Theta(1.0 / q_r_max_H,                            q_rtd_C, 7.0)
        E_E_H = x1 / x2 * p_rtd_H * 10 ** (-3)

        return L_H / 3.6 / E_E_H, E_E_H
