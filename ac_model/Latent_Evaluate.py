import numpy as np
import ac_model.main as ac
from ac_model.calc_refrigerant import get_f_p_sgas, get_f_H_gas_comp_in, get_f_H_gas_comp_out,\
                                      get_f_S_gas, get_f_H_liq, calc_e_ref_H_th
import archenv.archenv as ae

c_p_air = 1006.0                # 空気の比熱       [J/Kg・K]
rho_air = 1.2                   # 空気の密度       [kg/m3]
L_wtr   = 2500.8 - 2.3668 * 27  # 水の蒸発潜熱     [kJ/kg]
c_p_w   = 1.846                 # 水蒸気の定圧比熱 [J/Kg・K]
F       = 101325                # 大気圧           [Pa]
A_f_hex = 0.23559               # 室内機熱交換器の前面面積のうち熱交換に有効な面積 [m2]
A_e_hex = 6.396                 # 室内機熱交換器の表面積のうち熱交換に有効な面積 [m2]

# 表1 単位面積当たりの必要暖房能力及び冷房能力 [W/m2]
table_1 = [
    (73.91, 64.32, 62.65, 66.99, 72.64, 61.34, 64.55, 00.00),
    (37.61, 36.55, 42.34, 54.08, 61.69, 60.79, 72.53, 61.56)
]
A_f_hex_small_H = 0.2        #定格冷却能力が5.6kW未満の場合のA_f,hex
A_e_hex_small_H = 6.2        #定格冷却能力が5.6kW未満の場合のA_e,hex
A_f_hex_large_H = 0.3        #定格冷却能力が5.6kW以上の場合のA_f,hex
A_e_hex_large_H = 10.6       #定格冷却能力が5.6kW以上の場合のA_e,hex

class AirconModel(ac.AirconSpec):
    """
    DuctCentralModelモデルを適用するエアコンの派生クラス。
    """
    def __init__(self, spec):
        super().__init__(spec)

    def estimate_COP(self, mode, inputdata):
        q_hs_rtd_C  = self.spec['Q']['cooling']['rtd']        * 1000
        q_hs_mid_C  = self.spec['Q']['cooling']['mid']        * 1000
        q_hs_min_C  = self.spec['Q']['cooling']['min']        * 1000
        P_hs_rtd_C  = self.spec['P']['cooling']['rtd']        * 1000
        V_fan_rtd_C = self.spec['V_inner']['cooling']['rtd']  * 3600
        V_fan_mid_C = self.spec['V_inner']['cooling']['mid']  * 3600
        P_fan_rtd_C = self.spec['P_fan']['cooling']['rtd']    * 1000
        P_fan_mid_C = self.spec['P_fan']['cooling']['mid']    * 1000
        V_hs_dsgn_C = self.spec['V_inner']['cooling']['dsgn'] * 3600

        q_hs_rtd_H  = self.spec['Q']['heating']['rtd']        * 1000
        q_hs_mid_H  = self.spec['Q']['heating']['mid']        * 1000
        q_hs_min_H  = self.spec['Q']['heating']['min']        * 1000
        P_hs_rtd_H  = self.spec['P']['heating']['rtd']        * 1000
        V_fan_rtd_H = self.spec['V_inner']['heating']['rtd']  * 3600
        V_fan_mid_H = self.spec['V_inner']['heating']['mid']  * 3600
        P_fan_rtd_H = self.spec['P_fan']['heating']['rtd']    * 1000
        P_fan_mid_H = self.spec['P_fan']['heating']['mid']    * 1000
        V_hs_dsgn_H = self.spec['V_inner']['heating']['dsgn'] * 3600

        Theta          = inputdata['T_ex']
        h              = inputdata['X_ex'] / ae.absolute_humidity_from_e(ae.saturation_vapor_pressure(Theta)) * 100
        Theta_hs_in    = inputdata['T_in']
        X_hs_in        = inputdata['X_in']
        V_hs_supply    = inputdata['V_inner'] * 3600

        q_hs_H         = inputdata['Q_S']

        if mode == 'cooling':
            q_hs_CS = inputdata['Q_S']
            q_hs_CL = inputdata['Q_L'] 
            return self._estimate_cooling_COP(q_hs_rtd_C, q_hs_mid_C, q_hs_min_C, P_hs_rtd_C,
                                              V_fan_rtd_C, V_fan_mid_C, P_fan_rtd_C, P_fan_mid_C, V_hs_dsgn_C,
                                              Theta, h, Theta_hs_in, X_hs_in, q_hs_CS, q_hs_CL, V_hs_supply)

        elif mode == 'heating':
            q_hs_H  = inputdata['Q_S']
            return self._estimate_heating_COP(q_hs_rtd_H, q_hs_mid_H, q_hs_min_H, P_hs_rtd_H,
                                              V_fan_rtd_H, V_fan_mid_H, P_fan_rtd_H, P_fan_mid_H, V_hs_dsgn_H,
                                              Theta, h, Theta_hs_in, X_hs_in, q_hs_H, V_hs_supply)        
        else:
            return 0.0

    def _estimate_cooling_COP(self, q_hs_rtd_C, q_hs_mid_C, q_hs_min_C, P_hs_rtd_C,
                              V_fan_rtd_C, V_fan_mid_C, P_fan_rtd_C, P_fan_mid_C, V_hs_dsgn_C,
                              Theta, h, Theta_hs_in, X_hs_in, q_hs_CS, q_hs_CL, V_hs_supply): 
        
        Theta_hs_out = Theta_hs_in - q_hs_CS / (c_p_air * rho_air * V_hs_supply / 3600)
        print('Theta_hs_out', Theta_hs_out)

        # (4) 日付dの時刻tにおける1時間当たりの熱源機の平均冷房能力 [-]
        q_hs_C  = q_hs_CS + q_hs_CL
        print('q_hs_CS', q_hs_CS)
        print('q_hs_CL', q_hs_CL)

        #潜熱評価モデル（細井先生）の計算式
        A_f_hex = A_f_hex_small_H if q_hs_rtd_C < 5600 else A_f_hex_large_H
        A_e_hex = A_e_hex_small_H if q_hs_rtd_C < 5600 else A_e_hex_large_H

        #潜熱評価モデル（細井先生）の計算式
        x = q_hs_C / 1000
        if q_hs_C > 0:  E_E_fan_C = (1.4675 * x**3 - 8.5886 * x**2  + 20.217 * x  + 50.0) * 10 ** (-3)
        else:           E_E_fan_C = 0.0
        print('E_E_fan_C', E_E_fan_C)

        ########################################################################################################################
        #a = np.clip(V_hs_supply, 400, None)
        x = np.clip(V_hs_supply, 360, None) / (3600 * A_f_hex)
        # (36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)]
        #alpha_dash_c_hex_C = 0.050 * np.log((a / 3600) / A_f_hex) + 0.073

        #潜熱評価モデル（細井先生）の計算式
        alpha_dash_c_hex_C = 0.0631 * x + 0.0015
        # (36a) 冷房時の室内熱交換器表面の顕熱伝達率 [W/(m2・K)]
        alpha_c_hex_C = alpha_dash_c_hex_C * (c_p_air + c_p_w * X_hs_in)
        print('  alpha_c_hex_C', alpha_c_hex_C)
        print('  alpha_dash_c_hex_C', alpha_dash_c_hex_C)

        # (32)
        Theta_sur_f_hex_C = ((Theta_hs_in + Theta_hs_out) / 2) - (c_p_air * rho_air * V_hs_supply * (Theta_hs_in - Theta_hs_out) / 3600) / (A_e_hex * alpha_c_hex_C)
        print('  Theta_sur_f_hex_C', Theta_sur_f_hex_C)

        # (28) 冷房時の冷媒の蒸発温度 [℃]
        Theta_ref_evp_C = -50 if Theta_sur_f_hex_C < -50 else Theta_sur_f_hex_C
        print('  Theta_ref_evp_C', Theta_ref_evp_C)

        # (27) 冷房時の冷媒の凝縮温度 [℃]
        Theta_ref_cnd_C = np.maximum(Theta + 27.4 - 1.35 * Theta_ref_evp_C, Theta)
        if Theta_ref_cnd_C > 65:                    Theta_ref_cnd_C = 65
        if Theta_ref_cnd_C < Theta_ref_evp_C + 5.0: Theta_ref_cnd_C = Theta_ref_evp_C + 5.0
        print('  Theta_ref_cnd_C', Theta_ref_cnd_C)

        # (29) 冷房時の冷媒の過冷却度 [℃]
        Theta_ref_SC_C = np.maximum(0.772 * Theta_ref_cnd_C - 25.6, 0)
        print('  Theta_ref_SC_C', Theta_ref_SC_C)

        # (30) 冷房時の冷媒の過熱度 [℃]
        Theta_ref_SH_C = np.maximum(0.194 * Theta_ref_cnd_C - 3.86, 0)
        print('  Theta_ref_SH_C', Theta_ref_SH_C)

        # 4_8_a (1) ヒートポンプサイクルの理論暖房効率 [-]
        e_dash_th = calc_e_ref_H_th(Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C)
        print('  e_dash_th', e_dash_th)

        # (18) 冷房時のヒートポンプサイクルの理論効率 [-]
        e_th_C = e_dash_th - 1
        print('e_th_C', e_th_C)
        ########################################################################################################################

        #潜熱評価モデル（細井先生）の計算式
        x = q_hs_C / 1000
        e_r_C = -0.0316 * x**2 + 0.2944 * x

        print('e_r_C', e_r_C)

        # (8) 冷房時の熱源機の効率 [-]
        e_hs_C = e_th_C * e_r_C
        print('e_hs_C', e_hs_C)

        # (6) 冷房時の圧縮機の消費電力量（kWh/h）
        E_E_comp_C = (q_hs_C / e_hs_C) * 10 ** (-3) if q_hs_C > 0 else 0
        print('E_E_comp_C', E_E_comp_C)

        # (2)
        E_E_C_d_t = E_E_comp_C + E_E_fan_C

        return q_hs_C / e_hs_C * 10 ** (-3), E_E_C_d_t
    
    def _estimate_heating_COP(self, q_hs_rtd_H, q_hs_mid_H, q_hs_min_H, P_hs_rtd_H,
                              V_fan_rtd_H, V_fan_mid_H, P_fan_rtd_H, P_fan_mid_H, V_hs_dsgn_H,
                              Theta, h, Theta_hs_in, X_hs_in, q_hs_H, V_hs_supply):
        
        
        C_df_H = 0.77 if (Theta < 5.0) and (h >= 80.0) else 1.0

        q_hs_H = q_hs_H * (1 / C_df_H)
        print('q_hs_H', q_hs_H)

        #潜熱評価モデル（細井先生）の計算式
        A_f_hex = A_f_hex_small_H if q_hs_rtd_C < 5600 else A_f_hex_large_H
        A_e_hex = A_e_hex_small_H if q_hs_rtd_C < 5600 else A_e_hex_large_H

        #潜熱評価モデル（細井先生）の計算式
        x = q_hs_H / 1000
        if q_hs_H > 0:  E_E_fan_H = (1.4675 * x**3 - 8.5886 * x**2  + 20.217 * x  + 50.0) * 10 ** (-3)
        else:           E_E_fan_H = 0.0
        print('E_E_fan_H', E_E_fan_H)

        ########################################################################################################################

        # (35) 暖房時の室内熱交換器表面の顕熱伝達率 [-]
        alpha_c_hex_H = (-0.0017 * ((V_hs_supply / 3600) / A_f_hex) ** 2 + 0.044 * ((V_hs_supply / 3600) / A_f_hex) + 0.0271) * 10 ** 3
        print('  alpha_c_hex_H', alpha_c_hex_H)

        # (31) 暖房時の室内機熱交換器の表面温度 [℃]
        Theta_sur_f_hex_H = ((Theta_hs_in + Theta_hs_out) / 2) + (c_p_air * rho_air * V_hs_supply * (Theta_hs_out - Theta_hs_in) / 3600) / (A_e_hex * alpha_c_hex_H)
        print('  Theta_sur_f_hex_H', Theta_sur_f_hex_H)

        # (23) 暖房時の冷媒の凝縮温度 [℃]
        Theta_ref_cnd_H = 65 if Theta_sur_f_hex_H > 65 else Theta_sur_f_hex_H
        print('  Theta_ref_cnd_H', Theta_ref_cnd_H)

        # (24) 暖房時の冷媒の蒸発温度 [℃]
        Theta_ref_evp_H = Theta - (0.100 * Theta_ref_cnd_H + 2.95)
        if Theta_ref_evp_H < -50:                       Theta_ref_evp_H = -50
        if Theta_ref_evp_H > Theta_ref_cnd_H - 5.0:     Theta_ref_evp_H = Theta_ref_cnd_H - 5.0
        print('  Theta_ref_evp_H', Theta_ref_evp_H)

        # (25) 暖房時の冷媒の過冷却度 [℃]
        Theta_ref_SC_H = 0.245 * Theta_ref_cnd_H - 1.72
        print('  Theta_ref_SC_H', Theta_ref_SC_H)

        # (26) 暖房時の冷媒の過熱度 [℃]
        Theta_ref_SH_H = 4.49 - 0.036 * Theta_ref_cnd_H
        print('  Theta_ref_SH_H', Theta_ref_SH_H)

        # (17) 日付dの時刻tにおける暖房時のヒートポンプサイクルの理論効率 [-]
        e_th_H = calc_e_ref_H_th(Theta_ref_evp_H, Theta_ref_cnd_H, Theta_ref_SC_H, Theta_ref_SH_H)
        print('e_th_H', e_th_H)
        ########################################################################################################################

        #潜熱評価モデル（細井先生）の計算式
        x = q_hs_H / 1000
        e_r_H = -0.0316 * x**2 + 0.2944 * x
        print('e_r_H', e_r_H)

        # (7) 日付dの時刻tにおける暖房時の熱源機の効率 [-]
        e_hs_H = e_th_H * e_r_H
        print('e_hs_H', e_hs_H)

        # (5)
        E_E_comp_H = (q_hs_H / e_hs_H) * 10 ** (-3) if q_hs_H > 0 else 0
        print('E_E_comp_H', E_E_comp_H)

        # (1)
        E_E_H_d_t = E_E_comp_H + E_E_fan_H

        return q_hs_H / e_hs_H * 10 ** (-3), E_E_H_d_t
