import ac_model.main as ac
import ac_model.calc_refrigerant

class AirconModel_DuctCentral(ac.AirconSpec):
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

        # (38) 1時間当たりの送風機の消費電力量のうちの冷房設備への付加分 [kWh/h]
        if q_hs_C > 0:  E_E_fan_C = P_fan_rtd_C * V_hs_supply / V_hs_dsgn_C * 10 ** (-3)
        else:           E_E_fan_C = 0.0
        print('E_E_fan_C', E_E_fan_C)

        ########################################################################################################################
        a = np.clip(V_fan_mid_C, 400, None)
        # (36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)]
        alpha_dash_c_hex_C = 0.050 * np.log((a / 3600) / A_f_hex) + 0.073
        # (36a) 冷房時の室内熱交換器表面の顕熱伝達率 [W/(m2・K)]
        alpha_c_hex_C = alpha_dash_c_hex_C * (c_p_air + c_p_w * 0.010376)

        print('  alpha_c_hex_C(mid)',      alpha_c_hex_C)
        print('  alpha_dash_c_hex_C(mid)', alpha_dash_c_hex_C)

        # (34) 冷房時の室内機熱交換器の表面温度 [℃]
        def func(x):
            q_hs_CS = (27 - x) / (3600 / (2 * c_p_air * rho_air * V_fan_mid_C) + 1 / (A_e_hex * alpha_c_hex_C))

            T = x + 273.16 #(6) # 絶対温度（K）

            # 表1 式(5b)における係数の値
            a1 = -6096.9385
            a2 = 21.2409642
            a3 = -0.02711193
            a4 = 0.00001673952
            a5 = 2.433502
            b1 = -6024.5282
            b2 = 29.32707
            b3 = 0.010613863
            b4 = -0.000013198825
            b5 = -0.49382577

            #(5b)
            if x > 0:   k = a1 / T + a2 + a3 * T + a4 * T ** 2 + a5 * np.log(T)
            else:       k = b1 / T + b2 + b3 * T + b4 * T ** 2 + b5 * np.log(T)

            P_vs = np.exp(k)                            #(5a) 飽和水蒸気圧（Pa）
            X_surf_hex_C = 0.622 * (P_vs / (F - P_vs))  #(3) 飽和空気の絶対湿度（kg/kg(DA)）

            q_hs_CL = np.maximum((0.010376 - X_surf_hex_C) / (3600 / (2 * L_wtr * rho_air * V_fan_mid_C * 10 ** 3) + 1 / (L_wtr * A_e_hex * alpha_dash_c_hex_C * 10 ** 3)), 0)

            return q_hs_CS + q_hs_CL - q_hs_mid_C

        # x = fsolve(func, x0) は、点 x0 を開始点として方程式 fun(x) = 0 (ゼロの配列) の解を求めようとする
        Theta_sur_f_hex_C = optimize.bisect(func, -273.15, 99.96)
        print('  Theta_sur_f_hex_C(mid)', Theta_sur_f_hex_C)

        # (28) 冷房時の冷媒の蒸発温度（℃）
        Theta_ref_evp_C = -50 if Theta_sur_f_hex_C < -50 else Theta_sur_f_hex_C
        print('  Theta_ref_evp_C(mid)', Theta_ref_evp_C)

        # (27) 冷房時の冷媒の凝縮温度（℃）
        Theta_ref_cnd_C = np.maximum(35 + 27.4 - 1.35 * Theta_ref_evp_C, 35)
        if Theta_ref_cnd_C > 65:                    Theta_ref_cnd_C = 65
        if Theta_ref_cnd_C < Theta_ref_evp_C + 5.0: Theta_ref_cnd_C = Theta_ref_evp_C + 5.0
        print('  Theta_ref_cnd_C(mid)', Theta_ref_cnd_C)

        # (29) 冷房時の冷媒の過冷却度（℃）
        Theta_ref_SC_C = np.maximum(0.772 * Theta_ref_cnd_C - 25.6, 0)
        print('  Theta_ref_SC_C(mid)', Theta_ref_SC_C)

        # (30) 冷房時の冷媒の過熱度（℃）
        Theta_ref_SH_C = np.maximum(0.194 * Theta_ref_cnd_C - 3.86, 0)
        print('  Theta_ref_SH_C(mid)', Theta_ref_SH_C)

        # 4_8_a (1) ヒートポンプサイクルの理論暖房効率 [-]
        e_dash_th_mid_C = calc_e_ref_H_th(Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C)
        print('  e_dash_th_mid_C(mid)', e_dash_th_mid_C)

        # (22) 中間暖房能力運転時のヒートポンプサイクルの理論効率 [-]
        e_th_mid_C = e_dash_th_mid_C - 1

        print('e_th_mid_C', e_th_mid_C)
        ########################################################################################################################
        a = np.clip(V_fan_rtd_C, 400, None)
        # (36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)]
        alpha_dash_c_hex_C = 0.050 * np.log((a / 3600) / A_f_hex) + 0.073
        # (36a) 冷房時の室内熱交換器表面の顕熱伝達率 [W/(m2・K)]
        alpha_c_hex_C = alpha_dash_c_hex_C * (c_p_air + c_p_w * 0.010376)

        # (34) 冷房時の室内機熱交換器の表面温度 [℃]
        def func(x):
            q_hs_CS = (27 - x) / (3600 / (2 * c_p_air * rho_air * V_fan_rtd_C) + 1 / (A_e_hex * alpha_c_hex_C))

            T = x + 273.16 #(6) # 絶対温度（K）

            # 表1 式(5b)における係数の値
            a1 = -6096.9385
            a2 = 21.2409642
            a3 = -0.02711193
            a4 = 0.00001673952
            a5 = 2.433502
            b1 = -6024.5282
            b2 = 29.32707
            b3 = 0.010613863
            b4 = -0.000013198825
            b5 = -0.49382577

            #(5b)
            if x > 0:   k = a1 / T + a2 + a3 * T + a4 * T ** 2 + a5 * np.log(T)
            else:       k = b1 / T + b2 + b3 * T + b4 * T ** 2 + b5 * np.log(T)

            P_vs = np.exp(k)                            #(5a) 飽和水蒸気圧（Pa）
            X_surf_hex_C = 0.622 * (P_vs / (F - P_vs))  #(3) 飽和空気の絶対湿度（kg/kg(DA)）

            q_hs_CL = np.maximum((0.010376 - X_surf_hex_C) / (3600 / (2 * L_wtr * rho_air * V_fan_rtd_C * 10 ** 3) + 1 / (L_wtr * A_e_hex * alpha_dash_c_hex_C * 10 ** 3)), 0)

            return q_hs_CS + q_hs_CL - q_hs_rtd_C

        # x = fsolve(func, x0) は、点 x0 を開始点として方程式 fun(x) = 0 (ゼロの配列) の解を求めようとする
        Theta_sur_f_hex_C = optimize.bisect(func, -273.15, 99.96)

        # (28) 冷房時の冷媒の蒸発温度（℃）
        Theta_ref_evp_C = -50 if Theta_sur_f_hex_C < -50 else Theta_sur_f_hex_C

        # (27) 冷房時の冷媒の凝縮温度（℃）
        Theta_ref_cnd_C = np.maximum(35 + 27.4 - 1.35 * Theta_ref_evp_C, 35)
        if Theta_ref_cnd_C > 65:                    Theta_ref_cnd_C = 65
        if Theta_ref_cnd_C < Theta_ref_evp_C + 5.0: Theta_ref_cnd_C = Theta_ref_evp_C + 5.0

        # (29) 冷房時の冷媒の過冷却度（℃）
        Theta_ref_SC_C = np.maximum(0.772 * Theta_ref_cnd_C - 25.6, 0)

        # (30) 冷房時の冷媒の過熱度（℃）
        Theta_ref_SH_C = np.maximum(0.194 * Theta_ref_cnd_C - 3.86, 0)

        # 4_8_a (1) ヒートポンプサイクルの理論暖房効率 [-]
        e_dash_th_rtd_C = calc_e_ref_H_th(Theta_ref_evp_C, Theta_ref_cnd_C, Theta_ref_SC_C, Theta_ref_SH_C)

        # (21) 中間暖房能力運転時のヒートポンプサイクルの理論効率 [-]
        e_th_rtd_C = e_dash_th_rtd_C - 1
        print('e_th_rtd_C(rtd)', e_th_rtd_C)
        ########################################################################################################################
        a = np.clip(V_hs_supply, 400, None)
        # (36b) 冷房時の室内熱交換器表面の潜熱伝達率 [kg/(m2・s)]
        alpha_dash_c_hex_C = 0.050 * np.log((a / 3600) / A_f_hex) + 0.073
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
        # (12) 定格冷房能力運転時の熱源機の効率 [-]
        e_hs_rtd_C = q_hs_rtd_C / (P_hs_rtd_C - P_fan_rtd_C)
        e_r_rtd_C = e_hs_rtd_C / e_th_rtd_C
        e_r_rtd_C = np.clip(e_r_rtd_C, 0, 1.0)
        print('e_r_rtd_C', e_r_rtd_C)

        # (16) 最小冷房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
        e_r_min_C = e_r_rtd_C * 0.65
        print('e_r_min_C', e_r_min_C)

        # (14) 中間冷房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
        e_r_mid_C = e_r_rtd_C * 0.95
        #e_hs_mid_C = q_hs_mid_C / (P_hs_mid_C - P_fan_mid_C)   #定格能力試験と中間能力試験の値を入力する場合はこちら
        #e_r_mid_C = e_hs_mid_C / e_th_mid_C

        e_r_mid_C = np.clip(e_r_mid_C, 0, 1.0)
        print('e_r_mid_C', e_r_mid_C)

        # (10) 冷房時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
        if q_hs_C <= q_hs_min_C:
            e_r_C = e_r_min_C - (q_hs_min_C - q_hs_C) * (e_r_min_C / q_hs_min_C)

        if (q_hs_min_C < q_hs_C) & (q_hs_C <= q_hs_mid_C):
            e_r_C = e_r_mid_C - (q_hs_mid_C - q_hs_C) * ((e_r_mid_C - e_r_min_C) / (q_hs_mid_C - q_hs_min_C))

        if (q_hs_mid_C < q_hs_C) & (q_hs_C <= q_hs_rtd_C):
            e_r_C = e_r_rtd_C - (q_hs_rtd_C - q_hs_C) * ((e_r_rtd_C - e_r_mid_C) / (q_hs_rtd_C - q_hs_mid_C))

        if (q_hs_rtd_C < q_hs_C) and (e_r_rtd_C > 0.4):
            e_r_C = np.clip(e_r_rtd_C - (q_hs_C - q_hs_rtd_C) * (e_r_rtd_C / q_hs_rtd_C), 0.4, None)

        if (q_hs_rtd_C < q_hs_C) and (e_r_rtd_C <= 0.4):
            e_r_C = e_r_rtd_C
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

        #(37) 1時間当たりの送風機の消費電力量のうちの暖房設備への付加分 [kWh/h]
        if q_hs_H > 0:  E_E_fan_H = P_fan_rtd_H * V_hs_supply / V_hs_dsgn_H * 10 ** (-3)
        else:           E_E_fan_H = 0.0
        print('E_E_fan_H', E_E_fan_H)

        ########################################################################################################################

        # (35) 暖房時の室内熱交換器表面の顕熱伝達率 [-]
        alpha_c_hex_H = (-0.0017 * ((V_fan_mid_H / 3600) / A_f_hex) ** 2 + 0.044 * ((V_fan_mid_H / 3600) / A_f_hex) + 0.0271) * 10 ** 3
        print('  alpha_c_hex_H', alpha_c_hex_H)

        # (33) 暖房時の室内機熱交換器の表面温度 [℃]
        Theta_sur_f_hex_H = 20 + (q_hs_mid_H / (2 * c_p_air * rho_air * V_fan_mid_H)) * 3600 + (q_hs_mid_H / (A_e_hex * alpha_c_hex_H))
        print('  Theta_sur_f_hex_H', Theta_sur_f_hex_H)

        # (23) 暖房時の冷媒の凝縮温度 [℃]
        Theta_ref_cnd_H = 65 if Theta_sur_f_hex_H > 65 else Theta_sur_f_hex_H
        print('  Theta_ref_cnd_H', Theta_ref_cnd_H)

        # (24) 暖房時の冷媒の蒸発温度 [℃]
        Theta_ref_evp_H = 7 - (0.100 * Theta_ref_cnd_H + 2.95)
        if Theta_ref_evp_H < -50:                       Theta_ref_evp_H = -50
        if Theta_ref_evp_H > Theta_ref_cnd_H - 5.0:     Theta_ref_evp_H = Theta_ref_cnd_H - 5.0
        print('  Theta_ref_evp_H', Theta_ref_evp_H)

        # (25) 暖房時の冷媒の過冷却度 [℃]
        Theta_ref_SC_H = 0.245 * Theta_ref_cnd_H - 1.72
        print('  Theta_ref_SC_H', Theta_ref_SC_H)

        # (26) 暖房時の冷媒の過熱度 [℃]
        Theta_ref_SH_H = 4.49 - 0.036 * Theta_ref_cnd_H
        print('  Theta_ref_SH_H', Theta_ref_SH_H)

        # 4_8_a (1) ヒートポンプサイクルの理論暖房効率
        e_dash_th_mid_H = calc_e_ref_H_th(Theta_ref_evp_H, Theta_ref_cnd_H, Theta_ref_SC_H, Theta_ref_SH_H)
        print('  e_dash_th_mid_H', e_dash_th_mid_H)

        # (20) 中間暖房能力運転時のヒートポンプサイクルの理論効率 [-]
        e_th_mid_H = e_dash_th_mid_H
        print('e_th_mid_H', e_th_mid_H)
        ########################################################################################################################

        # (35) 暖房時の室内熱交換器表面の顕熱伝達率 [-]
        alpha_c_hex_H = (-0.0017 * ((V_fan_rtd_H / 3600) / A_f_hex) ** 2 + 0.044 * ((V_fan_rtd_H / 3600) / A_f_hex) + 0.0271) * 10 ** 3
        print('  alpha_c_hex_H', alpha_c_hex_H)

        # (33) 暖房時の室内機熱交換器の表面温度 [℃]
        Theta_sur_f_hex_H = 20 + (q_hs_rtd_H / (2 * c_p_air * rho_air * V_fan_rtd_H)) * 3600 + (q_hs_rtd_H / (A_e_hex * alpha_c_hex_H))
        print('  Theta_sur_f_hex_H', Theta_sur_f_hex_H)

        # (23) 暖房時の冷媒の凝縮温度 [℃]
        Theta_ref_cnd_H = 65 if Theta_sur_f_hex_H > 65 else Theta_sur_f_hex_H
        print('  Theta_ref_cnd_H', Theta_ref_cnd_H)

        # (24) 暖房時の冷媒の蒸発温度 [℃]
        Theta_ref_evp_H = 7 - (0.100 * Theta_ref_cnd_H + 2.95)
        if Theta_ref_evp_H < -50:                       Theta_ref_evp_H = -50
        if Theta_ref_evp_H > Theta_ref_cnd_H - 5.0:     Theta_ref_evp_H = Theta_ref_cnd_H - 5.0
        print('  Theta_ref_evp_H', Theta_ref_evp_H)

        # (25) 暖房時の冷媒の過冷却度 [℃]
        Theta_ref_SC_H = 0.245 * Theta_ref_cnd_H - 1.72
        print('  Theta_ref_SC_H', Theta_ref_SC_H)

        # (26) 暖房時の冷媒の過熱度 [℃]
        Theta_ref_SH_H = 4.49 - 0.036 * Theta_ref_cnd_H
        print('  Theta_ref_SH_H', Theta_ref_SH_H)

        # (19) 定格暖房能力運転時のヒートポンプサイクルサイクルの理論効率 [-]
        e_th_rtd_H = calc_e_ref_H_th(Theta_ref_evp_H, Theta_ref_cnd_H, Theta_ref_SC_H, Theta_ref_SH_H)
        print('e_th_rtd_H', e_th_rtd_H)
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

        # (11) 定格暖房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
        e_hs_rtd_H = q_hs_rtd_H / (P_hs_rtd_H - P_fan_rtd_H)
        e_r_rtd_H  = e_hs_rtd_H / e_th_rtd_H
        e_r_rtd_H  = np.clip(e_r_rtd_H, 0, 1.0)
        print('e_r_rtd_H', e_r_rtd_H)

        # (15) 最小暖房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
        e_r_min_H = e_r_rtd_H * 0.65
        print('e_r_min_H', e_r_min_H)

        # (13) 中間暖房能力運転時のヒートポンプサイクルの理論効率に対する熱源機の効率の比
        e_r_mid_H = e_r_rtd_H * 0.95
        #e_hs_mid_H = q_hs_mid_H / (P_hs_mid_H - P_fan_mid_H)    #定格能力試験と中間能力試験の値を入力する場合はこちら
        #e_r_mid_H = e_hs_mid_H / e_th_mid_H
        e_r_mid_H = np.clip(e_r_mid_H, 0, 1.0)
        print('e_r_mid_H', e_r_mid_H)

        # (9) 暖房時のヒートポンプサイクルの理論効率に対する熱源機の効率の比 [-]
        if q_hs_H <= q_hs_min_H:
            e_r_H = e_r_min_H - (q_hs_min_H - q_hs_H) * (e_r_min_H / q_hs_min_H)

        if (q_hs_min_H < q_hs_H) & (q_hs_H <= q_hs_mid_H):
            e_r_H = e_r_mid_H - (q_hs_mid_H - q_hs_H) * ((e_r_mid_H - e_r_min_H) / (q_hs_mid_H - q_hs_min_H))

        if (q_hs_mid_H < q_hs_H) & (q_hs_H <= q_hs_rtd_H):
            e_r_H = e_r_rtd_H - (q_hs_rtd_H - q_hs_H) * ((e_r_rtd_H - e_r_mid_H) / (q_hs_rtd_H - q_hs_mid_H))

        if (q_hs_rtd_H < q_hs_H) and e_r_rtd_H > 0.4:
            e_r_H = np.clip(e_r_rtd_H - (q_hs_H - q_hs_rtd_H) * (e_r_rtd_H / q_hs_rtd_H), 0.4, None)

        if (q_hs_rtd_H < q_hs_H) and e_r_rtd_H <= 0.4:
            e_r_H = e_r_rtd_H
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
