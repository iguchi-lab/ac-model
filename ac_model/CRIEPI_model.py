import ac_model.main as ac
import archenv.archenv as ae
import archenv.JIS as JIS

BF          = 0.2                    # バイパスファクター
KEYS_CRIEPI = ['min', 'rtd', 'max']  # CRIEPIモデルで使用するキー
ERROR_THRESHOLD = 1e-3               # 温度の許容誤差
MAX_TEMP        = 50                 # 最大温度探索範囲

def avoid_over_saturation(Td, X, threshold=ERROR_THRESHOLD, max_temp=MAX_TEMP):
    """
    飽和超過を回避するため、等エンタルピー線に沿って温度を調整する。
    """
    enthalpy = ae.total_enthalpy(Td, x=X)
    left, right = Td, max_temp

    while right - left > threshold:
        mid = (left + right) / 2
        X_new = (enthalpy - ae.const.SPECIFIC_HEAT_AIR * mid) / (
            ae.const.SPECIFIC_HEAT_WATER_VAPOR * mid + ae.const.LATENT_HEAT_VAPORIZATION
        )
        X_saturation = ae.absolute_humidity_from_e(ae.saturation_vapor_pressure(mid))

        if X_saturation < X_new:
            left = mid
        else:
            right = mid

    if ae.absolute_humidity_from_e(ae.saturation_vapor_pressure(left)) > X:
        return left
    else:
        return right

class AirconModel_CRIEPI(ac.AirconSpec):
    """
    CRIEPIモデルを適用するエアコンの派生クラス。
    """
    def __init__(self, spec):
        super().__init__(spec)
        self.coeffs = {mode: {} for mode in ac.Const.MODES}  # CRIEPIモデル用
        self.Pc     = {mode: {} for mode in ac.Const.MODES}  # CRIEPIモデル用
        self.prepare_CRIEPI_model()

    def _calculate_efficiency(self, mode, key):
        """
        熱効率を計算する。
        """
        Q       = self.spec['Q'][mode][key]
        P       = self.spec['P'][mode][key]
        V_inner = self.spec['V_inner'][mode]['rtd']
        V_outer = self.spec['V_outer'][mode]['rtd']

        if mode == 'cooling':
            return self._calc_cooling_efficiency(Q, P, V_inner, V_outer)
        else:
            return self._calc_heating_efficiency(Q, P, V_inner, V_outer)

    def _calc_cooling_efficiency(self, Q, P, V_inner, V_outer):
        """
        冷房時の熱効率を計算。
        """
        M_evp = (1 - BF) * ae.air_density(JIS.T_C_IN) * V_inner
        T_evp = JIS.T_C_IN - Q / (M_evp * ae.air_specific_heat(JIS.X_C_IN))
        T_evp = avoid_over_saturation(T_evp, JIS.X_C_IN)
        M_cnd = (1 - BF) * ae.air_density(JIS.T_C_EX) * V_outer
        T_cnd = JIS.T_C_EX + (Q + P) / (M_cnd * ae.air_specific_heat(JIS.X_C_EX))
        return (T_evp + 273.15) / (T_cnd - T_evp)

    def _calc_heating_efficiency(self, Q, P, V_inner, V_outer):
        """
        暖房時の熱効率を計算。
        """
        M_evp = (1 - BF) * ae.air_density(JIS.T_H_EX) * V_outer
        T_evp = JIS.T_H_EX - (Q - P) / (M_evp * ae.air_specific_heat(JIS.X_H_EX))
        T_evp = avoid_over_saturation(T_evp, JIS.X_H_EX)
        M_cnd = (1 - BF) * ae.air_density(JIS.T_H_IN) * V_inner
        T_cnd = JIS.T_H_IN + Q / (M_cnd * ae.air_specific_heat(JIS.X_H_IN))
        return (T_cnd + 273.15) / (T_cnd - T_evp)

    def _solve_coefficients(self, COP, Q, eta_th):
        """
        係数RとPcを求める。
        """
        X_matrix = np.array([[1 / eta_th['min'], 1 / Q['min']], [1 / eta_th['rtd'], 1 / Q['rtd']]])
        Y_matrix = np.array([[1 / COP['min']], [1 / COP['rtd']]])
        solution = np.linalg.solve(X_matrix, Y_matrix)
        R_minrtd, Pc = 1 / solution[0][0], solution[1][0]
        R_max = COP['max'] * Q['max'] / (Q['max'] - COP['max'] * Pc) / eta_th['max']
        coeffs = np.polyfit([Q['min'], Q['rtd'], Q['max']], [R_minrtd, R_minrtd, R_max], 2)
        return coeffs, Pc

    def prepare_CRIEPI_model(self):
        """
        CRIEPIモデルの準備を行い、係数を求める。
        """
        for mode in ac.Const.MODES:
            eta_th = {key: self._calculate_efficiency(mode, key) for key in KEYS_CRIEPI}
            COP    = self.COP[mode]
            Q      = self.spec['Q'][mode]
            self.coeffs[mode], self.Pc[mode] = self._solve_coefficients(COP, Q, eta_th)
        return self.coeffs, self.Pc

    def estimate_COP(self, mode, inputdata):
        """
        COPを推定し、消費電力も求める。
        """
        T_in    = inputdata['T_in']
        T_ex    = inputdata['T_ex']
        X_in    = inputdata['X_in']
        X_ex    = inputdata['X_ex']
        Q       = inputdata['Q']
        V_inner = inputdata['V_inner']
        V_outer = inputdata['V_outer']

        coeff_R = self.coeffs[mode][0] * Q ** 2 + self.coeffs[mode][1] * Q + self.coeffs[mode][2]
        Pc      = self.Pc[mode]

        if mode == 'cooling':
            return self._estimate_cooling_COP(T_in, X_in, T_ex, X_ex, Q, V_inner, V_outer, coeff_R, Pc)
        else:
            return self.e_stimate_heating_COP(T_in, X_in, T_ex, X_ex, Q, V_inner, V_outer, coeff_R, Pc)

    def _estimate_cooling_COP(self, T_in, X_in, T_ex, X_ex, Q, V_inner, V_outer,
                              coeff_R, Pc, max_iterations=100, tolerance=1e-3, COP=5):
        """
        暖房時のCOPを推定し、消費電力も求める。
        """
        M_evp = (1 - BF) * ae.air_density(T_in) * V_inner
        M_cnd = (1 - BF) * ae.air_density(T_ex) * V_outer
        T_evp = T_in - Q / (M_evp * ae.air_specific_heat(X_in))
        T_evp = avoid_over_saturation(T_evp, X_in)

        for _ in range(max_iterations):
            T_cnd = T_ex + (Q + Q / COP) / (M_cnd * ae.air_specific_heat(X_ex))
            Refrigeration_COP = coeff_R * (T_evp + 273.15) / (T_cnd - T_evp)
            calc_COP = Refrigeration_COP * Q / (Q + Pc * Refrigeration_COP)
            if math.isclose(calc_COP, COP, abs_tol=tolerance):
                break
            COP = calc_COP
        return COP, Q / COP

    def _estimate_heating_COP(self, T_in, X_in, T_ex, X_ex, Q, V_inner, V_outer,
                              coeff_R, Pc, max_iterations=100, tolerance=1e-3, COP=5):
        """
        冷房時のCOPを推定し、消費電力も求める。
        """
        M_evp = (1 - BF) * ae.air_density(T_ex) * V_outer
        M_cnd = (1 - BF) * ae.air_density(T_in) * V_inner
        T_cnd = T_in + Q / (M_cnd * ae.air_specific_heat(X_in))

        for _ in range(max_iterations):
            T_evp = T_ex - (Q - Q / COP) / (M_evp * ae.air_specific_heat(X_ex))
            T_evp = avoid_over_saturation(T_evp, X_ex)
            Refrigeration_COP = coeff_R * (T_cnd + 273.15) / (T_cnd - T_evp)
            calc_COP = Refrigeration_COP * Q / (Q + Pc * Refrigeration_COP)
            if math.isclose(calc_COP, COP, abs_tol=tolerance):
                break
            COP = calc_COP
        return COP, Q / COP
