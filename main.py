import math
import numpy as np
import pandas as pd

class Const:
    """
    定数を管理するクラス
    """
    MODES           = ['cooling', 'heating'] # 冷房、暖房
    ERROR_THRESHOLD = 1e-3                   # 温度の許容誤差
    MAX_TEMP        = 50                     # 最大温度探索範囲

class AirconSpec:
    """
    エアコンの仕様を管理する基底クラス。
    """
    def __init__(self, spec):
        self.spec = spec
        self.COP = {mode: {} for mode in Const.MODES}
        self._validate_spec()
        self._calculate_COP()

    def _validate_spec(self):
        """仕様データのバリデーションを行う。"""
        for key in self.spec:
            if not all(mode in self.spec.get(key, {}) for mode in Const.MODES):
                raise KeyError(f"{key} に {', '.join(Const.MODES)} がありません。")

    def _calculate_COP(self):
        """COP（成績係数）を計算する。"""
        if 'Q' in self.spec and 'P' in self.spec:
            for mode in Const.MODES:
                for key in set(self.spec['Q'][mode].keys()) & set(self.spec['P'][mode].keys()):
                    power = self.spec['P'][mode][key]
                    if power > 0:
                        self.COP[mode][key] = self.spec['Q'][mode][key] / power
