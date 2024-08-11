# https://kyoroid.github.io/algorithm/string/rolling_hash.html
class RollingHash():

    def __init__(self, S, b=3491, m=999999937):
        """任意の基数と法でハッシュを生成する"""
        n = len(S)
        self.prefix = prefix = [0] * (n+1)
        self.power = power = [1] * (n+1)
        self.b = b
        self.m = m
        for i in range(n):
            c = ord(S[i])
            prefix[i+1] = (prefix[i] * b + c) % m
            power[i+1] = (power[i] * b) % m

    def get(self, l, r):
        """S[l, r) のハッシュを求める"""
        return (self.prefix[r] - self.power[r-l] * self.prefix[l]) % self.m

    def concat(self, h1, h2, l2):
        """S1+S2 のハッシュを、それぞれのハッシュから求める"""
        return (self.power[l2] * h1 + h2) % self.m

    def lcp(self, l1, r1, l2, r2):
        """S[l1, r1) とS[l2, r2) の最大共通接頭辞を求める"""
        # LCPの最小値 (範囲内)
        low = 0
        # LCPの最大値 + 1 (範囲外)
        high = min(r1-l1, r2-l2) + 1
        while high - low > 1:
            mid = (high + low) // 2
            if self.get(l1, l1 + mid) == self.get(l2, l2 + mid):
                low = mid
            else:
                high = mid
        return low
