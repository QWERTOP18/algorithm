class Rerooting:
    def __init__(self, N):
        self.N = N
        self.G = [[] for _ in range(N)]  # グラフの初期化
        self.dp = [[] for _ in range(N)]  # dp[v][i]: vから出るi番目の有向辺に対応する部分木のDP
        self.ans = [None] * N  # ans[v]: 頂点vを根とする木の答え
        
        self.identity = -1  # 単位元
        self.merge = lambda dp_cum, d: max(dp_cum, d)  # 二項演算
        self.add_root = lambda d: d + 1  # 根を追加する演算
    
    def add_edge(self, a, b):
        self.G[a].append(b)
    
    def build(self):
        self.dfs(0, -1)  # DFS で木 DP を計算
        self.bfs(0, self.identity, -1)  # BFS で残りの部分木に対応する DP を計算

    def dfs(self, v, p):
        dp_cum = self.identity
        deg = len(self.G[v])
        self.dp[v] = [self.identity] * deg
        
        for i in range(deg):
            u = self.G[v][i]
            if u == p:
                continue
            self.dp[v][i] = self.dfs(u, v)
            dp_cum = self.merge(dp_cum, self.dp[v][i])
        
        return self.add_root(dp_cum)
    
    def bfs(self, v, dp_p, p):
        deg = len(self.G[v])
        
        # 前の DFS で計算した部分木の DP を保存
        for i in range(deg):
            if self.G[v][i] == p:
                self.dp[v][i] = dp_p
        
        dp_l = [self.identity] * (deg + 1)
        dp_r = [self.identity] * (deg + 1)
        
        # 左からの累積
        for i in range(deg):
            dp_l[i + 1] = self.merge(dp_l[i], self.dp[v][i])
        
        # 右からの累積
        for i in range(deg - 1, -1, -1):
            dp_r[i] = self.merge(dp_r[i + 1], self.dp[v][i])
        
        # 頂点 v の答え
        self.ans[v] = self.add_root(dp_l[deg])
        
        for i in range(deg):
            u = self.G[v][i]
            if u == p:
                continue
            new_dp_p = self.add_root(self.merge(dp_l[i], dp_r[i + 1]))
            self.bfs(u, new_dp_p, v)
    

if __name__ == "__main__":
    import sys
    input = sys.stdin.read
    data = input().split()
    
    N = int(data[0])
    reroot = Rerooting(N)
    
    index = 1
    for _ in range(N - 1):
        u = int(data[index]) - 1
        v = int(data[index + 1]) - 1
        reroot.add_edge(u, v)
        reroot.add_edge(v, u)
        index += 2
    
    reroot.build()
    
    for i in range(N):
        print(f"頂点{i + 1}: {reroot.ans[i]}")
