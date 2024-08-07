class Matrix():
    def __init__(self, mat, mod=None):
        self.mat = mat
        self.n = len(mat)
        self.m = len(mat[0])
        self.mod = mod
    def __mul__(self, other):
        ret = Matrix([[0]*other.m for _ in range(self.n)], self.mod)
        for i in range(self.n):
            for j in range(other.m):
                for k in range(self.m):
                    ret[i][j] += self.mat[i][k]*other.mat[k][j]
                    if self.mod is not None: ret[i][j] %= self.mod
        return ret
    def __add__(self, other):
        ret = Matrix([[self.mat[i][j] for j in range(self.m)] for i in range(self.n)], self.mod)
        for i in range(other.n):
            for j in range(other.m):
                ret[i][j] += other.mat[i][j]
                if self.mod is not None: ret[i][j] %= self.mod
        return ret
    def __sub__(self, other):
        ret = Matrix([[self.mat[i][j] for j in range(self.m)] for i in range(self.n)], self.mod)
        for i in range(other.n):
            for j in range(other.m):
                ret[i][j] -= other.mat[i][j]
                if self.mod is not None: ret[i][j] %= self.mod
        return ret
    def __pow__(self, scalar):
        a = Matrix([[self.mat[i][j] for j in range(self.m)] for i in range(self.n)], self.mod)
        ret = Matrix.e(self.n, self.mod)
        while scalar:
            if scalar&1:
                ret *= a
            a *= a
            scalar >>= 1
        return ret
    def scalar_mul(self, a):
        ret = Matrix([[self.mat[i][j] for j in range(self.m)] for i in range(self.n)], self.mod)
        for i in range(self.n):
            for j in range(self.m):
                ret[i][j] *= a
                if self.mod is not None: ret[i][j] %= self.mod
        return ret
    def __repr__(self) -> str:
        return self.mat.__repr__()
    def __getitem__(self, i):
        return self.mat[i]
    def __setitem__(self, i, x):
        self.mat[i] = x
    def __len__(self):
        return len(self.mat)
    def t(self):
        return Matrix([list(column) for column in zip(*self.mat)], self.mod)
    def turn(matrix):
        if type(matrix) != 'Matrix':
            return Matrix([list(column) for column in zip(*matrix)])
        return Matrix([list(column) for column in zip(*matrix.mat)], matrix.mod)
    def e(size, mod):
        return Matrix([[i == j for j in range(size)] for i in range(size)], mod)
    def zeros(h, w):
        return Matrix([[0]*w for _ in range(h)])
    # 行列Aの法mod上の掃き出し法/Aは拡大係数行列か？
    # mod=2の場合はBitsetを用いて定数倍(約wordsize=64倍)の高速化を行う
    def gauss_jordan(A, mod, is_extended = False):
        h, w = len(A), len(A[0])
        rank = 0
        if mod != 2:
            for j in range(w):
                if is_extended and j == w-1: break
                pivot = None
                for i in range(rank, h):
                    if A[i][j] != 0:
                        pivot = i
                        break
                if pivot is None: continue
                A[rank], A[pivot] = A[pivot], A[rank]
                inv = pow(A[rank][j], mod-2, mod)
                for k in range(j, w):
                    A[rank][k] = A[rank][k]*inv%mod
                for i in range(h):
                    if i == rank or A[i][j] == 0: continue
                    a = A[i][j]
                    for k in range(j, w):
                        A[i][k] = (A[i][k]-A[rank][k]*a)%mod
                rank += 1
        else:
            if type(A[0]) != Bitset:
                for i in range(h):
                    A[i] = Bitset(w, A[i])
            for j in range(w):
                if is_extended and j == w-1: break
                pivot = None
                for i in range(rank, h):
                    if A[i][j] != 0:
                        pivot = i
                        break
                if pivot is None: continue
                A[rank], A[pivot] = A[pivot], A[rank]
                for i in range(h):
                    if i == rank or A[i][j] == 0: continue
                    A[i] ^= A[rank]
                rank += 1
        return rank
    # 連立一次方程式 Ax=b を法mod上で解く
    def linear_equation(A, b, mod):
        h, w = len(A), len(A[0])
        M = [A[i]+[b[i]] for i in range(h)] # 拡大係数行列
        rank = Matrix.gauss_jordan(M, mod, True)
        # 解が無ければNoneをreturn
        # rank行以下の拡大部分にα≠0が残っていた場合、0=αという等式が残っていることになり矛盾
        for i in range(rank, h):
            if M[i][w]: return None, None
        ret = [M[i][w] for i in range(h)]
        return ret, rank
    # 行列Aの行列式を法mod上で求める
    def det(A, mod):
        n = len(A)
        tmp = [[a for a in row] for row in A]
        rank = 0
        ret = 1
        if mod != 2:
            for j in range(n):
                pivot = None
                for i in range(rank, n):
                    if tmp[i][j] != 0:
                        pivot = i
                        break
                if pivot is None: return 0
                if rank != pivot:
                    tmp[rank], tmp[pivot] = tmp[pivot], tmp[rank]
                    ret *= -1
                inv = pow(tmp[rank][j], mod-2, mod)
                ret = ret*tmp[rank][j]%mod
                for k in range(j, n):
                    tmp[rank][k] = tmp[rank][k]*inv%mod
                for i in range(n):
                    if i == rank or tmp[i][j] == 0: continue
                    a = tmp[i][j]
                    for k in range(j, n):
                        tmp[i][k] = (tmp[i][k]-tmp[rank][k]*a)%mod
                rank += 1
        else:
            if type(tmp[0]) != Bitset:
                for i in range(n):
                    tmp[i] = Bitset(n, tmp[i])
            for j in range(n):
                pivot = None
                for i in range(rank, n):
                    if tmp[i][j] != 0:
                        pivot = i
                        break
                if pivot is None: return 0
                tmp[rank], tmp[pivot] = tmp[pivot], tmp[rank]
                for i in range(n):
                    if i == rank or tmp[i][j] == 0: continue
                    tmp[i] ^= A[rank]
                rank += 1
        return ret
    # 行列Aの逆行列を法mod上で求める
    def inv(A, mod):
        n = len(A)
        tmp = [[a for a in row] for row in A]
        rank = 0
        ret = Matrix.e(n, mod)
        if mod != 2:
            for j in range(n):
                pivot = None
                for i in range(rank, n):
                    if tmp[i][j] != 0:
                        pivot = i
                        break
                if pivot is None: return None
                if rank != pivot:
                    tmp[rank], tmp[pivot] = tmp[pivot], tmp[rank]
                    ret[rank], ret[pivot] = ret[pivot], ret[rank]
                inv = pow(tmp[rank][j], mod-2, mod)
                for k in range(n):
                    tmp[rank][k] = tmp[rank][k]*inv%mod
                    ret[rank][k] = ret[rank][k]*inv%mod
                for i in range(n):
                    if i == rank or tmp[i][j] == 0: continue
                    a = tmp[i][j]
                    for k in range(n):
                        tmp[i][k] = (tmp[i][k]-tmp[rank][k]*a)%mod
                        ret[i][k] = (ret[i][k]-ret[rank][k]*a)%mod
                rank += 1
        else:
            if type(tmp[0]) != Bitset:
                for i in range(n):
                    tmp[i] = Bitset(n, tmp[i])
            for j in range(n):
                pivot = None
                for i in range(rank, n):
                    if tmp[i][j] != 0:
                        pivot = i
                        break
                if pivot is None: return None
                tmp[rank], tmp[pivot] = tmp[pivot], tmp[rank]
                ret[rank], ret[pivot] = ret[pivot], ret[rank]
                for i in range(n):
                    if i == rank or tmp[i][j] == 0: continue
                    tmp[i] ^= tmp[rank]
                    ret[i] ^= ret[rank]
                rank += 1
        return ret
