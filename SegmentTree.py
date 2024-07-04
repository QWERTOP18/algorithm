class SegmentTree:
    def __init__(self, size):
        self.size = size
        self.tree = [0] * (2 * size)
    
    def update(self, idx, value):
        idx += self.size
        self.tree[idx] = max(self.tree[idx], value)
        while idx > 1:
            idx //= 2
            self.tree[idx] = max(self.tree[2 * idx], self.tree[2 * idx + 1])
    
    def query(self, left, right):
        res = 0
        left += self.size
        right += self.size
        while left < right:
            if left % 2 == 1:
                res = max(res, self.tree[left])
                left += 1
            if right % 2 == 1:
                right -= 1
                res = max(res, self.tree[right])
            left //= 2
            right //= 2
        return res
