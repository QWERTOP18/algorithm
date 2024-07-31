W,N=map(int,input().split())
Spice=[]
for i in range(N):
  L,R,V=map(int,input().split())
  Spice.append((L,R,V))
  
INF=10**18
LN=1<<(W+1).bit_length()
D=[-INF]*LN*2



def f(x,y):
  return max(x,y)
  
def update(x,y):
  x+=LN
  D[x]=y
  x>>=1
  while x>0:
    D[x]=f(D[x*2],D[x*2+1])
    x>>=1
    
def query(L,R):
  L+=LN
  R+=LN
  ret=-INF
  while L<R:
    if L%2:ret=f(ret,D[L]);L+=1
    if R%2:R-=1;ret=f(ret,D[R])
    L>>=1;R>>=1
  return ret
  
update(0,0)

for l,r,v in Spice:
  for i in range(W,0,-1):
    t = max(query(i,i+1),query(max(0,i-r),min(W,i-l)+1)+v)
    update(i,t)
    
ans = query(W,W+1)
print(ans if ans>=0 else -1)
  
