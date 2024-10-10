def dijkstra(s):
  dist = [INF] * N
  dist[s] = 0
  que = []
  heappush(que,(0,s))
  while(que):
    cost,v = heappop(que)
    for vv,weight in G[v]:
      if dist[vv] > weight + cost:
        dist[vv] = weight + cost
        heappush(que,(dist[vv],vv))
        
  return dist
