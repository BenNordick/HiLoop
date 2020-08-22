class DisjointCache:
    def __init__(self, nodes):
        self._nodes = nodes
        self._cache = dict()
    def isdisjoint(self, holder1, holder2):
        if (self._nodes < 10):
            return holder1.value.isdisjoint(holder2.value)
        first, second = holder1.order(holder2)
        if (first in self._cache):
            subcache = self._cache[first]
        else:
            subcache = dict()
            self._cache[first] = subcache
        if (second in subcache):
            return subcache[second]
        else:
            result = holder1.value.isdisjoint(holder2.value)
            subcache[second] = result
            return result
