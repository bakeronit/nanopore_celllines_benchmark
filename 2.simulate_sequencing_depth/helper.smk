#wildcard_constraints:
#    sample = "|".join(samples),
#    depth = "|".join(["60x","45x","30x","15x"])

def my_dirty_combinator(*args):
    depth_t = args[1]
    depth_n = args[3]
    combs = [[depth_t[n],i] for n in range(4) for i in depth_n[n:n+3]]
    result = [[]]
    result = [ [t] + y for y in combs for t in args[0]]
    [r.insert(2,args[2][0]) for r in result]
    for w in result:
        yield tuple(w)



