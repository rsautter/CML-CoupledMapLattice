def doublingMap(x,par):
    if x < 0.5:
        return 2.0 * x
    else:
        return 2.0 * x - 1.0
    
def logisticMap(x, par):
    if(len(par)>0):
        a = par[0]
    else:
        a = 4.0
    return (a * x * (1.0 - x))
    
def somMap(x, par):
    a = par[0]
    b = par[1]
    gamma = 0.8 / (1.0 + a)
    if x < gamma:
        return (a * x + 0.2)
    elif x < 0.8:
        return (x - 0.8) / a + 1.0
    else:
        return (1.0 - x) / b