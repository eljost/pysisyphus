def canonical_order(L):
    inds = list()
    for i in range(L + 1):
        l = L - i
        for n in range(i + 1):
            m = i - n
            inds.append((l, m, n))
    return inds


def cca_order(l):
    inds = [
        (l, 0, 0),
    ]
    a = l
    b = c = 0
    for _ in range(((l + 1) * (l + 2) // 2) - 1):
        if c < l - a:
            b -= 1
            c += 1
        else:
            a -= 1
            c = 0
            b = l - a
        inds.append((a, b, c))
    return inds
