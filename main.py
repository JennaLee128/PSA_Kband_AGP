# pairwise sequence alignment ~ Kband ~ Affine gap penalty

match = 1
mismatch = -2
h = 2
g = 1

def p(x, y):
    return match if x == y else mismatch

def InsiderStrip(i, j, k, diff = 0):
    # judge whether the loc is in the k-band area
    return (-k <= j - i <= k + diff)

def Init(m, n, k, diff = 0):
    a = [[-float('Inf')] * (n + 1) for _ in range(m + 1)]
    b = [[-float('Inf')] * (n + 1) for _ in range(m + 1)]
    c = [[-float('Inf')] * (n + 1) for _ in range(m + 1)]
    for i in range(1, k+1):
        a[i][0] = -float('Inf')
        b[i][0] = -float('Inf')
        c[i][0] = -h - g * i
    for j in range(1, k+diff+1):
        a[0][j] = -float('Inf')
        b[0][j] = -h - g * j
        c[0][j] = -float('Inf')
    a[0][0] = 0
    b[0][0] = -h
    return a, b, c


def TraceBack(a, b, c, A, B):
    # Trace back the optimal path, and get the optimal alignment
    seq_A = ""
    seq_B = ""
    m = len(A)
    n = len(B)
    i = m
    j = n

    ret = max(a[m][n], b[m][n], c[m][n])
    flag = '_'
    if ret == a[m][n]:
        flag = 'a'
    elif ret == b[m][n]:
        flag = 'b'
    else:
        flag = 'c'
    while (i > 0 or j > 0):
        if i > 0 and j > 0:
            if flag == 'a':
                seq_A += A[i - 1]
                seq_B += B[j - 1]
                if ret == p(A[i - 1], B[j - 1]) + a[i - 1][j - 1]:
                    ret = a[i - 1][j - 1]
                    flag = 'a'
                elif ret == p(A[i - 1], B[j - 1]) + b[i - 1][j - 1]:
                    ret = b[i - 1][j - 1]
                    flag = 'b'
                elif ret == p(A[i - 1], B[j - 1]) + c[i - 1][j - 1]:
                    ret = c[i - 1][j - 1]
                    flag = 'c'
                i -= 1
                j -= 1
                continue
            elif flag == 'b':
                seq_A += '-'
                seq_B += B[j - 1]
                if ret == -h - g + a[i][j - 1]:
                    ret = a[i][j - 1]
                    flag = 'a'
                elif ret == -g + b[i][j - 1]:
                    ret = b[i][j - 1]
                    flag = 'b'
                elif ret == -h - g + c[i][j - 1]:
                    ret = c[i][j - 1]
                    flag = 'c'
                j -= 1
                continue
            else:
                seq_A += A[i - 1]
                seq_B += '-'
                if ret == -h - g + a[i - 1][j]:
                    ret = a[i - 1][j]
                    flag = 'a'
                elif ret == -h - g + b[i - 1][j]:
                    ret = b[i - 1][j]
                    flag = 'b'
                elif ret == -g + c[i - 1][j]:
                    ret = c[i - 1][j]
                    flag = 'c'
                i -= 1
                continue
        if i > 0:
            seq_A += A[i - 1]
            seq_B += '-'
            i -= 1
            continue
        if j > 0:
            seq_A += '-'
            seq_B += B[j - 1]
            j -= 1
            continue
        else:
            print(i, j)
            raise ValueError('i,j are Error')
    return seq_A[::-1], seq_B[::-1]


def PSA_Kband_AGP(A: str, B: str,  get_score=0):
    """
    Affine gap penalty ~ PSA ~ Kband
    """
    # len(A)=0 or len(B)=0
    if len(A) == len(B) == 0:
        return 0, '', ''
    elif len(A) == 0:
        return -2*len(B), '-'*len(B), B
    elif len(B) == 0:
        return -2*len(A), A, '-'*len(A)
    # n>=m
    # record the loc of A & B
    state_ex = 0
    if len(A) > len(B):
        A, B = B, A
        state_ex = 1
    m = len(A)
    n = len(B)
    diff = n - m
    k = 1
    old = -float('Inf')

    # to compute the optimal score
    while k <= m:
        a, b, c = Init(m, n, k, diff)
        for i in range(1, m+1):
            for d in range(-k, k + diff + 1):
                j = i + d
                if 1 <= j <= n:
                    # A[i] ~ B[j]
                    a[i][j] = p(A[i - 1], B[j - 1]) + max(a[i - 1][j - 1], b[i - 1][j - 1], c[i - 1][j - 1])
                    if InsiderStrip(i - 1, j, k, diff):
                        # A[j] ~ _
                        b[i][j] = max(-h - g + a[i][j - 1], -g + b[i][j - 1], -h - g + c[i][j - 1])
                    if InsiderStrip(i, j - 1, k, diff):
                        # _ ~ B[i]
                        c[i][j] = max(-h - g + a[i - 1][j], -h - g + b[i - 1][j], -g + c[i - 1][j])

        new = max(a[m][n], b[m][n], c[m][n])
        if old == new or (k * 2) > m:
            if get_score:
                return new, k
            else:
                break
        else:
            old = new
            k *= 2

    SeqA, SeqB = TraceBack(a, b, c, A, B)
    # exchange the loc of A & B
    if state_ex:
        SeqA, SeqB = SeqB, SeqA
    return new, SeqA, SeqB


if __name__ == "__main__":
    A = "ACTTTGCCATTGACCCCCCCCCCCCAATTTTGA"
    B = "ACGTGTTGCCATTCCAATTTTCATTA"
    new, SeqA, SeqB = PSA_Kband_AGP(A, B)
    print(new)
    print(SeqA)
    print(SeqB)