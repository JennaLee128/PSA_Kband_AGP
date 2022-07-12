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
    a = [[-float('Inf')] * (diff+2*k+1) for _ in range(m + 1)]
    b = [[-float('Inf')] * (diff+2*k+1) for _ in range(m + 1)]
    c = [[-float('Inf')] * (diff+2*k+1) for _ in range(m + 1)]
    for i in range(1, k+1):
        c[i][k - i] = -h - g * i
    for j in range(1, k+diff+1):
        b[0][j + k] = -h - g * j
    a[0][k] = 0
    return a, b, c


def TraceBack(a, b, c, A, B, k):
    # Trace back the optimal path, and get the optimal alignment
    seq_A = ""
    seq_B = ""
    m = len(A)
    n = len(B)
    diff = n - m
    i = m
    b_j = n
    j = diff + k

    ret = max(a[-1][-1-k], b[-1][-1-k], c[-1][-1-k])
    if ret == a[-1][-1-k]:
        flag = 'a'
    elif ret == b[-1][-1-k]:
        flag = 'b'
    else:
        flag = 'c'

    while (i > 0 or j > k):
        if flag == 'a' and i > 0 and j > 0:
            seq_A += A[i - 1]
            seq_B += B[b_j - 1]
            if ret == p(A[i - 1], B[b_j - 1]) + a[i - 1][j]:
                ret = a[i - 1][j]
                flag = 'a'
            elif ret == p(A[i - 1], B[b_j - 1]) + b[i - 1][j]:
                ret = b[i - 1][j]
                flag = 'b'
            elif ret == p(A[i - 1], B[b_j - 1]) + c[i - 1][j]:
                ret = c[i - 1][j]
                flag = 'c'
            i -= 1
            b_j -= 1
            continue
        elif flag == 'b' and j > 0:
            seq_A += '-'
            seq_B += B[b_j - 1]
            if ret == -h - g + a[i][j - 1]:
                ret = a[i][j - 1]
                flag = 'a'
            elif ret == -g + b[i][j - 1]:
                ret = b[i][j - 1]
                flag = 'b'
            elif ret == -h - g + c[i][j - 1]:
                ret = c[i][j - 1]
                flag = 'c'
            b_j -= 1
            j -= 1
            continue
        elif flag == 'c' and i > 0 and j + 1 <= 2 * k + diff:
            seq_A += A[i - 1]
            seq_B += '-'
            if ret == -h - g + a[i - 1][j + 1]:
                ret = a[i - 1][j + 1]
                flag = 'a'
            elif ret == -h - g + b[i - 1][j + 1]:
                ret = b[i - 1][j + 1]
                flag = 'b'
            elif ret == -g + c[i - 1][j + 1]:
                ret = c[i - 1][j + 1]
                flag = 'c'
            i -= 1
            j += 1
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
                if 1 <= d + i <= n:
                    j = d + k
                    # A[i] ~ B[j]
                    a[i][j] = p(A[i - 1], B[d + i - 1]) + max(a[i - 1][j], b[i - 1][j], c[i - 1][j])
                    if InsiderStrip(i - 1, d + i, k, diff):
                        # A[j] ~ _
                        # c[i][j] = max(-h - g + a[i - 1][j + 1], -g + c[i - 1][j + 1])
                        c[i][j] = max(-h - g + a[i - 1][j + 1], -h - g + b[i - 1][j + 1], -g + c[i - 1][j + 1])
                    if InsiderStrip(i, d + i - 1, k, diff):
                        # _ ~ B[i]
                        # b[i][j] = max(-h - g + a[i][j - 1], -g + b[i][j - 1])
                        b[i][j] = max(-h - g + a[i][j - 1], -g + b[i][j - 1], -h - g + c[i][j - 1])

        new = max(a[-1][-1-k], b[-1][-1-k], c[-1][-1-k])
        if old == new or (k * 2) > m:
            if get_score:
                return new, k
            else:
                break
        else:
            old = new
            k *= 2

    SeqA, SeqB = TraceBack(a, b, c, A, B, k)
    # exchange the loc of A & B
    if state_ex:
        SeqA, SeqB = SeqB, SeqA
    return new, SeqA, SeqB

def input_seq(f):
    X = []
    X_name = []
    for line in f:
        if not line.startswith('>'):
            X.append(line.replace('\n', ''))  # 去掉行尾的换行符
        else:
            X_name.append(line.replace('\n', ''))
    f.close()
    return ''.join(X), ''.join(X_name)

def output_seq(X, X_name, f):
    f.write(X_name +'\n')
    n = 0
    seq = []
    for i in X:
        seq.append(i)
        n += 1
        if n > 80:
            f.write(''.join(seq)+'\n')
            n = 0
            seq = []
    f.close()

if __name__ == "__main__":
    # input sequences
    f1 = open('C:/Users/LYZ/Desktop/1.fasta', 'r')
    f2 = open('C:/Users/LYZ/Desktop/2.fasta', 'r')
    A, A_name = input_seq(f1)
    B, B_name = input_seq(f2)

    print('序列A的长度为：', len(A))
    print('序列B的长度为：', len(B))
    new, SeqA, SeqB = PSA_Kband_AGP(A, B)

    # output sequences
    f3 = open('C:/Users/LYZ/Desktop/9.fasta', 'w')
    f4 = open('C:/Users/LYZ/Desktop/10.fasta', 'w')
    output_seq(SeqA, A_name, f3)
    output_seq(SeqB, B_name, f4)
    print('得分为：', new)
    print('比对结束！\n')