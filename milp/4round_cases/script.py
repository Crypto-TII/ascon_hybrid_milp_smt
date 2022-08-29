in_file = 'possible_configurations.txt'

# 3-round invfeasible configurations
ZZF = [[1, 4, 10], [1, 5, 9], [1, 6, 8], [1, 7, 7], [1, 8, 6], [1, 9, 5], [1, 10, 4], [1, 11, 3], [2, 4, 9], [2, 5, 8], [2, 6, 7], [2, 7, 6], [2, 8, 5], [2, 9, 4], [2, 10, 3], [3, 3, 9], [3, 4, 8], [3, 5, 7], [3, 6, 6], [3, 7, 5], [3, 8, 4], [3, 9, 3], [3, 10, 2], [4, 4, 7], [4, 5, 6], [4, 6, 5], [4, 7, 4], [5, 4, 6], [5, 5, 5], [5, 6, 4], [5, 7, 3], [6, 4, 5], [6, 5, 4], [7, 3, 5], [7, 4, 4], [9, 3, 3], [1, 4, 11], [1, 5, 10], [1, 6, 9], [1, 7, 8], [1, 8, 7], [1, 9, 6], [1, 10, 5], [1, 11, 4], [2, 5, 9], [2, 6, 8], [2, 7, 7], [2, 8, 6], [2, 9, 5], [2, 10, 4], [2, 11, 3], [3, 3, 10], [3, 4, 9], [3, 5, 8], [3, 6, 7], [3, 7, 6], [3, 8, 5], [3, 9, 4], [3, 10, 3], [4, 4, 8], [4, 5, 7], [4, 6, 6], [4, 7, 5], [4, 8, 4], [4, 9, 3], [4, 10, 2], [5, 4, 7], [5, 5, 6], [5, 6, 5], [5, 7, 4], [6, 4, 6], [6, 5, 5], [6, 6, 4], [6, 7, 3], [7, 3, 6], [7, 4, 5], [7, 5, 4], [8, 4, 4], [9, 3, 4], [10, 2, 4], [10, 3, 3], [1, 4, 12], [1, 5, 11], [1, 6, 10], [1, 7, 9], [1, 8, 8], [1, 9, 7], [1, 10, 6], [1, 11, 5], [2, 5, 10], [2, 6, 9], [2, 7, 8], [2, 8, 7], [2, 9, 6], [2, 10, 5], [2, 11, 4], [3, 3, 11], [3, 4, 10], [3, 5, 9], [3, 6, 8], [3, 7, 7], [3, 8, 6], [3, 9, 5], [3, 10, 4], [3, 11, 3], [4, 4, 9], [4, 5, 8], [4, 6, 7], [4, 7, 6], [4, 8, 5], [4, 9, 4], [4, 10, 3], [5, 4, 8], [5, 5, 7], [5, 6, 6], [5, 7, 5], [5, 8, 4], [5, 9, 3], [5, 10, 2], [6, 4, 7], [6, 5, 6], [6, 6, 5], [6, 7, 4], [7, 3, 7], [7, 4, 6], [7, 5, 5], [7, 6, 4], [7, 7, 3], [8, 4, 5], [8, 5, 4], [9, 3, 5], [9, 4, 4], [10, 2, 5], [10, 3, 4], [11, 3, 3], [1, 4, 13], [1, 5, 12], [1, 6, 11], [1, 7, 10], [1, 8, 9], [1, 9, 8], [1, 10, 7], [1, 11, 6], [2, 5, 11], [2, 6, 10], [2, 7, 9], [2, 8, 8], [2, 9, 7], [2, 10, 6], [2, 11, 5], [2, 12, 4], [2, 13, 3], [3, 3, 12], [3, 4, 11], [3, 5, 10], [3, 6, 9], [3, 7, 8], [3, 8, 7], [3, 9, 6], [3, 10, 5], [3, 11, 4], [4, 5, 9], [4, 6, 8], [4, 7, 7], [4, 8, 6], [4, 9, 5], [4, 10, 4], [4, 11, 3], [5, 4, 9], [5, 5, 8], [5, 6, 7], [5, 7, 6], [5, 8, 5], [5, 9, 4], [5, 10, 3], [6, 4, 8], [6, 5, 7], [6, 6, 6], [6, 7, 5], [6, 8, 4], [6, 9, 3], [6, 10, 2], [7, 3, 8], [7, 4, 7], [7, 5, 6], [7, 6, 5], [7, 7, 4], [8, 4, 6], [8, 5, 5], [8, 6, 4], [8, 7, 3], [9, 3, 6], [9, 4, 5], [9, 5, 4], [10, 2, 6], [10, 3, 5], [10, 4, 4], [11, 3, 4], [1, 4, 14], [1, 5, 13], [1, 6, 12], [1, 7, 11], [1, 8, 10], [1, 9, 9], [1, 10, 8], [1, 11, 7], [2, 5, 12], [2, 6, 11], [2, 7, 10], [2, 8, 9], [2, 9, 8], [2, 10, 7], [2, 11, 6], [2, 12, 5], [2, 13, 4], [3, 4, 12], [3, 5, 11], [3, 6, 10], [3, 7, 9], [3, 8, 8], [3, 9, 7], [3, 10, 6], [3, 11, 5], [3, 12, 4], [3, 13, 3], [4, 5, 10], [4, 6, 9], [4, 7, 8], [4, 8, 7], [4, 9, 6], [4, 10, 5], [4, 11, 4], [5, 4, 10], [5, 5, 9], [5, 6, 8], [5, 7, 7], [5, 8, 6], [5, 9, 5], [5, 10, 4], [5, 11, 3], [6, 4, 9], [6, 5, 8], [6, 6, 7], [6, 7, 6], [6, 8, 5], [6, 9, 4], [6, 10, 3], [7, 3, 9], [7, 4, 8], [7, 5, 7], [7, 6, 6], [7, 7, 5], [7, 8, 4], [7, 9, 3], [7, 10, 2], [8, 4, 7], [8, 5, 6], [8, 6, 5], [8, 7, 4], [9, 3, 7], [9, 4, 6], [9, 5, 5], [9, 6, 4], [9, 7, 3], [10, 2, 7], [10, 3, 6], [10, 4, 5], [10, 5, 4], [11, 3, 5], [11, 4, 4], [13, 3, 3], [1, 4, 15], [1, 6, 13], [1, 7, 12], [1, 8, 11], [1, 9, 10], [1, 10, 9], [1, 11, 8], [2, 6, 12], [2, 7, 11], [2, 8, 10], [2, 9, 9], [2, 10, 8], [2, 11, 7], [2, 12, 6], [2, 13, 5], [2, 14, 4], [2, 15, 3], [2, 16, 2], [3, 4, 13], [3, 5, 12], [3, 6, 11], [3, 7, 10], [3, 8, 9], [3, 9, 8], [3, 10, 7], [3, 11, 6], [3, 12, 5], [3, 13, 4], [4, 5, 11], [4, 6, 10], [4, 7, 9], [4, 8, 8], [4, 9, 7], [4, 10, 6], [4, 11, 5], [4, 12, 4], [4, 13, 3], [5, 4, 11], [5, 5, 10], [5, 6, 9], [5, 7, 8], [5, 8, 7], [5, 9, 6], [5, 10, 5], [5, 11, 4], [6, 4, 10], [6, 7, 7], [6, 8, 6], [6, 9, 5], [6, 10, 4], [6, 11, 3], [7, 3, 10], [7, 4, 9], [7, 5, 8], [7, 6, 7], [7, 7, 6], [7, 8, 5], [7, 9, 4], [7, 10, 3], [8, 4, 8], [8, 5, 7], [8, 6, 6], [8, 7, 5], [8, 8, 4], [8, 9, 3], [8, 10, 2], [9, 3, 8], [9, 4, 7], [9, 5, 6], [9, 6, 5], [9, 7, 4], [10, 2, 8], [10, 3, 7], [10, 4, 6], [10, 5, 5], [10, 6, 4], [10, 7, 3], [11, 3, 6], [11, 4, 5], [11, 5, 4], [12, 4, 4], [13, 3, 4], [1, 4, 16], [1, 6, 14], [1, 7, 13], [1, 8, 12], [1, 9, 11], [1, 10, 10], [1, 11, 9], [2, 7, 12], [2, 8, 11], [2, 9, 10], [2, 10, 9], [2, 11, 8], [2, 12, 7], [2, 13, 6], [2, 14, 5], [2, 15, 4], [2, 16, 3], [3, 4, 14], [3, 7, 11], [3, 8, 10], [3, 9, 9], [3, 10, 8], [3, 11, 7], [3, 12, 6], [3, 13, 5], [3, 14, 4], [3, 15, 3], [3, 16, 2], [1, 4, 17], [1, 6, 15], [1, 7, 14], [1, 8, 13], [1, 9, 12], [1, 10, 11], [1, 11, 10], [2, 9, 11], [2, 10, 10], [2, 11, 9], [2, 12, 8], [2, 13, 7], [2, 14, 6], [2, 15, 5], [2, 16, 4], [2, 17, 3], [2, 18, 2], [3, 4, 15], [3, 15, 4], [3, 16, 3], [4, 16, 2], [7, 3, 12], [9, 10, 3], [10, 3, 9], [10, 10, 2], [11, 3, 8], [12, 7, 3], [15, 3, 4], [16, 2, 4], [1, 4, 18], [1, 6, 16], [1, 7, 15], [1, 8, 14], [1, 9, 13], [1, 10, 12], [1, 11, 11], [2, 7, 14], [2, 8, 13], [2, 9, 12], [2, 10, 11], [2, 11, 10], [2, 12, 9], [2, 13, 8], [2, 14, 7], [2, 15, 6], [2, 16, 5], [2, 17, 4], [2, 18, 3], [3, 4, 16], [3, 15, 5], [3, 16, 4], [3, 17, 3], [3, 18, 2], [5, 16, 2], [7, 3, 13], [7, 13, 3], [9, 3, 11], [9, 11, 3], [10, 3, 10], [10, 10, 3], [11, 9, 3], [11, 10, 2], [12, 7, 4], [13, 3, 7], [13, 7, 3], [14, 5, 4], [15, 3, 5], [15, 4, 4], [16, 2, 5], [16, 3, 4], [17, 3, 3], [1, 4, 19], [1, 6, 17], [1, 7, 16], [1, 8, 15], [1, 9, 14], [1, 10, 13], [1, 11, 12], [2, 7, 15], [2, 8, 14], [2, 9, 13], [2, 10, 12], [2, 11, 11], [2, 12, 10], [2, 13, 9], [2, 14, 8], [2, 15, 7], [2, 16, 6], [2, 17, 5], [2, 18, 4], [2, 19, 3], [2, 20, 2], [3, 4, 17], [3, 15, 6], [3, 16, 5], [3, 17, 4], [3, 18, 3], [4, 17, 3], [4, 18, 2], [5, 16, 3], [6, 16, 2], [8, 13, 3], [10, 3, 11], [10, 11, 3], [11, 9, 4], [11, 10, 3], [12, 9, 3], [12, 10, 2], [13, 3, 8], [13, 7, 4], [14, 7, 3], [15, 3, 6], [15, 5, 4], [16, 2, 6], [16, 3, 5], [16, 4, 4], [17, 3, 4], [18, 2, 4], [18, 3, 3], [1, 4, 20], [1, 6, 18], [1, 7, 17], [1, 8, 16], [1, 9, 15], [1, 10, 14], [1, 11, 13], [2, 8, 15], [2, 9, 14], [2, 10, 13], [2, 11, 12], [2, 12, 11], [2, 13, 10], [2, 14, 9], [2, 15, 8], [2, 16, 7], [2, 17, 6], [2, 18, 5], [2, 19, 4], [2, 20, 3], [3, 4, 18], [3, 15, 7], [3, 16, 6], [3, 17, 5], [3, 18, 4], [3, 19, 3], [3, 20, 2], [4, 18, 3], [5, 18, 2], [6, 16, 3], [7, 15, 3], [7, 16, 2], [9, 13, 3], [10, 3, 12], [11, 11, 3], [12, 10, 3], [13, 3, 9], [13, 9, 3], [13, 10, 2], [14, 7, 4], [15, 3, 7], [15, 7, 3], [16, 2, 7], [16, 5, 4], [17, 3, 5], [17, 4, 4], [18, 2, 5], [18, 3, 4], [19, 3, 3], [1, 4, 21], [1, 6, 19], [1, 7, 18], [1, 8, 17], [1, 9, 16], [1, 10, 15], [1, 11, 14], [2, 8, 16], [2, 9, 15], [2, 10, 14], [2, 11, 13], [2, 12, 12], [2, 13, 11], [2, 14, 10], [2, 15, 9], [2, 16, 8], [2, 17, 7], [2, 18, 6], [2, 19, 5], [2, 20, 4], [2, 21, 3], [2, 22, 2], [3, 4, 19], [3, 17, 6], [3, 18, 5], [3, 19, 4], [3, 20, 3], [4, 20, 2], [5, 18, 3], [6, 17, 3], [6, 18, 2], [7, 16, 3], [8, 15, 3], [8, 16, 2], [10, 2, 14], [10, 3, 13], [10, 13, 3], [12, 11, 3], [13, 3, 10], [13, 10, 3], [14, 9, 3], [14, 10, 2], [15, 3, 8], [15, 7, 4], [16, 2, 8], [16, 7, 3], [17, 3, 6], [17, 5, 4], [18, 2, 6], [18, 3, 5], [18, 4, 4], [19, 3, 4], [20, 2, 4], [20, 3, 3], [1, 4, 22], [1, 6, 20], [1, 8, 18], [1, 9, 17], [1, 10, 16], [1, 11, 15], [2, 9, 16], [2, 10, 15], [2, 11, 14], [2, 12, 13], [2, 13, 12], [2, 14, 11], [2, 15, 10], [2, 16, 9], [2, 17, 8], [2, 18, 7], [2, 19, 6], [2, 20, 5], [2, 21, 4], [2, 22, 3], [3, 4, 20], [3, 17, 7], [3, 18, 6], [3, 19, 5], [3, 20, 4], [3, 21, 3], [3, 22, 2], [5, 20, 2], [6, 18, 3], [7, 17, 3], [7, 18, 2], [8, 16, 3], [9, 15, 3], [9, 16, 2], [10, 2, 15], [10, 3, 14], [11, 13, 3], [13, 11, 3], [14, 10, 3], [15, 3, 9], [15, 10, 2], [16, 2, 9], [16, 3, 8], [16, 7, 4], [17, 3, 7], [17, 7, 3], [18, 2, 7], [18, 3, 6], [18, 5, 4], [19, 3, 5], [19, 4, 4], [20, 2, 5], [20, 3, 4], [21, 3, 3], [1, 4, 23], [1, 6, 21], [1, 8, 19], [1, 9, 18], [1, 10, 17], [1, 11, 16], [2, 9, 17], [2, 10, 16], [2, 11, 15], [2, 12, 14], [2, 13, 13], [2, 14, 12], [2, 15, 11], [2, 16, 10], [2, 17, 9], [2, 18, 8], [2, 19, 7], [2, 20, 6], [2, 21, 5], [2, 22, 4], [3, 4, 21], [3, 18, 7], [3, 19, 6], [3, 20, 5], [3, 21, 4], [3, 22, 3], [4, 22, 2], [6, 20, 2], [8, 17, 3], [8, 18, 2], [9, 16, 3], [10, 2, 16], [10, 3, 15], [10, 15, 3], [10, 16, 2], [12, 13, 3], [15, 3, 10], [15, 10, 3], [16, 3, 9], [16, 10, 2], [17, 3, 8], [17, 7, 4], [18, 2, 8], [18, 7, 3], [19, 5, 4], [20, 2, 6], [21, 3, 4], [22, 2, 4], [22, 3, 3], [1, 4, 24], [1, 6, 22], [1, 8, 20], [1, 9, 19], [1, 10, 18], [1, 11, 17], [2, 10, 17], [2, 11, 16], [2, 12, 15], [2, 13, 14], [2, 14, 13], [2, 15, 12], [2, 16, 11], [2, 17, 10], [2, 18, 9], [2, 19, 8], [2, 20, 7], [2, 21, 6], [2, 22, 5], [3, 4, 22], [3, 20, 6], [3, 21, 5], [3, 22, 4], [3, 23, 3], [5, 22, 2], [7, 20, 2], [9, 17, 3], [9, 18, 2], [10, 2, 17], [10, 3, 16], [10, 16, 3], [11, 15, 3], [11, 16, 2], [13, 13, 3], [15, 3, 11], [16, 2, 11], [16, 3, 10], [16, 10, 3], [17, 10, 2], [18, 3, 8], [18, 7, 4], [19, 7, 3], [20, 2, 7], [20, 5, 4], [22, 2, 5], [22, 3, 4], [1, 4, 25], [1, 6, 23], [1, 8, 21], [1, 9, 20], [1, 10, 19], [1, 11, 18], [2, 10, 18], [2, 11, 17], [2, 12, 16], [2, 13, 15], [2, 14, 14], [2, 15, 13], [2, 16, 12], [2, 17, 11], [2, 18, 10], [2, 19, 9], [2, 20, 8], [2, 21, 7], [2, 22, 6], [3, 4, 23], [3, 21, 6], [3, 22, 5], [3, 23, 4], [3, 24, 3], [6, 22, 2], [8, 20, 2], [9, 18, 3], [10, 2, 18], [10, 3, 17], [10, 17, 3], [10, 18, 2], [11, 16, 3], [12, 15, 3], [12, 16, 2], [16, 3, 11], [17, 10, 3], [18, 3, 9], [18, 10, 2], [19, 7, 4], [20, 2, 8], [20, 7, 3], [21, 5, 4], [22, 2, 6], [23, 3, 4], [24, 3, 3]]

# 3-round all valid configurations up to 21 active Sboxes
ZZT_21 = [[2, 4, 10], [2, 4, 11], [2, 4, 12], [4, 4, 10], 
       [2, 4, 13], [3, 3, 13], [4, 4, 11], 
       [2, 4, 14], [2, 5, 13],
       [3, 3, 14], [4, 4, 12], [6, 5, 9], [6, 6, 8],
       [2, 4, 15], [2, 5, 14], [2, 6, 13], [3, 3, 15],
       [3, 5, 13], [3, 6, 12], [4, 4, 13], [4, 6, 11], [5, 4, 12], [5, 5, 11],
       [5, 6, 10], [5, 7, 9], [5, 9, 7], [6, 5, 10], [6, 6, 9], [7, 4, 10],
       [8, 4, 9], [9, 4, 8], [10, 2, 9]]


# Minimum number of active Sboxes in output after the linear operation
# For ex. LIN[1] = 4 means (1 + 1) active Sboxes at input and at least 4 active Sboxes at the output
LIN = [ 3, 4, 3, 4, 4, 4, 3, 4,
        3, 2, 3, 4, 3, 4, 3, 2,
        3, 2, 3, 2, 3, 2, 3, 3,
        3, 2, 3, 2, 3, 2, 1, 2,
        1, 2, 1, 2, 3, 2, 2, 2,
        2, 2, 2, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 2, 1, 1,
        1, 1, 1, 1, 2, 1, 1, 2]


#Step 1: Only add those (a, b, c) which are not in ZZF
ZZ1 = []
with open(in_file, 'r') as fc:
    LL = fc.readlines()
    for r in range(len(LL)-1):
        L = [int(i) for i in LL[r].split()]
        if([L[0], L[1], L[2]] not in ZZF and [L[1], L[2], L[3]] not in ZZF):
            ZZ1.append(L)
            



#Step 2: Filter based on ZZT_21        
ZZ2 = []
for i in range(len(ZZ1)):
    T = ZZ1[i]
    if(T[0] + T[1] + T[2] <= 21 and T[1] + T[2] + T[3] <= 21):
        if([T[0], T[1], T[2]] in ZZT_21 and [T[1], T[2], T[3]] in ZZT_21):
            ZZ2.append(T)
    elif(T[0] + T[1] + T[2] <= 21 and T[1] + T[2] + T[3] >= 22):
        if(([T[0], T[1], T[2]] in ZZT_21)):
            ZZ2.append(T)
    elif(T[0] + T[1] + T[2] >= 22 and T[1] + T[2] + T[3] <= 21):
        if(([T[1], T[2], T[3]] in ZZT_21)):
            ZZ2.append(T)
    elif(T[0] + T[1] + T[2] >= 22 and T[1] + T[2] + T[3] >= 22):
        ZZ2.append(T)       
    


# Step 3: Filter based on LIN
ZZ3 = []
for i in range(len(ZZ2)):
    T = ZZ2[i]
    if(T[1] >= LIN[T[0] - 1]):
        if(T[2] >= LIN[T[1] - 1]):
            if(T[3] >= LIN[T[2] - 1]):
                ZZ3.append(T)
    

# Step 4: Filter based on: 
    # (1) 2 active sboxes at input will give at least 45 active Sboxes for 4 rounds
    # (2) Sum of active sboxes in round (0,2) and (1, 3)
ZZ4 = []
for i in range(len(ZZ3)):
    T = ZZ3[i]
    if(T[0] == 2):
        continue
    elif(T[1] == 2 and T[3]<=9):
        continue
    elif((T[0] == 3 and T[2]<=11) or (T[1] == 3 and T[3]<=11)):
        continue
    else:
        ZZ4.append(T)

#Step 5: Filter based on maximum number of active Sboxes for a linear layer
ZZ5 = []
for i in range(len(ZZ4)):
    T = ZZ4[i]
    if(T[0] == 3 and T[1] >= 33):
        continue
    elif((T[1] == 2 and T[2] >= 23) or (T[1] == 3 and T[2] >= 33)):
        continue
    elif((T[2] == 2 and T[3] >= 23) or (T[2] == 3 and T[3] >= 33)):
        continue
    else:
        ZZ5.append(T)



ZZ5.sort()
for z in ZZ5:
    print z
print(len(ZZ5))

     
     