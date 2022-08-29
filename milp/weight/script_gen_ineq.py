# Copyright 2022 Rusydi H. Makarim and Raghvendra Rohit
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#Minimized Product of Sums:
W2 = "(A+C)(D'+J')(B+C'+H)(C'+E')(B'+C)(C'+H'+I)(A'+C')(D'+E)(B'+I')(C+E+H')(D'+H+I)(D'+H'+I')(C+D+F+J)(C'+I'+J')(E'+F'+J)(C+E+G)(E'+F+J')(C+F'+G'+J')(B'+F'+G'+J)(B'+F+G+J)(C'+F+G'+J')(E+F'+G+J')"

#Minimized Product of Sums:
W3 = "(A'+D'+E)(B'+C+D')(A+B'+E'+G)(A'+B+C)(B+D+E+F')(A+B'+D+E')(B+C+F+J)(A'+C+E)(A+C'+D+F)(A+D+E+H)(B'+D+H'+I+J)(B+D+G+H+J)(A+B+C'+D'+I)(A'+B+D+E'+J)(B'+E+H'+I+J')(A'+B'+D'+G')(B+C+E+F)(A'+B+D'+J')(B'+C'+D+E+H')(B'+E'+H'+I'+J')(B'+D'+G'+I')(B+C+E'+F'+J')(B+C'+D'+E)(A+D'+E'+F'+I')(B+D+F+G+H)(B'+D+E'+H+I'+J)(B'+C'+D+F'+G+J')(B'+C'+D+F+G+J)(B'+D+F+G'+H+J')(C'+D+F'+G'+H'+J)(B+C+D+G)(B'+E'+G+H+I+J')(A'+D'+H'+I')(A'+D'+H+I)(C'+D'+E+G')(C+D'+E+J)(C+E'+F'+G+I')(A'+E'+F'+G'+H+I)(C'+D+E'+F+G'+H'+J')(B+E'+F'+G'+H+J')(B+D+F'+G+H'+J')(C+D+F+G+I)(B+E+F+H)(C+D'+F+I')(E+F'+G'+H+J)(A+B'+H+I'+J')(D'+E+H'+I'+J)(C'+D'+G+H+I+J)(A'+C+F'+G'+I)(A'+C+F+G'+I')"

#Minimized Product of Sums:
W4 = "(C'+F+G+J)(D+F+G+H+J)(D+F+G'+H+J')(D+F'+G+H+J')(B+D)(D+F'+G'+H+J)(B+E')(A'+B+C'+F)(A'+D+E')(A+B+C)(B+F+J)(A+D+E)(A+D'+E'+I)(B+C+F'+J')(A'+E'+H'+I')(A'+B'+D'+E+G+I+J)(A+E+H'+I'+J)(A+B'+C'+D')(C'+D'+E')(B'+C'+E+G)(D+E+H)(C'+D+E)(D'+E'+H+I)(A+C+E+H'+I+J')(A+C+E+H+I'+J')(A'+B'+C+D'+E+G'+I+J')(A+E+H+I+J)(A'+B'+D'+E+G+I'+J')(A'+B'+C+D'+E+G'+I'+J)(A+C'+E+J)(C+D+E'+F'+G'+H'+J')(C+D+E'+F+G+H'+J')(C+D+E'+F+G'+H'+J)(C+D+E'+F'+G+H'+J)(C'+E'+F+G'+J')(C'+E'+F'+G+J')(C'+E'+F'+G'+J)"

def sp(X):
    Y = X.replace("(", "")
    Y = Y.split(")")
    return Y

def check(x, y):
    a = 0
    if x == y:
        a = 1
    elif x == y + str("'"):
        a = -1
    return a

def str_eq(X):
    A = [0]*10
    for x in X:
        A[0] |= check(x, 'A')
        A[1] |= check(x, 'B')
        A[2] |= check(x, 'C')
        A[3] |= check(x, 'D')
        A[4] |= check(x, 'E')
        A[5] |= check(x, 'F')
        A[6] |= check(x, 'G')
        A[7] |= check(x, 'H')
        A[8] |= check(x, 'I')
        A[9] |= check(x, 'J')
    return A

P = sp(W2)
for i in range(len(P)-1):
    temp = P[i].split("+")
    A = str_eq(temp)
    print("{}".format(A))

print(len(P)-1)
