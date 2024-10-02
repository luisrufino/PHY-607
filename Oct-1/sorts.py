import numpy as np
import time

a = [1,2,3,34,4,3]
def sort_list(A):
    for n in range(len(A) -1, 0, -1):
        for i in range(n):
            if A[i] > A[i+1]:
                A[i], A[i+1] = A[i+1], A[i]
    return A


def head_sort(A):
    N = len(A)
    tmps =[]
    for i in range(0,N,2):
        curr = A[i]
        next = A[i+1]
        if curr > next:
            tmp = [next, curr]
            tmps.append(tmp)
        if next > curr:
            tmp = [curr, next]
            tmps.append(tmp)
        else:
            tmp = [curr, next]
            tmps.append(tmp)




    return A
# print(f"OG list: {a}")
# sorted_a = sort_list(a)
# print(f"Sorted list: {sorted_a}")

sorted_a = head_sort(a)
print(f"Sorted list: {sorted_a}")

