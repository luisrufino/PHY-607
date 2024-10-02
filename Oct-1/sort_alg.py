import numpy as np
import matplotlib.pyplot as plt


def sort(list1):
    list2 = np.empty(len(list1))
    for i in range(len(list1)):
        for j in range(len(list1)):
            if list1[i] > list1[j]:
                list2[i] = (list1[j])
                list2[j] = list1[i]
    return list2


def bubble(list1):
    for n in range(len(list1) - 1, 0, -1):
        for i in range(n):
            if list1[i] > list1[i + 1]:
                list1[i], list1[i + 1] = list1[i + 1], list1[i]
    return list1


def heap_con(list1, n, i):
    max1 = i
    left = 2 * i + 1
    right = 2 * i + 2

    if left < n and list1[i] < list1[left]:
        max1 = left
    if right < n and list1[max1] < list1[right]:
        max1 = right
    if max1 != i:
        (list1[i], list1[max1]) = (list1[max1], list1[i])

        heap_con(list1, n, max1)


def heap_sort(list1):
    n = len(list1)

    for i in range(n // 2, -1, -1):
        heap_con(list1, n, i)

    for i in range(n - 1, 0, -1):
        (list1[i], list1[0]) = (list1[0], list1[i])
        heap_con(list1, i, 0)

    return list1


list1 = [11, 798, 32, 9, 3, 2, 6, 80, 4, 3, 1, 7, 5]
print(list1)
list2 = heap_sort(list1)
list3 = bubble(list1)
print(list2, list3)