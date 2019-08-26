import numpy as np

def binary_search(s,arr):
    '''
    Locate s in arr.
    s: int.
    arr: int[].
    '''
    l=0
    r=len(arr)-1
    while l<r:
        mid=(l+r)//2
        if arr[mid]<s:
            l = mid+1
            continue
        if arr[mid]>=s:
            r = mid-1
    return (l+r)//2


