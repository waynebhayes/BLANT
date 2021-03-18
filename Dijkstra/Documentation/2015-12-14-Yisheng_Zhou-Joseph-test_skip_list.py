#Yisheng Zhou, 12.11.15
from skip_list import *
import random

#initial tested list
t1=SkipList( ) #max priority queue
t2=SkipList(1) #min priority queue
number=1000   #test time:number of operations

#test length t1
true_length=0
i=0
while i<number:
    i+=1
    operation=random.randint(-1,1) #1:add new node; -1:remove node; 0:do nothing      
    if operation==1:
        t1.add((random.uniform(0,1),[]))
        true_length=true_length+1  
    elif operation==-1:
        popout=t1.pop()
        if (popout!=None):
            true_length=true_length-1
    assert(true_length==len(t1))
print("test max pq length done")

#test length t2
true_length=0
i=0
while i<number:
    i+=1
    operation=random.randint(-1,1) #1:add new node; -1:remove node; 0:do nothing      
    if operation==1:
        t2.add((random.uniform(0,1),[]))
        true_length=true_length+1  
    if operation==-1:
        popout=t2.pop()
        if (popout!=None):
            true_length=true_length-1
    assert(true_length==len(t2))
print("test min pq length done")

#test order t1
while random.randint(0,number)!=0:
    t1.add((random.uniform(0,1),[]))
i=0
list_length=len(t1)
while i<number:
    i=i+1
    index1=random.randint(0,list_length-1)
    index2=random.randint(0,list_length-1)
    if index1<index2:
        assert(t1[index1][0]>=t1[index2][0])
    if index1>index2:
        assert(t1[index1][0]<=t1[index2][0])
print("test max pq order done")

#test order t2
while random.randint(0,number)!=0:
    t2.add((random.uniform(0,1),[]))
i=0
list_length=len(t2)
while i<number:
    i=i+1
    index1=random.randint(0,list_length-1)
    index2=random.randint(0,list_length-1)
    if index1<index2:
        assert(t2[index1][0]<=t2[index2][0])
    if index1>index2:
        assert(t2[index1][0]>=t2[index2][0])
print("test min pq order done")
       
#test maxheight~=log2(length)
import math
t3=SkipList()
i=0
sum_different=0
count=0a
while i<number:
    i+=1
    operation=random.randint(-1,1) #1:add new node; -1:remove node; 0:do nothing      
    if operation==1:
        t3.add((random.uniform(0,1),[]))
    elif operation==-1:
        t3.pop()
    if len(t3)>1:
        sum_different=sum_different+abs(t3.maxHeight-math.log(t3.maxHeight)/math.log(len(t3)))
        count+=1
print("avg difference between log(n) and height:", sum_different/count)








