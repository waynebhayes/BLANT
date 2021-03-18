#Yisheng Zhou, 12.11.15
import random

class SkipNode:   
    def __init__(self, height = 0, value = None, info = None):
        self.next = [None]*height #list of pointers
        self.value = value #value for list sorting
        self.info = info #other information


class SkipList:
    def __init__(self,min_switch=0): #initial list
        self.head = SkipNode() #first node, no information
        self.len = 0 #list size
        self.maxHeight = 0 #max height
        self.switch=-1 #originally max priority queue
        if min_switch: #switch: change to min priority queue
            self.switch=1
        
    def __len__(self): #get the list size
        return self.len
    
    def find(self, value, update = None): #get node with the given value
        if update == None:
            update = self.updateList(value) #path to the given value
        if len(update) > 0:
            candidate = update[0].next[0] #bottom node of the entire route
            if candidate != None:
                if candidate.value == value:
                    return candidate
        return None #cannot find node with this value
    
    def __contains__(self, value): #check whether a value in this list
        return self.find(value, update=None) != None 

    def flip_coin(self): #determine the height of the new node
        count = 1
        while random.randint(0,1)!=0: #if toll head, have one more hierarchy
            count+=1
        return count

    def updateList(self, value): #get a route to the node
        update = [None]*self.maxHeight #save the route
        x = self.head #from the left headnode
        for i in reversed(range(self.maxHeight)): #from top to bottom
            while x.next[i] != None and x.next[i].value < value: #next_value>=given_value => walk down
                x = x.next[i] #go to next value
            update[i] = x #save the route
        return update
        
    def add(self,input_tuple):
        #self.switch determines min_heap(1) or max_heap(-1)
        value=input_tuple[0]*self.switch
        info=input_tuple[1]
        #flip coin to get the height of the new node
        node = SkipNode(self.flip_coin(), value , info)
        #update skip list if has new highest: maxHeight, head_node
        self.maxHeight = max(self.maxHeight, len(node.next))
        if len(self.head.next) < self.maxHeight:
            self.head.next.extend([None]*(self.maxHeight- len(self.head.next)))
        #find route to the given value
        update = self.updateList(value)
        #update every node in the given route
        for i in range(len(node.next)):
            #every rank point to the next
            node.next[i] = update[i].next[i]
            #update the left node
            update[i].next[i] = node
        self.len += 1

    def remove(self, value):
        #find the route to the given value
        update = self.updateList(value)
        #whether the value exist, the route can accelerate the find progress
        x = self.find(value, update)
        #if find the value
        if x != None:
            for i in reversed(range(len(x.next))):
                #update every rank, change the left node pointer
                update[i].next[i] = x.next[i]
                #check whether this node is the only node in this rank
                if self.head.next[i] == None:
                    #if it is, decrease the height of the skiplist
                    self.maxHeight -= 1
            #change list size
            self.len -= 1
        #1 means removed, 0 means cannot find the value
            return 1
        return 0
    
    def __str__(self):
        #return the list as a string
        result=[]
        x = self.head
        #append every node in the lowest rank
        while x.next[0] != None:
            x = x.next[0]
            result.append((x.value*self.switch,x.info))
        return str(result)

    def __getitem__(self,index):
        #return value by index
        #very slow, like find items by index in linked list
        #check index
        assert(index>=0)
        assert(index<self.len)
        #iteration find node
        x = self.head
        i = -1
        while True:
            if i==index:
                return (x.value*self.switch,x.info)
            if x.next[0] == None:
                break
            i=i+1
            x = x.next[0]
        return None

    def pop(self,domain=0.1):
        #check not empty
        if self.__len__()==0:
            return None
        #get a random value in the given random value range
        rand_value=random.uniform(self.head.next[0].value,self.head.next[0].value+domain)
        #if this value exist, directly remove it
        if self.__contains__(rand_value):
            temp=self.find(rand_value)
            self.remove(rand_value)
            return (temp.value*self.switch,temp.info)
        #if this value not in list, remove the nearest one
        else:
            update = self.updateList(rand_value)
            #get the left one & right one, compare them
            candidate1 = update[0]
            candidate2 = update[0].next[0]
            #if the right one not exist, return the left on
            if candidate2==None:
                self.remove(candidate1.value)
                return (candidate1.value*self.switch,candidate1.info)
            #otherwise compare them
            if abs(candidate1.value-rand_value)<=abs(candidate2.value-rand_value):
                self.remove(candidate1.value)
                return (candidate1.value*self.switch,candidate1.info)
            if abs(candidate2.value-rand_value)<abs(candidate1.value-rand_value):
                self.remove(candidate2.value)
                return (candidate2.value*self.switch,candidate2.info)
            return None

























        
