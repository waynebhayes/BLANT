import random
import numpy
#import fastrand


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
   
    '''
    def find_by_name(self, value, node_info, update = None): #get node with the given value
        if update == None:
            update = self.updateList(value) #path to the given value
        if len(update) > 0:
            candidate = update[0].next[0] #bottom node of the entire route
            if candidate != None:
                print("candinfo: " + str(candidate.info))
                print("candval: " + str(candidate.value))
                if candidate.value == value and candidate.info == node_info:
                    return candidate
    '''

    def print_list(self):
        for i in range(self.maxHeight-1, -1, -1):
            print("level:",i, end = " ")
            curr = self.head.next[i]
            while curr != None:
                print(curr.info, end = "->")
                curr = curr.next[i]
            print()


    def find_by_name(self, value, node_info, update = None): #get node with the given value
        if update == None:
            update = self.updateListInfo(value, node_info) #path to the given value
        if len(update) > 0:
            candidate = update[0].next[0] #bottom node of the entire route
            #found = None
            if candidate != None and candidate.value == value and candidate.info == node_info:
                return candidate
            return None


    def __contains__(self, value): #check whether a value in this list
        return self.find(value) != None

    def flip_coin(self): #determine the height of the new node
        count = 1
        while random.randint(0,1)!=0: #if toll head, have one more hierarchy
        #while numpy.random.randint(0,2)!=0: #if toll head, have one more hierarchy
        #while fastrand.pcg32bounded(2) != 0:
            count+=1
        return count

    def updateListInfo(self, value:float, info:tuple): #get a route to the node
        update = [None]*self.maxHeight #save the route
        x = self.head #from the left headnode
        for i in reversed(range(self.maxHeight)): #from top to bottom
            #while x.next[i] != None and x.next[i].value < value: #next_value>=given_value => walk down
            #while x.next[i] != None and x.next[i].value < value and x.next[i].info < info: #next_value>=given_value => walk down
            while x.next[i] != None and x.next[i].value < value: 
                x = x.next[i] 
            while x.next[i] != None and x.next[i].value == value and x.next[i].info < info: #next_value>=given_value => walk down
                x = x.next[i] #go to next value
            update[i] = x #save the route
            # print("Update list info [level",i,"]:", x.info, info, value)
        return update
        
    def updateList(self, value:float): #get a route to the node
        update = [None]*self.maxHeight #save the route
        x = self.head #from the left headnode
        for i in reversed(range(self.maxHeight)): #from top to bottom
            while x.next[i] != None and x.next[i].value < value: 
                x = x.next[i] 
            update[i] = x #save the route
        return update

    def add(self,input_tuple, debug=False):
        #self.switch determines min_heap(1) or max_heap(-1)
        value=input_tuple[0]*self.switch
        info=input_tuple[1]
        if debug:
            print("adding " + str(info) + ": " + str(value))
        #flip coin to get the height of the new node
        node = SkipNode(self.flip_coin(), value , info)
        #update skip list if has new highest: maxHeight, head_node
        self.maxHeight = max(self.maxHeight, len(node.next))
        if len(self.head.next) < self.maxHeight:
            self.head.next.extend([None]*(self.maxHeight- len(self.head.next)))
        #find route to the given value
        update = self.updateListInfo(value, info)
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
        # x.info means removed, 0 means cannot find the value
            return x.info
        return 0
    
    def remove_by_name(self, value, node_info, update=None):
        #find the route to the given value
        value = value*self.switch 
        #print("removing by " + str(node_info))
        #print("removing by " + str(value))
        if update == None:
            update = self.updateListInfo(value, node_info)
        #whether the value exist, the route can accelerate the find progress
        # print(update[0].info)
        x = self.find_by_name(value, node_info, update)
        #if find the value
        if x != None:
            for i in reversed(range(len(x.next))):
                update[i].next[i] = x.next[i]
                #check whether this node is the only node in this rank
                if self.head.next[i] == None:
                    #print("Decreasing height")
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

    def updateListEnd(self, value:float): #get a route to the node
        update = [None]*self.maxHeight #save the route
        x = self.head #from the left headnode
        for i in reversed(range(self.maxHeight)): #from top to bottom
            while x.next[i] != None and x.next[i].value <= value:
                x = x.next[i]
            update[i] = x #save the route
        return update

    def findEnd(self, value, update=None):
        if update == None:
            update = self.updateListEnd(value) #path to the given value
        if len(update) > 0:
            candidate = update[0]
            if candidate != None:
                if candidate.value == value:
                    return candidate
        return None #cannot find node with this value

    def removeRand(self, value, start=None, end=None):
        # find the start nodes and the end nodes that have such value
        if start == None:
            start = self.updateList(value)[0].next[0]
        if end == None:
            end = self.updateListEnd(value)[0]
        # print("nodes range for value", value, "is [", start.info, end.info, "]")
        # assert start.value == end.value, "start and end have different value!"+ str(start.value) + str(end.value)
        # random choose one info
        g1_info = random.randint(min(start.info[0], end.info[0]), max(start.info[0], end.info[0]))
        g2_info = random.randint(min(start.info[1], end.info[1]), max(start.info[1], end.info[1]))
        # print("Random node info:", (g1_info,g2_info))
        # get the nearest path to that node
        update = self.updateListInfo(value, (g1_info, g2_info))
        # TODO: low probability to choose the last one?
        # cand1 = update[0]
        # cand2 = update[0].next[0]
        # if cand1 and cand2 and g1_info + g2_info - cand1.info[1] -cand2.info[2] < cand1.info[1]

        cand1 = update[0] if update[0] and update[0].info != None else update[0].next[0]
        # print("candidate for removeRand", cand1.info, cand1.value)
        self.remove_by_name(cand1.value*self.switch, cand1.info)
        return (cand1.value, cand1.info)

    def pop(self, domain=0.1, fast=True):
        if fast:
            # random value, but the first one with that value
            return self.popfast(domain)
        else:
            return self.popUniform(domain)


    def popUniform(self, domain=0.1):
        # check not empty
        if self.__len__() == 0:
            raise IndexError()
        # get a random value in the given random value range
        rand_value = random.uniform(self.head.next[0].value, self.head.next[0].value + domain)
        cand1 = self.updateList(rand_value)[0]
        cand2 = self.updateListEnd(rand_value)[0]
        # print(cand1.value, cand2.value)
        if cand1 and cand2 and cand1.next[0] and cand1.next[0].value == cand2.value:
            ret = self.removeRand(cand1.value, start=cand1, end=cand2)
        elif cand2 and abs(cand2.value-rand_value) < abs(cand1.value-rand_value):
            ret = self.removeRand(cand2.value, start=cand2, end=None)
        else:
            ret = self.removeRand(cand1.value, start=None, end=cand1)
        return ret


    def popfast(self, domain=0.1):
        '''Pop out the first node with a random value'''
        #check not empty
        if self.__len__()==0:
            raise IndexError()
        #get a random value in the given random value range
        #rand_value=random.uniform(self.head.next[0].value,self.head.next[0].value+domain)
        rand_value=random.uniform(self.head.next[0].value, self.head.next[0].value+domain)
        #if this value exist, directly remove it
        if self.__contains__(rand_value):
            info = self.remove(rand_value)
            return (rand_value * self.switch, info)
        #if this value not in list, remove the nearest one
        else:
            update = self.updateList(rand_value)
            #get the left one & right one, compare them
            candidate1 = update[0]
            candidate2 = update[0].next[0]
            #if the right one not exist, return the left one
            if candidate2==None:
                info = self.remove(candidate1.value)
                return (candidate1.value*self.switch,info)
            #otherwise compare them
            if abs(candidate1.value-rand_value)<=abs(candidate2.value-rand_value):
                info = self.remove(candidate1.value)
                return (candidate1.value*self.switch,info)
            if abs(candidate2.value-rand_value)<abs(candidate1.value-rand_value):
                info = self.remove(candidate2.value)
                return (candidate2.value*self.switch,info)
            return None

    def opop(self, domain=0.1):
        '''Reservoir Sampling: slow!'''
        #that's the new one
        if self.__len__()==0:
            raise IndexError()
        tail=self.head.next[0].value+domain
        temp=self.head.next[0]
        choice=(temp.value,temp.info)
        n=0
        while temp and temp.value<tail :
            if temp.next[0]==None:
                break
            temp=temp.next[0]
            n+=1
            
            #j=random.randint(0,n)
            j=numpy.random.randint(0,n+1)
            if j==0:
                choice=(temp.value,temp.info)
        self.remove(choice[0])
        return (choice[0]*self.switch,choice[1])
