
// Node definition
struct SkipListNode
{
    SkipListNode[] next;
    size_t height;
    float value;
    void** info;
}

// Skip List definition
struct SkipList
{
    SkipListNode* head;
    size_t len;
    size_t maxHeight;
    int switch;
}

// Node initializer
SkipListNode* node_init(size_t height = 0, float value = NULL, void** info= NULL)
{
    struct SkipListNode self;
    self.next = (SkipListNode*)calloc(height, sizeof(SkipListNode*));
    self.height = height;
    self.value = value;
    self.info = info;
    return &self;
}

// Skip List Initializer
SkipList* list_init(bool min)
{
    struct SkipList self;
    self.head = node_init();
    self.len = 0;
    self.maxHeight = 0;
    self.switch = min ? -1 : 1;
    return &self;
}

// coin flip to determine height of new node 
// return height of new node
// TODO: better random generator (?)
int flip_coin()
{
    int height = 1;
    srand(time(NULL));
    while(rand()%2 != 0)
    {
        height++;
    }
    return height;
}

// find route to node, by value, through all levels 
// return array of last nodes at each level such that node.value < value at each level
SkipListNode** updateList(SkipList* self, float value)
{
    SkipListNode* path[];
    SkipListNode* x = self->head;
    for(int i = self->maxHeight-1; i>-1; i--)
    {
        while(x->next[i] != NULL && x->next[i]->value < value)
        {
            x = x->next[i];
        }
        path[i] = x;
    }
    return path;
}

// find a node by value 
// uses updateList (<) to find path
// return pointer to node
SkipListNode* find(SkipList* self, float value, SkipListNode** update = NULL)
{
    if(update == NULL)
    {
        update = updateList(self, value);
    }
    candidate = update[0]->next[0];
    if(candidate != NULL)
    {
        if(candidate->value == value)
        {
            return candidate;
        }
    }
    return NULL;
}

// check if value is in list
// use find()
// return true if found, false otherwise
bool __contains__(SkipList* self, float value)
{
    return (bool)(self.find(value) != NULL);
}


// insert a node with the given value and info
// uses updateListInfo
// return 0 on success, -1 on failure
int add(SkipList* self, float value, void** info)
{
    SkipList* node = node_init(flip_coin(), value, info);
    self->maxHeight = max(self->maxHeight, node->height);
    // if self->head->height < self->maxHeight, extend the head
    update = self.updateListInfo(value);
    for(int i = 0; i<node->height; i++)
    {
        node->next[i] = update[i]->next[i];
        update[i]->next[i] = node;
    }
    self->len++;
    return 0;
}

// remove node with value from the list
// uses updateList
// return 0 on failure, info on success
void** remove(self, value):
{
    x = find(self, value, update);
    if(x != NULL)
    {
        for(int i = x->height-1; i > -1; i--)
        {
            update[i]->next[i] = x->next[i];
            if(self->head->next[i] == NULL)
            {
                self->maxHeight--;
            }
        }
        self->len--;
        return x->info;
    }
    return 0;
}


// find route to node, by value, through all levels
// return array with last nodes at each level such that node.value <= value at each level
SkipListNode** updateListEnd(SkipList* self, float value)
{
    SkipListNode* path[];
    SkipListNode* x = self->head;
    for(int i = self->maxHeight-1; i>-1; i--)
    {
        while(x->next[i] != NULL && x->next[i]->value <= value)
        {
            x = x->next[i];
        }
        path[i] = x;
    }
    return path;
}

// find a node by value 
// uses updateListEnd (<=) to find path
// return pointer to node
SkipListNode* findEnd(SkipList* self, float value, SkipListNode** update=NULL)
{
    if(update == NULL)
    {
        update = updateListEnd(self, value);
    }
    candidate = update[0]->next[0];
    if(candidate != NULL)
    {
        if(candidate->value == value)
        {
            return candidate;
        }
    }
    return NULL;
}


// find route to node, by value, through all levels 
// return array with last node s.t. node.value<value and node.info<info at each level
SkipListNode** updateListInfo(SkipList* self, float value, void** tuple)
{
    SkipListNode* path[];
    SkipListNode* x = self->head;
    for(int i = self->maxHeight-1; i>-1; i--)
    {
        while(x->next[i] != NULL && x->next[i]->value < value)
        {
            x = x->next[i];
        }
        while(x->next[i] != NULL && x->next[i]->value == value && x->next->info<info)
        {
            x = x->next[i];
        }
        path[i] = x;
    }
    return path;
}

// remove node with value and info from the list
// return 1 on success, 0 if node not found
int remove_by_name(SkipList* self, float value, void** node_info, SkipListNode** update=None):
{
    float value = value*(self->switch);
    if(update == NULL)
    {
        update = updateListInfo(self, value, node_info);
    }
    SkipListNode* x = find_by_name(self, value, node_info, update);
    if(x != NULL)
    {
        for(int i = x->next->height-1; i > -1; i--)
        {
            update[i]->next[i] = x->next[i];
            if(self->head->next[i] == NULL)
            {
                self->maxHeight--;
            }
        }
        self->len--;
        return 1;
    }
    return 0;
}

// find a node by info 
// (return pointer to node)
SkipListNode* find_by_name(SkipList* self, float value, void** node_info, SkipListNode** update = NULL)
{
    if(update == NULL)
    {
        update = updateListInfo(self, value, node_info);
    }
    candidate = update[0]->next[0];
    if(candidate != NULL && candidate->value == value && candidate->info == node_info)
    {
        return candidate;
    }
    return NULL;
}

// print out full list
// return nothing
void print_list(SkipList* self)
{
    for(int i = self->maxHeight-1, i>-1, i--)
    {
        printf("level: %d", i);
        curr = self->head->next[i];
        while(curr != NULL)
        {
            printf("%d -> ", curr->value);
            // printf(curr->info);
            curr = curr->next[i];
        }
        printf("\n");
    }
}

// return node at index in list
// NULL if not found
SkipListNode* __getitem__(SkipList* self, int index):
{
    assert(index>=0);
    assert(index<self->len);
    SkipListNode* x = self->head;
    int i = -1;
    while(true)
    {
        if(i==index)
        {
            return x;
        }
        if(x->next[0] == NULL)
        {
            break;
        }
        i++;
        x = x->next[0];
    }
    return NULL;
}

// remove random node with given value
// uses updateList and updateListEnd to find upper and lower bounds
// return value and info of removed node 
// (just info for now until I figure out a good way to do this)
void** removeRand(SkipList* self, float value, SkipListNode** start=NULL, SkipListNode** end=NULL)
{
    if(start == NULL)
    {
        start = updateList(self, value)[0]->next[0];
    }
    if(end == NULL)
    {
        end = updateListEnd(self, value)[0]
    }
    srand(time(NULL));
    int max0 = max(start->info[0], end->info[0]);
    int min0 = min(start->info[0],end->info[0]);
    int diff0 = max0-min0;
    int max1 = max(start->info[1], end->info[1]);
    int min1 = min(start->info[1],end->info[1]);
    int diff1 = max1 - min1;
    int g1_info = min0 + rand()%(diff0);
    int g2_info = min1 + rand()%(diff1);
    int[] g_info = {g1_info, g2_info};

    SkipListNode** update = updateListInfo(self, value, g_info);

    SkipListNode* = (update[0] && update[0]->info != NULL) ? update[0] : update[0]->next[0];

    remove_by_name(self, (self->switch)*(cand1->value), cand1->info);
    return cand1->info;
}


// if fast, popfast(domain)
// else, popUniform(domain)
// return whatever the relevant function returns
// (going to make other pop functions return info for now)
void** pop(self, domain=0.1, fast=True):
{
    if(fast)
    {
        return self.popfast(domain);
    }
    else
    {
        return self.popUniform(domain);
    }
}

// remove a random node in the range (node[0].value, node[0].value+domain)
// use removeRand once the random value has been selected
// return the same as removeRand
void** popUniform(SkipList* self, float domain=0.1):
{
    assert(self->len > 0);
    srand(time(NULL));
    float rand_value = self->head->next[0]->value+ domain*(float)rand()/(float)(RAND_MAX);
    SkipListNode* cand1 = updateList(self,rand_value)[0];
    SkipListNode* cand2 = updateListEnd(self,rand_value)[0];

    if (cand1 && cand2 && cand1->next[0] && cand1->next[0]->value == cand2->value)
    {
        void** ret = removeRand(self, cand1->value, start = cand1, end= cand2);
    }
    else
    {
        if(cand2 && abs(cand2->value-rand_value) < abs(cand1->value-rand_value))
        {
            void** ret = removeRand(self, cand2->value, start = cand2, end=None);
        }
        else
        {
            void** ret = self.removeRand(cand1.value, start = None, end = cand1);
        }
    }
    return ret;
}

// pop first one that works for a randomly generated value
void** popfast(SkipList* self, float domain=0.1):
{
    assert(self->len > 0);
    srand(time(NULL));
    float rand_value = self->head->next[0]->value+ domain*(float)rand()/(float)(RAND_MAX);
    if(__contains__(self, rand_value))
    {
        void** info = remove(self, rand_value);
        return info;
    }
    else
    {
        SkipListNode** update = updateList(self, rand_value);
        SkipList* candidate1 = update[0];
        SkipList* candidate2 = update[0]->next[0];

        if(candidate2 == NULL)
        {
            void** info = remove(self, candidate1->value);
            return info;
        }
        if(abs(candidate1->value-rand_value)<=abs(candidate2->value-rand_value))
        {
            void** info = remove(self, candidate1->value);
            return info;
        }
        if(abs(candidate2.value-rand_value)<abs(candidate1.value-rand_value))
        {
            void** info = remove(self, candidate2->value);
            return info;
        }
        return NULL;
    }
}

// another removal algorithm
void** opop(SkipList* self, float domain=0.1):
{
    assert(self->len > 0);
    float tail = self->head->next[0]->value + domain;
    SkipListNode* temp = self->head->next[0];
    float choice_value = temp->value;
    void** choice_info = temp->info;

    int n = 0;
    
    srand(time(NULL));
    while(temp && temp->value < tail)
    {
        if(temp->next[0] == NULL)
        {
            break;
        }
        temp = temp->next[0];
        n++;

        int j = rand()%(n+2);
        if(j==0)
        {
            choice_value = temp->value;
            choice_info = temp->info;
        }
    }
    remove(self, choice_value);
    return choice_info;
}

// TODO: figure out a good way to do this in C
// return list as a string
// treat as regular linked list:
// only worry about bottom level
char* __str__(SkipList* self):
{
    return NULL;
}


int main()
{
    return 0;
}