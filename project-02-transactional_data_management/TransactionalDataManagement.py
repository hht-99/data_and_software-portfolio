
#We can store the base state and each transaction’s temporary changes.

#store = {}             # base committed state
#transactions = []      # stack of temporary change sets

#Each BEGIN pushes a new dictionary on the stack.
#Each SET adds key–value pairs to the current transaction (top of stack).
#Each COMMIT merges top transaction into the one below (or the base).
#Each ROLLBACK pops the last transaction off the stack.

#Example:

class KeyValueStore:
    def __init__(self):
        self.store = {}                 #base database where commited values stored
        self.transactions = []          #stack (list of changes to make)

#each transaction adds temp layer of changes on top
#the bottom layer is permanent, each BEGIN adds a new temp layer on top
#rollbacks delete the top layer
#commit merges from the bottom up, meaning from the oldest to the latest

    def begin(self):                    #this starts a transaction by creating an empty dict and adding to the stack
        self.transactions.append({})    #new dict will hold any changes until commit/rollaback (multiple BEGINS = nested transactions)

    def set(self, key, value):                  #adds or updates a value, when inside a transaction writes to the top-most layer - not base
        if self.transactions:                   #outside of a transaction this writes directly to base (commited store)
            self.transactions[-1][key] = value
        else:
            self.store[key] = value

    def get(self, key):                             #retrieves a value              
        for layer in reversed(self.transactions):   #looks through the transactions, newest to oldest, checking for the key specified   
            if key in layer:                        #if found, reurns the most recent version
                return layer[key]                   #if not found in transaction, returns the base (commited store)
        return self.store.get(key, None)            #allows you to see your latest change even if not commited

    def rollback(self):                 #Undoes the last transaction, removing its changes
        if self.transactions:           #if no transaction exists then prints a warning
            self.transactions.pop()     #crucially doesn't effect commited store
        else:
            print("NO TRANSACTION")

    def commit(self):                           #applies changes permanently, takes all transaction layers and merges them in base (self.store)
        while self.transactions:                #this takes the oldest transaction from the stack and merges it with the base
            changes = self.transactions.pop(0)  #then it removes the oldest and moves onto the next until it reaches the latest
            self.store.update(changes)

#Possible to have partial commits that merge inner layers with the outer layers (the latest transaction with the previous)
#This doesn't effect the base or update it

    def partialcommit(self):
        if not self.transactions:
            print("NO TRANSACTION")
            return
        if len(self.transactions) == 1:
            self.store.update(self.transactions.pop())
        else:
            top = self.transactions.pop()
            self.transactions[-1].update(top)

#Practice Questions:

#Prompt 1:
    #Create a KeyValueStore class that supports:
    #set(key, value)
    #get(key)

    #You should be able to:
        #store = KeyValueStore()
        #store.set("x", 10)
        #print(store.get("x"))

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def set(self, key, value):
        self.store[key] = value
    def get(self, key):
        if self.store:
            return self.store[key]

store = KeyValueStore()
store.set("x", 10)
print(store.get("x"))

#Prompt 2:
    #Extend your store to support begin(), commit(), and rollback() methods.
    #When a transaction is active, changes should be temporary until committed.

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    #adds empty dict to the transactions stack
        self.transactions.append({})
    def set(self, key, value):          #modifies the latest dict in the stack
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
    def get(self, key):
        if self.store:
            return self.store[key]  
    def commit(self):
        while self.transactions:
            changes = self.transactions.pop(0)
            self.store.update(changes)

store = KeyValueStore()
store.set("x", 5)
store.begin()
store.set("x", 10)
store.rollback()
print(store.get("x")) 

#Prompt 3:
    #Update your implementation so multiple begin() calls create nested transaction layers.
    #Each rollback() should undo only the most recent layer.
    #Each commit() should merge all layers into the base store.

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:                        
                return layer[key]                   
        return self.store.get(key, None) 
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
            self.store.update(changes)
            print("Transactions Value", self.transactions, "Updated Stored Value", self.store)

store = KeyValueStore()
store.set("x", 1)
store.begin()
store.set("x", 2)
store.begin()
store.set("x", 3)
store.rollback()
print(store.get("x"))
store.commit()
print(store.get("x"))

#Prompt 4
    #Modify get(key) so it returns "NULL" if the key doesn’t exist in any layer.

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:                        
                return layer[key]
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Key and Value stored but not present in any transaction layer \nStored Value =")
        else:
            print("If Null, Key not stored and not present in any layer")                 
        return key_value
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
            self.store.update(changes)
            print("Transactions Value", self.transactions, "Updated Stored Value", self.store)


store = KeyValueStore()
store.set("y", 1)
store.set("x", 3)
print(store.get("y"))

store = KeyValueStore()
store.begin()
store.set("x", 3)
print(store.get("y"))

#Prompt 5 and 6:
    #Add a delete(key) command that removes a key from the active transaction or base store.
    #If rolled back, the deletion should be undone.
    #Ensure that calling rollback() or commit() when there’s no active transaction prints: NO TRANSACTION

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def deletekey(self, key):
        if self.transactions:
            if key in self.transactions[-1]:
                print("Specified key found in top layer, deleting specified key")
                self.transactions[-1][key] = None
            else:
                print("Specified key not found in top layer")
                if key in self.store:
                    print("Specified key found in base store,\nwill delete specified key from base store after commit")
                    self.transactions.append({})
                    self.transactions[-1][key] = None
                else:
                    print("Specified key not found in base store")
        else:
            print("No active transaction")
            if key in self.store:
                print("Deleting key from base store")
                del self.store[key]
            else:
                print("Key not found in base store")
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
        else:
            print("No active transaction")
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:
                if layer[key] == None:    
                    print("Key set to None, commit to delete")
                    return layer[key]
                print("Key found in active transaction, Value =")
                return layer[key]
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Key and Value stored but not present in any active transaction layer \nStored Value =")
        else:
            print("If Null, Key not stored and not present in any layer")                 
        return key_value
    def getbasevalue(self, key):                                       
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Stored Value =")
        else:
            print("If Null, Key not stored")                 
        return key_value
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
                for key, value in changes.items():
                    if value is None:
                        print(key, "deleted")
                        self.store.pop(key, None)
                    else:
                        self.store[key] = value
                        print("All transactions committed")
                        print("Updated Stored Value", self.store)
        else:
            print("No active transaction")





store = KeyValueStore()
store.set("a", 100)
store.begin()
store.deletekey("a")
print(store.get("a"))  # → NULL
store.rollback()
print(store.get("a"))  # → 100

store = KeyValueStore()
store.rollback()
store.commit()





#Prompt 7:
    #Change your commit() behaviour so it commits only the most recent transaction layer into the one below, not all layers at once.
    #I'm going to add a partial commit function instead so you can choose

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def deletekey(self, key):
        if self.transactions:
            if key in self.transactions[-1]:
                print("Specified key found in top layer, deleting specified key")
                self.transactions[-1][key] = None
            else:
                print("Specified key not found in top layer")
                if key in self.store:
                    print("Specified key found in base store,\nwill delete specified key from base store after commit")
                    self.transactions.append({})
                    self.transactions[-1][key] = None
                else:
                    print("Specified key not found in base store")
        else:
            print("No active transaction")
            if key in self.store:
                print("Deleting key from base store")
                del self.store[key]
            else:
                print("Key not found in base store")
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
        else:
            print("No active transaction")
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:
                if layer[key] == None:    
                    print("Key set to None, commit to delete")
                    return layer[key]
                print("Key found in active transaction, Value =")
                return layer[key]
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Key and Value stored but not present in any active transaction layer \nStored Value =")
        else:
            print("If Null, Key not stored and not present in any layer")                 
        return key_value
    def getbasevalue(self, key):                                       
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Stored Value =")
        else:
            print("If Null, Key not stored")                 
        return key_value
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
                for key, value in changes.items():
                    if value is None:
                        print(key, "deleted")
                        self.store.pop(key, None)
                    else:
                        self.store[key] = value
                        print("All transactions committed")
                        print("Updated Stored Value", self.store)
        else:
            print("No active transaction")        
    def partialcommit(self):
        if len(self.transactions) < 2:
            print("Not enough nested transactions to perform partial commit")
            return
        if len(self.transactions) < 1:
            print("No active transaction")
            return
        top_layer = self.transactions.pop(-1)
        next_layer = self.transactions[-1]
        for key, value in top_layer.items():
                if value is None:
                    next_layer[key] = None
                else:
                    next_layer[key] = value
        print("Top layer merged with previous layer successfully")


store = KeyValueStore()
store.begin()
store.set("x", 10)
store.begin()
store.set("y", 20)

store = KeyValueStore()
store.set("x", 1)
store.begin()
store.set("x", 2)
store.begin()
store.set("y", 99)
store.deletekey("x")
store.partialcommit()
store.rollback()
print(store.store)

#Prompt 8:
    #Add a log() function that prints the current stack of transactions (from oldest to newest), showing the pending changes at each level.

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def deletekey(self, key):
        if self.transactions:
            if key in self.transactions[-1]:
                print("Specified key found in top layer, deleting specified key")
                self.transactions[-1][key] = None
            else:
                print("Specified key not found in top layer")
                if key in self.store:
                    print("Specified key found in base store,\nwill delete specified key from base store after commit")
                    self.transactions.append({})
                    self.transactions[-1][key] = None
                else:
                    print("Specified key not found in base store")
        else:
            print("No active transaction")
            if key in self.store:
                print("Deleting key from base store")
                del self.store[key]
            else:
                print("Key not found in base store")
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
        else:
            print("No active transaction")
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:
                if layer[key] == None:    
                    print("Key set to None, commit to delete")
                    return layer[key]
                print("Key found in active transaction, Value =")
                return layer[key]
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Key and Value stored but not present in any active transaction layer \nStored Value =")
        else:
            print("If Null, Key not stored and not present in any layer")                 
        return key_value
    def getbasevalue(self, key):                                       
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Stored Value =")
        else:
            print("If Null, Key not stored")                 
        return key_value
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
                for key, value in changes.items():
                    if value is None:
                        print(key, "deleted")
                        self.store.pop(key, None)
                    else:
                        self.store[key] = value
                        print("All transactions committed")
                        print("Updated Stored Value", self.store)
        else:
            print("No active transaction")        
    def partialcommit(self):
        if len(self.transactions) < 2:
            print("Not enough nested transactions to perform partial commit")
            return
        if len(self.transactions) < 1:
            print("No active transaction")
            return
        top_layer = self.transactions.pop(-1)
        next_layer = self.transactions[-1]
        for key, value in top_layer.items():
                if value is None:
                    next_layer[key] = None
                else:
                    next_layer[key] = value
        print("Top layer merged with previous layer successfully")
    def log(self):                                       
        print("Log:")
        if self.store:
            print("\nBase store:", self.store)
        else:
            print("No stored keys or values")
        if not self.transactions:
            print("No active transactions")
            return
        print("\nActice transactions:")
        for i, layer in enumerate(self.transactions, start=1):
            print(f"    Level {i}: {layer}")

store = KeyValueStore()
store.set("a", 10)
store.begin()
store.set("b", 20)
store.begin()
store.set("a", 15)
store.deletekey("b")
store.log()

#Prompt 9:
    #Simulate two independent key-value stores (store1, store2).
    #Demonstrate that a transaction in store1 doesn’t affect store2 — i.e. each has its own stack and base store.

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def deletekey(self, key):
        if self.transactions:
            if key in self.transactions[-1]:
                print("Specified key found in top layer, deleting specified key")
                self.transactions[-1][key] = None
            else:
                print("Specified key not found in top layer")
                if key in self.store:
                    print("Specified key found in base store,\nwill delete specified key from base store after commit")
                    self.transactions.append({})
                    self.transactions[-1][key] = None
                else:
                    print("Specified key not found in base store")
        else:
            print("No active transaction")
            if key in self.store:
                print("Deleting key from base store")
                del self.store[key]
            else:
                print("Key not found in base store")
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
        else:
            print("No active transaction")
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:
                if layer[key] == None:    
                    print("Key set to None, commit to delete")
                    return layer[key]
                print("Key found in active transaction, Value =")
                return layer[key]
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Key and Value stored but not present in any active transaction layer \nStored Value =")
        else:
            print("If Null, Key not stored and not present in any layer")                 
        return key_value
    def getbasevalue(self, key):                                       
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Stored Value =")
        else:
            print("If Null, Key not stored")                 
        return key_value
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
                for key, value in changes.items():
                    if value is None:
                        print(key, "deleted")
                        self.store.pop(key, None)
                    else:
                        self.store[key] = value
                        print("All transactions committed")
                        print("Updated Stored Value", self.store)
        else:
            print("No active transaction")        
    def partialcommit(self):
        if len(self.transactions) < 2:
            print("Not enough nested transactions to perform partial commit")
            return
        if len(self.transactions) < 1:
            print("No active transaction")
            return
        top_layer = self.transactions.pop(-1)
        next_layer = self.transactions[-1]
        for key, value in top_layer.items():
                if value is None:
                    next_layer[key] = None
                else:
                    next_layer[key] = value
        print("Top layer merged with previous layer successfully")
    def log(self):                                       
        print("Log:")
        if self.store:
            print("\nBase store:", self.store)
        else:
            print("No stored keys or values")
        if not self.transactions:
            print("No active transactions")
            return
        print("\nActice transactions:")
        for i, layer in enumerate(self.transactions, start=1):
            print(f"    Level {i}: {layer}")

store1 = KeyValueStore()
store2 = KeyValueStore()

store1.set("x", 5)
store2.set("x", 99)
store1.begin()
store1.set("x", 10)
store1.get("x")     #should be 10
store2.get("x")     #should be 99
store1.rollback()
store1.get("x")     #shoudl be 5


#Prompt 10:
    #Design your store so that every transaction keeps a checksum or version number that increments on each commit.
    #Print the version number with every get() call.

class KeyValueStore:
    def __init__(self):
        self.store = {}
        self.transactions = []
        self.version = 0
    def begin(self):                    
        self.transactions.append({})
    def set(self, key, value):          
        if self.transactions:
            self.transactions[-1][key] = value
        else:
            self.store[key] = value
    def deletekey(self, key):
        if self.transactions:
            if key in self.transactions[-1]:
                print("Specified key found in top layer, deleting specified key")
                self.transactions[-1][key] = None
            else:
                print("Specified key not found in top layer")
                if key in self.store:
                    print("Specified key found in base store,\nwill delete specified key from base store after commit")
                    self.transactions.append({})
                    self.transactions[-1][key] = None
                else:
                    print("Specified key not found in base store")
        else:
            print("No active transaction")
            if key in self.store:
                print("Deleting key from base store")
                del self.store[key]
            else:
                print("Key not found in base store")
    def rollback(self):
        if self.transactions:
            self.transactions.pop(-1)
        else:
            print("No active transaction")
    def get(self, key):                                       
        for layer in reversed(self.transactions):      
            if key in layer:
                if layer[key] == None:    
                    print("Key set to None, commit to delete")
                    return layer[key], f"Version: {self.version}"
                print("Key found in active transaction, Value =")
                return layer[key], f"Version: {self.version}"
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Key and Value stored but not present in any active transaction layer \nStored Value =")
        else:
            print("If Null, Key not stored and not present in any layer")                 
        return key_value, f"Version: {self.version}"
    def getbasevalue(self, key):                                       
        key_value = self.store.get(key, "Null")
        if key_value != "Null":
            print("Stored Value =")
        else:
            print("If Null, Key not stored")                 
        return key_value, f"Version: {self.version}"
    def commit(self):
        changes = self.store
        if self.transactions:
            while self.transactions:
                changes = self.transactions.pop(0)
                for key, value in changes.items():
                    if value is None:
                        print(key, "deleted")
                        self.store.pop(key, None)
                    else:
                        self.store[key] = value
                        self.version += 1
                        print(f"All transactions committed, store version now: {self.version}")
                        print("Updated Stored Value", self.store)
        else:
            print("No active transaction")        
    def partialcommit(self):
        if len(self.transactions) < 2:
            print("Not enough nested transactions to perform partial commit")
            return
        if len(self.transactions) < 1:
            print("No active transaction")
            return
        top_layer = self.transactions.pop(-1)
        next_layer = self.transactions[-1]
        for key, value in top_layer.items():
                if value is None:
                    next_layer[key] = None
                else:
                    next_layer[key] = value
        print("Top layer merged with previous layer successfully")
    def log(self):                                       
        print("Log:")
        if self.store:
            print("\nBase store:", self.store)
        else:
            print("No stored keys or values")
        if not self.transactions:
            print("No active transactions")
            return
        print("\nActice transactions:")
        for i, layer in enumerate(self.transactions, start=1):
            print(f"    Level {i}: {layer}")

store = KeyValueStore()

store.set("a", 5)
store.begin()
store.set("a", 10)
store.commit()

store.get("a")
store.set("b", 7)
store.begin()
store.set("b", 20)
store.commit()
store.get("b")