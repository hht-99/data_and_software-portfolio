
#Full script written to deal with transactional data management
    #Layered transaction handling
    #Rollback and commit support
    #Tombstone-based deletion
    #Version tracking
    #Nested partial commits
    #Logging of internal state


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
        print("\nActive transactions:")
        for i, layer in enumerate(self.transactions, start=1):
            print(f"    Level {i}: {layer}")